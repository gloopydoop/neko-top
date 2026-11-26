program rom_projection_pairing_diagnostic
  !> No-KSP pairing probe for the actual ROM continuity and correction maps.
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use fluid_pnpn, only: fluid_pnpn_t
  use field, only: field_t
  use field_math, only: field_rzero, field_col3
  use operators, only: grad, ortho
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2, glsc2, glsc3, glmin, glmax
  use num_types, only: rp
  use comm, only: pe_rank
  use rom_projection_operator, only: rom_projection_ax_t
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
  use utils, only: neko_error
  implicit none

  integer, parameter :: n_probe = 4
  type(simulation_t) :: simulation
  type(json_file) :: parameters
  type(fluid_pnpn_t), pointer :: fluid
  type(rom_projection_ax_t) :: schur
  type(field_t) :: q, sq, cq, pzero, pwork
  type(field_t) :: gx, gy, gz, fx, fy, fz, uzero, vzero, wzero
  type(field_t) :: vu, vv, vw, cv
  character(len=1024) :: case_file, csv_file
  integer :: argc, n, iq, iv, i, unit
  logical :: found
  real(kind=rp) :: xmin, xmax, ymin, ymax, zmin, zmax, x, y, z
  real(kind=rp) :: qnorm, sqnorm, cnorm, gnorm, vnorm
  real(kind=rp) :: lhs, rhs, pair_scale, pair_error
  real(kind=rp) :: schur_quad, correction_mass_sq, schur_mass_error
  real(kind=rp) :: cg_operator_error, cg_operator_relative
  real(kind=rp) :: rho_min, rho_max
  type(c_ptr) :: event

  argc = command_argument_count()
  if (argc < 1) call neko_error( &
       'Usage: ROM_rugby_ball_projection_pairing_diagnostic CASE [CSV]')
  call get_command_argument(1, case_file)
  if (argc >= 2) then
     call get_command_argument(2, csv_file)
  else
     csv_file = 'projection_pairing_diagnostic.csv'
  end if

  call neko_init()
  call neko_top_register_types()
  parameters = json_read_file(trim(case_file))
  call parameters%update('unsteady', .false., found)
  call user_setup(simulation%neko_case%user)
  call simulation%init(parameters)
  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('Projection pairing diagnostic requires Pn/Pn fluid')
  end select
  n = fluid%u%dof%size()

  call q%init(fluid%dm_Xh, 'pairing_q')
  call sq%init(fluid%dm_Xh, 'pairing_s_q')
  call cq%init(fluid%dm_Xh, 'pairing_c_gq')
  call pzero%init(fluid%dm_Xh, 'pairing_pzero')
  call pwork%init(fluid%dm_Xh, 'pairing_pwork')
  call gx%init(fluid%dm_Xh, 'pairing_gx')
  call gy%init(fluid%dm_Xh, 'pairing_gy')
  call gz%init(fluid%dm_Xh, 'pairing_gz')
  call fx%init(fluid%dm_Xh, 'pairing_fx')
  call fy%init(fluid%dm_Xh, 'pairing_fy')
  call fz%init(fluid%dm_Xh, 'pairing_fz')
  call uzero%init(fluid%dm_Xh, 'pairing_uzero')
  call vzero%init(fluid%dm_Xh, 'pairing_vzero')
  call wzero%init(fluid%dm_Xh, 'pairing_wzero')
  call vu%init(fluid%dm_Xh, 'pairing_vu')
  call vv%init(fluid%dm_Xh, 'pairing_vv')
  call vw%init(fluid%dm_Xh, 'pairing_vw')
  call cv%init(fluid%dm_Xh, 'pairing_cv')
  call field_rzero(uzero); call field_rzero(vzero); call field_rzero(wzero)
  call field_rzero(pzero)
  xmin = glmin(fluid%u%dof%x%x, n); xmax = glmax(fluid%u%dof%x%x, n)
  ymin = glmin(fluid%u%dof%y%x, n); ymax = glmax(fluid%u%dof%y%x, n)
  zmin = glmin(fluid%u%dof%z%x, n); zmax = glmax(fluid%u%dof%z%x, n)
  rho_min = glmin(fluid%rho%x, n); rho_max = glmax(fluid%rho%x, n)

  call schur%init(fluid, q, pzero, pwork, gx, gy, gz, fx, fy, fz, &
       uzero, vzero, wzero, mass_weighted_velocity_average = .true.)

  if (pe_rank == 0) then
     open(newunit=unit, file=trim(csv_file), status='replace', action='write')
     write(unit, '(A)') 'probe,velocity_probe,lhs_q_Cv,rhs_Gq_v,' // &
          'pairing_abs_error,pairing_norm_product_error,schur_q_Sq,' // &
          'correction_mass_sq,schur_mass_abs_error,schur_mass_relative_error,' // &
          'schur_vs_CG_abs_error,schur_vs_CG_relative_error,' // &
          'q_norm,Schur_norm,Cv_norm,Gq_norm,v_norm'
  end if

  ! Four restricted pressure probes and three independent no-slip velocities.
  do iq = 1, n_probe
     call make_pressure_probe(iq)
     call continuous_scalar(q)
     call fluid%bcs_prs_projector%apply(q%x, n)

     ! Actual S = C G application, with the same ROM-local mass-weighted
     ! gather/average and velocity BC path as rom_project_velocity_column.
     call schur%compute(sq%x, q%x, fluid%c_Xh, fluid%msh, fluid%Xh)
     call fluid%gs_Xh%op(sq%x, n, GS_OP_ADD)
     call fluid%bcs_prs_projector%apply(sq%x, n)
     call q_mass_norms(qnorm, sqnorm)
     schur_quad = glsc3(q%x, fluid%c_Xh%mult, sq%x, n)

     ! Independently materialize Gq using the correction-map ingredients.
     call grad(gx%x, gy%x, gz%x, q%x, fluid%c_Xh)
     call col2(gx%x, fluid%c_Xh%B, n)
     call col2(gy%x, fluid%c_Xh%B, n)
     call col2(gz%x, fluid%c_Xh%B, n)
     call fluid%gs_Xh%op(gx%x, gy%x, gz%x, n, GS_OP_ADD)
     call col2(gx%x, fluid%c_Xh%Binv, n)
     call col2(gy%x, fluid%c_Xh%Binv, n)
     call col2(gz%x, fluid%c_Xh%Binv, n)
     call fluid%bcs_vel_projector%apply(gx%x, gy%x, gz%x, n)
     gnorm = velocity_mass_norm(gx, gy, gz)
     correction_mass_sq = gnorm * gnorm
     schur_mass_error = abs(schur_quad - correction_mass_sq)
     call apply_continuity(gx, gy, gz, cq)
     call pressure_mass_norm(cq, cnorm)
     cg_operator_error = sqrt(max(0.0_rp, glsc3(sq%x-cq%x, &
          fluid%c_Xh%mult, sq%x-cq%x, n)))
     cg_operator_relative = cg_operator_error / max(sqnorm, cnorm, tiny(1.0_rp))

     do iv = 1, 3
        call make_velocity_probe(iv)
        call apply_continuity(vu, vv, vw, cv)
        call pressure_mass_norm(cv, cnorm)
        lhs = glsc3(q%x, fluid%c_Xh%mult, cv%x, n)
        rhs = glsc3(gx%x, fluid%c_Xh%B, vu%x, n) + &
             glsc3(gy%x, fluid%c_Xh%B, vv%x, n) + &
             glsc3(gz%x, fluid%c_Xh%B, vw%x, n)
        pair_scale = sqrt(max(0.0_rp, &
             glsc3(q%x, fluid%c_Xh%mult, q%x, n) * &
             glsc3(cv%x, fluid%c_Xh%mult, cv%x, n)))
        pair_error = abs(lhs - rhs) / max(pair_scale, tiny(1.0_rp))
        vnorm = velocity_mass_norm(vu, vv, vw)
        if (pe_rank == 0) then
           write(unit, '(I0,",",I0,15(",",ES24.16))') iq, iv, lhs, rhs, &
                abs(lhs - rhs), pair_error, schur_quad, correction_mass_sq, &
                schur_mass_error, schur_mass_error / &
                max(abs(schur_quad), correction_mass_sq, tiny(1.0_rp)), &
                cg_operator_error, cg_operator_relative, &
                qnorm, sqnorm, cnorm, gnorm, vnorm
        end if
     end do
  end do

  if (pe_rank == 0) close(unit)
  call schur%free()
  call q%free(); call sq%free(); call cq%free(); call pzero%free()
  call pwork%free(); call gx%free(); call gy%free(); call gz%free()
  call fx%free(); call fy%free(); call fz%free(); call uzero%free()
  call vzero%free(); call wzero%free(); call vu%free(); call vv%free()
  call vw%free(); call cv%free()
  if (pe_rank == 0) then
     write(*, '(A,2(1X,ES12.4))') 'rho min/max = ', rho_min, rho_max
     write(*, '(A,A)') 'pairing CSV: ', trim(csv_file)
  end if
  call neko_finalize()

contains

  subroutine continuous_scalar(a)
    type(field_t), intent(inout) :: a
    call fluid%gs_Xh%op(a%x, n, GS_OP_ADD)
    call col2(a%x, fluid%c_Xh%mult, n)
  end subroutine continuous_scalar

  subroutine continuous_vector(a, b, c)
    type(field_t), intent(inout) :: a, b, c
    call fluid%gs_Xh%op(a%x, b%x, c%x, n, GS_OP_ADD)
    call col2(a%x, fluid%c_Xh%mult, n)
    call col2(b%x, fluid%c_Xh%mult, n)
    call col2(c%x, fluid%c_Xh%mult, n)
  end subroutine continuous_vector

  subroutine make_pressure_probe(k)
    integer, intent(in) :: k
    real(kind=rp) :: a, b, c
    integer :: j
    do j = 1, n
       a = normalize(fluid%u%dof%x%x(j,1,1,1), xmin, xmax)
       b = normalize(fluid%u%dof%y%x(j,1,1,1), ymin, ymax)
       c = normalize(fluid%u%dof%z%x(j,1,1,1), zmin, zmax)
       select case (k)
       case (1); q%x(j,1,1,1) = 1.0_rp
       case (2); q%x(j,1,1,1) = a
       case (3); q%x(j,1,1,1) = a*b
       case default; q%x(j,1,1,1) = a*a + 0.7_rp*b*b + 0.3_rp*c
       end select
    end do
  end subroutine make_pressure_probe

  subroutine make_velocity_probe(k)
    integer, intent(in) :: k
    real(kind=rp) :: a, b, c, bubble
    integer :: j
    do j = 1, n
       a = normalize(fluid%u%dof%x%x(j,1,1,1), xmin, xmax)
       b = normalize(fluid%u%dof%y%x(j,1,1,1), ymin, ymax)
       c = normalize(fluid%u%dof%z%x(j,1,1,1), zmin, zmax)
       bubble = (a*(1.0_rp-a)*b*(1.0_rp-b)*c*(1.0_rp-c))**2
       select case (k)
       case (1)
          vu%x(j,1,1,1) = bubble
          vv%x(j,1,1,1) = bubble*(a-b)
          vw%x(j,1,1,1) = bubble*(b-c)
       case (2)
          vu%x(j,1,1,1) = bubble*(1.0_rp+a)
          vv%x(j,1,1,1) = bubble*(1.0_rp+b*b)
          vw%x(j,1,1,1) = bubble*(1.0_rp+c*c)
       case default
          vu%x(j,1,1,1) = bubble*(a+b)
          vv%x(j,1,1,1) = bubble*(b+c)
          vw%x(j,1,1,1) = bubble*(c+a)
       end select
    end do
    call continuous_vector(vu, vv, vw)
    call fluid%bcs_vel_projector%apply(vu%x, vv%x, vw%x, n)
  end subroutine make_velocity_probe

  subroutine apply_continuity(u, v, w, out)
    type(field_t), intent(inout) :: u, v, w, out
    event = C_NULL_PTR
    call field_col3(fx, u, fluid%rho)
    call field_col3(fy, v, fluid%rho)
    call field_col3(fz, w, fluid%rho)
    call col2(fx%x, fluid%c_Xh%B, n)
    call col2(fy%x, fluid%c_Xh%B, n)
    call col2(fz%x, fluid%c_Xh%B, n)
    call field_rzero(pzero)
    call fluid%prs_res%compute(pzero, pwork, u, v, w, uzero, vzero, wzero, &
         fx, fy, fz, fluid%c_Xh, fluid%gs_Xh, fluid%bc_prs_surface, &
         fluid%bc_sym_surface, fluid%Ax_prs, 0.0_rp, 1.0_rp, uzero, &
         fluid%rho, event)
    if (.not. fluid%prs_dirichlet) call ortho(pwork%x, fluid%glb_n_points, n)
    call fluid%gs_Xh%op(pwork%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(pwork%x, n)
    out%x = pwork%x
  end subroutine apply_continuity

  function normalize(value, low, high) result(res)
    real(kind=rp), intent(in) :: value, low, high
    real(kind=rp) :: res
    res = (value-low) / max(high-low, tiny(1.0_rp))
  end function normalize

  subroutine pressure_mass_norm(a, norm)
    type(field_t), intent(in) :: a
    real(kind=rp), intent(out) :: norm
    norm = sqrt(max(0.0_rp, glsc3(a%x, fluid%c_Xh%mult, a%x, n)))
  end subroutine pressure_mass_norm

  function velocity_mass_norm(a, b, c) result(norm)
    type(field_t), intent(in) :: a, b, c
    real(kind=rp) :: norm
    norm = sqrt(max(0.0_rp, glsc3(a%x, fluid%c_Xh%B, a%x, n) + &
         glsc3(b%x, fluid%c_Xh%B, b%x, n) + &
         glsc3(c%x, fluid%c_Xh%B, c%x, n)))
  end function velocity_mass_norm

  subroutine q_mass_norms(q_value, sq_value)
    real(kind=rp), intent(out) :: q_value, sq_value
    q_value = sqrt(max(0.0_rp, glsc3(q%x, fluid%c_Xh%mult, q%x, n)))
    sq_value = sqrt(max(0.0_rp, glsc3(sq%x, fluid%c_Xh%mult, sq%x, n)))
  end subroutine q_mass_norms

end program rom_projection_pairing_diagnostic
