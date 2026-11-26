program rom_curl_div_compatibility_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use fluid_pnpn, only: fluid_pnpn_t
  use field, only: field_t
  use field_math, only: field_rzero, field_col3
  use operators, only: grad, div, ortho
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2, glsc3, glmin, glmax
  use num_types, only: rp
  use comm, only: pe_rank
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
  use utils, only: neko_error
  implicit none

  type(simulation_t) :: simulation
  type(json_file) :: parameters
  type(fluid_pnpn_t), pointer :: fluid
  type(field_t) :: psi, u, v, w, ux, uy, uz, divergence
  type(field_t) :: p, p_res, zero_u, zero_v, zero_w, fx, fy, fz
  character(len=1024) :: case_file, output_file
  integer :: argc, n, i, unit
  logical :: found
  real(kind=rp) :: strong_local, strong_gs, strong_bc
  real(kind=rp) :: weak_gs, weak_bc, bc_gs
  real(kind=rp) :: velocity_gs, velocity_bc, x, y, z
  real(kind=rp) :: xmin, xmax, ymin, ymax, zmin, zmax
  type(c_ptr) :: event

  argc = command_argument_count()
  if (argc < 1) then
     call neko_error('Usage: curl_div_compatibility_diagnostic CASE [CSV]')
  end if
  call get_command_argument(1, case_file)
  if (argc >= 2) then
     call get_command_argument(2, output_file)
  else
     output_file = 'curl_div_compatibility.csv'
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
     call neko_error('Curl/div compatibility diagnostic requires Pn/Pn fluid')
  end select

  n = fluid%u%dof%size()
  call psi%init(fluid%dm_Xh, 'curl_div_psi')
  call u%init(fluid%dm_Xh, 'curl_div_u')
  call v%init(fluid%dm_Xh, 'curl_div_v')
  call w%init(fluid%dm_Xh, 'curl_div_w')
  call ux%init(fluid%dm_Xh, 'curl_div_ux')
  call uy%init(fluid%dm_Xh, 'curl_div_uy')
  call uz%init(fluid%dm_Xh, 'curl_div_uz')
  call divergence%init(fluid%dm_Xh, 'curl_div_divergence')
  call p%init(fluid%dm_Xh, 'curl_div_p')
  call p_res%init(fluid%dm_Xh, 'curl_div_p_res')
  call zero_u%init(fluid%dm_Xh, 'curl_div_zero_u')
  call zero_v%init(fluid%dm_Xh, 'curl_div_zero_v')
  call zero_w%init(fluid%dm_Xh, 'curl_div_zero_w')
  call fx%init(fluid%dm_Xh, 'curl_div_fx')
  call fy%init(fluid%dm_Xh, 'curl_div_fy')
  call fz%init(fluid%dm_Xh, 'curl_div_fz')
  call field_rzero(zero_u); call field_rzero(zero_v); call field_rzero(zero_w)

  ! A smooth, nontrivial scalar potential. The 2D curl is represented by
  ! the rotated SEM gradient (d_y psi, -d_x psi, 0). In 3D this only probes
  ! the corresponding z-directed curl component; it does not establish a
  ! general vector-potential identity.
  xmin = glmin(fluid%u%dof%x%x, n); xmax = glmax(fluid%u%dof%x%x, n)
  ymin = glmin(fluid%u%dof%y%x, n); ymax = glmax(fluid%u%dof%y%x, n)
  zmin = glmin(fluid%u%dof%z%x, n); zmax = glmax(fluid%u%dof%z%x, n)
  do i = 1, n
     x = fluid%u%dof%x%x(i, 1, 1, 1)
     y = fluid%u%dof%y%x(i, 1, 1, 1)
     z = fluid%u%dof%z%x(i, 1, 1, 1)
     x = (x - xmin) / max(xmax - xmin, tiny(1.0_rp))
     y = (y - ymin) / max(ymax - ymin, tiny(1.0_rp))
     z = (z - zmin) / max(zmax - zmin, tiny(1.0_rp))
     psi%x(i, 1, 1, 1) = (x * (1.0_rp - x))**2 * &
          (y * (1.0_rp - y))**2 * (z * (1.0_rp - z))**2
  end do
  call continuous_scalar(psi)
  call grad(ux%x, uy%x, uz%x, psi%x, fluid%c_Xh)
  u%x = uy%x
  v%x = -ux%x
  w%x = 0.0_rp

  ! Measure element-local derivative output before gather-scatter.
  call measure_strong(strong_local, velocity_gs)
  call continuous_vector(u, v, w)
  call measure_strong(strong_gs, velocity_gs)
  call measure_weak(weak_gs)

  ux%x = u%x
  uy%x = v%x
  uz%x = w%x
  call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
  ux%x = ux%x - u%x
  uy%x = uy%x - v%x
  uz%x = uz%x - w%x
  bc_gs = sqrt(max(0.0_rp, glsc3(ux%x, fluid%c_Xh%B, ux%x, n) + &
       glsc3(uy%x, fluid%c_Xh%B, uy%x, n) + &
       glsc3(uz%x, fluid%c_Xh%B, uz%x, n)))
  velocity_bc = velocity_norm(u, v, w)
  call measure_strong(strong_bc, velocity_bc)
  call measure_weak(weak_bc)

  if (pe_rank == 0) then
     open(newunit=unit, file=trim(output_file), status='replace', action='write')
     write(unit, '(A)') 'stage,strong_divergence,velocity_mass_norm,' // &
          'relative_strong_divergence,native_weak_continuity,' // &
          'weak_over_velocity_norm,bc_removed_mass_norm,bc_over_velocity_norm'
     write(unit, '(A,7(",",ES24.16))') 'element_local', strong_local, &
          velocity_gs, strong_local / max(velocity_gs, tiny(1.0_rp)), &
          -1.0_rp, -1.0_rp, -1.0_rp, -1.0_rp
     write(unit, '(A,7(",",ES24.16))') 'gather_scatter', strong_gs, &
          velocity_gs, strong_gs / max(velocity_gs, tiny(1.0_rp)), &
          weak_gs, weak_gs / max(velocity_gs, tiny(1.0_rp)), -1.0_rp, -1.0_rp
     write(unit, '(A,7(",",ES24.16))') 'after_velocity_bc', strong_bc, &
          velocity_bc, strong_bc / max(velocity_bc, tiny(1.0_rp)), &
          weak_bc, weak_bc / max(velocity_bc, tiny(1.0_rp)), bc_gs, &
          bc_gs / max(velocity_gs, tiny(1.0_rp))
     close(unit)
     write(*, '(A,ES12.4)') 'curl-div local strong/velocity = ', &
          strong_local / max(velocity_gs, tiny(1.0_rp))
     write(*, '(A,ES12.4)') 'curl-div assembled strong/velocity = ', &
          strong_gs / max(velocity_gs, tiny(1.0_rp))
     write(*, '(A,ES12.4)') 'curl-div assembled native weak/velocity = ', &
          weak_gs / max(velocity_gs, tiny(1.0_rp))
     write(*, '(A,ES12.4)') 'after BC strong/velocity = ', &
          strong_bc / max(velocity_bc, tiny(1.0_rp))
     write(*, '(A,ES12.4)') 'after BC native weak/velocity = ', &
          weak_bc / max(velocity_bc, tiny(1.0_rp))
     write(*, '(A,ES12.4)') 'BC removed mass norm/velocity = ', &
          bc_gs / max(velocity_gs, tiny(1.0_rp))
     write(*, '(A,A)') 'CSV evidence: ', trim(output_file)
  end if

  call psi%free(); call u%free(); call v%free(); call w%free()
  call ux%free(); call uy%free(); call uz%free(); call divergence%free()
  call p%free(); call p_res%free(); call zero_u%free(); call zero_v%free()
  call zero_w%free(); call fx%free(); call fy%free(); call fz%free()
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

  real(kind=rp) function velocity_norm(a, b, c)
    type(field_t), intent(in) :: a, b, c
    velocity_norm = sqrt(max(0.0_rp, glsc3(a%x, fluid%c_Xh%B, a%x, n) + &
         glsc3(b%x, fluid%c_Xh%B, b%x, n) + &
         glsc3(c%x, fluid%c_Xh%B, c%x, n)))
  end function velocity_norm

  subroutine measure_strong(value, norm_u)
    real(kind=rp), intent(out) :: value, norm_u
    call div(divergence%x, u%x, v%x, w%x, fluid%c_Xh)
    value = sqrt(max(0.0_rp, glsc3(divergence%x, fluid%c_Xh%B, &
         divergence%x, n)))
    norm_u = velocity_norm(u, v, w)
  end subroutine measure_strong

  subroutine measure_weak(value)
    real(kind=rp), intent(out) :: value
    event = C_NULL_PTR
    call field_col3(fx, u, fluid%rho)
    call field_col3(fy, v, fluid%rho)
    call field_col3(fz, w, fluid%rho)
    call col2(fx%x, fluid%c_Xh%B, n)
    call col2(fy%x, fluid%c_Xh%B, n)
    call col2(fz%x, fluid%c_Xh%B, n)
    call field_rzero(p)
    call field_rzero(p_res)
    call fluid%prs_res%compute(p, p_res, u, v, w, zero_u, zero_v, zero_w, &
         fx, fy, fz, fluid%c_Xh, fluid%gs_Xh, fluid%bc_prs_surface, &
         fluid%bc_sym_surface, fluid%Ax_prs, 0.0_rp, 1.0_rp, zero_u, &
         fluid%rho, event)
    if (.not. fluid%prs_dirichlet) call ortho(p_res%x, fluid%glb_n_points, n)
    call fluid%gs_Xh%op(p_res%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(p_res%x, n)
    value = sqrt(max(0.0_rp, glsc3(p_res%x, fluid%c_Xh%mult, p_res%x, n)))
  end subroutine measure_weak

end program rom_curl_div_compatibility_diagnostic
