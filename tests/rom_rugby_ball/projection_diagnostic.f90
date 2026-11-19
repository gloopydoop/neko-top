program rom_projection_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use simulation_m, only: simulation_t
  use file, only: file_t
  use fld_file_data, only: fld_file_data_t
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use fluid_pnpn, only: fluid_pnpn_t
  use field, only: field_t
  use field_math, only: field_col3, field_rzero
  use operators, only: ortho
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2, glsc3
  use krylov, only: ksp_monitor_t
  use cg, only: cg_t
  use precon, only: pc_t
  use identity, only: ident_t
  use rom_projection_operator, only: rom_projection_ax_t
  use rom_basis_handoff, only: rom_project_basis
  use comm, only: pe_rank
  use num_types, only: rp
  use utils, only: neko_error
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR

  implicit none

  type(simulation_t) :: simulation
  type(fluid_pnpn_t), pointer :: fluid
  type(file_t) :: modes_file
  type(fld_file_data_t) :: modes_data
  type(json_file) :: parameters
  real(kind=rp), allocatable :: basis(:, :)
  character(len=1024) :: case_file, modes_path
  logical :: found, run_production
  character(len=64) :: production_arg
  integer :: argc, n, n_modes, field_size, mode

  argc = command_argument_count()
  if (argc < 2) then
     call neko_error('Usage: rom_projection_diagnostic CASE MODES.fld [production]')
  end if
  call get_command_argument(1, case_file)
  call get_command_argument(2, modes_path)
  run_production = .false.
  if (argc >= 3) then
     call get_command_argument(3, production_arg)
     run_production = trim(production_arg) == 'production'
  end if

  call neko_init()
  call neko_top_register_types()

  parameters = json_read_file(trim(case_file))
  ! The diagnostic needs only the initialized spatial operators. Disable the
  ! unsteady state-recovery peer so this executable cannot enter a DNS/adjoint
  ! lifecycle or wait for Python.
  call parameters%update('unsteady', .false., found)
  call simulation%init(parameters)
  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('Projection diagnostic requires a Pn/Pn fluid')
  end select

  call modes_data%init()
  call modes_file%init(trim(modes_path))
  call modes_file%set_counter(0)
  call modes_file%read(modes_data)
  if (modes_data%n_scalars < 3 .or. mod(modes_data%n_scalars, 3) /= 0) then
     call neko_error('Projection diagnostic mode file must contain 3*N scalar fields')
  end if

  n = simulation%fluid%u%dof%size()
  field_size = modes_data%s(1)%size()
  if (field_size /= n) then
     call neko_error('Projection diagnostic mode file does not match the case partition')
  end if
  n_modes = modes_data%n_scalars / 3
  allocate(basis(3 * n, n_modes))
  do mode = 1, n_modes
     basis(1:n, mode) = modes_data%s(3 * (mode - 1) + 1)%x
     basis(n + 1:2 * n, mode) = modes_data%s(3 * (mode - 1) + 2)%x
     basis(2 * n + 1:3 * n, mode) = modes_data%s(3 * (mode - 1) + 3)%x
  end do
  call diagnose_actual_modes(fluid, basis, n, n_modes, &
       'rom_schur_mode_diagnostic.csv')
  call modes_data%free()
  call modes_file%free()

  if (run_production) then
     ! This command is the strict scientific acceptance diagnostic.  Keep it
     ! fail-closed even though the runnable engineering case opts into the
     ! explicit warning fallback.
     call parameters%update('rom.strict_admissibility', .true., found)
     call rom_project_basis(simulation, basis, n_modes)
     if (n_modes < 1) call neko_error('Projection diagnostic returned no modes')
  end if

  if (allocated(basis)) deallocate(basis)
  call simulation%free()
  call neko_finalize()

contains

  !> Compare constrained and native pressure solves on actual POD modes.
  subroutine diagnose_actual_modes(fluid, basis, n, n_modes, diagnostics_file)
    type(fluid_pnpn_t), intent(inout), target :: fluid
    real(kind=rp), intent(in) :: basis(:, :)
    integer, intent(in) :: n, n_modes
    character(len=*), intent(in) :: diagnostics_file
    type(field_t) :: p_work, p_rhs, p_op, q_op, phi
    type(field_t) :: u, v, w, u_zero, v_zero, w_zero
    type(field_t) :: f_x, f_y, f_z, grad_x, grad_y, grad_z
    type(rom_projection_ax_t) :: projection_ax
    type(ksp_monitor_t) :: result_hsmg, result_identity, result_native
    type(ksp_monitor_t) :: result_cg
    type(cg_t) :: cg_solver
    type(ident_t), target :: identity_pc
    class(pc_t), pointer :: saved_pc
    real(kind=rp), allocatable :: rhs(:), residual_hsmg(:), residual_identity(:)
    real(kind=rp), allocatable :: residual_native(:)
    real(kind=rp) :: saved_abs, saved_rel
    integer :: saved_max, mode, unit, niter_limit
    real(kind=rp) :: rhs_norm, res_hsmg, res_identity, res_cg, res_native
    type(c_ptr) :: event

    call p_work%init(fluid%dm_Xh, 'actual_mode_diag_p_work')
    call p_rhs%init(fluid%dm_Xh, 'actual_mode_diag_p_rhs')
    call p_op%init(fluid%dm_Xh, 'actual_mode_diag_p_operator')
    call q_op%init(fluid%dm_Xh, 'actual_mode_diag_q_operator')
    call phi%init(fluid%dm_Xh, 'actual_mode_diag_phi')
    call u%init(fluid%dm_Xh, 'actual_mode_diag_u')
    call v%init(fluid%dm_Xh, 'actual_mode_diag_v')
    call w%init(fluid%dm_Xh, 'actual_mode_diag_w')
    call u_zero%init(fluid%dm_Xh, 'actual_mode_diag_u_zero')
    call v_zero%init(fluid%dm_Xh, 'actual_mode_diag_v_zero')
    call w_zero%init(fluid%dm_Xh, 'actual_mode_diag_w_zero')
    call f_x%init(fluid%dm_Xh, 'actual_mode_diag_f_x')
    call f_y%init(fluid%dm_Xh, 'actual_mode_diag_f_y')
    call f_z%init(fluid%dm_Xh, 'actual_mode_diag_f_z')
    call grad_x%init(fluid%dm_Xh, 'actual_mode_diag_grad_x')
    call grad_y%init(fluid%dm_Xh, 'actual_mode_diag_grad_y')
    call grad_z%init(fluid%dm_Xh, 'actual_mode_diag_grad_z')
    call field_rzero(u_zero); call field_rzero(v_zero); call field_rzero(w_zero)
    call projection_ax%init(fluid, q_op, p_work, p_op, grad_x, grad_y, grad_z, &
         f_x, f_y, f_z, u_zero, v_zero, w_zero)
    allocate(rhs(n), residual_hsmg(n), residual_identity(n), residual_native(n))
    saved_abs = fluid%ksp_prs%abs_tol
    saved_rel = fluid%ksp_prs%rel_tol
    saved_max = fluid%ksp_prs%max_iter
    niter_limit = max(saved_max, 200)
    if (pe_rank == 0) open(newunit=unit, file=trim(diagnostics_file), &
         status='replace', action='write')
    if (pe_rank == 0) write(unit, '(A)') &
         'mode,rhs_norm,hsmg_residual,hsmg_relative,hsmg_iterations,hsmg_converged,' // &
         'identity_residual,identity_relative,identity_iterations,identity_converged,' // &
         'cg_residual,cg_relative,cg_iterations,cg_converged,' // &
         'native_residual,native_relative,native_iterations,native_converged'

    call cg_solver%init(n, niter_limit, rel_tol=0.0_rp, abs_tol=1.0e-12_rp)

    do mode = 1, n_modes
       u%x = reshape(basis(1:n, mode), shape(u%x))
       v%x = reshape(basis(n + 1:2*n, mode), shape(v%x))
       w%x = reshape(basis(2*n + 1:3*n, mode), shape(w%x))
       call fluid%gs_Xh%op(u%x, v%x, w%x, n, GS_OP_ADD)
       call col2(u%x, fluid%c_Xh%mult, n)
       call col2(v%x, fluid%c_Xh%mult, n)
       call col2(w%x, fluid%c_Xh%mult, n)
       call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
       call assemble_rhs(fluid, p_work, p_rhs, u, v, w, u_zero, v_zero, w_zero, &
            f_x, f_y, f_z, n, event)
       rhs = reshape(p_rhs%x, [n])
       rhs_norm = pressure_norm(rhs, fluid, n)

       saved_pc => fluid%ksp_prs%M
       fluid%ksp_prs%abs_tol = 1.0e-12_rp
       fluid%ksp_prs%rel_tol = 0.0_rp
       fluid%ksp_prs%max_iter = niter_limit
       call field_rzero(phi)
       call fluid%pc_prs%update()
       result_hsmg = fluid%ksp_prs%solve(projection_ax, phi, rhs, n, fluid%c_Xh, &
            fluid%bcs_prs_projector, fluid%gs_Xh, niter=niter_limit)
       call true_residual(projection_ax, phi, rhs, p_op, fluid, n, res_hsmg)

       fluid%ksp_prs%M => identity_pc
       call field_rzero(phi)
       result_identity = fluid%ksp_prs%solve(projection_ax, phi, rhs, n, fluid%c_Xh, &
            fluid%bcs_prs_projector, fluid%gs_Xh, niter=niter_limit)
       call true_residual(projection_ax, phi, rhs, p_op, fluid, n, res_identity)

       call field_rzero(phi)
       result_cg = cg_solver%solve(projection_ax, phi, rhs, n, fluid%c_Xh, &
            fluid%bcs_prs_projector, fluid%gs_Xh, niter=niter_limit)
       call true_residual(projection_ax, phi, rhs, p_op, fluid, n, res_cg)

       fluid%ksp_prs%M => saved_pc
       call field_rzero(phi)
       result_native = fluid%ksp_prs%solve(fluid%Ax_prs, phi, rhs, n, fluid%c_Xh, &
            fluid%bcs_prs_projector, fluid%gs_Xh, niter=niter_limit)
       call native_true_residual(phi, rhs, p_op, fluid, n, res_native)
       if (pe_rank == 0) write(unit, '(I0,17(",",ES24.16))') mode, rhs_norm, &
            res_hsmg, res_hsmg/max(rhs_norm,tiny(1.0_rp)), &
            real(result_hsmg%iter, rp), merge(1.0_rp,0.0_rp,result_hsmg%converged), &
            res_identity, res_identity/max(rhs_norm,tiny(1.0_rp)), &
            real(result_identity%iter, rp), merge(1.0_rp,0.0_rp,result_identity%converged), &
            res_cg, res_cg/max(rhs_norm,tiny(1.0_rp)), &
            real(result_cg%iter, rp), merge(1.0_rp,0.0_rp,result_cg%converged), &
            res_native, res_native/max(rhs_norm,tiny(1.0_rp)), &
            real(result_native%iter, rp), merge(1.0_rp,0.0_rp,result_native%converged)
       if (pe_rank == 0) write(*, '(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') &
            'actual mode ', mode, ' RHS/HSMG/identity residual = ', rhs_norm, ' / ', &
            res_hsmg, ' / ', res_identity
    end do
    call cg_solver%free()
    if (pe_rank == 0) close(unit)
    fluid%ksp_prs%abs_tol = saved_abs
    fluid%ksp_prs%rel_tol = saved_rel
    fluid%ksp_prs%max_iter = saved_max
    call p_work%free(); call p_rhs%free(); call p_op%free(); call q_op%free(); call phi%free()
    call u%free(); call v%free(); call w%free()
    call u_zero%free(); call v_zero%free(); call w_zero%free()
    call f_x%free(); call f_y%free(); call f_z%free()
    call grad_x%free(); call grad_y%free(); call grad_z%free()
    deallocate(rhs, residual_hsmg, residual_identity, residual_native)
  end subroutine diagnose_actual_modes

  !> Assemble the native weak-continuity right-hand side.
  subroutine assemble_rhs(fluid, p, p_res, u, v, w, u_zero, v_zero, w_zero, &
       f_x, f_y, f_z, n, event)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(field_t), intent(inout) :: p, p_res, u, v, w, f_x, f_y, f_z
    type(field_t), intent(in) :: u_zero, v_zero, w_zero
    integer, intent(in) :: n
    type(c_ptr), intent(inout) :: event
    call field_col3(f_x, u, fluid%rho); call field_col3(f_y, v, fluid%rho)
    call field_col3(f_z, w, fluid%rho)
    call col2(f_x%x, fluid%c_Xh%B, n); call col2(f_y%x, fluid%c_Xh%B, n)
    call col2(f_z%x, fluid%c_Xh%B, n); call field_rzero(p)
    call fluid%prs_res%compute(p, p_res, u, v, w, u_zero, v_zero, w_zero, &
         f_x, f_y, f_z, fluid%c_Xh, fluid%gs_Xh, fluid%bc_prs_surface, &
         fluid%bc_sym_surface, fluid%Ax_prs, 0.0_rp, 1.0_rp, u_zero, fluid%rho, event)
    if (.not. fluid%prs_dirichlet) call ortho(p_res%x, fluid%glb_n_points, n)
    call fluid%gs_Xh%op(p_res%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(p_res%x, n)
  end subroutine assemble_rhs

  !> Measure a true residual for a matrix-free pressure operator.
  subroutine true_residual(operator, solution, rhs, work, fluid, n, value)
    type(rom_projection_ax_t), intent(inout) :: operator
    type(field_t), intent(in) :: solution
    real(kind=rp), intent(in) :: rhs(:)
    type(field_t), intent(inout) :: work
    type(fluid_pnpn_t), intent(inout) :: fluid
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: value
    call operator%compute(work%x, solution%x, fluid%c_Xh, fluid%msh, fluid%Xh)
    call fluid%gs_Xh%op(work%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(work%x, n)
    work%x = work%x - reshape(rhs, shape(work%x))
    call fluid%bcs_prs_projector%apply(work%x, n)
    value = pressure_norm(reshape(work%x, [n]), fluid, n)
  end subroutine true_residual

  !> Measure a true residual for Neko's native pressure operator.
  subroutine native_true_residual(solution, rhs, work, fluid, n, value)
    type(field_t), intent(in) :: solution
    real(kind=rp), intent(in) :: rhs(:)
    type(field_t), intent(inout) :: work
    type(fluid_pnpn_t), intent(inout) :: fluid
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: value
    call fluid%Ax_prs%compute(work%x, solution%x, fluid%c_Xh, fluid%msh, fluid%Xh)
    call fluid%gs_Xh%op(work%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(work%x, n)
    work%x = work%x - reshape(rhs, shape(work%x))
    call fluid%bcs_prs_projector%apply(work%x, n)
    value = pressure_norm(reshape(work%x, [n]), fluid, n)
  end subroutine native_true_residual

  !> Compute the native pressure-space right-hand-side norm.
  function pressure_norm(a, fluid, n) result(value)
    real(kind=rp), intent(in) :: a(:)
    type(fluid_pnpn_t), intent(in) :: fluid
    integer, intent(in) :: n
    real(kind=rp) :: value
    value = sqrt(max(0.0_rp, glsc3(a, fluid%c_Xh%mult, a, n)))
  end function pressure_norm

end program rom_projection_diagnostic
