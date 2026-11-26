program rom_projection_algebra_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use fluid_pnpn, only: fluid_pnpn_t
  use field, only: field_t
  use field_math, only: field_add2s2, field_col3, field_rzero, field_copy
  use operators, only: grad, div, ortho
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2, glsc3
  use krylov, only: ksp_monitor_t
  use num_types, only: rp
  use comm, only: pe_rank
  use logger, only: neko_log
  use utils, only: neko_error
  use rom_projection_operator, only: rom_projection_ax_t
  use rom_basis_handoff, only: rom_project_velocity_column
  use rom_basis_orthogonalization, only: rom_mass_orthonormalize
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR

  implicit none

  type(simulation_t) :: simulation
  type(json_file) :: parameters
  character(len=1024) :: case_file
  character(len=1024) :: diagnostics_file
  type(field_t) :: phi, phi_solved, ax_phi, ax_work, p_work, p_rhs
  type(field_t) :: p_op, q_op
  type(field_t) :: u, v, w, u_zero, v_zero, w_zero
  type(field_t) :: f_x, f_y, f_z, grad_x, grad_y, grad_z
  type(field_t) :: divergence, bc_u, bc_v, bc_w
  type(field_t) :: c_before, c_after
  type(field_t) :: psi, psi_two, sol_u, sol_v, sol_w
  type(field_t) :: sol_two_u, sol_two_v, sol_two_w
  type(field_t) :: mixed_out_u, mixed_out_v, mixed_out_w
  type(ksp_monitor_t) :: ksp_result
  type(rom_projection_ax_t) :: projection_ax
  type(fluid_pnpn_t), pointer :: fluid
  real(kind=rp) :: weak_before, weak_after, strong_before, strong_after
  real(kind=rp) :: bc_before, bc_after, mass_before, mass_after
  real(kind=rp) :: pressure_residual, pressure_relative
  real(kind=rp) :: constrained_residual, constrained_relative
  real(kind=rp) :: identity_plus, identity_minus, identity_error
  real(kind=rp) :: schur_error, schur_relative
  real(kind=rp) :: identity_relative, correction_sign
  real(kind=rp) :: pressure_rhs_norm, correction_norm
  real(kind=rp) :: saved_abs_tol, saved_rel_tol
  real(kind=rp) :: x, mass_change, weak_after_normalized
  real(kind=rp) :: pure_weak_after, pure_weak_after_normalized
  real(kind=rp) :: mixed_weak_before, mixed_weak_after
  real(kind=rp) :: mixed_strong_before, mixed_strong_after
  real(kind=rp) :: mixed_bc_before, mixed_bc_after
  real(kind=rp) :: mixed_mass_before, mixed_mass_after
  real(kind=rp) :: mixed_mass_change, mixed_weak_normalized
  real(kind=rp) :: native_column_weak, native_column_strong
  real(kind=rp) :: native_column_bc, native_column_mass
  real(kind=rp) :: native_column_linear_residual
  integer :: native_column_iterations
  logical :: native_column_converged
  real(kind=rp) :: mixed_preservation_error, mixed_idempotence
  real(kind=rp) :: correction_orthogonality, basis_gram_error
  real(kind=rp) :: mixed_rhs_norm, mixed_schur_residual, mixed_schur_relative
  real(kind=rp) :: sol_norm, correction_mass_norm, difference_mass_norm
  real(kind=rp) :: schur_symmetry_max, native_symmetry_max
  real(kind=rp) :: schur_quadratic_min, native_quadratic_min
  real(kind=rp) :: schur_probe_rank_indicator, native_probe_rank_indicator
  real(kind=rp) :: schur_constant_norm, native_constant_norm
  real(kind=rp) :: rhs_compatibility, adjoint_lhs, adjoint_rhs
  real(kind=rp) :: adjoint_relative_error, probe_dot_ij, probe_dot_ji
  integer, parameter :: n_probe = 5
  integer, parameter :: n_pressure_cases = 6
  integer, parameter :: pressure_case_max(n_pressure_cases) = [50, 100, 200, 400, 800, 1600]
  real(kind=rp), allocatable :: pressure_probes(:, :), schur_probes(:, :)
  real(kind=rp), allocatable :: native_probes(:, :), probe_s_norm(:)
  real(kind=rp), allocatable :: probe_a_norm(:), probe_s_a_relative(:)
  real(kind=rp), allocatable :: probe_s_quadratic(:), probe_a_quadratic(:)
  real(kind=rp) :: pressure_case_linear(n_pressure_cases)
  real(kind=rp) :: pressure_case_weak(n_pressure_cases)
  integer :: pressure_case_iterations(n_pressure_cases)
  logical :: pressure_case_converged(n_pressure_cases)
  real(kind=rp), allocatable :: mass_weights(:), mixed_basis(:, :)
  real(kind=rp), allocatable :: basis_norm(:), basis_cross(:)
  integer, allocatable :: mixed_final_index(:)
  integer :: mixed_accepted
  integer :: argc, n, i, unit
  integer :: saved_max_iter
  logical :: found, mixed_basis_ok
  type(c_ptr) :: event

  argc = command_argument_count()
  if (argc < 1) then
     call neko_error('Usage: ROM_rugby_ball_projection_algebra_diagnostic CASE [CSV]')
  end if
  call get_command_argument(1, case_file)
  if (argc >= 2) then
     call get_command_argument(2, diagnostics_file)
  else
     diagnostics_file = 'projection_algebra_diagnostic.csv'
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
     call neko_error('Projection algebra diagnostic requires a Pn/Pn fluid')
  end select

  n = fluid%u%dof%size()
  call phi%init(fluid%dm_Xh, 'projection_algebra_phi')
  call phi_solved%init(fluid%dm_Xh, 'projection_algebra_phi_solved')
  call ax_phi%init(fluid%dm_Xh, 'projection_algebra_ax_phi')
  call ax_work%init(fluid%dm_Xh, 'projection_algebra_ax_work')
  call p_work%init(fluid%dm_Xh, 'projection_algebra_p_work')
  call p_rhs%init(fluid%dm_Xh, 'projection_algebra_p_rhs')
  call p_op%init(fluid%dm_Xh, 'projection_algebra_p_operator')
  call q_op%init(fluid%dm_Xh, 'projection_algebra_q_operator')
  call u%init(fluid%dm_Xh, 'projection_algebra_u')
  call v%init(fluid%dm_Xh, 'projection_algebra_v')
  call w%init(fluid%dm_Xh, 'projection_algebra_w')
  call u_zero%init(fluid%dm_Xh, 'projection_algebra_u_zero')
  call v_zero%init(fluid%dm_Xh, 'projection_algebra_v_zero')
  call w_zero%init(fluid%dm_Xh, 'projection_algebra_w_zero')
  call f_x%init(fluid%dm_Xh, 'projection_algebra_f_x')
  call f_y%init(fluid%dm_Xh, 'projection_algebra_f_y')
  call f_z%init(fluid%dm_Xh, 'projection_algebra_f_z')
  call grad_x%init(fluid%dm_Xh, 'projection_algebra_grad_x')
  call grad_y%init(fluid%dm_Xh, 'projection_algebra_grad_y')
  call grad_z%init(fluid%dm_Xh, 'projection_algebra_grad_z')
  call divergence%init(fluid%dm_Xh, 'projection_algebra_divergence')
  call bc_u%init(fluid%dm_Xh, 'projection_algebra_bc_u')
  call bc_v%init(fluid%dm_Xh, 'projection_algebra_bc_v')
  call bc_w%init(fluid%dm_Xh, 'projection_algebra_bc_w')
  call c_before%init(fluid%dm_Xh, 'projection_algebra_c_before')
  call c_after%init(fluid%dm_Xh, 'projection_algebra_c_after')
  call psi%init(fluid%dm_Xh, 'projection_algebra_streamfunction')
  call psi_two%init(fluid%dm_Xh, 'projection_algebra_streamfunction_two')
  call sol_u%init(fluid%dm_Xh, 'projection_algebra_solenoidal_u')
  call sol_v%init(fluid%dm_Xh, 'projection_algebra_solenoidal_v')
  call sol_w%init(fluid%dm_Xh, 'projection_algebra_solenoidal_w')
  call sol_two_u%init(fluid%dm_Xh, 'projection_algebra_solenoidal_two_u')
  call sol_two_v%init(fluid%dm_Xh, 'projection_algebra_solenoidal_two_v')
  call sol_two_w%init(fluid%dm_Xh, 'projection_algebra_solenoidal_two_w')
  call mixed_out_u%init(fluid%dm_Xh, 'projection_algebra_mixed_out_u')
  call mixed_out_v%init(fluid%dm_Xh, 'projection_algebra_mixed_out_v')
  call mixed_out_w%init(fluid%dm_Xh, 'projection_algebra_mixed_out_w')

  call field_rzero(u_zero)
  call field_rzero(v_zero)
  call field_rzero(w_zero)
  call field_rzero(p_work)
  call field_rzero(phi)
  call field_rzero(psi)
  call field_rzero(psi_two)

  ! phi is zero, with zero normal derivative, on x=0 and x=2.  The other
  ! velocity components are zero, so the manufactured velocity satisfies the
  ! homogeneous forward velocity restriction before it is projected.
  do i = 1, n
     x = fluid%u%dof%x%x(i, 1, 1, 1)
     phi%x(i, 1, 1, 1) = x * x * (2.0_rp - x) * (2.0_rp - x)
  end do
  call continuous_scalar(fluid, phi, n)
  call fluid%bcs_prs_projector%apply(phi%x, n)

  ! This is the SEM strong derivative used by Neko's pressure Poisson
  ! Schur complement.  It is deliberately not replaced by an independently
  ! assembled strong-form Laplacian.
  call grad(u%x, v%x, w%x, phi%x, fluid%c_Xh)
  call continuous_vector(fluid, u, v, w, n)
  call rom_bc_residual(fluid, u, v, w, bc_u, bc_v, bc_w, n, bc_before)
  call rom_strong_divergence(fluid, u, v, w, divergence, n, strong_before)
  mass_before = rom_mass_norm(fluid, u, v, w, n)

  call assemble_weak_continuity(fluid, p_work, p_rhs, u, v, w, u_zero, &
       v_zero, w_zero, f_x, f_y, f_z, c_before, n, weak_before, event)

  ! Exercise the ROM candidate-column wrapper itself.  Its pressure solve uses
  ! Neko's native Ax_prs; the following legacy algebra block remains below as
  ! an independent comparison of the older constrained Schur construction.
  saved_abs_tol = fluid%ksp_prs%abs_tol
  saved_rel_tol = fluid%ksp_prs%rel_tol
  saved_max_iter = fluid%ksp_prs%max_iter
  ! Use a diagnostic linear tolerance below the production value so that the
  ! 50/100/200 caps distinguish Krylov truncation from the final continuity
  ! residual.  The wrapper and this test restore the production settings.
  fluid%ksp_prs%abs_tol = 1.0e-12_rp
  fluid%ksp_prs%rel_tol = 0.0_rp
  do i = 1, n_pressure_cases
     ! Hold the incoming manufactured field fixed while varying only the
     ! candidate projection Krylov cap.  The wrapper restores each setting;
  ! assign the next cap explicitly before the next independent trial.
     fluid%ksp_prs%max_iter = pressure_case_max(i)
     call field_copy(sol_u, u); call field_copy(sol_v, v); call field_copy(sol_w, w)
     call rom_project_velocity_column(fluid, sol_u, sol_v, sol_w, &
          native_column_weak, linear_residual = native_column_linear_residual, &
          iterations = native_column_iterations, converged = native_column_converged)
     pressure_case_weak(i) = native_column_weak
     pressure_case_linear(i) = native_column_linear_residual
     pressure_case_iterations(i) = native_column_iterations
     pressure_case_converged(i) = native_column_converged
  end do
  fluid%ksp_prs%abs_tol = saved_abs_tol
  fluid%ksp_prs%rel_tol = saved_rel_tol
  fluid%ksp_prs%max_iter = saved_max_iter
  ! Keep the 200-iteration trial as the representative final-column result.
  native_column_weak = pressure_case_weak(n_pressure_cases)
  native_column_linear_residual = pressure_case_linear(n_pressure_cases)
  native_column_iterations = pressure_case_iterations(n_pressure_cases)
  native_column_converged = pressure_case_converged(n_pressure_cases)
  call assemble_weak_continuity(fluid, p_work, p_rhs, sol_u, sol_v, sol_w, &
       u_zero, v_zero, w_zero, f_x, f_y, f_z, c_after, n, &
       native_column_weak, event)
  call rom_bc_residual(fluid, sol_u, sol_v, sol_w, bc_u, bc_v, bc_w, n, &
       native_column_bc)
  call rom_strong_divergence(fluid, sol_u, sol_v, sol_w, divergence, n, &
       native_column_strong)
  native_column_mass = rom_mass_norm(fluid, sol_u, sol_v, sol_w, n)

  call projection_ax%init(fluid, q_op, p_work, p_op, grad_x, grad_y, grad_z, &
       f_x, f_y, f_z, u_zero, v_zero, w_zero)

  call fluid%Ax_prs%compute(ax_phi%x, phi%x, fluid%c_Xh, fluid%msh, fluid%Xh)
  call fluid%gs_Xh%op(ax_phi%x, n, GS_OP_ADD)
  call fluid%bcs_prs_projector%apply(ax_phi%x, n)
  call projection_ax%compute(p_op%x, phi%x, fluid%c_Xh, fluid%msh, fluid%Xh)
  call fluid%gs_Xh%op(p_op%x, n, GS_OP_ADD)
  call fluid%bcs_prs_projector%apply(p_op%x, n)
  p_op%x = p_op%x - ax_phi%x
  schur_error = rom_pressure_norm(p_op, fluid, n)
  schur_relative = schur_error / max(rom_pressure_norm(ax_phi, fluid, n), &
       tiny(1.0_rp))
  c_before%x = c_before%x - ax_phi%x
  identity_plus = rom_pressure_norm(c_before, fluid, n)
  c_before%x = c_before%x + 2.0_rp * ax_phi%x
  identity_minus = rom_pressure_norm(c_before, fluid, n)
  identity_error = min(identity_plus, identity_minus)
  if (identity_plus <= identity_minus) then
     correction_sign = 1.0_rp
  else
     correction_sign = -1.0_rp
  end if
  call assemble_weak_continuity(fluid, p_work, p_rhs, u, v, w, u_zero, &
       v_zero, w_zero, f_x, f_y, f_z, c_before, n, weak_before, event)
  call fluid%Ax_prs%compute(ax_phi%x, phi%x, fluid%c_Xh, fluid%msh, fluid%Xh)
  call fluid%gs_Xh%op(ax_phi%x, n, GS_OP_ADD)
  call fluid%bcs_prs_projector%apply(ax_phi%x, n)
  identity_relative = identity_error / max(rom_pressure_norm(ax_phi, fluid, n), &
       tiny(1.0_rp))

  ! Solve A_p q = sign(C(u)).  The sign is measured above; for a consistent
  ! Neko Pn/Pn pair it must be +1.  The temporary tolerances are diagnostic
  ! settings and are restored before this executable exits.
  p_rhs%x = c_before%x
  call scale_pressure(fluid, p_rhs, correction_sign, n)
  pressure_rhs_norm = rom_pressure_norm(p_rhs, fluid, n)
  call field_rzero(phi_solved)
  saved_abs_tol = fluid%ksp_prs%abs_tol
  saved_rel_tol = fluid%ksp_prs%rel_tol
  saved_max_iter = fluid%ksp_prs%max_iter
  fluid%ksp_prs%abs_tol = 1.0e-12_rp
  fluid%ksp_prs%rel_tol = 0.0_rp
  fluid%ksp_prs%max_iter = max(saved_max_iter, 100)
  call fluid%pc_prs%update()
  ksp_result = fluid%ksp_prs%solve(fluid%Ax_prs, phi_solved, p_rhs%x, n, &
       fluid%c_Xh, fluid%bcs_prs_projector, fluid%gs_Xh, &
       niter = fluid%ksp_prs%max_iter)
  fluid%ksp_prs%abs_tol = saved_abs_tol
  fluid%ksp_prs%rel_tol = saved_rel_tol
  fluid%ksp_prs%max_iter = saved_max_iter

  call projection_ax%compute(p_op%x, phi_solved%x, fluid%c_Xh, fluid%msh, &
       fluid%Xh)
  call fluid%gs_Xh%op(p_op%x, n, GS_OP_ADD)
  call fluid%bcs_prs_projector%apply(p_op%x, n)
  call field_add2s2(p_op, p_rhs, -1.0_rp)
  call fluid%bcs_prs_projector%apply(p_op%x, n)
  constrained_residual = rom_pressure_norm(p_op, fluid, n)
  constrained_relative = constrained_residual / &
       max(pressure_rhs_norm, tiny(1.0_rp))
  call fluid%Ax_prs%compute(ax_phi%x, phi_solved%x, fluid%c_Xh, fluid%msh, &
       fluid%Xh)
  call fluid%gs_Xh%op(ax_phi%x, n, GS_OP_ADD)
  call fluid%bcs_prs_projector%apply(ax_phi%x, n)
  ax_work%x = ax_phi%x
  call field_add2s2(ax_work, p_rhs, -1.0_rp)
  call fluid%bcs_prs_projector%apply(ax_work%x, n)
  pressure_residual = rom_pressure_norm(ax_work, fluid, n)
  pressure_relative = pressure_residual / max(pressure_rhs_norm, tiny(1.0_rp))

  call grad(grad_x%x, grad_y%x, grad_z%x, phi_solved%x, fluid%c_Xh)
  call fluid%gs_Xh%op(grad_x%x, grad_y%x, grad_z%x, n, GS_OP_ADD)
  call col2(grad_x%x, fluid%c_Xh%mult, n)
  call col2(grad_y%x, fluid%c_Xh%mult, n)
  call col2(grad_z%x, fluid%c_Xh%mult, n)
  call fluid%bcs_vel_projector%apply(grad_x%x, grad_y%x, grad_z%x, n)
  correction_norm = sqrt(max(0.0_rp, &
       glsc3(grad_x%x, fluid%c_Xh%B, grad_x%x, n) + &
       glsc3(grad_y%x, fluid%c_Xh%B, grad_y%x, n) + &
       glsc3(grad_z%x, fluid%c_Xh%B, grad_z%x, n)))
  call field_add2s2(u, grad_x, -1.0_rp)
  call field_add2s2(v, grad_y, -1.0_rp)
  call field_add2s2(w, grad_z, -1.0_rp)
  call continuous_vector(fluid, u, v, w, n)
  call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)

  call rom_bc_residual(fluid, u, v, w, bc_u, bc_v, bc_w, n, bc_after)
  call rom_strong_divergence(fluid, u, v, w, divergence, n, strong_after)
  mass_after = rom_mass_norm(fluid, u, v, w, n)
  call assemble_weak_continuity(fluid, p_work, p_rhs, u, v, w, u_zero, &
       v_zero, w_zero, f_x, f_y, f_z, c_after, n, weak_after, event)
  mass_change = abs(mass_after - mass_before) / max(mass_before, tiny(1.0_rp))
  weak_after_normalized = weak_after / max(weak_before, tiny(1.0_rp))
  pure_weak_after = weak_after
  pure_weak_after_normalized = weak_after_normalized

  ! ------------------------------------------------------------------------- !
  ! Mixed manufactured test.  The stream functions vanish with their first
  ! derivatives on x=0,2 and y=0,1, so curl(psi) is an admissible solenoidal
  ! velocity.  Adding grad_h(phi) gives a non-admissible component that the
  ! pressure projection must remove without destroying the solenoidal part.
  call build_solenoidal_mode(fluid, psi, sol_u, sol_v, sol_w, n, 1)
  call build_solenoidal_mode(fluid, psi_two, sol_two_u, sol_two_v, sol_two_w, n, 2)
  call grad(grad_x%x, grad_y%x, grad_z%x, phi%x, fluid%c_Xh)
  call continuous_vector(fluid, grad_x, grad_y, grad_z, n)
  call fluid%bcs_vel_projector%apply(grad_x%x, grad_y%x, grad_z%x, n)
  u%x = sol_u%x + grad_x%x
  v%x = sol_v%x + grad_y%x
  w%x = sol_w%x + grad_z%x
  call continuous_vector(fluid, u, v, w, n)
  call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
  call rom_bc_residual(fluid, u, v, w, bc_u, bc_v, bc_w, n, mixed_bc_before)
  call rom_strong_divergence(fluid, u, v, w, divergence, n, mixed_strong_before)
  mixed_mass_before = rom_mass_norm(fluid, u, v, w, n)
  call assemble_weak_continuity(fluid, p_work, p_rhs, u, v, w, u_zero, &
       v_zero, w_zero, f_x, f_y, f_z, c_before, n, mixed_weak_before, event)

  call project_mixed_mode(fluid, projection_ax, u, v, w, mixed_out_u, &
       mixed_out_v, mixed_out_w, mixed_rhs_norm, mixed_schur_residual, &
       mixed_schur_relative)
  call rom_bc_residual(fluid, mixed_out_u, mixed_out_v, mixed_out_w, bc_u, &
       bc_v, bc_w, n, mixed_bc_after)
  call rom_strong_divergence(fluid, mixed_out_u, mixed_out_v, mixed_out_w, &
       divergence, n, mixed_strong_after)
  mixed_mass_after = rom_mass_norm(fluid, mixed_out_u, mixed_out_v, mixed_out_w, n)
  call assemble_weak_continuity(fluid, p_work, p_rhs, mixed_out_u, mixed_out_v, &
       mixed_out_w, u_zero, v_zero, w_zero, f_x, f_y, f_z, c_after, n, &
       mixed_weak_after, event)
  mixed_mass_change = abs(mixed_mass_after - mixed_mass_before) / &
       max(mixed_mass_before, tiny(1.0_rp))
  mixed_weak_normalized = mixed_weak_after / max(mixed_weak_before, tiny(1.0_rp))

  ! Compare with the known admissible component.  A residual-only test would
  ! not detect destruction of a divergence-free component.
  bc_u%x = mixed_out_u%x - sol_u%x
  bc_v%x = mixed_out_v%x - sol_v%x
  bc_w%x = mixed_out_w%x - sol_w%x
  difference_mass_norm = sqrt(max(0.0_rp, &
       glsc3(bc_u%x, fluid%c_Xh%B, bc_u%x, n) + &
       glsc3(bc_v%x, fluid%c_Xh%B, bc_v%x, n) + &
       glsc3(bc_w%x, fluid%c_Xh%B, bc_w%x, n)))
  sol_norm = rom_mass_norm(fluid, sol_u, sol_v, sol_w, n)
  mixed_preservation_error = difference_mass_norm / max(sol_norm, tiny(1.0_rp))

  ! For an M-orthogonal Helmholtz projection, the correction is orthogonal to
  ! every admissible velocity.  Test it against the manufactured solenoidal
  ! component explicitly.
  bc_u%x = u%x - mixed_out_u%x
  bc_v%x = v%x - mixed_out_v%x
  bc_w%x = w%x - mixed_out_w%x
  correction_mass_norm = sqrt(max(0.0_rp, &
       glsc3(bc_u%x, fluid%c_Xh%B, bc_u%x, n) + &
       glsc3(bc_v%x, fluid%c_Xh%B, bc_v%x, n) + &
       glsc3(bc_w%x, fluid%c_Xh%B, bc_w%x, n)))
  correction_orthogonality = abs(rom_field_inner(fluid, bc_u, bc_v, bc_w, &
       sol_u, sol_v, sol_w, n)) / max(correction_mass_norm * sol_norm, &
       tiny(1.0_rp))

  ! Idempotence: a second projection of an already projected column must not
  ! change it.  The second solve is still isolated and never runs DNS/adjoint.
  call project_mixed_mode(fluid, projection_ax, mixed_out_u, mixed_out_v, &
       mixed_out_w, u, v, w, mixed_rhs_norm, mixed_schur_residual, &
       mixed_schur_relative)
  bc_u%x = u%x - mixed_out_u%x
  bc_v%x = v%x - mixed_out_v%x
  bc_w%x = w%x - mixed_out_w%x
  difference_mass_norm = sqrt(max(0.0_rp, &
       glsc3(bc_u%x, fluid%c_Xh%B, bc_u%x, n) + &
       glsc3(bc_v%x, fluid%c_Xh%B, bc_v%x, n) + &
       glsc3(bc_w%x, fluid%c_Xh%B, bc_w%x, n)))
  mixed_idempotence = difference_mass_norm / max(mixed_mass_after, tiny(1.0_rp))

  ! A second independent stream function checks the final common-basis gate,
  ! not just the single-column projection.
  call grad(grad_x%x, grad_y%x, grad_z%x, phi%x, fluid%c_Xh)
  call continuous_vector(fluid, grad_x, grad_y, grad_z, n)
  call fluid%bcs_vel_projector%apply(grad_x%x, grad_y%x, grad_z%x, n)
  u%x = sol_two_u%x + grad_x%x
  v%x = sol_two_v%x + grad_y%x
  w%x = sol_two_w%x + grad_z%x
  call continuous_vector(fluid, u, v, w, n)
  call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
  call project_mixed_mode(fluid, projection_ax, u, v, w, mixed_out_u, &
       mixed_out_v, mixed_out_w, mixed_rhs_norm, mixed_schur_residual, &
       mixed_schur_relative)
  allocate(mass_weights(n), mixed_basis(3 * n, 2), basis_norm(2), &
       basis_cross(2), mixed_final_index(2))
  mass_weights = reshape(fluid%c_Xh%B, [n])
  mixed_basis(1:n, 1) = reshape(sol_u%x, [n])
  mixed_basis(n + 1:2 * n, 1) = reshape(sol_v%x, [n])
  mixed_basis(2 * n + 1:3 * n, 1) = reshape(sol_w%x, [n])
  mixed_basis(1:n, 2) = reshape(mixed_out_u%x, [n])
  mixed_basis(n + 1:2 * n, 2) = reshape(mixed_out_v%x, [n])
  mixed_basis(2 * n + 1:3 * n, 2) = reshape(mixed_out_w%x, [n])
  call rom_mass_orthonormalize(mixed_basis, n, mass_weights, 1.0e-8_rp, &
       mixed_accepted, mixed_final_index, basis_norm, basis_cross, mixed_basis_ok)
  if (mixed_basis_ok .and. mixed_accepted == 2) then
     basis_gram_error = max(abs(basis_norm(1) - 1.0_rp), &
          abs(basis_norm(2) - 1.0_rp), basis_cross(2))
  else
     basis_gram_error = huge(1.0_rp)
  end if

  ! ------------------------------------------------------------------------- !
  ! Independent constrained-Schur algebra probes.  These use several
  ! pressure directions, including a constant/nullspace candidate, rather
  ! than treating one Krylov solve as evidence that S is the right operator.
  allocate(pressure_probes(n, n_probe), schur_probes(n, n_probe), &
       native_probes(n, n_probe), probe_s_norm(n_probe), probe_a_norm(n_probe), &
       probe_s_a_relative(n_probe), probe_s_quadratic(n_probe), &
       probe_a_quadratic(n_probe))
  do mixed_accepted = 1, n_probe
     do i = 1, n
        x = fluid%u%dof%x%x(i, 1, 1, 1)
        if (mixed_accepted == 1) then
           q_op%x(i, 1, 1, 1) = 1.0_rp
        else if (mixed_accepted == 2) then
           q_op%x(i, 1, 1, 1) = x
        else if (mixed_accepted == 3) then
           q_op%x(i, 1, 1, 1) = fluid%u%dof%y%x(i, 1, 1, 1)
        else if (mixed_accepted == 4) then
           q_op%x(i, 1, 1, 1) = x * fluid%u%dof%y%x(i, 1, 1, 1)
        else
           q_op%x(i, 1, 1, 1) = x * x + fluid%u%dof%y%x(i, 1, 1, 1) ** 2
        end if
     end do
     call continuous_scalar(fluid, q_op, n)
     call fluid%bcs_prs_projector%apply(q_op%x, n)
     pressure_probes(:, mixed_accepted) = reshape(q_op%x, [n])

     call projection_ax%compute(p_op%x, q_op%x, fluid%c_Xh, fluid%msh, fluid%Xh)
     call fluid%gs_Xh%op(p_op%x, n, GS_OP_ADD)
     call fluid%bcs_prs_projector%apply(p_op%x, n)
     schur_probes(:, mixed_accepted) = reshape(p_op%x, [n])
     call fluid%Ax_prs%compute(ax_phi%x, q_op%x, fluid%c_Xh, fluid%msh, fluid%Xh)
     call fluid%gs_Xh%op(ax_phi%x, n, GS_OP_ADD)
     call fluid%bcs_prs_projector%apply(ax_phi%x, n)
     native_probes(:, mixed_accepted) = reshape(ax_phi%x, [n])

     probe_s_norm(mixed_accepted) = rom_pressure_norm(p_op, fluid, n)
     probe_a_norm(mixed_accepted) = rom_pressure_norm(ax_phi, fluid, n)
     probe_s_a_relative(mixed_accepted) = sqrt(max(0.0_rp, &
          glsc3(p_op%x - ax_phi%x, fluid%c_Xh%mult, p_op%x - ax_phi%x, n))) / &
          max(probe_a_norm(mixed_accepted), tiny(1.0_rp))
     probe_s_quadratic(mixed_accepted) = glsc3(q_op%x, fluid%c_Xh%mult, &
          p_op%x, n) / max(glsc3(q_op%x, fluid%c_Xh%mult, q_op%x, n), &
          tiny(1.0_rp))
     probe_a_quadratic(mixed_accepted) = glsc3(q_op%x, fluid%c_Xh%mult, &
          ax_phi%x, n) / max(glsc3(q_op%x, fluid%c_Xh%mult, q_op%x, n), &
          tiny(1.0_rp))
  end do
  schur_symmetry_max = 0.0_rp
  native_symmetry_max = 0.0_rp
  do i = 1, n_probe
     do mixed_accepted = 1, n_probe
        probe_dot_ij = glsc3(pressure_probes(:, i), fluid%c_Xh%mult, &
             schur_probes(:, mixed_accepted), n)
        probe_dot_ji = glsc3(pressure_probes(:, mixed_accepted), fluid%c_Xh%mult, &
             schur_probes(:, i), n)
        schur_symmetry_max = max(schur_symmetry_max, abs(probe_dot_ij - &
             probe_dot_ji) / max(sqrt(abs(glsc3(pressure_probes(:, i), &
             fluid%c_Xh%mult, pressure_probes(:, i), n) * &
             glsc3(pressure_probes(:, mixed_accepted), fluid%c_Xh%mult, &
             pressure_probes(:, mixed_accepted), n))), tiny(1.0_rp)))
        probe_dot_ij = glsc3(pressure_probes(:, i), fluid%c_Xh%mult, &
             native_probes(:, mixed_accepted), n)
        probe_dot_ji = glsc3(pressure_probes(:, mixed_accepted), fluid%c_Xh%mult, &
             native_probes(:, i), n)
        native_symmetry_max = max(native_symmetry_max, abs(probe_dot_ij - &
             probe_dot_ji) / max(sqrt(abs(glsc3(pressure_probes(:, i), &
             fluid%c_Xh%mult, pressure_probes(:, i), n) * &
             glsc3(pressure_probes(:, mixed_accepted), fluid%c_Xh%mult, &
             pressure_probes(:, mixed_accepted), n))), tiny(1.0_rp)))
     end do
  end do
  schur_quadratic_min = minval(probe_s_quadratic)
  native_quadratic_min = minval(probe_a_quadratic)
  schur_probe_rank_indicator = minval(abs(probe_s_quadratic)) / &
       maxval(abs(probe_s_quadratic))
  native_probe_rank_indicator = minval(abs(probe_a_quadratic)) / &
       maxval(abs(probe_a_quadratic))

  ! Distinguish an actual nullspace from a pressure boundary restriction.
  q_op%x = 1.0_rp
  call continuous_scalar(fluid, q_op, n)
  call fluid%bcs_prs_projector%apply(q_op%x, n)
  call projection_ax%compute(p_op%x, q_op%x, fluid%c_Xh, fluid%msh, fluid%Xh)
  call fluid%gs_Xh%op(p_op%x, n, GS_OP_ADD)
  call fluid%bcs_prs_projector%apply(p_op%x, n)
  schur_constant_norm = rom_pressure_norm(p_op, fluid, n)
  call fluid%Ax_prs%compute(ax_phi%x, q_op%x, fluid%c_Xh, fluid%msh, fluid%Xh)
  call fluid%gs_Xh%op(ax_phi%x, n, GS_OP_ADD)
  call fluid%bcs_prs_projector%apply(ax_phi%x, n)
  native_constant_norm = rom_pressure_norm(ax_phi, fluid, n)
  if (fluid%prs_dirichlet) then
     rhs_compatibility = 0.0_rp
  else
     rhs_compatibility = abs(glsc3(c_before%x, fluid%c_Xh%mult, &
          q_op%x, n)) / max(sqrt(glsc3(q_op%x, fluid%c_Xh%mult, q_op%x, n)), &
          tiny(1.0_rp))
  end if

  ! Verify the mass-adjoint relationship directly.  q is a pressure probe,
  ! while the test velocity is a homogeneous gradient field, so this is not a
  ! zero-on-both-sides solenoidal test.
  call grad(u%x, v%x, w%x, phi%x, fluid%c_Xh)
  call continuous_vector(fluid, u, v, w, n)
  call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
  q_op%x = reshape(pressure_probes(:, 4), shape(q_op%x))
  call projection_ax%compute(p_op%x, q_op%x, fluid%c_Xh, fluid%msh, fluid%Xh)
  call assemble_weak_continuity(fluid, p_work, p_rhs, u, v, w, u_zero, &
       v_zero, w_zero, f_x, f_y, f_z, c_after, n, weak_after, event)
  adjoint_lhs = glsc3(q_op%x, fluid%c_Xh%mult, c_after%x, n)
  adjoint_rhs = glsc3(grad_x%x, fluid%c_Xh%B, u%x, n) + &
       glsc3(grad_y%x, fluid%c_Xh%B, v%x, n) + &
       glsc3(grad_z%x, fluid%c_Xh%B, w%x, n)
  adjoint_relative_error = abs(adjoint_lhs - adjoint_rhs) / &
       max(abs(adjoint_lhs), abs(adjoint_rhs), tiny(1.0_rp))

  if (pe_rank == 0) then
     open(newunit = unit, file = trim(diagnostics_file), status = 'replace', &
          action = 'write')
     write(unit, '(A)') 'quantity,value'
     write(unit, '(A,",",ES24.16)') 'weak_continuity_before', weak_before
     write(unit, '(A,",",ES24.16)') 'weak_continuity_after', pure_weak_after
     write(unit, '(A,",",ES24.16)') 'weak_continuity_after_normalized', &
          pure_weak_after_normalized
     write(unit, '(A,",",ES24.16)') 'strong_divergence_before', strong_before
     write(unit, '(A,",",ES24.16)') 'strong_divergence_after', strong_after
     write(unit, '(A,",",ES24.16)') 'bc_residual_before', bc_before
     write(unit, '(A,",",ES24.16)') 'bc_residual_after', bc_after
     write(unit, '(A,",",ES24.16)') 'mass_norm_before', mass_before
     write(unit, '(A,",",ES24.16)') 'mass_norm_after', mass_after
     write(unit, '(A,",",ES24.16)') 'mass_norm_relative_change', mass_change
     write(unit, '(A,",",ES24.16)') 'identity_C_minus_A', identity_plus
     write(unit, '(A,",",ES24.16)') 'identity_C_plus_A', identity_minus
     write(unit, '(A,",",ES24.16)') 'identity_relative_error', identity_relative
     write(unit, '(A,",",ES24.16)') 'constrained_schur_minus_native_A', schur_error
     write(unit, '(A,",",ES24.16)') 'constrained_schur_relative_error', schur_relative
     write(unit, '(A,",",ES24.16)') 'identity_sign', correction_sign
     write(unit, '(A,",",ES24.16)') 'pressure_rhs_norm', pressure_rhs_norm
     write(unit, '(A,",",ES24.16)') 'constrained_schur_residual', constrained_residual
     write(unit, '(A,",",ES24.16)') 'constrained_schur_relative_residual', constrained_relative
     write(unit, '(A,",",ES24.16)') 'pressure_poisson_residual', pressure_residual
     write(unit, '(A,",",ES24.16)') 'pressure_relative_residual', pressure_relative
     write(unit, '(A,",",I0)') 'pressure_iterations', ksp_result%iter
     write(unit, '(A,",",L1)') 'pressure_converged', ksp_result%converged
     write(unit, '(A,",",ES24.16)') 'pressure_ksp_start', ksp_result%res_start
     write(unit, '(A,",",ES24.16)') 'pressure_ksp_final', ksp_result%res_final
     write(unit, '(A,",",ES24.16)') 'correction_mass_norm', correction_norm
     write(unit, '(A,",",ES24.16)') 'native_rom_column_weak_continuity', &
          native_column_weak
     write(unit, '(A,",",ES24.16)') 'native_rom_column_strong_divergence', &
          native_column_strong
     write(unit, '(A,",",ES24.16)') 'native_rom_column_bc_residual', &
          native_column_bc
     write(unit, '(A,",",ES24.16)') 'native_rom_column_mass_norm', &
          native_column_mass
     write(unit, '(A,",",ES24.16)') 'native_rom_column_linear_residual', &
          native_column_linear_residual
     write(unit, '(A,",",I0)') 'native_rom_column_iterations', &
          native_column_iterations
     write(unit, '(A,",",L1)') 'native_rom_column_converged', &
          native_column_converged
     do i = 1, n_pressure_cases
        write(unit, '(A,I0,A,",",ES24.16)') &
             'projection_cap_', pressure_case_max(i), '_linear_residual', &
             pressure_case_linear(i)
        write(unit, '(A,I0,A,",",ES24.16)') &
             'projection_cap_', pressure_case_max(i), '_final_weak_residual', &
             pressure_case_weak(i)
        write(unit, '(A,I0,A,",",I0)') &
             'projection_cap_', pressure_case_max(i), '_iterations', &
             pressure_case_iterations(i)
        write(unit, '(A,I0,A,",",L1)') &
             'projection_cap_', pressure_case_max(i), '_converged', &
             pressure_case_converged(i)
     end do
     write(unit, '(A,",",ES24.16)') 'native_rom_column_linear_residual', &
          native_column_linear_residual
     write(unit, '(A,",",I0)') 'native_rom_column_iterations', &
          native_column_iterations
     write(unit, '(A,",",L1)') 'native_rom_column_converged', &
          native_column_converged
     write(unit, '(A,",",ES24.16)') 'mixed_weak_continuity_before', mixed_weak_before
     write(unit, '(A,",",ES24.16)') 'mixed_weak_continuity_after', mixed_weak_after
     write(unit, '(A,",",ES24.16)') 'mixed_weak_continuity_after_normalized', &
          mixed_weak_normalized
     write(unit, '(A,",",ES24.16)') 'mixed_strong_divergence_before', mixed_strong_before
     write(unit, '(A,",",ES24.16)') 'mixed_strong_divergence_after', mixed_strong_after
     write(unit, '(A,",",ES24.16)') 'mixed_bc_residual_before', mixed_bc_before
     write(unit, '(A,",",ES24.16)') 'mixed_bc_residual_after', mixed_bc_after
     write(unit, '(A,",",ES24.16)') 'mixed_mass_norm_before', mixed_mass_before
     write(unit, '(A,",",ES24.16)') 'mixed_mass_norm_after', mixed_mass_after
     write(unit, '(A,",",ES24.16)') 'mixed_mass_relative_change', mixed_mass_change
     write(unit, '(A,",",ES24.16)') 'mixed_preservation_error', mixed_preservation_error
     write(unit, '(A,",",ES24.16)') 'mixed_idempotence_error', mixed_idempotence
     write(unit, '(A,",",ES24.16)') 'mixed_correction_M_orthogonality', correction_orthogonality
     write(unit, '(A,",",ES24.16)') 'mixed_constrained_schur_rhs', mixed_rhs_norm
     write(unit, '(A,",",ES24.16)') 'mixed_constrained_schur_residual', mixed_schur_residual
     write(unit, '(A,",",ES24.16)') 'mixed_constrained_schur_relative', mixed_schur_relative
     write(unit, '(A,",",ES24.16)') 'mixed_final_basis_gram_error', basis_gram_error
     write(unit, '(A,",",ES24.16)') 'schur_probe_symmetry_max', schur_symmetry_max
     write(unit, '(A,",",ES24.16)') 'native_probe_symmetry_max', native_symmetry_max
     write(unit, '(A,",",ES24.16)') 'schur_probe_quadratic_min', schur_quadratic_min
     write(unit, '(A,",",ES24.16)') 'native_probe_quadratic_min', native_quadratic_min
     write(unit, '(A,",",ES24.16)') 'schur_probe_rank_indicator', schur_probe_rank_indicator
     write(unit, '(A,",",ES24.16)') 'native_probe_rank_indicator', native_probe_rank_indicator
     write(unit, '(A,",",ES24.16)') 'schur_constant_probe_norm', schur_constant_norm
     write(unit, '(A,",",ES24.16)') 'native_constant_probe_norm', native_constant_norm
     write(unit, '(A,",",ES24.16)') 'rhs_range_compatibility_residual', rhs_compatibility
     write(unit, '(A,",",ES24.16)') 'divergence_gradient_adjoint_lhs', adjoint_lhs
     write(unit, '(A,",",ES24.16)') 'divergence_gradient_adjoint_rhs', adjoint_rhs
     write(unit, '(A,",",ES24.16)') 'divergence_gradient_adjoint_relative_error', &
          adjoint_relative_error
     do i = 1, n_probe
        write(unit, '(A,I0,A,ES24.16)') 'schur_probe_', i, '_norm', probe_s_norm(i)
        write(unit, '(A,I0,A,ES24.16)') 'native_probe_', i, '_norm', probe_a_norm(i)
        write(unit, '(A,I0,A,ES24.16)') 'schur_native_probe_', i, '_relative', &
             probe_s_a_relative(i)
        write(unit, '(A,I0,A,ES24.16)') 'schur_probe_', i, '_quadratic', &
             probe_s_quadratic(i)
        write(unit, '(A,I0,A,ES24.16)') 'native_probe_', i, '_quadratic', &
             probe_a_quadratic(i)
     end do
     close(unit)

     write(*, '(A)') 'Projection algebra diagnostic (manufactured u = grad_h(phi))'
     write(*, '(A,ES12.4,A,ES12.4,A,ES12.4)') &
          'weak continuity before/after/normalised = ', weak_before, ' / ', &
          pure_weak_after, ' / ', pure_weak_after_normalized
     write(*, '(A,ES12.4,A,ES12.4)') &
          'strong divergence before/after = ', strong_before, ' / ', strong_after
     write(*, '(A,ES12.4,A,ES12.4)') &
          'homogeneous BC residual before/after = ', bc_before, ' / ', bc_after
     write(*, '(A,ES12.4,A,ES12.4,A,ES12.4)') &
          'mass norm before/after/relative change = ', mass_before, ' / ', &
          mass_after, ' / ', mass_change
     write(*, '(A,ES12.4,A,ES12.4,A,ES12.4)') &
          'C-A / C+A / relative identity error = ', identity_plus, ' / ', &
          identity_minus, ' / ', identity_relative
     write(*, '(A,ES12.4,A,ES12.4)') &
          'constrained Schur minus native A / relative = ', schur_error, ' / ', &
          schur_relative
     write(*, '(A,ES12.4,A,ES12.4,A,ES12.4)') &
          'constrained Schur rhs/residual/relative = ', pressure_rhs_norm, ' / ', &
          constrained_residual, ' / ', constrained_relative
     write(*, '(A,ES12.4,A,ES12.4,A,ES12.4)') &
          'pressure Poisson rhs/residual/relative = ', pressure_rhs_norm, ' / ', &
          pressure_residual, ' / ', pressure_relative
     write(*, '(A,I0,A,L1)') 'pressure iterations/converged = ', &
          ksp_result%iter, ' / ', ksp_result%converged
     do i = 1, n_pressure_cases
        write(*, '(A,I0,A,ES12.4,A,ES12.4,A,I0,A,L1)') &
             'projection cap ', pressure_case_max(i), &
             ' linear/final weak/iterations/converged = ', pressure_case_linear(i), &
             ' / ', pressure_case_weak(i), ' / ', pressure_case_iterations(i), &
             ' / ', pressure_case_converged(i)
     end do
     write(*, '(A,A)') 'CSV evidence: ', trim(diagnostics_file)
     write(*, '(A,ES12.4,A,ES12.4,A,ES12.4)') &
          'mixed weak continuity before/after/normalised = ', mixed_weak_before, ' / ', &
          mixed_weak_after, ' / ', mixed_weak_normalized
     write(*, '(A,ES12.4,A,ES12.4)') &
          'mixed strong divergence before/after = ', mixed_strong_before, ' / ', &
          mixed_strong_after
     write(*, '(A,ES12.4,A,ES12.4)') &
          'mixed homogeneous BC before/after = ', mixed_bc_before, ' / ', mixed_bc_after
     write(*, '(A,ES12.4,A,ES12.4,A,ES12.4)') &
          'mixed mass norm before/after/relative change = ', mixed_mass_before, ' / ', &
          mixed_mass_after, ' / ', mixed_mass_change
     write(*, '(A,ES12.4)') 'mixed admissible-component preservation error = ', &
          mixed_preservation_error
     write(*, '(A,ES12.4)') 'mixed idempotence error = ', mixed_idempotence
     write(*, '(A,ES12.4)') 'mixed correction M-orthogonality = ', &
          correction_orthogonality
     write(*, '(A,ES12.4)') 'mixed final-basis Gram error = ', basis_gram_error
     write(*, '(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
          'native ROM column weak/strong/BC/mass = ', native_column_weak, ' / ', &
          native_column_strong, ' / ', native_column_bc, ' / ', native_column_mass
     write(*, '(A,ES12.4,A,I0,A,L1)') &
          'native ROM column pressure linear residual/iterations/converged = ', &
          native_column_linear_residual, ' / ', native_column_iterations, ' / ', &
          native_column_converged
     write(*, '(A,ES12.4,A,ES12.4)') 'Schur/native probe symmetry maxima = ', &
          schur_symmetry_max, ' / ', native_symmetry_max
     write(*, '(A,ES12.4,A,ES12.4)') 'Schur/native minimum quadratic forms = ', &
          schur_quadratic_min, ' / ', native_quadratic_min
     write(*, '(A,ES12.4,A,ES12.4)') 'Schur/native probe rank indicators = ', &
          schur_probe_rank_indicator, ' / ', native_probe_rank_indicator
     write(*, '(A,ES12.4,A,ES12.4)') 'Schur/native constant probe norms = ', &
          schur_constant_norm, ' / ', native_constant_norm
     write(*, '(A,ES12.4)') 'pressure RHS range compatibility residual = ', &
          rhs_compatibility
     write(*, '(A,ES12.4,A,ES12.4,A,ES12.4)') &
          'weak-divergence/gradient adjoint lhs/rhs/relative = ', adjoint_lhs, ' / ', &
          adjoint_rhs, ' / ', adjoint_relative_error
     do i = 1, n_probe
        write(*, '(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') 'probe ', i, &
             ' S/A norm/relative = ', probe_s_norm(i), ' / ', probe_a_norm(i), &
             ' / ', probe_s_a_relative(i)
     end do
     if (fluid%prs_dirichlet) then
        write(*, '(A)') 'pressure nullspace treatment = Dirichlet projector (no mean nullspace)'
     else
        write(*, '(A)') 'pressure nullspace treatment = global mean removal'
     end if
  end if

  call phi%free(); call phi_solved%free(); call ax_phi%free(); call ax_work%free()
  call p_work%free(); call p_rhs%free(); call p_op%free(); call q_op%free()
  call u%free(); call v%free(); call w%free()
  call u_zero%free(); call v_zero%free(); call w_zero%free()
  call f_x%free(); call f_y%free(); call f_z%free()
  call grad_x%free(); call grad_y%free(); call grad_z%free()
  call divergence%free(); call bc_u%free(); call bc_v%free(); call bc_w%free()
  call c_before%free(); call c_after%free()
  call psi%free(); call psi_two%free()
  call sol_u%free(); call sol_v%free(); call sol_w%free()
  call sol_two_u%free(); call sol_two_v%free(); call sol_two_w%free()
  call mixed_out_u%free(); call mixed_out_v%free(); call mixed_out_w%free()
  if (allocated(mass_weights)) deallocate(mass_weights)
  if (allocated(mixed_basis)) deallocate(mixed_basis)
  if (allocated(basis_norm)) deallocate(basis_norm)
  if (allocated(basis_cross)) deallocate(basis_cross)
  if (allocated(mixed_final_index)) deallocate(mixed_final_index)
  if (allocated(pressure_probes)) deallocate(pressure_probes)
  if (allocated(schur_probes)) deallocate(schur_probes)
  if (allocated(native_probes)) deallocate(native_probes)
  if (allocated(probe_s_norm)) deallocate(probe_s_norm)
  if (allocated(probe_a_norm)) deallocate(probe_a_norm)
  if (allocated(probe_s_a_relative)) deallocate(probe_s_a_relative)
  if (allocated(probe_s_quadratic)) deallocate(probe_s_quadratic)
  if (allocated(probe_a_quadratic)) deallocate(probe_a_quadratic)
  call simulation%free()
  call neko_finalize()

contains

  !> Gather and average one scalar field into Neko's continuous space.
  subroutine continuous_scalar(fluid, a, n)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(field_t), intent(inout) :: a
    integer, intent(in) :: n
    call fluid%gs_Xh%op(a%x, n, GS_OP_ADD)
    call col2(a%x, fluid%c_Xh%mult, n)
  end subroutine continuous_scalar

  !> Gather and average a vector field into Neko's continuous space.
  subroutine continuous_vector(fluid, a, b, c, n)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(field_t), intent(inout) :: a, b, c
    integer, intent(in) :: n
    call fluid%gs_Xh%op(a%x, b%x, c%x, n, GS_OP_ADD)
    call col2(a%x, fluid%c_Xh%mult, n)
    call col2(b%x, fluid%c_Xh%mult, n)
    call col2(c%x, fluid%c_Xh%mult, n)
  end subroutine continuous_vector

  !> Assemble the native weak continuity residual for a test velocity.
  subroutine assemble_weak_continuity(fluid, p, p_res, u, v, w, u_e, v_e, &
       w_e, f_x, f_y, f_z, c, n, norm, event)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(field_t), intent(inout) :: p, p_res, u, v, w, c
    type(field_t), intent(in) :: u_e, v_e, w_e
    type(field_t), intent(inout) :: f_x, f_y, f_z
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: norm
    type(c_ptr), intent(inout) :: event

    call field_col3(f_x, u, fluid%rho)
    call field_col3(f_y, v, fluid%rho)
    call field_col3(f_z, w, fluid%rho)
    call col2(f_x%x, fluid%c_Xh%B, n)
    call col2(f_y%x, fluid%c_Xh%B, n)
    call col2(f_z%x, fluid%c_Xh%B, n)
    call field_rzero(p)
    call fluid%prs_res%compute(p, p_res, u, v, w, u_e, v_e, w_e, f_x, f_y, f_z, &
         fluid%c_Xh, fluid%gs_Xh, fluid%bc_prs_surface, fluid%bc_sym_surface, &
         fluid%Ax_prs, 0.0_rp, 1.0_rp, u_e, fluid%rho, event)
    if (.not. fluid%prs_dirichlet) then
       call ortho(p_res%x, fluid%glb_n_points, n)
    end if
    call fluid%gs_Xh%op(p_res%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(p_res%x, n)
    c%x = p_res%x
    norm = rom_pressure_norm(c, fluid, n)
  end subroutine assemble_weak_continuity

  !> Scale a pressure field in the diagnostic coefficient space.
  subroutine scale_pressure(fluid, p, scale, n)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(field_t), intent(inout) :: p
    real(kind=rp), intent(in) :: scale
    integer, intent(in) :: n
    p%x = scale * p%x
    call fluid%bcs_prs_projector%apply(p%x, n)
  end subroutine scale_pressure

  !> Compute a pressure-space SEM mass norm.
  function rom_pressure_norm(a, fluid, n) result(value)
    type(field_t), intent(in) :: a
    type(fluid_pnpn_t), intent(in) :: fluid
    integer, intent(in) :: n
    real(kind=rp) :: value
    value = sqrt(max(0.0_rp, glsc3(a%x, fluid%c_Xh%mult, a%x, n)))
  end function rom_pressure_norm

  !> Compute a three-component velocity mass norm.
  function rom_mass_norm(fluid, u, v, w, n) result(value)
    type(fluid_pnpn_t), intent(in) :: fluid
    type(field_t), intent(in) :: u, v, w
    integer, intent(in) :: n
    real(kind=rp) :: value
    value = sqrt(max(0.0_rp, &
         glsc3(u%x, fluid%c_Xh%B, u%x, n) + &
         glsc3(v%x, fluid%c_Xh%B, v%x, n) + &
         glsc3(w%x, fluid%c_Xh%B, w%x, n)))
  end function rom_mass_norm

  !> Measure strong SEM divergence for a manufactured velocity.
  subroutine rom_strong_divergence(fluid, u, v, w, divergence, n, value)
    type(fluid_pnpn_t), intent(in) :: fluid
    type(field_t), intent(in) :: u, v, w
    type(field_t), intent(inout) :: divergence
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: value
    call div(divergence%x, u%x, v%x, w%x, fluid%c_Xh)
    value = sqrt(max(0.0_rp, glsc3(divergence%x, fluid%c_Xh%B, &
         divergence%x, n)))
  end subroutine rom_strong_divergence

  !> Measure the homogeneous velocity-BC residual of a test field.
  subroutine rom_bc_residual(fluid, u, v, w, bu, bv, bw, n, value)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(field_t), intent(in) :: u, v, w
    type(field_t), intent(inout) :: bu, bv, bw
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: value
    bu%x = u%x; bv%x = v%x; bw%x = w%x
    call fluid%bcs_vel_projector%apply(bu%x, bv%x, bw%x, n)
    bu%x = bu%x - u%x; bv%x = bv%x - v%x; bw%x = bw%x - w%x
    value = sqrt(max(0.0_rp, &
         glsc3(bu%x, fluid%c_Xh%B, bu%x, n) + &
         glsc3(bv%x, fluid%c_Xh%B, bv%x, n) + &
         glsc3(bw%x, fluid%c_Xh%B, bw%x, n)))
  end subroutine rom_bc_residual

  !> Build one stream-function curl used by the manufactured test.
  subroutine build_solenoidal_mode(fluid, psi, u, v, w, n, variant)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(field_t), intent(inout) :: psi, u, v, w
    integer, intent(in) :: n, variant
    integer :: i
    real(kind=rp) :: x, y

    do i = 1, n
       x = fluid%u%dof%x%x(i, 1, 1, 1)
       y = fluid%u%dof%y%x(i, 1, 1, 1)
       if (variant == 1) then
          psi%x(i, 1, 1, 1) = x * x * (2.0_rp - x) * (2.0_rp - x) * &
               y * y * (1.0_rp - y) * (1.0_rp - y)
       else
          psi%x(i, 1, 1, 1) = x * x * x * (2.0_rp - x) * (2.0_rp - x) * &
               y * y * (1.0_rp - y) * (1.0_rp - y) * (1.0_rp + 0.25_rp * y)
       end if
    end do
    call continuous_scalar(fluid, psi, n)
    call grad(grad_x%x, grad_y%x, grad_z%x, psi%x, fluid%c_Xh)
    u%x = grad_y%x
    v%x = -grad_x%x
    w%x = 0.0_rp
    call continuous_vector(fluid, u, v, w, n)
    call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
  end subroutine build_solenoidal_mode

  !> Project and measure a mixed manufactured pressure/solenoidal mode.
  subroutine project_mixed_mode(fluid, projection_ax, u_in, v_in, w_in, &
       u_out, v_out, w_out, rhs_norm, residual, relative)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(rom_projection_ax_t), intent(inout) :: projection_ax
    type(field_t), intent(inout) :: u_in, v_in, w_in
    type(field_t), intent(inout) :: u_out, v_out, w_out
    real(kind=rp), intent(out) :: rhs_norm, residual, relative
    type(ksp_monitor_t) :: local_ksp
    real(kind=rp) :: saved_abs, saved_rel
    integer :: saved_max, n_local
    type(c_ptr) :: local_event

    n_local = fluid%u%dof%size()
    local_event = C_NULL_PTR
    call assemble_weak_continuity(fluid, p_work, p_rhs, u_in, v_in, w_in, &
         u_zero, v_zero, w_zero, f_x, f_y, f_z, c_before, n_local, rhs_norm, &
         local_event)
    call field_rzero(phi_solved)
    saved_abs = fluid%ksp_prs%abs_tol
    saved_rel = fluid%ksp_prs%rel_tol
    saved_max = fluid%ksp_prs%max_iter
    fluid%ksp_prs%abs_tol = 1.0e-12_rp
    fluid%ksp_prs%rel_tol = 0.0_rp
    fluid%ksp_prs%max_iter = max(saved_max, 100)
    call fluid%pc_prs%update()
    local_ksp = fluid%ksp_prs%solve(fluid%Ax_prs, phi_solved, c_before%x, &
         n_local, fluid%c_Xh, fluid%bcs_prs_projector, fluid%gs_Xh, &
         niter = fluid%ksp_prs%max_iter)
    fluid%ksp_prs%abs_tol = saved_abs
    fluid%ksp_prs%rel_tol = saved_rel
    fluid%ksp_prs%max_iter = saved_max

    call projection_ax%compute(p_op%x, phi_solved%x, fluid%c_Xh, fluid%msh, &
         fluid%Xh)
    call fluid%gs_Xh%op(p_op%x, n_local, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(p_op%x, n_local)
    p_op%x = p_op%x - c_before%x
    call fluid%bcs_prs_projector%apply(p_op%x, n_local)
    residual = rom_pressure_norm(p_op, fluid, n_local)
    relative = residual / max(rhs_norm, tiny(1.0_rp))

    call grad(grad_x%x, grad_y%x, grad_z%x, phi_solved%x, fluid%c_Xh)
    call continuous_vector(fluid, grad_x, grad_y, grad_z, n_local)
    call fluid%bcs_vel_projector%apply(grad_x%x, grad_y%x, grad_z%x, n_local)
    u_out%x = u_in%x - grad_x%x
    v_out%x = v_in%x - grad_y%x
    w_out%x = w_in%x - grad_z%x
    call continuous_vector(fluid, u_out, v_out, w_out, n_local)
    call fluid%bcs_vel_projector%apply(u_out%x, v_out%x, w_out%x, n_local)
  end subroutine project_mixed_mode

  !> Compute a three-component SEM mass inner product.
  function rom_field_inner(fluid, au, av, aw, bu, bv, bw, n_local) result(value)
    type(fluid_pnpn_t), intent(in) :: fluid
    type(field_t), intent(in) :: au, av, aw, bu, bv, bw
    integer, intent(in) :: n_local
    real(kind=rp) :: value
    value = glsc3(au%x, fluid%c_Xh%B, bu%x, n_local) + &
         glsc3(av%x, fluid%c_Xh%B, bv%x, n_local) + &
         glsc3(aw%x, fluid%c_Xh%B, bw%x, n_local)
  end function rom_field_inner

end program rom_projection_algebra_diagnostic
