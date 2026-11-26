program rom_projection_transpose_diagnostic
  !> Numerical dot-product, idempotence and boundary check for the ROM-local
  !! candidate projection and its mass-adjoint.
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use fluid_pnpn, only: fluid_pnpn_t
  use field, only: field_t
  use field_math, only: field_rzero
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2, glsc2, glsc3
  use num_types, only: rp
  use utils, only: neko_error
  use rom_basis_handoff, only: rom_project_velocity_column, &
       rom_project_velocity_cotangent, rom_bc_residual_norm
  use, intrinsic :: iso_c_binding, only: c_ptr
  implicit none

  type(simulation_t) :: simulation
  type(json_file) :: parameters
  type(fluid_pnpn_t), pointer :: fluid
  type(field_t) :: x_u, x_v, x_w, px_u, px_v, px_w
  type(field_t) :: y_u, y_v, y_w, pty_u, pty_v, pty_w
  type(field_t) :: bc_u, bc_v, bc_w
  character(len=1024) :: case_file
  real(kind=rp) :: lhs, rhs, dot_error, idempotence
  real(kind=rp) :: weak_residual(3), linear_residual(3), true_relative(3)
  real(kind=rp) :: bc_residual, x_norm
  integer :: iterations(3), n, i, argc, saved_max
  logical :: converged(3)
  logical :: found

  argc = command_argument_count()
  if (argc < 1) call neko_error( &
       'Usage: ROM_rugby_ball_projection_transpose_diagnostic CASE')
  call get_command_argument(1, case_file)

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
     call neko_error('Projection transpose diagnostic requires a Pn/Pn fluid')
  end select

  n = fluid%u%dof%size()
  call x_u%init(fluid%dm_Xh, 'projection_transpose_x_u')
  call x_v%init(fluid%dm_Xh, 'projection_transpose_x_v')
  call x_w%init(fluid%dm_Xh, 'projection_transpose_x_w')
  call px_u%init(fluid%dm_Xh, 'projection_transpose_px_u')
  call px_v%init(fluid%dm_Xh, 'projection_transpose_px_v')
  call px_w%init(fluid%dm_Xh, 'projection_transpose_px_w')
  call y_u%init(fluid%dm_Xh, 'projection_transpose_y_u')
  call y_v%init(fluid%dm_Xh, 'projection_transpose_y_v')
  call y_w%init(fluid%dm_Xh, 'projection_transpose_y_w')
  call pty_u%init(fluid%dm_Xh, 'projection_transpose_pty_u')
  call pty_v%init(fluid%dm_Xh, 'projection_transpose_pty_v')
  call pty_w%init(fluid%dm_Xh, 'projection_transpose_pty_w')
  call bc_u%init(fluid%dm_Xh, 'projection_transpose_bc_u')
  call bc_v%init(fluid%dm_Xh, 'projection_transpose_bc_v')
  call bc_w%init(fluid%dm_Xh, 'projection_transpose_bc_w')
  call field_rzero(x_w); call field_rzero(y_w)
  do i = 1, n
     x_u%x(i,1,1,1) = 0.3_rp + fluid%u%dof%x%x(i,1,1,1) ** 2
     x_v%x(i,1,1,1) = -0.2_rp + fluid%u%dof%y%x(i,1,1,1)
     y_u%x(i,1,1,1) = 0.4_rp + fluid%u%dof%y%x(i,1,1,1) ** 2
     y_v%x(i,1,1,1) = 0.1_rp + fluid%u%dof%x%x(i,1,1,1)
  end do
  call continuous_vector(fluid, x_u, x_v, x_w, n)
  call continuous_vector(fluid, y_u, y_v, y_w, n)
  call fluid%bcs_vel_projector%apply(x_u%x, x_v%x, x_w%x, n)
  call fluid%bcs_vel_projector%apply(y_u%x, y_v%x, y_w%x, n)

  saved_max = fluid%ksp_prs%max_iter
  fluid%ksp_prs%abs_tol = 1.0e-12_rp
  fluid%ksp_prs%rel_tol = 0.0_rp
  fluid%ksp_prs%max_iter = 1600
  call field_copy(px_u, x_u); call field_copy(px_v, x_v); call field_copy(px_w, x_w)
  call rom_project_velocity_column(fluid, px_u, px_v, px_w, weak_residual(1), &
       linear_residual = linear_residual(1), iterations = iterations(1), &
       converged = converged(1), true_relative_residual = true_relative(1))
  ! The cotangent routine returns Euclidean P^T G. Supply M*y and test the
  ! Euclidean transpose identity with consistent coefficient pairings.
  call col2(y_u%x, fluid%c_Xh%B, n)
  call col2(y_v%x, fluid%c_Xh%B, n)
  call col2(y_w%x, fluid%c_Xh%B, n)
  call field_copy(pty_u, y_u); call field_copy(pty_v, y_v); call field_copy(pty_w, y_w)
  call rom_project_velocity_cotangent(fluid, pty_u, pty_v, pty_w, weak_residual(2), &
       linear_residual = linear_residual(2), iterations = iterations(2), &
       converged = converged(2), true_relative_residual = true_relative(2))
  lhs = glsc2(y_u%x, px_u%x, n) + glsc2(y_v%x, px_v%x, n) + &
       glsc2(y_w%x, px_w%x, n)
  rhs = glsc2(pty_u%x, x_u%x, n) + glsc2(pty_v%x, x_v%x, n) + &
       glsc2(pty_w%x, x_w%x, n)
  dot_error = abs(lhs-rhs) / max(abs(lhs), abs(rhs), tiny(1.0_rp))

  call field_copy(x_u, px_u); call field_copy(x_v, px_v); call field_copy(x_w, px_w)
  call rom_project_velocity_column(fluid, x_u, x_v, x_w, weak_residual(3), &
       linear_residual = linear_residual(3), iterations = iterations(3), &
       converged = converged(3), true_relative_residual = true_relative(3))
  idempotence = sqrt(max(0.0_rp, &
       glsc3(x_u%x-px_u%x, fluid%c_Xh%B, x_u%x-px_u%x, n) + &
       glsc3(x_v%x-px_v%x, fluid%c_Xh%B, x_v%x-px_v%x, n) + &
       glsc3(x_w%x-px_w%x, fluid%c_Xh%B, x_w%x-px_w%x, n))) / &
       max(sqrt(max(0.0_rp, glsc3(px_u%x, fluid%c_Xh%B, px_u%x, n) + &
       glsc3(px_v%x, fluid%c_Xh%B, px_v%x, n) + &
       glsc3(px_w%x, fluid%c_Xh%B, px_w%x, n))), tiny(1.0_rp))
  call rom_bc_residual_norm(fluid, px_u, px_v, px_w, bc_u, bc_v, bc_w, n, &
       bc_residual)
  x_norm = sqrt(max(0.0_rp, &
       glsc3(x_u%x, fluid%c_Xh%B, x_u%x, n) + &
       glsc3(x_v%x, fluid%c_Xh%B, x_v%x, n) + &
       glsc3(x_w%x, fluid%c_Xh%B, x_w%x, n)))
  fluid%ksp_prs%max_iter = saved_max

  do i = 1, 3
     write(*,'(A,I0,A,L1,A,I0)') 'projection solve ', i, ' converged=', &
          converged(i), ' iterations=', iterations(i)
     write(*,'(A,I0,A,ES12.4)') 'projection solve ', i, ' KSP residual=', &
          linear_residual(i)
     write(*,'(A,I0,A,ES12.4)') 'projection solve ', i, &
          ' true relative Schur residual=', true_relative(i)
     write(*,'(A,I0,A,ES12.4)') 'projection solve ', i, &
          ' weak continuity=', weak_residual(i)
  end do
  write(*,'(A,ES12.4)') 'Euclidean transpose relative error for M-adjoint pairing = ', dot_error
  write(*,'(A,ES12.4)') 'candidate projection idempotence error = ', idempotence
  write(*,'(A,ES12.4)') 'homogeneous BC residual / input mass norm = ', &
       bc_residual / max(x_norm, tiny(1.0_rp))
  if (.not. all(converged)) call neko_error( &
       'projection transpose diagnostic requires all Schur solves to converge')
  if (maxval(true_relative) > 1.0e-8_rp .or. maxval(weak_residual) > 1.0e-8_rp) &
       call neko_error('projection transpose diagnostic Schur/continuity gate failed')
  if (dot_error > 1.0e-6_rp .or. idempotence > 1.0e-6_rp .or. &
       bc_residual / max(x_norm, tiny(1.0_rp)) > 1.0e-8_rp) call neko_error( &
       'projection transpose diagnostic transpose/idempotence/BC gate failed')

  call x_u%free(); call x_v%free(); call x_w%free()
  call px_u%free(); call px_v%free(); call px_w%free()
  call y_u%free(); call y_v%free(); call y_w%free()
  call pty_u%free(); call pty_v%free(); call pty_w%free()
  call bc_u%free(); call bc_v%free(); call bc_w%free()
  call neko_finalize()

contains

  subroutine continuous_vector(fluid, u, v, w, n)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(field_t), intent(inout) :: u, v, w
    integer, intent(in) :: n
    call col2(u%x, fluid%c_Xh%B, n)
    call col2(v%x, fluid%c_Xh%B, n)
    call col2(w%x, fluid%c_Xh%B, n)
    call fluid%gs_Xh%op(u%x, v%x, w%x, n, GS_OP_ADD)
    call col2(u%x, fluid%c_Xh%Binv, n)
    call col2(v%x, fluid%c_Xh%Binv, n)
    call col2(w%x, fluid%c_Xh%Binv, n)
  end subroutine continuous_vector

  subroutine field_copy(dst, src)
    type(field_t), intent(inout) :: dst
    type(field_t), intent(in) :: src
    dst%x = src%x
  end subroutine field_copy

end program rom_projection_transpose_diagnostic
