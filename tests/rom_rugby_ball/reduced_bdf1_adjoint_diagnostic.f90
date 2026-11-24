!> Finite-difference test of the exact reduced BDF1 discrete adjoint.
!!
!! The parameter perturbs both the constant and implicit reduced operators.
!! The objective is a time integral of a quadratic reduced state functional.
!! The test therefore checks the actual IMEX-BDF1 transition, not a
!! continuous-time surrogate.
program rom_reduced_bdf1_adjoint_diagnostic
  use neko, only: neko_init, neko_finalize
  use num_types, only: rp
  use matrix, only: matrix_t
  use vector, only: vector_t
  use rom_model, only: rom_model_t
  use utils, only: neko_error
  implicit none

  integer, parameter :: n = 2, n_steps = 6
  real(kind=rp), parameter :: dt = 0.05_rp, parameter_value = 0.4_rp
  real(kind=rp), parameter :: epsilon = 1.0e-6_rp
  type(rom_model_t) :: model
  type(matrix_t) :: mass, linear, implicit_linear, state_trajectory, adjoint_trajectory
  type(matrix_t) :: objective_hessian, running_gradient
  type(matrix_t), allocatable :: quadratic(:)
  type(vector_t) :: constant, initial_state, terminal_gradient
  type(vector_t) :: lambda, next_adjoint, state
  real(kind=rp) :: objective, objective_plus, objective_minus, predicted
  real(kind=rp) :: finite_difference, relative_error
  integer :: i, step

  call neko_init()
  allocate(quadratic(n))
  call mass%init(n,n); call linear%init(n,n); call implicit_linear%init(n,n)
  call objective_hessian%init(n,n)
  call constant%init(n); call initial_state%init(n); call terminal_gradient%init(n)
  call running_gradient%init(n,n_steps)
  do i = 1, n
     call quadratic(i)%init(n,n)
     quadratic(i)%x = 0.0_rp
  end do
  mass%x = reshape([1.2_rp, 0.05_rp, 0.1_rp, 0.9_rp], [n,n])
  linear%x = reshape([0.11_rp, -0.08_rp, 0.04_rp, 0.07_rp], [n,n])
  objective_hessian%x = reshape([1.0_rp, 0.2_rp, 0.2_rp, 0.7_rp], [n,n])
  initial_state%x = [0.12_rp, -0.21_rp]

  call evaluate(parameter_value, objective, state_trajectory)
  call build_objective_derivatives(state_trajectory, running_gradient, terminal_gradient)
  call model_at(parameter_value, model)
  call adjoint_trajectory%init(n, n_steps+1)
  call model%discrete_adjoint_bdf1(state_trajectory, dt, running_gradient, &
       terminal_gradient, adjoint_trajectory)

  predicted = 0.0_rp
  call lambda%init(n); call next_adjoint%init(n); call state%init(n)
  do step = 1, n_steps
     next_adjoint%x = adjoint_trajectory%x(:,step+1)
     state%x = state_trajectory%x(:,step+1)
     call model%solve_step_transpose(dt, next_adjoint, lambda)
     predicted = predicted + dt * dot_product(lambda%x, [0.31_rp, -0.17_rp] + &
          matmul(reshape([-0.10_rp, 0.02_rp, 0.03_rp, -0.05_rp], [n,n]), state%x))
  end do
  call lambda%free(); call next_adjoint%free(); call state%free()
  call model%free()
  call adjoint_trajectory%free()
  call state_trajectory%free()

  call evaluate(parameter_value + epsilon, objective_plus, state_trajectory)
  call state_trajectory%free()
  call evaluate(parameter_value - epsilon, objective_minus, state_trajectory)
  call state_trajectory%free()
  finite_difference = (objective_plus - objective_minus) / (2.0_rp * epsilon)
  relative_error = abs(finite_difference - predicted) / max(1.0_rp, abs(finite_difference))
  write(*,'(A,ES16.8,A,ES16.8,A,ES12.4)') &
       'Reduced BDF1 adjoint FD/predicted/relative error = ', finite_difference, &
       ' / ', predicted, ' / ', relative_error
  if (relative_error > 2.0e-9_rp) call neko_error( &
       'reduced BDF1 discrete adjoint finite-difference gate failed')

  call running_gradient%free(); call terminal_gradient%free()
  call objective_hessian%free(); call mass%free(); call linear%free(); call implicit_linear%free()
  call constant%free(); call initial_state%free()
  do i = 1, n; call quadratic(i)%free(); end do
  deallocate(quadratic)
  call neko_finalize()

contains

  !> Build a parameterized reduced model and integrate its BDF1 trajectory.
  subroutine evaluate(parameter, value, states)
    real(kind=rp), intent(in) :: parameter
    real(kind=rp), intent(out) :: value
    type(matrix_t), intent(inout) :: states
    type(rom_model_t) :: local_model
    type(vector_t) :: local_constant
    type(matrix_t) :: local_implicit
    integer :: k
    call local_constant%init(n)
    call local_implicit%init(n,n)
    local_constant%x = [0.20_rp, -0.12_rp] + parameter * [0.31_rp, -0.17_rp]
    local_implicit%x = reshape([-0.30_rp, 0.02_rp, 0.03_rp, -0.20_rp], [n,n]) + &
         parameter * reshape([-0.10_rp, 0.02_rp, 0.03_rp, -0.05_rp], [n,n])
    call model_at_with_implicit(local_model, local_constant, local_implicit)
    call states%init(n,n_steps+1)
    call local_model%forward_bdf1(initial_state, dt, states)
    value = 0.0_rp
    do k = 1, n_steps
       value = value + 0.5_rp * dt * dot_product(states%x(:,k+1), &
            matmul(objective_hessian%x, states%x(:,k+1)))
    end do
    call local_model%free(); call local_constant%free(); call local_implicit%free()
  end subroutine evaluate

  !> Initialize the reference model used by the adjoint check.
  subroutine model_at(parameter, local_model)
    real(kind=rp), intent(in) :: parameter
    type(rom_model_t), intent(inout) :: local_model
    type(vector_t) :: local_constant
    type(matrix_t) :: local_implicit
    call local_constant%init(n); call local_implicit%init(n,n)
    local_constant%x = [0.20_rp, -0.12_rp] + parameter * [0.31_rp, -0.17_rp]
    local_implicit%x = reshape([-0.30_rp, 0.02_rp, 0.03_rp, -0.20_rp], [n,n]) + &
         parameter * reshape([-0.10_rp, 0.02_rp, 0.03_rp, -0.05_rp], [n,n])
    call model_at_with_implicit(local_model, local_constant, local_implicit)
    call local_constant%free(); call local_implicit%free()
  end subroutine model_at

  !> Initialize a model from the fixed and parameterized reduced operators.
  subroutine model_at_with_implicit(local_model, local_constant, local_implicit)
    type(rom_model_t), intent(inout) :: local_model
    type(vector_t), intent(in) :: local_constant
    type(matrix_t), intent(in) :: local_implicit
    call local_model%init(mass, local_constant, linear, quadratic, local_implicit)
  end subroutine model_at_with_implicit

  !> Assemble the running and terminal state derivatives for the discrete sum.
  subroutine build_objective_derivatives(states, running, terminal)
    type(matrix_t), intent(in) :: states
    type(matrix_t), intent(inout) :: running
    type(vector_t), intent(inout) :: terminal
    integer :: k
    running%x = 0.0_rp
    terminal%x = dt * matmul(objective_hessian%x, states%x(:,n_steps+1))
    do k = 1, n_steps-1
       running%x(:,k+1) = matmul(objective_hessian%x, states%x(:,k+1))
    end do
  end subroutine build_objective_derivatives

end program rom_reduced_bdf1_adjoint_diagnostic
