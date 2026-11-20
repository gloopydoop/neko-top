program rom_model_smoke
  use num_types, only: rp
  use neko, only: neko_init, neko_finalize
  use vector, only: vector_t
  use matrix, only: matrix_t
  use rom_model, only: rom_model_t
  implicit none

  integer, parameter :: n = 3, n_steps = 5
  integer :: i, j, k
  real(kind=rp), parameter :: dt = 0.07_rp
  real(kind=rp), parameter :: eps = 2.0e-6_rp
  real(kind=rp), parameter :: direct_design = 0.37_rp
  real(kind=rp) :: eta_plus, eta_minus, fd, predicted
  type(rom_model_t) :: model, model_plus, model_minus
  type(matrix_t) :: mass, linear, implicit_linear, states, states_plus, states_minus
  type(matrix_t) :: running_gradient, adjoint
  type(matrix_t) :: quadratic(n)
  type(vector_t) :: constant, initial, terminal_gradient, direction
  type(vector_t) :: constant_plus, constant_minus
  type(matrix_t) :: linear_plus, linear_minus
  type(vector_t) :: design_constant
  type(vector_t) :: design_adjoint
  type(matrix_t) :: design_linear
  real(kind=rp) :: objective_plus, objective_minus

  call neko_init()

  call mass%init(n, n)
  call linear%init(n, n)
  call implicit_linear%init(n, n)
  call constant%init(n)
  do i = 1, n
     mass%x(i, i) = 1.0_rp + 0.2_rp * real(i, rp)
     linear%x(i, i) = -0.11_rp * real(i, rp)
     implicit_linear%x(i, i) = -0.04_rp * real(i, rp)
     constant%x(i) = 0.03_rp * real(i, rp)
  end do
  linear%x(1, 2) = 0.08_rp
  linear%x(2, 3) = -0.05_rp
  implicit_linear%x(1, 3) = 0.02_rp
  do i = 1, n
     call quadratic(i)%init(n, n)
     do j = 1, n
        quadratic(i)%x(j, j) = 0.01_rp * real(i + j, rp)
        if (j < n) quadratic(i)%x(j, j + 1) = -0.004_rp * real(i, rp)
     end do
  end do
  call model%init(mass, constant, linear, quadratic, implicit_linear)

  call initial%init(n)
  initial%x = [0.22_rp, -0.13_rp, 0.31_rp]
  call terminal_gradient%init(n)
  terminal_gradient%x = [-0.17_rp, 0.23_rp, 0.11_rp]
  call running_gradient%init(n, n_steps)
  do k = 1, n_steps
     running_gradient%x(:, k) = [0.03_rp * real(k, rp), &
          -0.02_rp * real(k + 1, rp), 0.015_rp * real(2 * k, rp)]
  end do
  call states%init(n, n_steps + 1)
  call model%forward_bdf1(initial, dt, states)
  call model%discrete_adjoint_bdf1(states, dt, running_gradient, &
       terminal_gradient, adjoint)

  call direction%init(n)
  direction%x = [0.13_rp, -0.21_rp, 0.09_rp]
  call states_plus%init(n, n_steps + 1)
  call states_minus%init(n, n_steps + 1)
  call model%forward_bdf1(initial, dt, states_plus)
  states_plus%x(:, 1) = initial%x + eps * direction%x
  call model%forward_bdf1(vector_from_column(states_plus, 1), dt, states_plus)
  call model%forward_bdf1(initial, dt, states_minus)
  states_minus%x(:, 1) = initial%x - eps * direction%x
  call model%forward_bdf1(vector_from_column(states_minus, 1), dt, states_minus)
  fd = (trajectory_objective(states_plus, running_gradient, terminal_gradient, dt) - &
       trajectory_objective(states_minus, running_gradient, terminal_gradient, dt)) / (2.0_rp * eps)
  predicted = dot_product(adjoint%x(:, 1), direction%x)
  call assert_close('BDF1 initial-state adjoint', predicted, fd, 2.0e-7_rp)

  call design_constant%init(n)
  call design_adjoint%init(n)
  design_constant%x = [0.09_rp, -0.04_rp, 0.07_rp]
  call design_linear%init(n, n)
  design_linear%x = reshape([0.02_rp, -0.01_rp, 0.00_rp, &
       0.00_rp, 0.03_rp, -0.02_rp, 0.01_rp, 0.00_rp, 0.015_rp], [n, n])
  call constant_plus%init(n); call constant_minus%init(n)
  call linear_plus%init(n, n); call linear_minus%init(n, n)
  eta_plus = eps
  eta_minus = -eps
  constant_plus%x = constant%x + eta_plus * design_constant%x
  constant_minus%x = constant%x + eta_minus * design_constant%x
  linear_plus%x = linear%x + eta_plus * design_linear%x
  linear_minus%x = linear%x + eta_minus * design_linear%x
  call model_plus%init(mass, constant_plus, linear_plus, quadratic, implicit_linear)
  call model_minus%init(mass, constant_minus, linear_minus, quadratic, implicit_linear)
  call model_plus%forward_bdf1(initial, dt, states_plus)
  call model_minus%forward_bdf1(initial, dt, states_minus)
  objective_plus = trajectory_objective(states_plus, running_gradient, terminal_gradient, dt) + &
       direct_design * eta_plus
  objective_minus = trajectory_objective(states_minus, running_gradient, terminal_gradient, dt) + &
       direct_design * eta_minus
  fd = (objective_plus - objective_minus) / (2.0_rp * eps)
  predicted = direct_design
  do k = 1, n_steps
     direction%x = design_constant%x + matmul(design_linear%x, states%x(:, k))
     call model%solve_step_transpose(dt, adjoint_column(adjoint, k + 1), design_adjoint)
     predicted = predicted + dt * dot_product(design_adjoint%x, direction%x)
  end do
  call assert_close('BDF1 design adjoint', predicted, fd, 3.0e-7_rp)

  call model_minus%free(); call model_plus%free(); call model%free()
  call mass%free(); call linear%free(); call implicit_linear%free()
  call constant%free(); call constant_plus%free(); call constant_minus%free()
  call initial%free(); call terminal_gradient%free(); call direction%free()
  call design_constant%free(); call design_linear%free()
  call design_adjoint%free()
  call states%free(); call states_plus%free(); call states_minus%free()
  call running_gradient%free(); call adjoint%free()
  do i = 1, n
     call quadratic(i)%free()
  end do
  call neko_finalize()
  write(*, '(a)') 'rom_model_smoke: PASS'

contains

  !> Copy one matrix column into a reduced vector container.
  function vector_from_column(matrix, column) result(vector)
    type(matrix_t), intent(in) :: matrix
    integer, intent(in) :: column
    type(vector_t) :: vector
    call vector%init(matrix%get_nrows())
    vector%x = matrix%x(:, column)
  end function vector_from_column

  !> Return one reduced adjoint column as a vector container.
  function adjoint_column(matrix, column) result(vector)
    type(matrix_t), intent(in) :: matrix
    integer, intent(in) :: column
    type(vector_t) :: vector
    call vector%init(matrix%get_nrows())
    vector%x = matrix%x(:, column)
  end function adjoint_column

  !> Evaluate the scalar manufactured trajectory objective.
  function trajectory_objective(states, running, terminal, timestep) result(value)
    type(matrix_t), intent(in) :: states, running
    type(vector_t), intent(in) :: terminal
    real(kind=rp), intent(in) :: timestep
    real(kind=rp) :: value
    integer :: step
    value = dot_product(terminal%x, states%x(:, states%get_ncols()))
    do step = 1, running%get_ncols()
       value = value + timestep * dot_product(running%x(:, step), states%x(:, step))
    end do
  end function trajectory_objective

  !> Abort the smoke test if two scalar diagnostics differ too much.
  subroutine assert_close(label, expected, actual, tolerance)
    character(len=*), intent(in) :: label
    real(kind=rp), intent(in) :: expected, actual, tolerance
    real(kind=rp) :: error
    error = abs(expected - actual)
    if (error > tolerance * max(1.0_rp, abs(expected), abs(actual))) then
       write(*, '(a,2(es16.8,1x),a,es16.8)') trim(label)//': FAIL expected/actual=', &
            expected, actual, ' error=', error
       error stop 1
    end if
    write(*, '(a,2(es16.8,1x),a,es16.8)') trim(label)//': PASS expected/actual=', &
         expected, actual, ' error=', error
  end subroutine assert_close

end program rom_model_smoke
