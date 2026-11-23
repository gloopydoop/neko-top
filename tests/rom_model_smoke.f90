program rom_model_smoke
  use num_types, only: rp
  use neko, only: neko_init, neko_finalize
  use vector, only: vector_t
  use matrix, only: matrix_t
  use rom_model, only: rom_model_t
  implicit none

  integer, parameter :: n = 3, e = 4, n_steps = 5
  integer :: i, j, k
  real(kind=rp), parameter :: dt = 0.07_rp
  real(kind=rp), parameter :: eps = 2.0e-6_rp
  real(kind=rp), parameter :: direct_design = 0.37_rp
  real(kind=rp) :: eta_plus, eta_minus, fd, predicted
  type(rom_model_t) :: model, model_plus, model_minus
  type(rom_model_t) :: enriched_model, enriched_reference
  type(matrix_t) :: mass, linear, implicit_linear, states, states_plus, states_minus
  type(matrix_t) :: running_gradient, adjoint
  type(matrix_t) :: quadratic(n)
  type(matrix_t) :: enriched_mass, enriched_implicit, enriched_linear
  type(matrix_t) :: enriched_states, enriched_running, enriched_adjoint
  type(matrix_t) :: enriched_reference_adjoint
  type(matrix_t) :: enriched_quadratic(e), time_jacobian(n_steps)
  type(vector_t) :: constant, initial, terminal_gradient, direction
  type(vector_t) :: constant_plus, constant_minus
  type(matrix_t) :: linear_plus, linear_minus
  type(vector_t) :: design_constant
  type(vector_t) :: design_adjoint
  type(matrix_t) :: design_linear
  real(kind=rp) :: objective_plus, objective_minus
  type(vector_t) :: enriched_initial, enriched_terminal, enriched_zero_constant
  integer :: enriched_i

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
  ! Exercise the nonsymmetric cross-mass path used by Petrov--Galerkin.
  mass%x(1,2) = 0.06_rp
  mass%x(2,1) = -0.025_rp
  mass%x(2,3) = 0.04_rp
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

  ! Verify the independent, larger-dimensional time-varying adjoint map.
  ! With a constant Jacobian and identity mass, it must agree with the
  ! existing BDF1 discrete adjoint of the corresponding enriched model.
  call enriched_mass%init(e, e)
  call enriched_implicit%init(e, e)
  call enriched_linear%init(e, e)
  call enriched_zero_constant%init(e)
  enriched_zero_constant%x = 0.0_rp
  do enriched_i = 1, e
     enriched_mass%x(enriched_i, enriched_i) = 1.0_rp
     enriched_implicit%x(enriched_i, enriched_i) = -0.025_rp * real(enriched_i, rp)
     enriched_linear%x(enriched_i, enriched_i) = -0.08_rp * real(enriched_i, rp)
     if (enriched_i < e) enriched_linear%x(enriched_i, enriched_i + 1) = 0.03_rp
     call enriched_quadratic(enriched_i)%init(e, e)
     enriched_quadratic(enriched_i)%x = 0.0_rp
  end do
  ! The time Jacobian is supplied as the residual derivative, as required by
  ! rom_model_discrete_adjoint_bdf1_timevarying.
  do enriched_i = 1, n_steps
     call time_jacobian(enriched_i)%init(e, e)
     time_jacobian(enriched_i)%x = enriched_linear%x
  end do
  call enriched_model%init_linearized_adjoint(enriched_mass, enriched_implicit)
  call enriched_reference%init(enriched_mass, enriched_zero_constant, enriched_linear, &
       enriched_quadratic, enriched_implicit)
  call enriched_initial%init(e)
  enriched_initial%x = [0.11_rp, -0.07_rp, 0.19_rp, 0.05_rp]
  call enriched_terminal%init(e)
  enriched_terminal%x = [-0.04_rp, 0.12_rp, 0.06_rp, -0.09_rp]
  call enriched_running%init(e, n_steps)
  do enriched_i = 1, n_steps
     enriched_running%x(:, enriched_i) = [0.01_rp * real(enriched_i, rp), &
          -0.015_rp, 0.02_rp, 0.005_rp * real(enriched_i, rp)]
  end do
  call enriched_states%init(e, n_steps + 1)
  call enriched_reference%forward_bdf1(enriched_initial, dt, enriched_states)
  call enriched_reference%discrete_adjoint_bdf1(enriched_states, dt, enriched_running, &
       enriched_terminal, enriched_reference_adjoint)
  call enriched_model%discrete_adjoint_bdf1_timevarying(dt, enriched_running, &
       enriched_terminal, time_jacobian, enriched_adjoint)
  if (enriched_adjoint%get_nrows() /= e .or. &
       enriched_adjoint%get_ncols() /= n_steps + 1) then
     error stop 'time-varying enriched adjoint has the wrong dimension'
  end if
  call assert_matrix_close('time-varying enriched BDF1 adjoint', &
       enriched_reference_adjoint, enriched_adjoint, 3.0e-7_rp)

  call model_minus%free(); call model_plus%free(); call model%free()
  call enriched_reference%free(); call enriched_model%free()
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
  call enriched_mass%free(); call enriched_implicit%free(); call enriched_linear%free()
  call enriched_zero_constant%free(); call enriched_initial%free(); call enriched_terminal%free()
  call enriched_states%free(); call enriched_running%free(); call enriched_adjoint%free()
  call enriched_reference_adjoint%free()
  do enriched_i = 1, e
     call enriched_quadratic(enriched_i)%free()
  end do
  do enriched_i = 1, n_steps
     call time_jacobian(enriched_i)%free()
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

  !> Abort the smoke test if two reduced matrix diagnostics differ too much.
  subroutine assert_matrix_close(label, expected, actual, tolerance)
    character(len=*), intent(in) :: label
    type(matrix_t), intent(in) :: expected, actual
    real(kind=rp), intent(in) :: tolerance
    real(kind=rp) :: error, scale
    if (any(expected%get_dims() /= actual%get_dims())) then
       error stop trim(label)//': matrix dimensions differ'
    end if
    error = maxval(abs(expected%x - actual%x))
    scale = max(1.0_rp, maxval(abs(expected%x)), maxval(abs(actual%x)))
    if (error > tolerance * scale) then
       write(*, '(a,es16.8)') trim(label)//': FAIL error=', error
       error stop 1
    end if
    write(*, '(a,es16.8)') trim(label)//': PASS error=', error
  end subroutine assert_matrix_close

end program rom_model_smoke
