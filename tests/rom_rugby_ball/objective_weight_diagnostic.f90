!> Focused ROM objective-weight and accumulator diagnostic.
!!
!! This test constructs the public reduced objective state directly.  It does
!! not require a mesh or a POD file, so it isolates the time accumulation and
!! weight handling in rom_problem_t from basis and full-order errors.
program rom_objective_weight_diagnostic
  use neko, only: neko_init, neko_finalize
  use num_types, only: rp
  use vector, only: vector_t
  use rom_problem, only: rom_problem_t
  use utils, only: neko_error
  implicit none

  type(rom_problem_t) :: problem
  type(vector_t) :: density
  real(kind=rp) :: objective, expected, raw_v, raw_b, raw_s
  real(kind=rp) :: running_norm
  integer :: i, j

  call neko_init()

  problem%n_steps = 2
  problem%dt = 0.5_rp
  problem%time_start = 0.0_rp
  problem%time_end = 1.0_rp
  problem%forward_ready = .true.
  problem%has_viscous = .true.
  problem%has_brinkman = .true.
  problem%has_scalar_mixing = .true.
  problem%viscous_start = 0.0_rp
  problem%viscous_end = 1.0_rp
  problem%brinkman_start = 0.0_rp
  problem%brinkman_end = 1.0_rp
  problem%scalar_mixing_start = 0.0_rp
  problem%scalar_mixing_end = 1.0_rp
  problem%viscous_domain = 2.0_rp
  problem%brinkman_domain = 2.0_rp
  problem%scalar_mixing_domain = 2.0_rp

  call problem%states%init(2, 3)
  problem%states%x = reshape([0.0_rp, 0.0_rp, 0.4_rp, -0.2_rp, &
       0.7_rp, 0.3_rp], [2, 3])
  call problem%running_gradient%init(2, 2)
  call problem%terminal_gradient%init(2)
  call problem%cost_sum_a%init(2)
  call problem%cost_sum_aa%init(2, 2)
  call problem%weights%init(2); problem%weights%x = 1.0_rp
  call problem%brinkman_mask%init(2); problem%brinkman_mask%x = 1.0_rp
  call problem%lambda_sum%init(2); problem%lambda_sum%x = 0.0_rp
  call problem%correlation%init(2, 2); problem%correlation%x = 0.0_rp
  allocate(problem%velocity_modes(3,2), problem%test_velocity_modes(3,2))
  do i = 1, 3
     call problem%lift_velocity(i)%init(2); problem%lift_velocity(i)%x = 0.0_rp
     do j = 1, 2
        call problem%velocity_modes(i,j)%init(2)
        call problem%test_velocity_modes(i,j)%init(2)
        problem%velocity_modes(i,j)%x = 0.1_rp * real(i+j, rp)
        problem%test_velocity_modes(i,j)%x = problem%velocity_modes(i,j)%x
     end do
  end do

  call problem%viscous_constant%init(1)
  call problem%viscous_linear%init(2)
  call problem%viscous_quadratic%init(2, 2)
  problem%viscous_constant%x(1) = 0.8_rp
  problem%viscous_linear%x = [0.2_rp, -0.1_rp]
  problem%viscous_quadratic%x = reshape([1.2_rp, 0.1_rp, 0.1_rp, 0.9_rp], [2, 2])

  problem%candidate_objective_constant = 0.6_rp
  call problem%candidate_objective_lift%init(2)
  call problem%candidate_objective_brinkman%init(2, 2)
  problem%candidate_objective_lift%x = [0.15_rp, -0.2_rp]
  problem%candidate_objective_brinkman%x = reshape([0.7_rp, 0.05_rp, 0.05_rp, 0.5_rp], [2, 2])

  call problem%scalar_mixing_constant%init(1)
  call problem%scalar_mixing_linear%init(2)
  call problem%scalar_mixing_quadratic%init(2, 2)
  problem%scalar_mixing_constant%x(1) = 0.4_rp
  problem%scalar_mixing_linear%x = [-0.1_rp, 0.25_rp]
  problem%scalar_mixing_quadratic%x = reshape([0.8_rp, -0.02_rp, -0.02_rp, 0.6_rp], [2, 2])

  ! All three objectives are configured, but all have zero optimization
  ! weight.  Raw physical values must still be accumulated and the weighted
  ! total and state forcing must be exactly zero.
  problem%viscous_weight = 0.0_rp
  problem%brinkman_weight = 0.0_rp
  problem%scalar_mixing_weight = 0.0_rp
  call problem%compute_objectives(objective)
  raw_v = problem%viscous_value
  raw_b = problem%brinkman_value
  raw_s = problem%scalar_mixing_value
  if (min(raw_v, raw_b, raw_s) <= 0.0_rp) call neko_error( &
       'zero-weight objective did not retain a raw physical value')
  if (abs(objective) > 1.0e-14_rp) call neko_error( &
       'zero-weight objectives contaminated total ROM objective')
  running_norm = maxval(abs(problem%running_gradient%x))
  running_norm = max(running_norm, maxval(abs(problem%terminal_gradient%x)))
  if (running_norm > 1.0e-14_rp) call neko_error( &
       'zero-weight objectives contaminated ROM state forcing')
  call problem%compute_density(density)
  if (maxval(abs(density%x)) > 1.0e-14_rp) call neko_error( &
       'zero-weight Brinkman objective contaminated ROM design density')

  ! The weighted total must be a linear combination of the raw, normalized
  ! component values.  Repeating compute_objectives also checks that all
  ! objective accumulators are reset between candidate designs.
  problem%viscous_weight = 120.0_rp
  problem%brinkman_weight = 0.0_rp
  problem%scalar_mixing_weight = 0.0_rp
  call problem%compute_objectives(objective)
  if (abs(objective - 120.0_rp*raw_v) > 1.0e-13_rp) call neko_error( &
       'viscous-only objective weight failed')

  problem%viscous_weight = 0.0_rp
  problem%brinkman_weight = 15.0_rp
  problem%scalar_mixing_weight = 0.0_rp
  call problem%compute_objectives(objective)
  if (abs(objective - 15.0_rp*raw_b) > 1.0e-13_rp) call neko_error( &
       'Brinkman-only objective weight failed')

  problem%viscous_weight = 0.0_rp
  problem%brinkman_weight = 0.0_rp
  problem%scalar_mixing_weight = 15.0_rp
  call problem%compute_objectives(objective)
  if (abs(objective - 15.0_rp*raw_s) > 1.0e-13_rp) call neko_error( &
       'scalar-only objective weight failed')

  problem%viscous_weight = 120.0_rp
  problem%brinkman_weight = 15.0_rp
  problem%scalar_mixing_weight = 15.0_rp
  call problem%compute_objectives(objective)
  expected = 120.0_rp * raw_v + 15.0_rp * raw_b + 15.0_rp * raw_s
  if (abs(objective - expected) > 1.0e-13_rp * max(1.0_rp, abs(expected))) &
       call neko_error('weighted ROM objective failed linear superposition')

  problem%viscous_weight = 120.0_rp
  problem%brinkman_weight = 150.0_rp
  problem%scalar_mixing_weight = 150.0_rp
  call problem%compute_objectives(objective)
  expected = 120.0_rp * raw_v + 150.0_rp * raw_b + 150.0_rp * raw_s
  if (abs(objective - expected) > 1.0e-13_rp * max(1.0_rp, abs(expected))) &
       call neko_error('high Brinkman/scalar objective weights failed')

  problem%brinkman_weight = -15.0_rp
  problem%scalar_mixing_weight = 15.0_rp
  call problem%compute_objectives(objective)
  expected = 120.0_rp * raw_v - 15.0_rp * raw_b + 15.0_rp * raw_s
  if (abs(objective - expected) > 1.0e-13_rp * max(1.0_rp, abs(expected))) &
       call neko_error('signed ROM objective weight failed linear superposition')

  if (abs(problem%viscous_value - raw_v) > 1.0e-14_rp .or. &
       abs(problem%brinkman_value - raw_b) > 1.0e-14_rp .or. &
       abs(problem%scalar_mixing_value - raw_s) > 1.0e-14_rp) call neko_error( &
       'raw ROM objective accumulator was retained across evaluations')

  write(*,'(A,3(ES14.6,1X),A,ES14.6)') &
       'ROM zero-weight raw values (viscous brinkman scalar) / signed total = ', &
       raw_v, raw_b, raw_s, ' / ', objective
  call problem%free()
  call density%free()
  call neko_finalize()
end program rom_objective_weight_diagnostic
