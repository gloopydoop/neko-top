!> @file scalar_objective_equivalence_diagnostic.f90
!! @brief Compare the ROM scalar-mixing quadratic form with Neko-TOP.
program rom_scalar_objective_equivalence_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use fluid_pnpn, only: fluid_pnpn_t
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use json_utils, only: json_get
  use design, only: design_t, design_factory
  use scalar_mixing_objective, only: scalar_mixing_objective_t
  use point_zone, only: point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  use rom_array_test_adapters, only: test_assemble_scalar_mixing, &
       test_scalar_mixing_value
  use comm, only: NEKO_COMM, MPI_REAL_PRECISION, pe_rank
  use mpi_f08, only: MPI_Allreduce, MPI_MAX, MPI_SUM, MPI_IN_PLACE
  use num_types, only: rp
  use vector, only: vector_t
  use utils, only: neko_error
  implicit none

  type(simulation_t) :: simulation
  type(fluid_pnpn_t), pointer :: fluid
  class(design_t), allocatable :: design
  type(scalar_mixing_objective_t) :: scalar_objective
  type(json_file) :: parameters, design_parameters
  class(point_zone_t), pointer :: objective_zone
  real(kind=rp), allocatable :: weights(:), mask(:), lift(:), modes(:, :)
  real(kind=rp), allocatable :: state(:), direction(:), phi_base(:), phi_plus(:), phi_minus(:)
  real(kind=rp), allocatable :: x_base(:), x_plus(:), x_minus(:), design_direction(:)
  real(kind=rp), allocatable :: linear(:), quadratic(:, :)
  real(kind=rp) :: constant, domain, rom_value, neko_value
  real(kind=rp) :: neko_plus, neko_minus, fd_state, analytic_state, state_error
  real(kind=rp) :: design_plus, design_minus, fd_design, global_error
  real(kind=rp) :: reference, eps_state, eps_design
  real(kind=rp) :: xcoord, ycoord
  integer :: argc, n, r, i
  logical :: found
  character(len=1024) :: case_file
  type(vector_t) :: x_vector

  argc = command_argument_count()
  if (argc < 1) call neko_error( &
       'Usage: ROM_rugby_ball_scalar_objective_equivalence_diagnostic CASE')
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
     call neko_error('Scalar objective diagnostic requires a Pn/Pn fluid')
  end select
  if (.not. associated(simulation%scalars)) call neko_error( &
       'Scalar objective diagnostic requires a passive scalar')

  call json_get(parameters, 'optimization.design', design_parameters)
  call design_factory(design, design_parameters, simulation)

  call scalar_objective%init_from_attributes(design, simulation, 1.0_rp, &
       'Scalar Mixing diagnostic', 'objective_domain', 0.5_rp, 's', 0.0_rp, huge(0.0_rp))

  n = fluid%u%dof%size()
  r = 3
  reference = 0.5_rp
  allocate(weights(n), mask(n), lift(n), modes(n,r), state(r), direction(r))
  allocate(phi_base(n), phi_plus(n), phi_minus(n))
  allocate(x_base(n), x_plus(n), x_minus(n))
  allocate(linear(r), quadratic(r,r))
  weights = reshape(fluid%c_Xh%B, [n])
  mask = 0.0_rp
  objective_zone => neko_point_zone_registry%get_point_zone('objective_domain')
  mask(objective_zone%mask%get()) = 1.0_rp
  domain = sum(weights * mask)
  call MPI_Allreduce(MPI_IN_PLACE, domain, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)

  ! Use a nontrivial spatial lift and three independent scalar modes.  The
  ! physical field passed to Neko and the reduced reconstruction below are
  ! exactly the same field, so this test isolates objective algebra rather
  ! than POD quality or time integration.
  do i = 1, n
     xcoord = fluid%u%dof%x%x(i,1,1,1)
     ycoord = fluid%u%dof%y%x(i,1,1,1)
     lift(i) = 0.13_rp + 0.07_rp * sin(0.71_rp*xcoord) * cos(0.37_rp*ycoord)
     modes(i,1) = 0.31_rp * cos(0.41_rp*xcoord) + 0.04_rp * ycoord
     modes(i,2) = 0.17_rp * sin(0.53_rp*xcoord + 0.19_rp*ycoord)
     modes(i,3) = 0.09_rp * cos(0.29_rp*xcoord - 0.47_rp*ycoord)
     x_base(i) = 0.20_rp + 0.60_rp * (0.5_rp + 0.5_rp * &
          sin(0.23_rp*xcoord) * cos(0.31_rp*ycoord))
  end do
  state = [0.27_rp, -0.18_rp, 0.11_rp]
  phi_base = lift + matmul(modes, state)
  simulation%scalars%scalar_fields(1)%scalar%s%x = reshape(phi_base, &
       shape(simulation%scalars%scalar_fields(1)%scalar%s%x))

  call test_assemble_scalar_mixing(weights, mask, lift, modes, reference, &
       constant, linear, quadratic)
  rom_value = test_scalar_mixing_value(state, constant, linear, quadratic, domain)
  call scalar_objective%update_value(design)
  neko_value = scalar_objective%get_value()
  global_error = abs(rom_value - neko_value)
  call MPI_Allreduce(MPI_IN_PLACE, global_error, 1, MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM)
  if (pe_rank == 0) write(*,'(A,ES14.6,A,ES14.6,A,ES14.6)') &
       'Scalar objective ROM/Neko/absolute error = ', rom_value, ' / ', neko_value, &
       ' / ', global_error
  if (global_error > 5.0e-12_rp) call neko_error( &
       'Scalar objective value equivalence gate failed')

  ! Compare the reduced analytic derivative with centred finite differences of
  ! the original objective evaluated on the identical reconstructed fields.
  direction = [0.19_rp, -0.23_rp, 0.17_rp]
  eps_state = 1.0e-6_rp
  phi_plus = lift + matmul(modes, state + eps_state*direction)
  phi_minus = lift + matmul(modes, state - eps_state*direction)
  simulation%scalars%scalar_fields(1)%scalar%s%x = reshape(phi_plus, &
       shape(simulation%scalars%scalar_fields(1)%scalar%s%x))
  call scalar_objective%update_value(design)
  neko_plus = scalar_objective%get_value()
  simulation%scalars%scalar_fields(1)%scalar%s%x = reshape(phi_minus, &
       shape(simulation%scalars%scalar_fields(1)%scalar%s%x))
  call scalar_objective%update_value(design)
  neko_minus = scalar_objective%get_value()
  fd_state = (neko_plus - neko_minus) / (2.0_rp*eps_state)
  analytic_state = dot_product(direction, linear + matmul(quadratic, state)) / domain
  state_error = abs(fd_state - analytic_state)
  if (pe_rank == 0) write(*,'(A,ES14.6,A,ES14.6,A,ES14.6)') &
       'Scalar reduced derivative FD/ROM/absolute error = ', fd_state, ' / ', &
       analytic_state, ' / ', state_error
  if (state_error > 5.0e-7_rp * max(1.0_rp, abs(fd_state))) call neko_error( &
       'Scalar objective reduced derivative gate failed')

  ! The native scalar-mixing objective has no explicit dependence on the
  ! Brinkman field: update_sensitivity is intentionally empty.  Verify that
  ! changing the mapped design while holding the physical scalar fixed gives
  ! the same zero direct derivative, rather than silently introducing a ROM
  ! design term.
  call x_vector%init(n)
  x_vector%x = x_base
  call design%update_design(x_vector)
  eps_design = 1.0e-6_rp
  allocate(design_direction(n))
  do i = 1, n
     design_direction(i) = 0.11_rp + 0.03_rp * sin(real(i, rp))
  end do
  x_plus = x_base + eps_design * design_direction
  x_minus = x_base - eps_design * design_direction
  x_vector%x = x_plus
  call design%update_design(x_vector)
  call scalar_objective%update_value(design)
  design_plus = scalar_objective%get_value()
  x_vector%x = x_minus
  call design%update_design(x_vector)
  call scalar_objective%update_value(design)
  design_minus = scalar_objective%get_value()
  fd_design = (design_plus - design_minus) / (2.0_rp*eps_design)
  if (pe_rank == 0) write(*,'(A,ES14.6)') &
       'Scalar direct mapped-design derivative FD/analytic = ', fd_design
  if (abs(fd_design) > 5.0e-12_rp) call neko_error( &
       'Scalar objective unexpectedly depends directly on design')

  x_vector%x = x_base
  call design%update_design(x_vector)
  simulation%scalars%scalar_fields(1)%scalar%s%x = reshape(phi_base, &
       shape(simulation%scalars%scalar_fields(1)%scalar%s%x))
  call scalar_objective%free()
  call x_vector%free()
  call design%free()
  deallocate(weights, mask, lift, modes, state, direction, phi_base, phi_plus, phi_minus)
  deallocate(x_base, x_plus, x_minus, design_direction, linear, quadratic)
  call simulation%free()
  call neko_finalize()

contains

end program rom_scalar_objective_equivalence_diagnostic
