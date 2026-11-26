!> One-pass CPU smoke of the production algebraic adaptive sensitivity path.
!!
!! The input basis is the coupled physical ROM basis exported as a raw
!! column-major matrix. This exercises the production forward, adjoint and
!! sensitivity entry points once; it is not a derivative verification.
program adaptive_sensitivity_smoke_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use fluid_pnpn, only: fluid_pnpn_t
  use design, only: design_t, design_factory
  use problem, only: problem_t
  use rom_optimizer, only: rom_optimizer_t
  use json_module, only: json_file
  use json_utils, only: json_get
  use json_utils_ext, only: json_read_file
  use vector, only: vector_t
  use num_types, only: rp
  use comm, only: NEKO_COMM, pe_rank
  use mpi_f08, only: MPI_Comm_size
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use utils, only: neko_error
  implicit none

  integer, parameter :: basis_rank = 8, basis_fields = 4
  type(simulation_t) :: simulation
  type(problem_t) :: problem
  type(rom_optimizer_t) :: optimizer
  class(design_t), allocatable :: design
  type(fluid_pnpn_t), pointer :: fluid
  type(json_file) :: parameters, design_parameters
  type(vector_t) :: density
  real(kind=rp), allocatable :: basis(:, :), final_state(:)
  real(kind=rp), allocatable :: constraints(:)
  real(kind=rp) :: objective
  integer :: argc, n, nranks, ierr, unit, file_bytes, expected_bytes
  logical :: found
  character(len=1024) :: case_file, basis_file

  argc = command_argument_count()
  if (argc /= 2) call neko_error( &
       'Usage: adaptive_sensitivity_smoke_diagnostic CASE TRIAL_BASIS.bin')
  call get_command_argument(1, case_file)
  call get_command_argument(2, basis_file)

  call neko_init()
  call neko_top_register_types()
  parameters = json_read_file(trim(case_file))
  ! Keep the case's physical and optimization parameters. The diagnostic
  ! selects the same-space adaptive path required by initialize_from_basis
  ! and enables the sensitivity contraction under test in memory only.
  call parameters%update('unsteady', .false., found)
  if (.not. found) call neko_error('Could not disable the unsteady driver')
  call parameters%update('optimization.solver.rom.basis.forward', &
       'forward', found)
  if (.not. found) call neko_error('Could not select a forward trial basis')
  call parameters%update('optimization.solver.rom.basis.adjoint', &
       'forward', found)
  if (.not. found) call neko_error('Could not select a forward test basis')
  call parameters%update('optimization.solver.rom.basis_transform', &
       'algebraic', found)
  if (.not. found) call neko_error('Could not select algebraic mode adaptation')
  call parameters%update('optimization.solver.rom.differentiate_basis', &
       .true., found)
  if (.not. found) call neko_error('Could not enable basis differentiation')
  call parameters%update('optimization.solver.rom.enable_output', .false., found)
  call parameters%update('optimization.solver.rom.write_trajectory_fields', &
       .false., found)
  if (.not. found) call neko_error('Could not disable trajectory output for the smoke')
  call parameters%update('optimization.solver.rom.export_reduced_system', &
       .false., found)
  if (.not. found) call neko_error('Could not disable reduced-system export')
  call parameters%update('optimization.solver.enable_output', .false., found)
  call user_setup(simulation%neko_case%user)
  call simulation%init(parameters)

  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('Adaptive sensitivity smoke requires Pn/Pn')
  end select
  call json_get(parameters, 'optimization.design', design_parameters)
  call design_factory(design, design_parameters, simulation)
  call problem%init(parameters, design, simulation)
  call optimizer%init_from_json(parameters, problem, design, simulation)
  ! The supplied-basis diagnostic path does not enter the usual optimizer
  ! lifecycle that binds the ROM time grid, so bind it from the unchanged case.
  optimizer%rom_problem%time_start = &
       real(simulation%neko_case%time%start_time, rp)
  optimizer%rom_problem%time_end = real(simulation%neko_case%time%end_time, rp)
  optimizer%rom_problem%dt = real(simulation%neko_case%time%dt, rp)
  optimizer%rom_problem%n_steps = nint((optimizer%rom_problem%time_end - &
       optimizer%rom_problem%time_start) / optimizer%rom_problem%dt)
  if (optimizer%rom_problem%dt <= 0.0_rp .or. &
       optimizer%rom_problem%n_steps < 1) &
       call neko_error('Smoke case has an invalid ROM time grid')
  if (.not. optimizer%rom_problem%experimental_algebraic .or. &
       .not. optimizer%rom_problem%differentiate_basis .or. &
       optimizer%rom_problem%sensitivity_source /= 'exact_algebraic') &
       call neko_error('Smoke setup did not select exact algebraic sensitivity')

  call MPI_Comm_size(NEKO_COMM, nranks, ierr)
  if (nranks /= 1) call neko_error( &
       'Archived basis smoke input is partitioned for one MPI rank')
  n = fluid%u%dof%size()
  expected_bytes = basis_fields * n * basis_rank * storage_size(0.0_rp) / 8
  inquire(file=trim(basis_file), size=file_bytes, iostat=ierr)
  if (ierr /= 0 .or. file_bytes /= expected_bytes) &
       call neko_error('Trial basis file does not match case mesh/layout/rank')
  allocate(basis(basis_fields*n,basis_rank))
  open(newunit=unit, file=trim(basis_file), status='old', action='read', &
       access='stream', form='unformatted', iostat=ierr)
  if (ierr /= 0) call neko_error('Could not open trial basis binary')
  read(unit, iostat=ierr) basis
  close(unit)
  if (ierr /= 0) call neko_error('Could not read trial basis binary')

  call optimizer%initialize_from_basis(basis, basis_rank, basis_fields, simulation)
  allocate(constraints(1))
  call optimizer%rom_problem%compute(design, objective, constraints, final_state)
  if (optimizer%rom_problem%scalar_transport_cache_hit /= 1.0_rp) &
       call neko_error('Adaptive smoke did not use the scalar transport cache')
  call density%init(n)
  call optimizer%rom_problem%compute_sensitivity(density)
  if (.not. ieee_is_finite(objective) .or. &
       .not. all(ieee_is_finite(final_state)) .or. &
       .not. all(ieee_is_finite(density%x))) &
       call neko_error('Adaptive sensitivity smoke produced non-finite output')
  if (pe_rank == 0) write(*,'(A,I0,A,ES16.8)') &
       'Adaptive exact sensitivity smoke completed, modes=', basis_rank, &
       ', objective=', objective

  call density%free()
  deallocate(final_state, constraints, basis)
  call optimizer%free()
  call problem%free()
  deallocate(design)
  call neko_finalize()
end program adaptive_sensitivity_smoke_diagnostic
