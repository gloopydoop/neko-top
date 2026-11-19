!> MMA optimizer for a fixed-basis reduced Neko velocity problem.
!!
!! The reduced lifecycle and physical objective/sensitivity evaluation are
!! provided by rom_problem_t.  This type owns only MMA, full-order refresh,
!! validation diagnostics and optimizer logging.  A reduced candidate is never
!! rolled back by the ROM layer: validation reports scientific status but does
!! not veto the engineering optimization trajectory.
module rom_optimizer
  use optimizer, only: optimizer_t
  use mma, only: mma_t
  use problem, only: problem_t
  use design, only: design_t
  use brinkman_design, only: brinkman_design_t
  use simulation_m, only: simulation_t
  use rom_problem, only: rom_problem_t
  use rom_validation, only: rom_validation_metrics_t, rom_validation_accept
  use num_types, only: rp
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use vector, only: vector_t
  use matrix, only: matrix_t
  use device, only: DEVICE_TO_HOST, HOST_TO_DEVICE
  use neko_config, only: NEKO_BCKND_DEVICE
  use logger, only: neko_log
  use utils, only: neko_error
  use comm, only: pe_rank, NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only: MPI_Allreduce, MPI_MAX, MPI_SUM
  implicit none
  private

  public :: rom_optimizer_t

  type, extends(optimizer_t) :: rom_optimizer_t
     type(mma_t) :: mma
     type(rom_problem_t) :: rom_problem
     type(vector_t) :: base_x
     integer :: inner_iterations = 1
     integer :: inner_counter = 0
     real(kind=rp) :: tolerance = 1.0e-3_rp
     real(kind=rp) :: validation_tolerance = 5.0e-2_rp
     real(kind=rp) :: validation_gradient_tolerance = 5.0e-2_rp
     logical :: validate_full_order = .true.
     real(kind=rp) :: mma_scale = 1.0_rp
     logical :: mma_auto_scale = .false.
     logical :: enable_output = .true.
     logical :: supplied_basis_for_step = .false.
     logical :: validation_log_initialized = .false.
     integer :: rom_log_unit = -1
     logical :: rom_log_initialized = .false.
     character(len=16) :: last_validation_status = 'initial'
   contains
     procedure, pass(this) :: init_from_json => rom_optimizer_init_from_json
     procedure, pass(this) :: free => rom_optimizer_free
     procedure, pass(this) :: initialize => rom_optimizer_initialize
     procedure, pass(this), public :: initialize_from_basis => &
          rom_optimizer_initialize_from_basis
     procedure, pass(this), public :: evaluate_fixed_basis => &
          rom_optimizer_evaluate_fixed_basis
     procedure, pass(this), public :: step_from_basis => &
          rom_optimizer_step_from_basis
     procedure, pass(this) :: step => rom_optimizer_step
     procedure, pass(this) :: validate => rom_optimizer_validate
     procedure, pass(this) :: write => rom_optimizer_write
     procedure, pass(this) :: save_checkpoint_components => &
          rom_optimizer_save_checkpoint_components
     procedure, pass(this) :: load_checkpoint_components => &
          rom_optimizer_load_checkpoint_components
  end type rom_optimizer_t

contains

  !> Read base optimizer, MMA, ROM and full-order validation settings.
  subroutine rom_optimizer_init_from_json(this, parameters, problem, design, &
       simulation)
    class(rom_optimizer_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(in) :: design
    type(simulation_t), optional, intent(in) :: simulation
    type(json_file) :: solver_parameters, rom_parameters
    integer :: max_iterations

    call this%free()
    call json_get(parameters, 'optimization.solver', solver_parameters)
    call this%read_base_settings(solver_parameters)
    call json_get_or_default(solver_parameters, 'max_iterations', &
         max_iterations, 20)
    call json_get_or_default(solver_parameters, 'tolerance', this%tolerance, &
         1.0e-3_rp)
    call json_get_or_default(solver_parameters, 'enable_output', &
         this%enable_output, .true.)

    call this%rom_problem%init_from_json(parameters, problem)
    this%inner_iterations = this%rom_problem%inner_iterations
    if ('rom' .in. parameters) then
       call json_get(parameters, 'rom', rom_parameters)
    else
       rom_parameters = solver_parameters
    end if
    call json_get_or_default(rom_parameters, 'validate_full_order', &
         this%validate_full_order, .true.)
    call json_get_or_default(rom_parameters, 'validation_tolerance', &
         this%validation_tolerance, 5.0e-2_rp)
    call json_get_or_default(rom_parameters, 'validation_gradient_tolerance', &
         this%validation_gradient_tolerance, 5.0e-2_rp)
    if (this%inner_iterations < 1) call neko_error('invalid rom_mma inner_iterations')
    this%rom_problem%enable_output = this%enable_output

    call this%base_x%init(design%size())
    call design%get_values(this%base_x)
    call this%mma%init(this%base_x, design%size(), problem%get_n_constraints(), &
         solver_parameters, this%mma_scale, this%mma_auto_scale)
    if (this%enable_output) call rom_optimizer_init_log(this)
    call this%init_base('ROM_MMA', max_iterations)
  end subroutine rom_optimizer_init_from_json

  !> Release MMA, ROM problem and optimizer-owned work arrays.
  subroutine rom_optimizer_free(this)
    class(rom_optimizer_t), intent(inout) :: this
    if (this%rom_log_initialized) then
       close(this%rom_log_unit)
    end if
    this%rom_log_unit = -1
    this%rom_log_initialized = .false.
    call this%rom_problem%free()
    call this%base_x%free()
    call this%mma%free()
    call this%free_base()
    this%supplied_basis_for_step = .false.
    this%validation_log_initialized = .false.
    this%last_validation_status = 'initial'
  end subroutine rom_optimizer_free

  !> Initialise the ROM problem against the live Neko simulation.
  subroutine rom_optimizer_initialize(this, problem, design, simulation)
    class(rom_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), optional, intent(inout) :: simulation
    if (.not. present(simulation)) call neko_error('rom_mma requires a simulation')
    call this%rom_problem%initialize(problem, design, simulation)
  end subroutine rom_optimizer_initialize

  !> Initialise the ROM problem from an externally supplied common basis.
  subroutine rom_optimizer_initialize_from_basis(this, basis, n_modes, &
       n_fields, simulation)
    class(rom_optimizer_t), intent(inout) :: this
    real(kind=rp), intent(in) :: basis(:, :)
    integer, intent(in) :: n_modes, n_fields
    type(simulation_t), intent(inout) :: simulation
    call this%rom_problem%initialize_from_basis(basis, n_modes, n_fields, simulation)
  end subroutine rom_optimizer_initialize_from_basis

  !> Evaluate a fixed-basis ROM and return physical objective sensitivities.
  subroutine rom_optimizer_evaluate_fixed_basis(this, design, objective, &
       gradient, constraint, final_state, chi_gradient)
    class(rom_optimizer_t), intent(inout) :: this
    class(design_t), intent(inout) :: design
    real(kind=rp), intent(out) :: objective, gradient(:), constraint
    real(kind=rp), allocatable, intent(out) :: final_state(:)
    real(kind=rp), intent(out), optional :: chi_gradient(:)
    type(vector_t) :: density, projected_gradient
    real(kind=rp), allocatable :: constraints(:)
    integer :: n, m

    n = design%size()
    m = 1
    if (size(gradient) /= n) call neko_error('ROM diagnostic gradient size mismatch')
    call density%init(n)
    call projected_gradient%init(n)
    allocate(constraints(m))
    call this%rom_problem%compute(design, objective, constraints, final_state)
    constraint = constraints(1)
    call this%rom_problem%compute_sensitivity(density%x)
    if (present(chi_gradient)) then
       if (size(chi_gradient) /= n) call neko_error( &
            'ROM diagnostic chi-gradient size mismatch')
       chi_gradient = density%x
    end if
    call this%rom_problem%project_density(design, density, projected_gradient)
    gradient = projected_gradient%x
    call density%free()
    call projected_gradient%free()
    deallocate(constraints)
  end subroutine rom_optimizer_evaluate_fixed_basis

  !> Run one MMA step using a basis supplied by the caller.
  function rom_optimizer_step_from_basis(this, iter, problem, design, &
       simulation, basis, n_modes, n_fields) result(converged)
    class(rom_optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter, n_modes, n_fields
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), intent(inout) :: simulation
    real(kind=rp), intent(in) :: basis(:, :)
    logical :: converged
    call this%rom_problem%initialize_from_basis(basis, n_modes, n_fields, simulation)
    this%supplied_basis_for_step = .true.
    converged = this%step(iter, problem, design, simulation)
    this%supplied_basis_for_step = .false.
  end function rom_optimizer_step_from_basis

  !> Run ROM MMA and full-order validation diagnostics for one iteration.
  function rom_optimizer_step(this, iter, problem, design, simulation) result(converged)
    class(rom_optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), optional, intent(inout) :: simulation
    type(vector_t) :: x, xold, fval, df0dx, dx, density, old_design
    type(matrix_t) :: dfdx
    real(kind=rp), allocatable :: rom_state(:)
    real(kind=rp) :: rom_objective, actual_objective
    real(kind=rp), allocatable :: predicted_constraints(:), actual_constraints(:), &
         predicted_gradient(:), actual_gradient(:), predicted_constraint_gradient(:, :), &
         actual_constraint_gradient(:, :)
    real(kind=rp) :: local_max, global_max, local_norm, global_norm
    real(kind=rp) :: objective_gradient_error, constraint_gradient_error
    type(rom_validation_metrics_t) :: validation_metrics
    logical :: converged, validation_checked, validation_passed
    integer :: inner, n, m, i

    if (.not. present(simulation)) call neko_error('rom_mma step requires simulation')
    call simulation%reset()
    if (.not. this%supplied_basis_for_step) then
       call this%rom_problem%refresh(simulation)
    end if
    if (this%enable_output) call this%rom_problem%initialize_outputs(design, simulation, iter)

    n = design%size(); m = problem%get_n_constraints()
    call x%init(n); call xold%init(n); call old_design%init(n)
    call fval%init(m); call df0dx%init(n); call dx%init(n); call density%init(n)
    call dfdx%init(m, n)
    call design%get_values(x)
    if (NEKO_BCKND_DEVICE .eq. 1) call x%copy_from(DEVICE_TO_HOST, sync = .true.)
    old_design = x
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call x%copy_from(HOST_TO_DEVICE, sync = .true.)
       call old_design%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if

    ! Evaluate the current design once before the first MMA update.  The
    ! resulting ROM forward/adjoint and gradients are consumed by the first
    ! inner update and then replaced only after the design has changed.
    call this%rom_problem%compute(design, rom_objective, final_state=rom_state)
    call this%rom_problem%compute_sensitivity(density%x)
    call this%rom_problem%compute_volume_constraint(design, fval%x(1), dfdx)
    call this%rom_problem%project_density(design, density, df0dx)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call fval%copy_from(HOST_TO_DEVICE, sync = .true.)
       call dfdx%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if

    do inner = 1, this%inner_iterations
       this%inner_counter = this%inner_counter + 1
       xold = x
       ! mma%update applies its configured xmin/xmax and mma.move_limit when
       ! constructing the subproblem. Do not alter that MMA solution with a
       ! second ROM-specific clip relative to the outer design.
       call this%mma%update(this%inner_counter, x, df0dx, fval, dfdx)
       call design%update_design(x)

       ! Evaluate the updated design exactly once.  This evaluation supplies
       ! the KKT check, output fields, and the cached state/gradient used by
       ! the next inner update; the next loop entrance does not repeat it.
       call this%rom_problem%compute(design, rom_objective, final_state=rom_state)
       call this%rom_problem%compute_sensitivity(density%x)
       call this%rom_problem%compute_volume_constraint(design, fval%x(1), dfdx)
       call this%rom_problem%project_density(design, density, df0dx)
       if (NEKO_BCKND_DEVICE .eq. 1) call fval%copy_from(HOST_TO_DEVICE, sync = .true.)
       if (NEKO_BCKND_DEVICE .eq. 1) call dfdx%copy_from(HOST_TO_DEVICE, sync = .true.)
       call this%mma%KKT(x, df0dx, fval, dfdx)

       if (NEKO_BCKND_DEVICE .eq. 1) call x%copy_from(DEVICE_TO_HOST, sync = .true.)
       local_max = maxval(abs(x%x - xold%x))
       call MPI_Allreduce(local_max, global_max, 1, MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM)
       local_norm = sum((x%x - xold%x) ** 2)
       call MPI_Allreduce(local_norm, global_norm, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
       this%max_design_change = global_max
       this%norm2_design_change = sqrt(global_norm)
       if (NEKO_BCKND_DEVICE .eq. 1) call fval%copy_from(DEVICE_TO_HOST, sync = .true.)
       call rom_optimizer_write_inner(this, iter, inner, rom_objective, fval%x)
       if (this%enable_output .and. allocated(rom_state)) &
            call this%rom_problem%write_outputs(design, simulation, inner, rom_state)
       if (pe_rank == 0) write(*,'(A,I0,A,I0,A,ES12.4,A,ES12.4)') &
            'ROM outer ', iter, ', inner ', inner, ': objective=', rom_objective, &
            ', max design change=', global_max
       if (this%get_stop_design_change() >= 0.0_rp .and. &
            global_max < this%get_stop_design_change()) exit
    end do

    if (allocated(predicted_constraints)) deallocate(predicted_constraints)
    allocate(predicted_constraints(m))
    predicted_constraints = fval%x
    allocate(predicted_gradient(n), actual_gradient(n), predicted_constraint_gradient(m, n), &
         actual_constraint_gradient(m, n))
    predicted_gradient = df0dx%x
    predicted_constraint_gradient = dfdx%x
    call design%get_values(x)
    if (NEKO_BCKND_DEVICE .eq. 1) call x%copy_from(DEVICE_TO_HOST, sync = .true.)
    local_max = maxval(abs(x%x - old_design%x))
    call MPI_Allreduce(local_max, global_max, 1, MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM)
    local_norm = sum((x%x - old_design%x) ** 2)
    call MPI_Allreduce(local_norm, global_norm, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
    this%max_design_change = global_max
    this%norm2_design_change = sqrt(global_norm)

    call simulation%set_design_iteration(iter)
    call problem%compute(design, simulation)
    call problem%get_objective_value(actual_objective)
    allocate(actual_constraints(m))
    call problem%get_constraint_values(fval)
    call fval%copy_from(DEVICE_TO_HOST, sync = .true.)
    actual_constraints = fval%x
    ! Compute the full-order sensitivity before deciding whether the ROM
    ! candidate is trustworthy.  Brinkman sensitivities are fields, so the
    ! exact SEM diagonal conversion is required to compare directional
    ! derivatives with the ROM/MMA convention.
    call problem%compute_sensitivity(design, simulation)
    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(df0dx)
       call des%convert_to_directional_derivative(df0dx)
       call problem%get_constraint_sensitivities(dfdx)
       call des%convert_to_directional_derivative(dfdx)
    class default
       call problem%get_objective_sensitivities(df0dx)
       call problem%get_constraint_sensitivities(dfdx)
    end select
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call df0dx%copy_from(DEVICE_TO_HOST, sync = .true.)
       call dfdx%copy_from(DEVICE_TO_HOST, sync = .true.)
    end if
    actual_gradient = df0dx%x
    actual_constraint_gradient = dfdx%x

    local_norm = sum((actual_gradient - predicted_gradient) ** 2)
    call MPI_Allreduce(local_norm, global_norm, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
    local_norm = sum(actual_gradient ** 2)
    call MPI_Allreduce(local_norm, global_max, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
    objective_gradient_error = sqrt(global_norm) / max(1.0_rp, sqrt(global_max))
    constraint_gradient_error = 0.0_rp
    do i = 1, m
       local_norm = sum((actual_constraint_gradient(i,:) - predicted_constraint_gradient(i,:)) ** 2)
       call MPI_Allreduce(local_norm, global_norm, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
       local_norm = sum(actual_constraint_gradient(i,:) ** 2)
       call MPI_Allreduce(local_norm, global_max, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
       constraint_gradient_error = max(constraint_gradient_error, &
            sqrt(global_norm) / max(1.0_rp, sqrt(global_max)))
    end do
    validation_metrics%objective = abs(actual_objective - rom_objective) / &
         max(1.0_rp, abs(actual_objective))
    validation_metrics%constraint = maxval(abs(actual_constraints - predicted_constraints) / &
         max(1.0_rp, abs(actual_constraints)))
    validation_metrics%objective_gradient = objective_gradient_error
    validation_metrics%constraint_gradient = constraint_gradient_error
    validation_checked = this%validate_full_order
    validation_passed = .true.
    if (validation_checked) then
       validation_passed = rom_validation_accept(validation_metrics, &
            this%validation_tolerance, this%validation_gradient_tolerance)
    end if
    if (.not. validation_checked) then
       this%last_validation_status = 'unchecked'
    else if (validation_passed) then
       this%last_validation_status = 'validated'
    else
       ! Validation remains a strict scientific diagnostic, but it must not
       ! discard the candidate from the engineering optimization trajectory.
       if (pe_rank == 0) call neko_log%warning( &
            'ROM candidate failed full-order validation; retaining candidate')
       this%last_validation_status = 'kept_unvalidated'
    end if
    call write_validation_record(this, iter, rom_objective, actual_objective, &
         predicted_constraints, actual_constraints, validation_metrics, &
         validation_checked, validation_passed)
    ! Full-order validation is reported independently of MMA convergence.  A
    ! failed validation must not restore old_design or prevent the updated
    ! design from being written and used by the next outer iteration.
    converged = this%mma%get_residumax() < this%tolerance

    call x%free(); call xold%free(); call old_design%free(); call fval%free()
    call df0dx%free(); call dx%free(); call density%free(); call dfdx%free()
    if (allocated(rom_state)) deallocate(rom_state)
    if (allocated(predicted_constraints)) deallocate(predicted_constraints)
    if (allocated(actual_constraints)) deallocate(actual_constraints)
    if (allocated(predicted_gradient)) deallocate(predicted_gradient)
    if (allocated(actual_gradient)) deallocate(actual_gradient)
    if (allocated(predicted_constraint_gradient)) deallocate(predicted_constraint_gradient)
    if (allocated(actual_constraint_gradient)) deallocate(actual_constraint_gradient)
  end function rom_optimizer_step


  !> Check that the current full-order constraints are feasible.
  subroutine rom_optimizer_validate(this, problem, design)
    class(rom_optimizer_t), intent(inout) :: this
    class(problem_t), intent(in) :: problem
    class(design_t), intent(in) :: design
    type(vector_t) :: values
    call values%init(problem%get_n_constraints())
    call problem%get_constraint_values(values)
    call values%copy_from(DEVICE_TO_HOST, sync=.true.)
    if (any(values%x > 0.0_rp)) call neko_error('ROM MMA validation failed: constraints are not satisfied')
    call values%free()
  end subroutine rom_optimizer_validate

  !> Initialise the ROM-specific mixed outer/inner optimization history.
  !!
  !! The generic optimizer log is intentionally not used here: its callback is
  !! invoked only once before the outer loop and once after each outer step.
  !! ROM inner iterations are evaluated inside `rom_optimizer_step`, so they
  !! need an explicitly ordered record stream of their own.
  subroutine rom_optimizer_init_log(this)
    class(rom_optimizer_t), intent(inout) :: this
    integer :: ios

    if (pe_rank /= 0) return
    open(newunit=this%rom_log_unit, file='optimization_data.csv', &
         status='replace', action='write', iostat=ios)
    if (ios /= 0) call neko_error('could not open ROM optimization_data.csv')
    write(this%rom_log_unit, '(A)') &
         'phase,outer_iteration,inner_iteration,source,total_objective,'// &
         'viscous_dissipation,viscous_dissipation.weight,'// &
         'brinkman_dissipation,brinkman_dissipation.weight,'// &
         'volume_constraint,volume_constraint.volume,max_design_change,'// &
         'norm2_design_change,status'
    this%rom_log_initialized = .true.
  end subroutine rom_optimizer_init_log

  !> Write one full-order outer record using the live Neko problem log.
  !! @param this ROM optimizer.
  !! @param iter Outer optimizer iteration.
  !! @param problem Full-order Neko problem supplying the values.
  subroutine rom_optimizer_write_outer(this, iter, problem)
    class(rom_optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    class(problem_t), intent(inout) :: problem
    real(kind=rp), allocatable :: values(:)
    real(kind=rp) :: total, viscous, viscous_weight, brinkman, brinkman_weight
    real(kind=rp) :: constraint, volume
    integer :: n

    if (.not. this%rom_log_initialized) return
    n = problem%get_log_size(.true.)
    allocate(values(n))
    call problem%get_log_values(values, .true.)

    total = values(1)
    viscous = 0.0_rp; viscous_weight = 0.0_rp
    brinkman = 0.0_rp; brinkman_weight = 0.0_rp
    constraint = 0.0_rp; volume = 0.0_rp
    if (n >= 3) then
       viscous = values(2); viscous_weight = values(3)
    end if
    if (n >= 5) then
       brinkman = values(4); brinkman_weight = values(5)
    end if
    if (n >= 7) then
       constraint = values(6); volume = values(7)
    end if
    call rom_optimizer_write_record(this, 'outer', iter, 0, 'full_order', &
         total, viscous, viscous_weight, brinkman, brinkman_weight, &
         constraint, volume, this%last_validation_status)
    deallocate(values)
  end subroutine rom_optimizer_write_outer

  !> Write one reduced inner-MMA record after its ROM objective evaluation.
  !! @param this ROM optimizer.
  !! @param outer Outer optimizer iteration containing this inner step.
  !! @param inner Inner reduced MMA iteration.
  !! @param objective Weighted ROM objective value.
  !! @param constraints Reduced constraint values.
  subroutine rom_optimizer_write_inner(this, outer, inner, objective, constraints)
    class(rom_optimizer_t), intent(inout) :: this
    integer, intent(in) :: outer, inner
    real(kind=rp), intent(in) :: objective
    real(kind=rp), intent(in) :: constraints(:)
    real(kind=rp) :: constraint, volume

    if (.not. this%rom_log_initialized) return
    ! Inner records must describe the reduced candidate that was just
    ! evaluated.  The full-order problem is deliberately not queried here:
    ! its state and objective are refreshed only at the outer boundary.
    constraint = constraints(1)
    if (this%rom_problem%volume_is_max) then
       volume = this%rom_problem%volume_limit + constraint
    else
       volume = this%rom_problem%volume_limit - constraint
    end if
    call rom_optimizer_write_record(this, 'inner', outer, inner, 'rom', &
         objective, this%rom_problem%viscous_value, &
         this%rom_problem%viscous_weight, this%rom_problem%brinkman_value, &
         this%rom_problem%brinkman_weight, constraint, volume, &
         'not_validated')
  end subroutine rom_optimizer_write_inner

  !> Append one ordered record to the ROM optimization history.
  !! @param this ROM optimizer.
  !! @param phase `outer` for a full-order refresh or `inner` for ROM MMA.
  !! @param outer Outer optimizer iteration.
  !! @param inner Inner reduced iteration, or zero for an outer record.
  !! @param source Objective source (`full_order` or `rom`).
  !! @param total Weighted total objective.
  !! @param viscous Unweighted viscous objective value.
  !! @param viscous_weight Viscous objective weight.
  !! @param brinkman Unweighted Brinkman objective value.
  !! @param brinkman_weight Brinkman objective weight.
  !! @param constraint Volume constraint value.
  !! @param volume Normalised design volume.
  !! @param status Validation status for the outer record.
  subroutine rom_optimizer_write_record(this, phase, outer, inner, source, &
       total, viscous, viscous_weight, brinkman, brinkman_weight, constraint, &
       volume, status)
    class(rom_optimizer_t), intent(inout) :: this
    character(len=*), intent(in) :: phase, source, status
    integer, intent(in) :: outer, inner
    real(kind=rp), intent(in) :: total, viscous, viscous_weight, brinkman
    real(kind=rp), intent(in) :: brinkman_weight, constraint, volume
    integer :: ios

    if (.not. this%rom_log_initialized) return
    write(this%rom_log_unit, &
         '(A,",",I0,",",I0,",",A,",",ES24.16,",",ES24.16,",",ES24.16,",",'// &
         'ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",'// &
         'ES24.16,",",A)', iostat=ios) phase, outer, inner, source, total, &
         viscous, viscous_weight, brinkman, brinkman_weight, constraint, volume, &
         this%max_design_change, this%norm2_design_change, status
    if (ios /= 0) call neko_error('could not write ROM optimization_data.csv')
  end subroutine rom_optimizer_write_record

  !> Write the ordered ROM/full-order history for one outer callback.
  subroutine rom_optimizer_write(this, iter, problem)
    class(rom_optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    class(problem_t), intent(inout) :: problem
    if (this%enable_output) call rom_optimizer_write_outer(this, iter, problem)
  end subroutine rom_optimizer_write

  !> Append ROM/full-order value and gradient errors to the validation CSV.
  subroutine write_validation_record(this, iter, predicted_objective, &
       actual_objective, predicted_constraints, actual_constraints, metrics, &
       validation_checked, validation_passed)
    class(rom_optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    real(kind=rp), intent(in) :: predicted_objective, actual_objective
    real(kind=rp), intent(in) :: predicted_constraints(:), actual_constraints(:)
    type(rom_validation_metrics_t), intent(in) :: metrics
    logical, intent(in) :: validation_checked, validation_passed
    integer :: unit, ios
    logical :: exists
    if (pe_rank /= 0) return
    inquire(file='rom_validation.csv', exist=exists)
    open(newunit=unit, file='rom_validation.csv', status='unknown', &
         position='append', action='write', iostat=ios)
    if (ios /= 0) call neko_error('could not open rom_validation.csv')
    if (.not. exists) write(unit,'(A)') &
         'outer_iteration,predicted_objective,actual_objective,'// &
         'predicted_constraint,actual_constraint,objective_error,'// &
         'constraint_error,objective_gradient_error,constraint_gradient_error,'//&
         'validation_checked,validation_passed,candidate_kept'
    write(unit,'(I0,8(",",ES24.16),",",L1,",",L1,",",L1)') iter, &
         predicted_objective, &
         actual_objective, predicted_constraints(1), actual_constraints(1), &
         metrics%objective, metrics%constraint, metrics%objective_gradient, &
         metrics%constraint_gradient, validation_checked, validation_passed, .true.
    close(unit)
    this%validation_log_initialized = .true.
  end subroutine write_validation_record

  !> Save the MMA state used by the ROM optimizer checkpoint.
  subroutine rom_optimizer_save_checkpoint_components(this, filename, overwrite)
    class(rom_optimizer_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: overwrite
    call this%mma%save_checkpoint(filename, overwrite)
  end subroutine rom_optimizer_save_checkpoint_components

  !> Restore the MMA state used by the ROM optimizer checkpoint.
  subroutine rom_optimizer_load_checkpoint_components(this, filename)
    class(rom_optimizer_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    call this%mma%load_checkpoint(filename)
  end subroutine rom_optimizer_load_checkpoint_components

end module rom_optimizer
