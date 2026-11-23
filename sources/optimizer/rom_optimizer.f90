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
  use constraint, only: constraint_t
  use dummy_constraint, only: dummy_constraint_t
  use simulation_m, only: simulation_t
  use rom_problem, only: rom_problem_t
  use rom_diagnostic_audit, only: rom_diagnostic_audit_t
  use rom_validation, only: rom_validation_metrics_t, rom_validation_accept
  use num_types, only: rp
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use vector, only: vector_t
  use matrix, only: matrix_t
  use vector_math, only: vector_copy, vector_sub3, vector_absval, vector_glmax, &
       vector_glsc2, vector_glsubnorm
  use matrix_math, only: matrix_copy, matrix_glsc2, matrix_glsubnorm
  use device, only: DEVICE_TO_HOST, HOST_TO_DEVICE
  use neko_config, only: NEKO_BCKND_DEVICE
  use logger, only: neko_log
  use utils, only: neko_error
  use comm, only: pe_rank
  use mpi_f08, only: MPI_WTIME
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: rom_optimizer_t

  type, extends(optimizer_t) :: rom_optimizer_t
     type(mma_t) :: mma
     type(rom_problem_t) :: rom_problem
     type(vector_t) :: base_x
     logical :: unconstrained_problem = .false.
     integer :: inner_iterations = 1
     integer :: inner_counter = 0
     integer :: outer_max_iterations = 0
     real(kind=rp) :: tolerance = 1.0e-3_rp
     real(kind=rp) :: validation_tolerance = 5.0e-2_rp
     real(kind=rp) :: validation_gradient_tolerance = 5.0e-2_rp
     logical :: validate_full_order = .true.
     real(kind=rp) :: mma_scale = 1.0_rp
     logical :: mma_auto_scale = .false.
     logical :: diagnostic_audit = .false.
     real(kind=rp) :: mma_asyinit = 0.2_rp
     real(kind=rp) :: mma_asyincr = 1.2_rp
     real(kind=rp) :: mma_asydecr = 0.7_rp
     real(kind=rp) :: mma_move_limit = 0.2_rp
     real(kind=rp) :: mma_xmin = 0.0_rp
     real(kind=rp) :: mma_xmax = 1.0_rp
     type(rom_diagnostic_audit_t) :: audit
     logical :: enable_output = .true.
     logical :: supplied_basis_for_step = .false.
     logical :: validation_log_initialized = .false.
     integer :: rom_log_unit = -1
     logical :: rom_log_initialized = .false.
     character(len=16) :: last_validation_status = 'initial'
     real(kind=rp) :: rom_forward_seconds = 0.0_rp
     real(kind=rp) :: rom_adjoint_seconds = 0.0_rp
     real(kind=rp) :: mma_seconds = 0.0_rp
     real(kind=rp) :: full_order_seconds = 0.0_rp
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
    class(constraint_t), allocatable :: dummy_con
    integer :: max_iterations

    call this%free()
    call json_get(parameters, 'optimization.solver', solver_parameters)
    call this%read_base_settings(solver_parameters)
    call json_get_or_default(solver_parameters, 'max_iterations', &
         max_iterations, 20)
    this%outer_max_iterations = max_iterations
    call json_get_or_default(solver_parameters, 'tolerance', this%tolerance, &
         1.0e-3_rp)
    call json_get_or_default(solver_parameters, 'enable_output', &
         this%enable_output, .true.)

    ! MMA requires at least one constraint row.  Match mma_optimizer's
    ! unconstrained-problem behavior by adding a zero-gradient dummy row to
    ! the full-order problem when the user supplied no constraints.  The ROM
    ! problem still parses the original JSON and does not fabricate a volume
    ! constraint.
    this%unconstrained_problem = problem%get_n_constraints() == 0
    if (this%unconstrained_problem) then
       call neko_log%message('Unconstrained ROM problem detected. Adding a dummy MMA constraint.')
       allocate(dummy_constraint_t::dummy_con)
       select type (con => dummy_con)
       type is (dummy_constraint_t)
          call con%init_from_attributes(design)
       end select
       call problem%add_constraint(dummy_con)
       if (allocated(dummy_con)) deallocate(dummy_con)
    end if

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
    call json_get_or_default(rom_parameters, 'diagnostic_audit', &
         this%diagnostic_audit, .false.)
    call json_get_or_default(solver_parameters, 'mma.asyinit', &
         this%mma_asyinit, 0.2_rp)
    call json_get_or_default(solver_parameters, 'mma.asyincr', &
         this%mma_asyincr, 1.2_rp)
    call json_get_or_default(solver_parameters, 'mma.asydecr', &
         this%mma_asydecr, 0.7_rp)
    call json_get_or_default(solver_parameters, 'mma.move_limit', &
         this%mma_move_limit, 0.2_rp)
    call json_get_or_default(solver_parameters, 'mma.xmin', this%mma_xmin, 0.0_rp)
    call json_get_or_default(solver_parameters, 'mma.xmax', this%mma_xmax, 1.0_rp)
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
    this%unconstrained_problem = .false.
    call this%audit%free()
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
    call this%rom_problem%compute_sensitivity(density)
    if (present(chi_gradient)) then
       if (size(chi_gradient) /= n) call neko_error( &
            'ROM diagnostic chi-gradient size mismatch')
       if (NEKO_BCKND_DEVICE .eq. 1) call density%copy_from(DEVICE_TO_HOST, sync=.true.)
       chi_gradient = density%x
    end if
    call this%rom_problem%project_density(design, density, projected_gradient)
    if (NEKO_BCKND_DEVICE .eq. 1) call projected_gradient%copy_from( &
         DEVICE_TO_HOST, sync=.true.)
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
    type(vector_t) :: predicted_gradient
    type(vector_t) :: exact_density, projected_density
    type(vector_t) :: exact_gradient, projected_gradient
    type(vector_t) :: exact_density_before, projected_density_before
    type(vector_t) :: exact_gradient_before, projected_gradient_before
    type(vector_t) :: exact_gradient_reference, projected_gradient_reference
    type(vector_t) :: selected_gradient_before, selected_density_before
    type(vector_t) :: first_update
    type(vector_t) :: fom_initial_raw, fom_initial_gradient
    type(vector_t) :: fom_candidate_raw, fom_candidate_gradient
    type(vector_t) :: fom_initial_constraints, fom_candidate_constraints
    type(matrix_t) :: fom_initial_constraint_gradient
    type(matrix_t) :: fom_candidate_constraint_gradient
    type(matrix_t) :: fom_initial_constraint_raw
    type(matrix_t) :: fom_candidate_constraint_raw
    type(matrix_t) :: dfdx, predicted_constraint_gradient
    real(kind=rp), allocatable :: rom_state(:)
    real(kind=rp) :: rom_objective, actual_objective
    real(kind=rp), allocatable :: predicted_constraints(:), actual_constraints(:)
    real(kind=rp) :: global_max, global_norm
    real(kind=rp) :: objective_gradient_error, constraint_gradient_error
    real(kind=rp) :: exact_gradient_error, projected_gradient_error
    real(kind=rp) :: reference_objective, rom_reference_objective
    real(kind=rp) :: predicted_improvement, actual_improvement
    real(kind=rp) :: rom_scalar_before, rom_scalar_after
    real(kind=rp) :: rom_viscous_before, rom_viscous_after
    real(kind=rp) :: rom_brinkman_before, rom_brinkman_after
    real(kind=rp) :: rom_objective_before
    real(kind=rp) :: fom_initial_scalar, fom_initial_viscous, fom_initial_brinkman
    real(kind=rp) :: fom_candidate_scalar, fom_candidate_viscous, fom_candidate_brinkman
    real(kind=rp) :: selected_gradient_device_error
    real(kind=rp) :: selected_gradient_device_error_before
    real(kind=rp) :: initial_full_objective
    real(kind=rp), allocatable :: initial_constraints(:)
    real(kind=rp), allocatable :: constraint_before(:)
    real(kind=rp), allocatable :: reference_log_values(:)
    integer :: reference_log_size
    logical :: reference_available
    type(rom_validation_metrics_t) :: validation_metrics
    logical :: converged, validation_checked, validation_passed
    integer :: inner, n, m
    real(kind=rp) :: timer_start
    logical :: audit_active
    character(len=512) :: audit_design_file

    if (.not. present(simulation)) call neko_error('rom_mma step requires simulation')
    reference_objective = 0.0_rp
    reference_available = .false.
    reference_log_size = problem%get_log_size(.true.)
    if (reference_log_size > 0) then
       allocate(reference_log_values(reference_log_size))
       call problem%get_log_values(reference_log_values, .true.)
       reference_objective = reference_log_values(1)
       reference_available = .true.
    end if
    this%rom_forward_seconds = 0.0_rp
    this%rom_adjoint_seconds = 0.0_rp
    this%mma_seconds = 0.0_rp
    this%full_order_seconds = 0.0_rp
    call simulation%reset()
    if (.not. this%supplied_basis_for_step) then
       call this%rom_problem%refresh(simulation)
    end if
    if (this%enable_output) call this%rom_problem%initialize_outputs(design, simulation, iter)

    ! The streamed/POD basis is not available during optimizer initialization.
    ! Create the audit only after refresh (or an externally supplied basis), so
    ! its metadata and vector partition describe the actual final basis.
    if (this%diagnostic_audit .and. .not. this%audit%initialized) then
       call this%audit%init('rom_diagnostic_audit', this%rom_problem, this%mma, &
            this%outer_max_iterations, this%inner_iterations, this%tolerance, &
            this%validation_tolerance, this%validation_gradient_tolerance, &
            this%mma_scale, this%mma_auto_scale, this%mma_asyinit, &
            this%mma_asyincr, this%mma_asydecr, this%mma_move_limit, &
            this%mma_xmin, this%mma_xmax)
    end if

    n = design%size(); m = problem%get_n_constraints()
    call x%init(n); call xold%init(n); call old_design%init(n)
    call fval%init(m); call df0dx%init(n); call dx%init(n); call density%init(n)
    call dfdx%init(m, n)
    call predicted_gradient%init(n)
    call predicted_constraint_gradient%init(m, n)
    call first_update%init(n)
    first_update = 0.0_rp
    call design%get_values(x)
    call vector_copy(old_design, x)

    audit_active = this%diagnostic_audit .and. this%rom_problem%forward_constrained_enriched
    if (audit_active) then
       call exact_density%init(n)
       call projected_density%init(n)
       call exact_gradient%init(n)
       call projected_gradient%init(n)
       call exact_density_before%init(n)
       call projected_density_before%init(n)
       call exact_gradient_before%init(n)
       call projected_gradient_before%init(n)
       call exact_gradient_reference%init(n)
       call projected_gradient_reference%init(n)
       call selected_gradient_before%init(n)
       call selected_density_before%init(n)
       call fom_initial_raw%init(n)
       call fom_initial_gradient%init(n)
       call fom_candidate_raw%init(n)
       call fom_candidate_gradient%init(n)
       call fom_initial_constraints%init(m)
       call fom_candidate_constraints%init(m)
       call fom_initial_constraint_gradient%init(m, n)
       call fom_candidate_constraint_gradient%init(m, n)
       call fom_initial_constraint_raw%init(m, n)
       call fom_candidate_constraint_raw%init(m, n)
       allocate(initial_constraints(m), constraint_before(m))

       ! rom_problem_initialize has just completed the matching full-order
       ! reference lifecycle.  Capture its objective, constraints and mapped
       ! gradient before the ROM changes any design state.
       call problem%get_objective_value(initial_full_objective)
       call problem%get_constraint_values(fom_initial_constraints)
       call fom_initial_constraints%copy_from(DEVICE_TO_HOST, sync=.true.)
       initial_constraints = fom_initial_constraints%x
       call problem%get_constraint_sensitivities(fom_initial_constraint_raw)
       call matrix_copy(fom_initial_constraint_gradient, fom_initial_constraint_raw)
       select type (des_initial => design)
       type is (brinkman_design_t)
          call des_initial%get_sensitivity(fom_initial_raw)
          call vector_copy(fom_initial_gradient, fom_initial_raw)
          call des_initial%project_sensitivity(fom_initial_gradient)
          call des_initial%project_sensitivity(fom_initial_constraint_gradient)
       class default
          call problem%get_objective_sensitivities(fom_initial_gradient)
          call vector_copy(fom_initial_raw, fom_initial_gradient)
       end select
       call fom_initial_gradient%copy_from(DEVICE_TO_HOST, sync=.true.)
       call fom_initial_raw%copy_from(DEVICE_TO_HOST, sync=.true.)
       call fom_initial_constraint_gradient%copy_from(DEVICE_TO_HOST, sync=.true.)
       call fom_initial_constraint_raw%copy_from(DEVICE_TO_HOST, sync=.true.)
       fom_initial_scalar = 0.0_rp
       fom_initial_viscous = 0.0_rp
       fom_initial_brinkman = 0.0_rp
       if (reference_available) then
          fom_initial_scalar = rom_optimizer_log_value(reference_log_values, &
               problem%get_log_header(.true.), 'scalar_mixing')
          fom_initial_viscous = rom_optimizer_log_value(reference_log_values, &
               problem%get_log_header(.true.), 'viscous_dissipation')
          fom_initial_brinkman = rom_optimizer_log_value(reference_log_values, &
               problem%get_log_header(.true.), 'brinkman_dissipation')
       end if
    end if

    ! Evaluate the current design once before the first MMA update.  The
    ! resulting ROM forward/adjoint and gradients are consumed by the first
    ! inner update and then replaced only after the design has changed.
    call this%rom_problem%compute(design, rom_objective, final_state=rom_state)
    rom_reference_objective = rom_objective
    call rom_optimizer_check_finite_scalar('initial ROM objective', rom_objective)
    if (allocated(rom_state)) call rom_optimizer_check_finite_vector( &
         'initial ROM state', rom_state)
    call this%rom_problem%compute_sensitivity(density)
    call rom_optimizer_check_finite_vector('initial ROM raw density', density%x)
    if (this%unconstrained_problem) then
       fval = -1.0_rp
       dfdx = 0.0_rp
    else
       call this%rom_problem%compute_volume_constraint(design, fval%x(1), dfdx)
    end if
    call this%rom_problem%project_density(design, density, df0dx)
    call rom_optimizer_check_finite_vector('initial ROM density', density%x)
    call rom_optimizer_check_finite_vector('initial ROM constraint', fval%x)
    call rom_optimizer_check_finite_vector('initial ROM objective sensitivity', df0dx%x)
    call rom_optimizer_check_finite_matrix('initial ROM constraint sensitivity', dfdx%x)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call fval%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if
    if (audit_active) then
       call this%rom_problem%compute_density(exact_density, 'exact_galerkin')
       call this%rom_problem%compute_density(projected_density, 'projected_adjoint')
       call this%rom_problem%project_density(design, exact_density, exact_gradient)
       call this%rom_problem%project_density(design, projected_density, projected_gradient)
       call vector_copy(exact_gradient_reference, exact_gradient)
       call vector_copy(projected_gradient_reference, projected_gradient)
       call rom_optimizer_check_finite_vector('initial exact Galerkin gradient', exact_gradient%x)
       call rom_optimizer_check_finite_vector('initial projected-adjoint gradient', projected_gradient%x)
       if (trim(this%rom_problem%sensitivity_source) == 'exact_galerkin') then
          selected_gradient_device_error = vector_glsubnorm(df0dx, exact_gradient) / &
               max(1.0_rp, sqrt(vector_glsc2(exact_gradient, exact_gradient)))
       else
          selected_gradient_device_error = vector_glsubnorm(df0dx, projected_gradient) / &
               max(1.0_rp, sqrt(vector_glsc2(projected_gradient, projected_gradient)))
       end if
       write(audit_design_file, '("vectors/outer_",I0,"_inner_0_design.bin")') iter
       call this%audit%write_evaluation(iter, 0, x, density, df0dx, exact_density, &
            exact_gradient, projected_density, projected_gradient, fval%x, dfdx)
       call this%audit%write_full_evaluation(iter, 0, x, fom_initial_raw, &
            fom_initial_gradient, initial_constraints, fom_initial_constraint_gradient)
       call this%audit%write_gradient_comparison('reference', iter, 0, exact_gradient, &
            projected_gradient, fom_initial_gradient)
       call this%audit%write_objective_comparison('reference', iter, 0, &
            trim(audit_design_file), &
            rom_objective, this%rom_problem%scalar_mixing_value, &
            this%rom_problem%viscous_value, this%rom_problem%brinkman_value, &
            initial_full_objective, fom_initial_scalar, fom_initial_viscous, &
            fom_initial_brinkman)
       call this%mma%KKT(x, df0dx, fval, dfdx)
       constraint_before = fval%x
       call this%audit%write_step(iter, 0, x, x, exact_gradient, projected_gradient, &
            df0dx, rom_objective, rom_objective, this%rom_problem%scalar_mixing_value, &
            this%rom_problem%scalar_mixing_value, this%rom_problem%viscous_value, &
            this%rom_problem%viscous_value, this%rom_problem%brinkman_value, &
            this%rom_problem%brinkman_value, constraint_before, constraint_before, &
            this%mma%get_residumax(), this%mma%get_residunorm(), &
            selected_gradient_device_error)
    end if
    ! Stamp the complete initial ROM state after its adjoint sensitivity has
    ! been formed.  This produces the inner-zero forward, design and
    ! sensitivity fields in the same output lifecycle as later inner steps.
    if (this%enable_output .and. allocated(rom_state)) then
       call this%rom_problem%write_outputs(design, simulation, 0, rom_state)
    end if

    do inner = 1, this%inner_iterations
       this%inner_counter = this%inner_counter + 1
       call vector_copy(xold, x)
       if (audit_active) then
          call vector_copy(selected_gradient_before, df0dx)
          call vector_copy(selected_density_before, density)
          call vector_copy(exact_density_before, exact_density)
          call vector_copy(projected_density_before, projected_density)
          call vector_copy(exact_gradient_before, exact_gradient)
          call vector_copy(projected_gradient_before, projected_gradient)
          rom_scalar_before = this%rom_problem%scalar_mixing_value
          rom_viscous_before = this%rom_problem%viscous_value
          rom_brinkman_before = this%rom_problem%brinkman_value
          rom_objective_before = rom_objective
          selected_gradient_device_error_before = selected_gradient_device_error
          constraint_before = fval%x
       end if
       ! mma%update applies its configured xmin/xmax and mma.move_limit when
       ! constructing the subproblem. Do not alter that MMA solution with a
       ! second ROM-specific clip relative to the outer design.
       timer_start = real(MPI_WTIME(), rp)
       call this%mma%update(this%inner_counter, x, df0dx, fval, dfdx)
       this%mma_seconds = real(MPI_WTIME(), rp) - timer_start
       call rom_optimizer_check_finite_vector('MMA design candidate', x%x)
       call design%update_design(x)

       ! Evaluate the updated design exactly once.  This evaluation supplies
       ! the KKT check, output fields, and the cached state/gradient used by
       ! the next inner update; the next loop entrance does not repeat it.
       timer_start = real(MPI_WTIME(), rp)
       call this%rom_problem%compute(design, rom_objective, final_state=rom_state)
       this%rom_forward_seconds = real(MPI_WTIME(), rp) - timer_start
       call rom_optimizer_check_finite_scalar('ROM objective', rom_objective)
       if (allocated(rom_state)) call rom_optimizer_check_finite_vector( &
            'ROM state', rom_state)
       timer_start = real(MPI_WTIME(), rp)
       call this%rom_problem%compute_sensitivity(density)
       this%rom_adjoint_seconds = real(MPI_WTIME(), rp) - timer_start
       call rom_optimizer_check_finite_vector('ROM raw density', density%x)
       if (this%unconstrained_problem) then
          fval = -1.0_rp
          dfdx = 0.0_rp
       else
          call this%rom_problem%compute_volume_constraint(design, fval%x(1), dfdx)
       end if
       call this%rom_problem%project_density(design, density, df0dx)
       call rom_optimizer_check_finite_vector('ROM density', density%x)
       call rom_optimizer_check_finite_vector('ROM constraint', fval%x)
       call rom_optimizer_check_finite_vector('ROM objective sensitivity', df0dx%x)
       call rom_optimizer_check_finite_matrix('ROM constraint sensitivity', dfdx%x)
       if (NEKO_BCKND_DEVICE .eq. 1) call fval%copy_from(HOST_TO_DEVICE, sync = .true.)
       call this%mma%KKT(x, df0dx, fval, dfdx)

       call vector_sub3(dx, x, xold)
       if (audit_active .and. inner == 1) call vector_copy(first_update, dx)
       call vector_absval(dx)
       global_max = vector_glmax(dx)
       global_norm = sqrt(vector_glsc2(dx, dx))
       this%max_design_change = global_max
       this%norm2_design_change = global_norm
       if (NEKO_BCKND_DEVICE .eq. 1) call fval%copy_from(DEVICE_TO_HOST, sync = .true.)
       if (audit_active) then
          call this%rom_problem%compute_density(exact_density, 'exact_galerkin')
          call this%rom_problem%compute_density(projected_density, 'projected_adjoint')
          call this%rom_problem%project_density(design, exact_density, exact_gradient)
          call this%rom_problem%project_density(design, projected_density, projected_gradient)
          if (trim(this%rom_problem%sensitivity_source) == 'exact_galerkin') then
             selected_gradient_device_error = vector_glsubnorm(df0dx, exact_gradient) / &
                  max(1.0_rp, sqrt(vector_glsc2(exact_gradient, exact_gradient)))
          else
             selected_gradient_device_error = vector_glsubnorm(df0dx, projected_gradient) / &
                  max(1.0_rp, sqrt(vector_glsc2(projected_gradient, projected_gradient)))
          end if
          call this%audit%write_evaluation(iter, inner, x, density, df0dx, exact_density, &
               exact_gradient, projected_density, projected_gradient, fval%x, dfdx)
          call this%audit%write_step(iter, inner, xold, x, exact_gradient_before, &
               projected_gradient_before, selected_gradient_before, &
               rom_objective_before, rom_objective, rom_scalar_before, &
               this%rom_problem%scalar_mixing_value, rom_viscous_before, &
               this%rom_problem%viscous_value, rom_brinkman_before, &
               this%rom_problem%brinkman_value, constraint_before, fval%x, &
               this%mma%get_residumax(), this%mma%get_residunorm(), &
               selected_gradient_device_error_before)
          selected_gradient_device_error_before = selected_gradient_device_error
       end if
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
    call vector_copy(predicted_gradient, df0dx)
    call matrix_copy(predicted_constraint_gradient, dfdx)
    exact_gradient_error = -1.0_rp
    projected_gradient_error = -1.0_rp
    if (this%rom_problem%forward_constrained_enriched) then
       if (.not. audit_active) then
          call exact_density%init(n)
          call projected_density%init(n)
          call exact_gradient%init(n)
          call projected_gradient%init(n)
       end if
       call this%rom_problem%compute_density(exact_density, 'exact_galerkin')
       call this%rom_problem%compute_density(projected_density, 'projected_adjoint')
       call rom_optimizer_check_finite_vector('exact Galerkin raw density', &
            exact_density%x)
       call rom_optimizer_check_finite_vector('projected-adjoint raw density', &
            projected_density%x)
       call this%rom_problem%project_density(design, exact_density, exact_gradient)
       call this%rom_problem%project_density(design, projected_density, projected_gradient)
       call rom_optimizer_check_finite_vector('exact Galerkin gradient', &
            exact_gradient%x)
       call rom_optimizer_check_finite_vector('projected-adjoint gradient', &
            projected_gradient%x)
    end if
    call design%get_values(x)
    call vector_sub3(dx, x, old_design)
    call vector_absval(dx)
    global_max = vector_glmax(dx)
    global_norm = sqrt(vector_glsc2(dx, dx))
    this%max_design_change = global_max
    this%norm2_design_change = global_norm

    timer_start = real(MPI_WTIME(), rp)
    call simulation%set_design_iteration(iter)
    call problem%compute(design, simulation)
    call problem%get_objective_value(actual_objective)
    allocate(actual_constraints(m))
    call problem%get_constraint_values(fval)
    call fval%copy_from(DEVICE_TO_HOST, sync = .true.)
    actual_constraints = fval%x
    ! Compute the full-order sensitivity before deciding whether the ROM
    ! candidate is trustworthy.  Use the same constant-average SEM weighting
    ! as the production mma_optimizer; otherwise the validation would compare
    ! the ROM's production convention against a different full-order gradient.
    if (this%enable_output) call simulation%write_forward(iter)
    call problem%compute_sensitivity(design, simulation)
    ! Keep the standard outer-loop physical outputs alongside the ROM fields.
    ! This writes the full-order state already computed for validation; it does
    ! not trigger another forward or adjoint lifecycle.
    if (this%enable_output) call simulation%write_adjoint(iter)
    select type (des => design)
    type is (brinkman_design_t)
       if (audit_active) then
          call des%get_sensitivity(fom_candidate_raw)
          call vector_copy(fom_candidate_gradient, fom_candidate_raw)
          call problem%get_constraint_sensitivities(fom_candidate_constraint_raw)
          call matrix_copy(fom_candidate_constraint_gradient, fom_candidate_constraint_raw)
       end if
       call des%get_sensitivity(df0dx)
       call des%project_sensitivity(df0dx)
       call problem%get_constraint_sensitivities(dfdx)
       call des%project_sensitivity(dfdx)
       if (audit_active) then
          call vector_copy(fom_candidate_gradient, df0dx)
          call matrix_copy(fom_candidate_constraint_gradient, dfdx)
       end if
    class default
       call problem%get_objective_sensitivities(df0dx)
       call problem%get_constraint_sensitivities(dfdx)
       if (audit_active) then
          call vector_copy(fom_candidate_raw, df0dx)
          call vector_copy(fom_candidate_gradient, df0dx)
          call matrix_copy(fom_candidate_constraint_raw, dfdx)
          call matrix_copy(fom_candidate_constraint_gradient, dfdx)
       end if
    end select
    this%full_order_seconds = real(MPI_WTIME(), rp) - timer_start
    if (this%rom_problem%forward_constrained_enriched) then
       exact_gradient_error = vector_glsubnorm(df0dx, exact_gradient) / &
            max(1.0_rp, sqrt(vector_glsc2(df0dx, df0dx)))
       projected_gradient_error = vector_glsubnorm(df0dx, projected_gradient) / &
            max(1.0_rp, sqrt(vector_glsc2(df0dx, df0dx)))
    end if
    if (audit_active) then
       call fom_candidate_raw%copy_from(DEVICE_TO_HOST, sync=.true.)
       call fom_candidate_gradient%copy_from(DEVICE_TO_HOST, sync=.true.)
       call fom_candidate_constraint_raw%copy_from(DEVICE_TO_HOST, sync=.true.)
       call fom_candidate_constraint_gradient%copy_from(DEVICE_TO_HOST, sync=.true.)
       call problem%get_log_values(reference_log_values, .true.)
       fom_candidate_scalar = rom_optimizer_log_value(reference_log_values, &
            problem%get_log_header(.true.), 'scalar_mixing')
       fom_candidate_viscous = rom_optimizer_log_value(reference_log_values, &
            problem%get_log_header(.true.), 'viscous_dissipation')
       fom_candidate_brinkman = rom_optimizer_log_value(reference_log_values, &
            problem%get_log_header(.true.), 'brinkman_dissipation')
       call fom_candidate_constraints%copy_from(DEVICE_TO_HOST, sync=.true.)
       fom_candidate_constraints%x = actual_constraints
       call this%audit%write_full_evaluation(iter, this%inner_iterations, x, &
            fom_candidate_raw, fom_candidate_gradient, actual_constraints, &
            fom_candidate_constraint_gradient)
       call this%audit%write_gradient_comparison('candidate', iter, &
            this%inner_iterations, exact_gradient, projected_gradient, &
            fom_candidate_gradient)
       write(audit_design_file, '("vectors/outer_",I0,"_inner_",I0,"_design.bin")') &
            iter, this%inner_iterations
       call this%audit%write_objective_comparison('candidate', iter, &
            this%inner_iterations, trim(audit_design_file), rom_objective, &
            this%rom_problem%scalar_mixing_value, this%rom_problem%viscous_value, &
            this%rom_problem%brinkman_value, actual_objective, fom_candidate_scalar, &
            fom_candidate_viscous, fom_candidate_brinkman)
       call this%audit%finite_difference(iter, this%rom_problem, design, old_design, &
            x, exact_gradient_reference, projected_gradient_reference, first_update)
    end if
    objective_gradient_error = vector_glsubnorm(df0dx, predicted_gradient) / &
         max(1.0_rp, sqrt(vector_glsc2(df0dx, df0dx)))
    constraint_gradient_error = matrix_glsubnorm(dfdx, predicted_constraint_gradient) / &
         max(1.0_rp, sqrt(matrix_glsc2(dfdx, dfdx)))
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
    predicted_improvement = rom_reference_objective - rom_objective
    actual_improvement = 0.0_rp
    if (reference_available) actual_improvement = reference_objective - actual_objective
    call write_gradient_diagnostic(this, iter, rom_reference_objective, &
         reference_objective, rom_objective, actual_objective, predicted_improvement, &
         actual_improvement, exact_gradient_error, projected_gradient_error, &
         objective_gradient_error, reference_available, validation_checked, &
         validation_passed)
    ! Full-order validation is reported independently of MMA convergence.  A
    ! failed validation must not restore old_design or prevent the updated
    ! design from being written and used by the next outer iteration.
    converged = this%mma%get_residumax() < this%tolerance

    call x%free(); call xold%free(); call old_design%free(); call fval%free()
    call df0dx%free(); call dx%free(); call density%free(); call dfdx%free()
    call predicted_gradient%free(); call predicted_constraint_gradient%free()
    if (this%rom_problem%forward_constrained_enriched) then
       call exact_density%free(); call projected_density%free()
       call exact_gradient%free(); call projected_gradient%free()
       if (audit_active) then
          call exact_gradient_reference%free(); call projected_gradient_reference%free()
          call exact_density_before%free(); call projected_density_before%free()
          call exact_gradient_before%free(); call projected_gradient_before%free()
          call selected_gradient_before%free(); call selected_density_before%free()
          call fom_initial_raw%free(); call fom_initial_gradient%free()
          call fom_candidate_raw%free(); call fom_candidate_gradient%free()
          call fom_initial_constraints%free(); call fom_candidate_constraints%free()
          call fom_initial_constraint_gradient%free()
          call fom_candidate_constraint_gradient%free()
          call fom_initial_constraint_raw%free()
          call fom_candidate_constraint_raw%free()
       end if
    end if
    if (allocated(rom_state)) deallocate(rom_state)
    if (allocated(predicted_constraints)) deallocate(predicted_constraints)
    if (allocated(actual_constraints)) deallocate(actual_constraints)
    if (allocated(reference_log_values)) deallocate(reference_log_values)
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
         'scalar_mixing,scalar_mixing.weight,'// &
         'viscous_dissipation,viscous_dissipation.weight,'// &
         'brinkman_dissipation,brinkman_dissipation.weight,'// &
         'volume_constraint,volume_constraint.volume,max_design_change,'// &
         'norm2_design_change,rom_forward_seconds,rom_adjoint_seconds,'// &
         'mma_seconds,full_order_seconds,status'
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
    real(kind=rp) :: total, scalar_mixing, scalar_mixing_weight
    real(kind=rp) :: viscous, viscous_weight, brinkman, brinkman_weight
    real(kind=rp) :: constraint, volume
    integer :: n
    character(len=:), allocatable :: header

    if (.not. this%rom_log_initialized) return
    n = problem%get_log_size(.true.)
    allocate(values(n))
    call problem%get_log_values(values, .true.)

    total = values(1)
    header = problem%get_log_header(.true.)
    scalar_mixing = rom_optimizer_log_value(values, header, 'scalar_mixing')
    scalar_mixing_weight = rom_optimizer_log_value(values, header, 'scalar_mixing.weight')
    viscous = rom_optimizer_log_value(values, header, 'viscous_dissipation')
    viscous_weight = rom_optimizer_log_value(values, header, 'viscous_dissipation.weight')
    brinkman = rom_optimizer_log_value(values, header, 'brinkman_dissipation')
    brinkman_weight = rom_optimizer_log_value(values, header, 'brinkman_dissipation.weight')
    constraint = rom_optimizer_log_value(values, header, 'volume_constraint')
    volume = rom_optimizer_log_value(values, header, 'volume_constraint.volume')
    if (this%unconstrained_problem) then
       constraint = -1.0_rp
       volume = 0.0_rp
    end if
    call rom_optimizer_write_record(this, 'outer', iter, 0, 'full_order', &
         total, scalar_mixing, scalar_mixing_weight, viscous, viscous_weight, &
         brinkman, brinkman_weight, constraint, volume, 0.0_rp, 0.0_rp, &
         0.0_rp, this%full_order_seconds, this%last_validation_status)
    deallocate(values)
    deallocate(header)
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
    if (this%unconstrained_problem) then
       volume = 0.0_rp
    else if (this%rom_problem%volume_is_max) then
       volume = this%rom_problem%volume_limit + constraint
    else
       volume = this%rom_problem%volume_limit - constraint
    end if
    call rom_optimizer_write_record(this, 'inner', outer, inner, 'rom', &
         objective, this%rom_problem%scalar_mixing_value, &
         this%rom_problem%scalar_mixing_weight, this%rom_problem%viscous_value, &
         this%rom_problem%viscous_weight, this%rom_problem%brinkman_value, &
         this%rom_problem%brinkman_weight, constraint, volume, &
         this%rom_forward_seconds, this%rom_adjoint_seconds, this%mma_seconds, &
         0.0_rp, 'not_validated')
  end subroutine rom_optimizer_write_inner

  !> Append one ordered record to the ROM optimization history.
  !! @param this ROM optimizer.
  !! @param phase `outer` for a full-order refresh or `inner` for ROM MMA.
  !! @param outer Outer optimizer iteration.
  !! @param inner Inner reduced iteration, or zero for an outer record.
  !! @param source Objective source (`full_order` or `rom`).
  !! @param total Weighted total objective.
  !! @param scalar_mixing Unweighted scalar-mixing objective value.
  !! @param scalar_mixing_weight Scalar-mixing objective weight.
  !! @param viscous Unweighted viscous objective value.
  !! @param viscous_weight Viscous objective weight.
  !! @param brinkman Unweighted Brinkman objective value.
  !! @param brinkman_weight Brinkman objective weight.
  !! @param constraint Volume constraint value.
  !! @param volume Normalised design volume.
  !! @param status Validation status for the outer record.
  subroutine rom_optimizer_write_record(this, phase, outer, inner, source, &
       total, scalar_mixing, scalar_mixing_weight, viscous, viscous_weight, &
       brinkman, brinkman_weight, constraint, volume, rom_forward_seconds, &
       rom_adjoint_seconds, mma_seconds, full_order_seconds, status)
    class(rom_optimizer_t), intent(inout) :: this
    character(len=*), intent(in) :: phase, source, status
    integer, intent(in) :: outer, inner
    real(kind=rp), intent(in) :: total, scalar_mixing, scalar_mixing_weight
    real(kind=rp), intent(in) :: viscous, viscous_weight, brinkman, brinkman_weight
    real(kind=rp), intent(in) :: constraint, volume, rom_forward_seconds
    real(kind=rp), intent(in) :: rom_adjoint_seconds, mma_seconds, full_order_seconds
    integer :: ios

    if (.not. this%rom_log_initialized) return
    write(this%rom_log_unit, &
         '(A,",",I0,",",I0,",",A,",",ES24.16,",",ES24.16,",",ES24.16,",",'// &
         'ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",'// &
         'ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",'// &
         'ES24.16,",",ES24.16,",",A)', iostat=ios) phase, outer, inner, &
         source, total, scalar_mixing, scalar_mixing_weight, viscous, viscous_weight, &
         brinkman, brinkman_weight, constraint, volume, this%max_design_change, &
         this%norm2_design_change, rom_forward_seconds, rom_adjoint_seconds, &
         mma_seconds, full_order_seconds, status
    if (ios /= 0) call neko_error('could not write ROM optimization_data.csv')
    ! Keep the mixed outer/inner history inspectable while the optimisation is
    ! running; otherwise Fortran may retain several inner records in its I/O
    ! buffer until the optimizer is freed.
    flush(this%rom_log_unit, iostat=ios)
    if (ios /= 0) call neko_error('could not flush ROM optimization_data.csv')
  end subroutine rom_optimizer_write_record

  !> Return one value from Neko's comma-separated objective/constraint log.
  !! The full-order problem owns the ordering; searching its header avoids
  !! assuming that a scalar objective is present or comes first.
  function rom_optimizer_log_value(values, header, label) result(value)
    real(kind=rp), intent(in) :: values(:)
    character(len=*), intent(in) :: header, label
    real(kind=rp) :: value
    character(len=len(header)) :: normal_header
    character(len=len(label)) :: normal_label
    integer :: position, column, i

    value = 0.0_rp
    normal_header = header
    normal_label = label
    do i = 1, len(normal_header)
       if (normal_header(i:i) >= 'A' .and. normal_header(i:i) <= 'Z') &
            normal_header(i:i) = achar(iachar(normal_header(i:i)) + 32)
       if (normal_header(i:i) == ' ') normal_header(i:i) = '_'
    end do
    do i = 1, len(normal_label)
       if (normal_label(i:i) >= 'A' .and. normal_label(i:i) <= 'Z') &
            normal_label(i:i) = achar(iachar(normal_label(i:i)) + 32)
       if (normal_label(i:i) == ' ') normal_label(i:i) = '_'
    end do
    position = index(normal_header, trim(normal_label))
    if (position <= 0) return
    column = 1
    do i = 1, position - 1
       if (header(i:i) == ',') column = column + 1
    end do
    if (column <= size(values)) value = values(column)
  end function rom_optimizer_log_value

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

  !> Append the selected and alternative enriched-adjoint gradient diagnostics.
  !!
  !! The first full-order value is the objective at the design entering the
  !! outer iteration.  The ROM candidate value is evaluated after the inner
  !! MMA loop, and the following full-order value is the same candidate after
  !! the normal validation lifecycle.  This keeps predicted and actual
  !! improvements on one explicitly labelled record without changing MMA.
  subroutine write_gradient_diagnostic(this, iter, rom_reference_objective, &
       reference_objective, rom_candidate_objective, actual_objective, &
       predicted_improvement, actual_improvement, exact_gradient_error, &
       projected_gradient_error, &
       selected_gradient_error, reference_available, validation_checked, &
       validation_passed)
    class(rom_optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    real(kind=rp), intent(in) :: rom_reference_objective, reference_objective
    real(kind=rp), intent(in) :: rom_candidate_objective, actual_objective
    real(kind=rp), intent(in) :: predicted_improvement
    real(kind=rp), intent(in) :: actual_improvement, exact_gradient_error
    real(kind=rp), intent(in) :: projected_gradient_error, selected_gradient_error
    logical, intent(in) :: reference_available, validation_checked
    logical, intent(in) :: validation_passed
    integer :: unit, ios
    logical :: exists

    if (pe_rank /= 0) return
    inquire(file='rom_gradient_diagnostics.csv', exist=exists)
    open(newunit=unit, file='rom_gradient_diagnostics.csv', status='unknown', &
         position='append', action='write', iostat=ios)
    if (ios /= 0) call neko_error('could not open rom_gradient_diagnostics.csv')
    if (.not. exists) write(unit, '(A)') &
         'outer_iteration,rom_reference_objective,full_order_reference_objective,'// &
         'rom_candidate_objective,full_order_candidate_objective,'// &
         'predicted_rom_improvement,actual_full_order_improvement,'// &
         'full_order_reference_available,exact_galerkin_gradient_error,'// &
         'projected_adjoint_gradient_error,selected_gradient_error,'// &
         'validation_checked,validation_passed,sensitivity_source'
    write(unit, '(I0,6(",",ES24.16),",",L1,3(",",ES24.16),",",L1,",",L1,",",A)') &
         iter, rom_reference_objective, reference_objective, rom_candidate_objective, &
         actual_objective, predicted_improvement, actual_improvement, &
         reference_available, exact_gradient_error, projected_gradient_error, &
         selected_gradient_error, validation_checked, validation_passed, &
         trim(this%rom_problem%sensitivity_source)
    close(unit)
  end subroutine write_gradient_diagnostic

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

  !> Stop the ROM optimization before MMA if a scalar is not finite.
  subroutine rom_optimizer_check_finite_scalar(label, value)
    character(len=*), intent(in) :: label
    real(kind=rp), intent(in) :: value
    if (.not. ieee_is_finite(value)) call neko_error( &
         'ROM failure: non-finite '//trim(label)//'; stopping before MMA')
  end subroutine rom_optimizer_check_finite_scalar

  !> Stop the ROM optimization before MMA if a vector contains a non-finite value.
  subroutine rom_optimizer_check_finite_vector(label, values)
    character(len=*), intent(in) :: label
    real(kind=rp), intent(in) :: values(:)
    if (.not. all(ieee_is_finite(values))) call neko_error( &
         'ROM failure: non-finite '//trim(label)//'; stopping before MMA')
  end subroutine rom_optimizer_check_finite_vector

  !> Stop the ROM optimization before MMA if a matrix contains a non-finite value.
  subroutine rom_optimizer_check_finite_matrix(label, values)
    character(len=*), intent(in) :: label
    real(kind=rp), intent(in) :: values(:, :)
    if (.not. all(ieee_is_finite(values))) call neko_error( &
         'ROM failure: non-finite '//trim(label)//'; stopping before MMA')
  end subroutine rom_optimizer_check_finite_matrix

end module rom_optimizer
