!> Reduced optimization problem for the Neko velocity equations.
!!
!! This type owns construction of the reduced operators, forward/adjoint
!! unsteady integration, objective evaluation and physical design sensitivity.
!! MMA and full-order trust decisions remain in rom_optimizer.
module rom_problem
  use problem, only: problem_t
  use design, only: design_t
  use brinkman_design, only: brinkman_design_t
  use simulation_m, only: simulation_t
  use field_output, only: field_output_t
  use field, only: field_t, field_wrapper_t
  use field_math, only: field_add2s2, field_copy
  use matrix, only: matrix_t
  use rom_model, only: rom_model_t
  use rom_lift, only: rom_lift_t
  use rom_assembly, only: rom_assemble_mass, rom_assemble_brinkman, &
       rom_assemble_lift_convection, rom_assemble_diffusion, &
       rom_assemble_convection
  use rom_objectives, only: rom_assemble_viscous_terms, &
       rom_assemble_brinkman_objective, rom_viscous_value, rom_brinkman_value
  use rom_basis_handoff, only: rom_receive_basis, rom_project_basis
  use operators, only: grad
  use num_types, only: rp, sp, dp
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default, json_get_or_lookup, &
       json_extract_item
  use vector, only: vector_t, vector_ptr_t
  use scratch_registry, only: neko_scratch_registry
  use vector_math, only: vector_glsc2, vector_glsc3, vector_col3, vector_add2, &
       vector_add2s2, vector_cmult, vector_vdot3
  use device, only: DEVICE_TO_HOST, HOST_TO_DEVICE
  use neko_config, only: NEKO_BCKND_DEVICE
  use point_zone, only: point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  use utils, only: neko_error, filename_path
  use comm, only: NEKO_COMM, MPI_REAL_PRECISION
  use mpi_f08, only: MPI_Allreduce, MPI_MAX
  implicit none
  private

  public :: rom_problem_t

  type :: rom_basis_storage_t
     integer :: n_rows = 0
     integer :: n_modes = 0
     integer :: n_fields = 0
     type(matrix_t) :: modes
  end type rom_basis_storage_t

  type, public :: rom_problem_t
     type(rom_model_t) :: rom
     type(rom_basis_storage_t) :: basis

     integer :: inner_iterations = 1
     integer :: rank_total_max = 1
     integer :: history_outer_iterations = 3
     real(kind=rp) :: scalar_metric_weight = 1.0_rp
     real(kind=rp) :: admissibility_tolerance = 1.0e-8_rp
     logical :: enable_output = .true.
     logical :: rom_ready = .false.
     logical :: write_trajectory_fields = .false.
     real(kind=rp) :: trajectory_output_interval = 0.5_rp

     integer :: n_steps = 0
     real(kind=rp) :: dt = 0.0_rp
     real(kind=rp) :: time_start = 0.0_rp
     real(kind=rp) :: time_end = 0.0_rp
     real(kind=rp) :: viscosity = 0.0_rp
     real(kind=rp) :: density = 1.0_rp

     logical :: has_viscous = .false.
     logical :: has_brinkman = .false.
     real(kind=rp) :: viscous_weight = 0.0_rp
     real(kind=rp) :: brinkman_weight = 0.0_rp
     real(kind=rp) :: viscous_start = 0.0_rp
     real(kind=rp) :: viscous_end = huge(0.0_rp)
     real(kind=rp) :: brinkman_start = 0.0_rp
     real(kind=rp) :: brinkman_end = huge(0.0_rp)
     character(len=80) :: viscous_mask_name = ''
     character(len=80) :: brinkman_mask_name = ''
     character(len=80) :: volume_mask_name = ''
     logical :: volume_is_max = .false.
     real(kind=rp) :: volume_limit = 0.0_rp
     real(kind=rp) :: volume_domain = 0.0_rp
     real(kind=rp) :: viscous_domain = 0.0_rp
     real(kind=rp) :: brinkman_domain = 0.0_rp

     type(vector_t) :: weights, viscous_mask, brinkman_mask, volume_mask
     type(vector_t), allocatable :: velocity_modes(:, :)
     type(vector_t), allocatable :: velocity_gradients(:, :, :)
     type(rom_lift_t) :: lift
     type(vector_t) :: lift_gradients(3,3)
     type(vector_t) :: lift_velocity(3)
     type(matrix_t) :: mass, fixed_linear, fixed_implicit
     type(vector_t) :: fixed_constant
     type(matrix_t), allocatable :: quadratic(:)
     type(vector_t) :: viscous_constant, viscous_linear
     type(matrix_t) :: viscous_quadratic
     type(vector_t) :: initial_state

     type(matrix_t) :: candidate_brinkman, candidate_objective_brinkman
     type(vector_t) :: candidate_brinkman_lift, candidate_objective_lift
     real(kind=rp) :: candidate_objective_constant = 0.0_rp
     real(kind=rp) :: candidate_volume = 0.0_rp
     type(matrix_t) :: states, adjoint, running_gradient, correlation, cost_sum_aa
     type(vector_t) :: terminal_gradient, lambda_sum, cost_sum_a
     real(kind=rp) :: cost_total = 0.0_rp
     real(kind=rp) :: objective_value = 0.0_rp
     real(kind=rp) :: viscous_value = 0.0_rp
     real(kind=rp) :: brinkman_value = 0.0_rp
     logical :: forward_ready = .false.
     logical :: objectives_ready = .false.
     logical :: backward_ready = .false.

     type(field_output_t) :: rom_forward_output
     type(field_output_t) :: rom_design_output
     type(field_output_t) :: rom_sensitivity_output
     type(field_output_t) :: rom_forward_trajectory_output
     type(field_output_t) :: rom_adjoint_trajectory_output
     type(field_wrapper_t) :: rom_velocity(3)
     type(field_wrapper_t) :: rom_adjoint_velocity(3)
     logical :: rom_outputs_initialized = .false.
   contains
     procedure, pass(this), public :: init_from_json => rom_problem_init_from_json
     procedure, pass(this), public :: free => rom_problem_free
     procedure, pass(this), public :: initialize => rom_problem_initialize
     procedure, pass(this), public :: initialize_from_basis => &
          rom_problem_initialize_from_basis
     procedure, pass(this), public :: compute => rom_problem_compute
     procedure, pass(this), public :: compute_sensitivity => &
          rom_problem_compute_sensitivity
     procedure, pass(this), public :: run_forward_unsteady => &
          rom_problem_run_forward_unsteady
     procedure, pass(this), public :: run_backward_unsteady => &
          rom_problem_run_backward_unsteady
     procedure, pass(this), public :: refresh => rom_problem_refresh
     procedure, pass(this), public :: project_density => rom_problem_project_density
     procedure, pass(this), public :: project_density_row => &
          rom_problem_project_density_row
     procedure, pass(this), public :: initialize_outputs => &
          rom_problem_initialize_outputs
     procedure, pass(this), public :: write_outputs => rom_problem_write_outputs
     procedure, pass(this), public :: compute_objectives => rom_problem_compute_objectives
     procedure, pass(this), public :: compute_density => rom_problem_compute_density
     procedure, pass(this), public :: compute_volume_constraint => &
          rom_problem_compute_volume_constraint
     procedure, pass(this), public :: verify_dirichlet_reconstruction => &
          rom_problem_verify_dirichlet_reconstruction
  end type rom_problem_t

contains

  !> Read ROM time-integration, rank, objective and admissibility settings.
  subroutine rom_problem_init_from_json(this, parameters, problem)
    class(rom_problem_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    class(problem_t), intent(in) :: problem
    type(json_file) :: rom_parameters
    logical :: dealias, variable_dt
    character(len=:), allocatable :: time_scheme, pressure_treatment
    integer :: time_order

    call this%free()
    if ('rom' .in. parameters) then
       call json_get(parameters, 'rom', rom_parameters)
    else
       call json_get(parameters, 'optimization.solver', rom_parameters)
    end if
    call json_get_or_default(rom_parameters, 'inner_iterations', &
         this%inner_iterations, 3)
    call json_get_or_default(rom_parameters, 'history_outer_iterations', &
         this%history_outer_iterations, 3)
    call json_get_or_default(rom_parameters, 'rank_total_max', &
         this%rank_total_max, 12)
    call json_get_or_default(rom_parameters, 'scalar_metric_weight', &
         this%scalar_metric_weight, 1.0_rp)
    call json_get_or_default(rom_parameters, 'admissibility_tolerance', &
         this%admissibility_tolerance, 1.0e-8_rp)
    call json_get_or_default(rom_parameters, 'enable_output', &
         this%enable_output, .true.)
    call json_get_or_default(rom_parameters, 'write_trajectory_fields', &
         this%write_trajectory_fields, .false.)
    call json_get_or_default(rom_parameters, 'trajectory_output_interval', &
         this%trajectory_output_interval, 0.5_rp)
    if (this%inner_iterations < 1) call neko_error( &
         'rom_mma requires inner_iterations >= 1')
    if (this%history_outer_iterations < 1) call neko_error( &
         'rom_mma requires history_outer_iterations >= 1')
    if (this%rank_total_max < 1) call neko_error( &
         'rom_mma requires rank_total_max >= 1')
    if (this%admissibility_tolerance <= 0.0_rp) call neko_error( &
         'rom_mma requires admissibility_tolerance > 0')
    if (this%write_trajectory_fields .and. &
         this%trajectory_output_interval <= 0.0_rp) call neko_error( &
         'rom_mma requires trajectory_output_interval > 0')
    call json_get_or_default(rom_parameters, 'time_scheme', time_scheme, &
         'bdf1')
    call json_get_or_default(rom_parameters, 'pressure_treatment', &
         pressure_treatment, 'solenoidal_galerkin')
    if (trim(time_scheme) /= 'bdf1') then
       call neko_error('rom_mma requires the explicitly documented time_scheme=bdf1')
    end if
    if (trim(pressure_treatment) /= 'solenoidal_galerkin') then
       call neko_error('rom_mma requires pressure_treatment=solenoidal_galerkin')
    end if
    call json_get(parameters, 'case.numerics.dealias', dealias)
    if (dealias) call neko_error('rom_mma requires case.numerics.dealias=false')
    call json_get_or_lookup(parameters, 'case.numerics.time_order', time_order)
    if (time_order /= 1) call neko_error('rom_mma currently requires case.numerics.time_order=1')
    call json_get_or_default(parameters, 'case.time.variable_timestep', &
         variable_dt, .false.)
    if (variable_dt) call neko_error('rom_mma does not support variable timestep yet')
    if (problem%get_n_constraints() /= 1) then
       call neko_error('rom_mma currently supports exactly one volume constraint')
    end if

    call read_rom_objectives(this, parameters)
    call read_rom_constraint(this, parameters)
  end subroutine rom_problem_init_from_json

  !> Initialise the ROM time grid and problem-dependent physical metadata.
  subroutine rom_problem_initialize(this, problem, design, simulation)
    class(rom_problem_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), intent(inout) :: simulation
    if (simulation%have_scalar) call neko_error('rom_mma scalar coupling is not implemented')
    ! Construct the lift before the full-order run starts.  State recovery has
    ! already used the same ROM-local construction for its initial stream.
    call this%lift%init(simulation%neko_case)
    if (.not. this%lift%enabled) call neko_error( &
         'ROM requires a spatial lift constructed from the case velocity BCs')
    call problem%compute(design, simulation)
    call problem%compute_sensitivity(design, simulation)
    this%n_steps = simulation%n_timesteps
    this%dt = real(simulation%neko_case%time%dt, rp)
    this%time_start = real(simulation%neko_case%time%start_time, rp)
    this%time_end = real(simulation%neko_case%time%end_time, rp)
    if (this%n_steps < 1 .or. this%dt <= 0.0_rp) call neko_error('invalid Neko ROM time grid')
    if (abs(this%time_start + this%dt * real(this%n_steps, rp) - this%time_end) &
         > 100.0_rp * epsilon(this%time_end) * max(1.0_rp, abs(this%time_end))) then
       call neko_error('rom_mma requires a fixed timestep ending exactly at case.time.end_time')
    end if
    call simulation%reset()
  end subroutine rom_problem_initialize

  !> Release ROM operators, trajectories, outputs and basis storage.
  subroutine rom_problem_free(this)
    class(rom_problem_t), intent(inout) :: this
    call this%rom%free()
    call this%basis%modes%free()
    call rom_problem_free_arrays(this)
    call rom_problem_free_outputs(this)
    call rom_problem_free_trajectories(this)
    call this%lift%free()
    this%rom_ready = .false.
  end subroutine rom_problem_free

  !> Release only forward and adjoint trajectory storage.
  subroutine rom_problem_free_trajectories(this)
    class(rom_problem_t), intent(inout) :: this
    call this%candidate_brinkman%free()
    call this%candidate_brinkman_lift%free()
    call this%candidate_objective_brinkman%free()
    call this%candidate_objective_lift%free()
    call this%states%free(); call this%adjoint%free()
    call this%running_gradient%free(); call this%terminal_gradient%free()
    call this%lambda_sum%free(); call this%correlation%free()
    call this%cost_sum_a%free(); call this%cost_sum_aa%free()
    this%forward_ready = .false.
    this%objectives_ready = .false.
    this%backward_ready = .false.
  end subroutine rom_problem_free_trajectories

  !> Load a common velocity basis and assemble all fixed reduced operators.
  subroutine rom_problem_initialize_from_basis(this, basis, n_modes, &
       n_fields, simulation)
    class(rom_problem_t), intent(inout) :: this
    real(kind=rp), intent(in) :: basis(:, :)
    integer, intent(in) :: n_modes, n_fields
    type(simulation_t), intent(inout) :: simulation
    integer :: expected_rows
    if (n_fields /= 3 .or. n_modes < 1 .or. size(basis, 2) /= n_modes) then
       call neko_error('ROM diagnostic basis must be a nonempty velocity basis')
    end if
    expected_rows = 3 * simulation%fluid%u%dof%size()
    if (size(basis, 1) /= expected_rows) then
       call neko_error('ROM diagnostic basis rows do not match the velocity fields')
    end if
    if (.not. this%lift%enabled) call this%lift%init(simulation%neko_case)
    if (.not. this%lift%enabled) call neko_error( &
         'ROM diagnostic basis requires a spatial velocity lift')
    call this%basis%modes%free()
    this%basis%n_rows = size(basis, 1)
    this%basis%n_modes = n_modes
    this%basis%n_fields = n_fields
    call this%basis%modes%init(size(basis,1), size(basis,2))
    this%basis%modes%x = basis
    call rom_problem_assemble_fixed(this, simulation)
    call rom_problem_setup_masks(this, simulation)
    call rom_problem_assemble_objective_operators(this)
  end subroutine rom_problem_initialize_from_basis

  !> Run the ROM forward/adjoint lifecycle and return objective data.
  !! The optional constraint output is retained for standalone ROM diagnostics.
  !! The production ROM-MMA path omits it and obtains the design-only volume
  !! constraint from the live full-order problem_t.
  subroutine rom_problem_compute(this, design, objective, constraints, final_state)
    class(rom_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design
    real(kind=rp), intent(out) :: objective
    real(kind=rp), intent(out), optional :: constraints(:)
    real(kind=rp), allocatable, intent(out) :: final_state(:)
    integer :: r
    if (present(constraints)) then
       call this%run_forward_unsteady(design, .true.)
    else
       call this%run_forward_unsteady(design, .false.)
    end if
    if (present(constraints)) then
       call this%compute_objectives(objective, constraints)
    else
       call this%compute_objectives(objective)
    end if
    call this%run_backward_unsteady()
    r = this%states%get_nrows()
    allocate(final_state(r))
    final_state = this%states%x(:, this%n_steps + 1)
  end subroutine rom_problem_compute

  !> Form the full-dimensional ROM objective design sensitivity.
  !! The optional constraint density is retained only for standalone reduced
  !! diagnostics; production ROM-MMA obtains it from problem_t.
  subroutine rom_problem_compute_sensitivity(this, density, constraint_density)
    class(rom_problem_t), intent(inout) :: this
    real(kind=rp), intent(out) :: density(:)
    real(kind=rp), intent(out), optional :: constraint_density(:, :)
    type(vector_t), pointer :: density_vector
    integer :: density_index
    if (.not. this%backward_ready) call neko_error( &
         'ROM sensitivity requested before run_backward_unsteady')
    call neko_scratch_registry%request_vector(density_vector, density_index, &
         this%weights%size(), .true.)
    call this%compute_density(density_vector)
    density = density_vector%x
    if (present(constraint_density)) then
       if (this%volume_is_max) then
          constraint_density(1,:) = this%volume_mask%x / this%volume_domain
       else
          constraint_density(1,:) = -this%volume_mask%x / this%volume_domain
       end if
    end if
    call neko_scratch_registry%relinquish_vector(density_index)
    nullify(density_vector)
  end subroutine rom_problem_compute_sensitivity

  !> Evaluate the design-only volume constraint and its MMA derivative.
  !!
  !! The full-order volume constraint is defined on the unfiltered design
  !! variable returned by `design%get_values`.  The mapped Brinkman amplitude
  !! is used by the flow operators, but it is not the volume constraint unless
  !! a separate constraint mapping is explicitly configured.  The ROM
  !! implementation currently supports the same unmapped volume constraint as
  !! the production ROM case.  The returned derivative is already multiplied
  !! by the SEM diagonal, matching the directional-derivative convention used
  !! by MMA.
  !! @param[in] this The ROM problem.
  !! @param[in] design The current unfiltered design.
  !! @param[out] value The volume constraint value.
  !! @param[inout] sensitivity The one-row MMA constraint derivative.
  subroutine rom_problem_compute_volume_constraint(this, design, value, &
       sensitivity)
    class(rom_problem_t), intent(in) :: this
    class(design_t), intent(in) :: design
    real(kind=rp), intent(out) :: value
    type(matrix_t), intent(inout) :: sensitivity
    type(vector_t), pointer :: raw_design, derivative
    integer :: raw_index, derivative_index, n
    real(kind=rp) :: volume

    n = this%weights%size()
    if (sensitivity%get_nrows() /= 1 .or. sensitivity%get_ncols() /= n) then
       call neko_error('ROM volume constraint sensitivity has wrong shape')
    end if

    select type (design)
    type is (brinkman_design_t)
       call neko_scratch_registry%request_vector(raw_design, raw_index, n, .true.)
       call design%get_values(raw_design)
    class default
       call neko_error('ROM volume constraint requires a Brinkman design')
    end select

    ! vector_glsc3 is already a global Neko contraction.
    volume = vector_glsc3(this%weights, this%volume_mask, raw_design)
    if (this%volume_is_max) then
       value = volume / this%volume_domain - this%volume_limit
    else
       value = this%volume_limit - volume / this%volume_domain
    end if

    call neko_scratch_registry%request_vector(derivative, derivative_index, &
         n, .true.)
    call vector_col3(derivative, this%weights, this%volume_mask)
    if (this%volume_is_max) then
       call vector_cmult(derivative, 1.0_rp / this%volume_domain)
    else
       call vector_cmult(derivative, -1.0_rp / this%volume_domain)
    end if
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call derivative%copy_from(DEVICE_TO_HOST, sync = .true.)
    end if
    sensitivity%x(1,:) = derivative%x
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call sensitivity%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if

    call neko_scratch_registry%relinquish_vector([raw_index, derivative_index])
    nullify(raw_design); nullify(derivative)
  end subroutine rom_problem_compute_volume_constraint

  !> Integrate the reduced forward equations through the configured time grid.
  subroutine rom_problem_run_forward_unsteady(this, design, compute_constraint)
    class(rom_problem_t), intent(inout) :: this
    type(vector_t), pointer :: candidate_constant
    type(matrix_t), pointer :: candidate_linear
    integer :: candidate_constant_index, candidate_linear_index
    class(design_t), intent(inout) :: design
    logical, intent(in), optional :: compute_constraint
    integer :: r
    logical :: need_constraint
    need_constraint = .false.
    if (present(compute_constraint)) need_constraint = compute_constraint
    call rom_problem_free_trajectories(this)
    r = this%mass%get_nrows()
    call this%candidate_brinkman%init(r,r)
    call this%candidate_brinkman_lift%init(r)
    call this%candidate_objective_brinkman%init(r,r)
    call this%candidate_objective_lift%init(r)
    if (need_constraint) then
       call rom_problem_assemble_candidate(this, design, this%candidate_brinkman, &
            this%candidate_brinkman_lift, this%candidate_objective_brinkman, &
            this%candidate_objective_lift, this%candidate_objective_constant, &
            this%candidate_volume)
    else
       call rom_problem_assemble_candidate(this, design, this%candidate_brinkman, &
            this%candidate_brinkman_lift, this%candidate_objective_brinkman, &
            this%candidate_objective_lift, this%candidate_objective_constant)
    end if
    call neko_scratch_registry%request_vector(candidate_constant, candidate_constant_index, &
         r, .true.)
    call neko_scratch_registry%request_matrix(candidate_linear, candidate_linear_index, &
         r, r, .true.)
    candidate_constant%x = this%fixed_constant%x - this%candidate_brinkman_lift%x
    candidate_linear%x = this%fixed_linear%x - this%candidate_brinkman%x
    call this%rom%init(this%mass, candidate_constant, candidate_linear, &
         this%quadratic, this%fixed_implicit)
    call this%states%init(r, this%n_steps + 1)
    call this%running_gradient%init(r, this%n_steps)
    call this%terminal_gradient%init(r)
    call this%cost_sum_a%init(r)
    call this%cost_sum_aa%init(r,r)
    call this%rom%forward_bdf1(this%initial_state, this%dt, this%states)
    call neko_scratch_registry%relinquish_vector(candidate_constant_index)
    call neko_scratch_registry%relinquish_matrix(candidate_linear_index)
    nullify(candidate_constant); nullify(candidate_linear)
    this%forward_ready = .true.
    if (this%write_trajectory_fields .and. this%rom_outputs_initialized) &
         call rom_problem_write_forward_trajectory(this)
  end subroutine rom_problem_run_forward_unsteady

  !> Evaluate time-integrated ROM objectives and their state derivatives.
  subroutine rom_problem_compute_objectives(this, objective, constraints)
    class(rom_problem_t), intent(inout) :: this
    real(kind=rp), intent(out) :: objective
    real(kind=rp), intent(out), optional :: constraints(:)
    type(vector_t), pointer :: state, state_gradient
    integer :: state_index, state_gradient_index
    integer :: step
    real(kind=rp) :: state_cost, wv, wb, viscous_time_weight, brinkman_time_weight
    real(kind=rp) :: viscous_volume, brinkman_volume, time_value
    if (.not. this%forward_ready) call neko_error( &
         'ROM objectives requested before run_forward_unsteady')
    this%running_gradient = 0.0_rp
    this%terminal_gradient = 0.0_rp
    this%cost_sum_a = 0.0_rp
    this%cost_sum_aa = 0.0_rp
    this%cost_total = 0.0_rp
    this%objective_value = 0.0_rp
    this%viscous_value = 0.0_rp
    this%brinkman_value = 0.0_rp
    viscous_volume = window_measure(this, this%viscous_start, this%viscous_end)
    brinkman_volume = window_measure(this, this%brinkman_start, this%brinkman_end)
    if (this%has_viscous .and. viscous_volume <= 0.0_rp) call neko_error( &
         'ROM viscous objective time window contains no timestep')
    if (this%has_brinkman .and. brinkman_volume <= 0.0_rp) call neko_error( &
         'ROM Brinkman objective time window contains no timestep')
    call neko_scratch_registry%request_vector(state, state_index, &
         this%states%get_nrows(), .true.)
    call neko_scratch_registry%request_vector(state_gradient, state_gradient_index, &
         this%states%get_nrows(), .true.)
    do step = 1, this%n_steps
       state%x = this%states%x(:,step+1)
       time_value = this%time_start + real(step, rp) * this%dt
       viscous_time_weight = 0.0_rp
       brinkman_time_weight = 0.0_rp
       if (this%has_viscous) viscous_time_weight = objective_time_weight( &
            this%viscous_start, this%viscous_end, time_value, viscous_volume, &
            this%dt)
       if (this%has_brinkman) brinkman_time_weight = objective_time_weight( &
            this%brinkman_start, this%brinkman_end, time_value, brinkman_volume, &
            this%dt)
       wv = viscous_time_weight * this%viscous_weight
       wb = brinkman_time_weight * this%brinkman_weight
       if (wv /= 0.0_rp) then
          state_cost = rom_viscous_value(state, &
               this%viscous_constant%x(1), this%viscous_linear, &
               this%viscous_quadratic, this%viscous_domain)
          this%viscous_value = this%viscous_value + this%dt * &
               viscous_time_weight * state_cost
          this%objective_value = this%objective_value + this%dt * wv * state_cost
          state_gradient%x = (this%viscous_linear%x + &
               matmul(this%viscous_quadratic%x, state%x)) / &
               max(this%viscous_domain, tiny(1.0_rp))
          if (step < this%n_steps) this%running_gradient%x(:,step+1) = &
               this%running_gradient%x(:,step+1) + wv * state_gradient%x
          if (step == this%n_steps) this%terminal_gradient%x = &
               this%terminal_gradient%x + this%dt * wv * state_gradient%x
       end if
       if (wb /= 0.0_rp) then
          state_cost = rom_brinkman_value(state, &
               this%candidate_objective_constant, this%candidate_objective_lift, &
               this%candidate_objective_brinkman, this%brinkman_domain)
          this%brinkman_value = this%brinkman_value + this%dt * &
               brinkman_time_weight * state_cost
          this%objective_value = this%objective_value + this%dt * wb * state_cost
          state_gradient%x = (this%candidate_objective_lift%x + &
               matmul(this%candidate_objective_brinkman%x, state%x)) / &
               this%brinkman_domain
          if (step < this%n_steps) this%running_gradient%x(:,step+1) = &
               this%running_gradient%x(:,step+1) + wb * state_gradient%x
          if (step == this%n_steps) this%terminal_gradient%x = &
               this%terminal_gradient%x + this%dt * wb * state_gradient%x
          this%cost_sum_a%x = this%cost_sum_a%x + this%dt * wb * state%x
          call add_outer_product(this%cost_sum_aa, state, state, this%dt * wb)
          this%cost_total = this%cost_total + this%dt * wb
       end if
    end do
    if (present(constraints)) then
       if (this%volume_is_max) then
          constraints(1) = this%candidate_volume / this%volume_domain - this%volume_limit
       else
          constraints(1) = this%volume_limit - this%candidate_volume / this%volume_domain
       end if
    end if
    objective = this%objective_value
    this%objectives_ready = .true.
    call neko_scratch_registry%relinquish_vector([state_index, state_gradient_index])
    nullify(state); nullify(state_gradient)
  end subroutine rom_problem_compute_objectives

  !> Integrate the reduced discrete adjoint backward through stored states.
  subroutine rom_problem_run_backward_unsteady(this)
    class(rom_problem_t), intent(inout) :: this
    type(vector_t), pointer :: lambda, next_adjoint, state
    integer :: lambda_index, next_adjoint_index, state_index
    integer :: step
    if (.not. this%objectives_ready) call neko_error( &
         'ROM adjoint requested before compute_objectives')
    call this%adjoint%free(); call this%lambda_sum%free(); call this%correlation%free()
    call this%adjoint%init(this%states%get_nrows(), this%n_steps + 1)
    call this%lambda_sum%init(this%states%get_nrows())
    call this%correlation%init(this%states%get_nrows(), this%states%get_nrows())
    call this%rom%discrete_adjoint_bdf1(this%states, this%dt, &
         this%running_gradient, this%terminal_gradient, this%adjoint)
    this%lambda_sum = 0.0_rp
    this%correlation = 0.0_rp
    do step = 1, this%n_steps
       if (step == 1) then
          call neko_scratch_registry%request_vector(lambda, lambda_index, &
               this%states%get_nrows(), .true.)
          call neko_scratch_registry%request_vector(next_adjoint, next_adjoint_index, &
               this%states%get_nrows(), .true.)
          call neko_scratch_registry%request_vector(state, state_index, &
               this%states%get_nrows(), .true.)
       end if
       next_adjoint%x = this%adjoint%x(:,step+1)
       state%x = this%states%x(:,step)
       call this%rom%solve_step_transpose(this%dt, next_adjoint, lambda)
       this%lambda_sum%x = this%lambda_sum%x + lambda%x
       call add_outer_product(this%correlation, lambda, state, 1.0_rp)
    end do
    this%backward_ready = .true.
    if (this%write_trajectory_fields .and. this%rom_outputs_initialized) &
         call rom_problem_write_adjoint_trajectory(this)
    call neko_scratch_registry%relinquish_vector([lambda_index, next_adjoint_index, state_index])
    nullify(lambda); nullify(next_adjoint); nullify(state)
  end subroutine rom_problem_run_backward_unsteady

  !> Assemble the mapped full-dimensional ROM design sensitivity density.
  subroutine rom_problem_compute_density(this, density)
    class(rom_problem_t), intent(in) :: this
    type(vector_t), intent(inout) :: density
    type(vector_t), pointer :: lift_square, direct, masked_direct, product, component_sum
    type(vector_ptr_t) :: mode_sum(3)
    integer :: lift_square_index, direct_index, masked_direct_index, product_index
    integer :: component_sum_index, mode_sum_index(3)
    integer :: j, c, r, n

    r = this%mass%get_nrows(); n = this%weights%size()
    if (density%size() /= n) call density%init(n)
    density = 0.0_rp
    mode_sum_index = 0
    call neko_scratch_registry%request_vector(lift_square, lift_square_index, n, .true.)
    call neko_scratch_registry%request_vector(direct, direct_index, n, .true.)
    call neko_scratch_registry%request_vector(masked_direct, masked_direct_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    call neko_scratch_registry%request_vector(component_sum, component_sum_index, n, .true.)
    call vector_vdot3(lift_square, this%lift_velocity(1), this%lift_velocity(2), &
         this%lift_velocity(3), this%lift_velocity(1), this%lift_velocity(2), &
         this%lift_velocity(3))
    direct = 0.0_rp
    if (this%has_brinkman .and. this%cost_total > 0.0_rp) then
       direct%x = this%cost_total * lift_square%x
       do c = 1, 3
          call neko_scratch_registry%request_vector(mode_sum(c)%ptr, mode_sum_index(c), &
               n, .true.)
          mode_sum(c)%ptr = 0.0_rp
          do j = 1, r
             call vector_add2s2(mode_sum(c)%ptr, this%velocity_modes(c,j), &
                  this%cost_sum_a%x(j))
          end do
          call vector_col3(component_sum, this%lift_velocity(c), mode_sum(c)%ptr)
          call vector_add2s2(direct, component_sum, 2.0_rp)
       end do
       do j = 1, r
          do c = 1, r
             component_sum = 0.0_rp
             call vector_col3(component_sum, this%velocity_modes(1,j), &
                  this%velocity_modes(1,c))
             call vector_col3(product, this%velocity_modes(2,j), this%velocity_modes(2,c))
             call vector_add2(component_sum, product)
             call vector_col3(product, this%velocity_modes(3,j), this%velocity_modes(3,c))
             call vector_add2(component_sum, product)
             call vector_add2s2(direct, component_sum, this%cost_sum_aa%x(j,c))
          end do
       end do
       call vector_col3(masked_direct, this%brinkman_mask, direct)
       call vector_cmult(masked_direct, 0.5_rp / &
            max(this%brinkman_domain, tiny(1.0_rp)))
       density = masked_direct
    end if
    do j = 1, r
       component_sum = 0.0_rp
       do c = 1, 3
          call vector_col3(product, this%velocity_modes(c,j), this%lift_velocity(c))
          call vector_add2(component_sum, product)
       end do
       call vector_add2s2(density, component_sum, -this%dt * this%lambda_sum%x(j))
       do c = 1, r
          component_sum = 0.0_rp
          do n = 1, 3
             call vector_col3(product, this%velocity_modes(n,j), this%velocity_modes(n,c))
             call vector_add2(component_sum, product)
          end do
          call vector_add2s2(density, component_sum, -this%dt * this%correlation%x(j,c))
       end do
    end do
    do c = 1, 3
       if (mode_sum_index(c) > 0) then
          call neko_scratch_registry%relinquish_vector(mode_sum_index(c))
          call mode_sum(c)%free()
       end if
    end do
    call neko_scratch_registry%relinquish_vector([lift_square_index, direct_index, &
         masked_direct_index, product_index, component_sum_index])
    nullify(lift_square); nullify(direct); nullify(masked_direct)
    nullify(product); nullify(component_sum)
  end subroutine rom_problem_compute_density

  !> Measure the configured objective time window on the ROM time grid.
  real(kind=rp) function window_measure(this, start_time, end_time) result(value)
    class(rom_problem_t), intent(in) :: this
    real(kind=rp), intent(in) :: start_time, end_time
    integer :: step
    real(kind=rp) :: t
    value = 0.0_rp
    do step = 1, this%n_steps
       t = this%time_start + real(step, rp) * this%dt
       if (t >= start_time - 1.0e-6_rp*this%dt .and. &
            t <= min(end_time, this%time_end) + 1.0e-6_rp*this%dt) &
            value = value + this%dt
    end do
  end function window_measure

  !> Return the normalised quadrature weight at one objective time.
  real(kind=rp) function objective_time_weight(start_time, end_time, &
       time_value, measure, dt) result(value)
    real(kind=rp), intent(in) :: start_time, end_time, time_value, measure, dt
    real(kind=rp) :: tolerance
    value = 0.0_rp
    tolerance = 1.0e-6_rp * dt
    if (measure > 0.0_rp .and. time_value >= start_time - tolerance .and. &
         time_value <= end_time + tolerance) value = 1.0_rp / measure
  end function objective_time_weight

  !> Parse the supported viscous and Brinkman objective configurations.
  subroutine read_rom_objectives(this, parameters)
    class(rom_problem_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(json_file) :: item
    character(len=:), allocatable :: typ, mask_name
    real(kind=rp) :: weight, start_time, end_time
    integer :: n_items, i

    this%has_viscous = .false.
    this%has_brinkman = .false.
    if (.not. parameters%valid_path('optimization.objectives')) then
       call neko_error('rom_mma requires viscous_dissipation and/or '&
            //'brinkman_dissipation objectives')
    end if
    call parameters%info('optimization.objectives', n_children = n_items)
    do i = 1, n_items
       call json_extract_item(parameters, 'optimization.objectives', i, item)
       call json_get(item, 'type', typ)
       call json_get_or_default(item, 'weight', weight, 1.0_rp)
       call json_get_or_default(item, 'mask_name', mask_name, '')
       call json_get_or_default(item, 'start_time', start_time, 0.0_rp)
       call json_get_or_default(item, 'end_time', end_time, huge(0.0_rp))
       select case (trim(typ))
       case ('viscous_dissipation')
          if (this%has_viscous) call neko_error('duplicate viscous objective')
          this%has_viscous = .true.
          this%viscous_weight = weight
          this%viscous_start = start_time
          this%viscous_end = end_time
          this%viscous_mask_name = trim(mask_name)
       case ('brinkman_dissipation')
          if (this%has_brinkman) call neko_error('duplicate Brinkman objective')
          this%has_brinkman = .true.
          this%brinkman_weight = weight
          this%brinkman_start = start_time
          this%brinkman_end = end_time
          this%brinkman_mask_name = trim(mask_name)
          block
             logical :: dealias_sensitivity, dealias_forcing
             call json_get_or_default(item, 'dealias_sensitivity', &
                  dealias_sensitivity, .false.)
             call json_get_or_default(item, 'dealias_forcing', &
                  dealias_forcing, dealias_sensitivity)
             if (dealias_sensitivity .or. dealias_forcing .or. &
                  dealias_sensitivity .neqv. dealias_forcing) then
                call neko_error('ROM Brinkman objective requires identical '&
                     //'native-GLL dealias_sensitivity=false and '&
                     //'dealias_forcing=false')
             end if
          end block
       case default
          call neko_error('unsupported objective in rom_mma: ' // trim(typ))
       end select
    end do
    if (.not. this%has_viscous .and. .not. this%has_brinkman) then
       call neko_error('rom_mma found no supported objective')
    end if
  end subroutine read_rom_objectives

  !> Parse the single supported volume constraint configuration.
  subroutine read_rom_constraint(this, parameters)
    class(rom_problem_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(json_file) :: item
    character(len=:), allocatable :: typ, mask_name
    logical :: is_max
    real(kind=rp) :: limit
    integer :: n_items

    if (.not. parameters%valid_path('optimization.constraints')) then
       call neko_error('rom_mma requires a volume constraint')
    end if
    call parameters%info('optimization.constraints', n_children = n_items)
    if (n_items /= 1) call neko_error('rom_mma requires one volume constraint')
    call json_extract_item(parameters, 'optimization.constraints', 1, item)
    call json_get(item, 'type', typ)
    if (trim(typ) /= 'volume') call neko_error('rom_mma only supports a volume constraint')
    call json_get_or_default(item, 'mask_name', mask_name, '')
    call json_get(item, 'limit', limit)
    call json_get_or_default(item, 'is_max', is_max, .false.)
    if ('mapping' .in. item) call neko_error( &
         'ROM volume constraint does not support a separate mapping')
    this%volume_mask_name = trim(mask_name)
    this%volume_limit = limit
    this%volume_is_max = is_max
  end subroutine read_rom_constraint

  !> Rebuild basis and fixed operators after a full-order refresh.
  subroutine rom_problem_refresh(this, simulation)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    call rom_problem_load_local_basis(this, simulation)
    call rom_problem_assemble_fixed(this, simulation)
    call rom_problem_setup_masks(this, simulation)
    call rom_problem_assemble_objective_operators(this)
  end subroutine rom_problem_refresh

  !> Receive and project the common local velocity basis from Python.
  subroutine rom_problem_load_local_basis(this, simulation)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    integer :: n_modes, n_fields
    real(kind=rp), allocatable :: local_basis(:, :)
    call rom_receive_basis(simulation%state_recover, simulation%fluid%u%x, &
         local_basis, n_modes, n_fields)
    if (n_fields /= 3) call neko_error('rom_mma requires a velocity-only common basis')
    call rom_project_basis(simulation, local_basis, n_modes)
    if (n_modes < 1) call neko_error('Python ROM peer returned an empty basis')
    this%basis%n_rows = size(local_basis, 1)
    this%basis%n_modes = n_modes
    this%basis%n_fields = n_fields
    if (this%basis%n_rows /= 3 * simulation%fluid%u%dof%size()) then
       call neko_error('ROM basis rows do not match local velocity fields')
    end if
    call this%basis%modes%free()
    call this%basis%modes%init(size(local_basis,1), size(local_basis,2))
    this%basis%modes%x = local_basis
    deallocate(local_basis)
  end subroutine rom_problem_load_local_basis

  !> Assemble mass, lift, convection and diffusion operators independent of design.
  subroutine rom_problem_assemble_fixed(this, simulation)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    integer :: n, r, i, j, c, d, offset
    type(matrix_t), pointer :: mass_local, linear_local, implicit_local
    type(matrix_t), allocatable :: quad_local(:)
    type(vector_t), pointer :: rhs_local, constant_local, initial_difference
    type(field_t), pointer :: field, gx, gy, gz
    integer :: mass_local_index, linear_local_index, implicit_local_index
    integer :: rhs_local_index, constant_local_index, initial_difference_index
    integer :: field_index, gx_index, gy_index, gz_index
    real(kind=rp) :: nu

    call rom_problem_free_arrays(this)
    n = simulation%fluid%u%dof%size()
    r = min(this%rank_total_max, this%basis%n_modes)
    call this%weights%init(n)
    this%weights%x = reshape(simulation%fluid%c_Xh%B, [n])
    allocate(this%velocity_modes(3,r), this%velocity_gradients(3,r,3))
    do i = 1, r
       do c = 1, 3
          do d = 1, 3
             call this%velocity_gradients(c,i,d)%init(n)
          end do
          call this%velocity_modes(c,i)%init(n)
       end do
    end do
    do c = 1, 3
       call this%lift_velocity(c)%free()
    end do
    call this%lift%copy_velocity(this%lift_velocity)
    call this%lift%sync_host()
    do c = 1, 3
       do d = 1, 3
          call this%lift_gradients(c,d)%init(n)
          this%lift_gradients(c,d)%x = reshape( &
               this%lift%gradient(c,d)%x, [n])
       end do
    end do
    do i = 1, r
       do c = 1, 3
          offset = (c - 1) * n
          this%velocity_modes(c, i)%x = this%basis%modes%x(offset + 1:offset + n, i)
       end do
    end do

    call neko_scratch_registry%request_field(field, field_index, &
         simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(gx, gx_index, simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(gy, gy_index, simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(gz, gz_index, simulation%fluid%u%dof, .true.)
    do i = 1, r
       do c = 1, 3
          field%x = reshape(this%velocity_modes(c,i)%x, shape(field%x))
          call grad(gx%x, gy%x, gz%x, field%x, simulation%fluid%c_Xh)
          this%velocity_gradients(c,i,1)%x = reshape(gx%x, [n])
          this%velocity_gradients(c,i,2)%x = reshape(gy%x, [n])
          this%velocity_gradients(c,i,3)%x = reshape(gz%x, [n])
       end do
    end do

    call neko_scratch_registry%request_matrix(mass_local, mass_local_index, r, r, .true.)
    call neko_scratch_registry%request_matrix(linear_local, linear_local_index, r, r, .true.)
    call neko_scratch_registry%request_matrix(implicit_local, implicit_local_index, &
         r, r, .true.)
    call neko_scratch_registry%request_vector(constant_local, constant_local_index, &
         r, .true.)
    allocate(quad_local(r))
    call rom_assemble_mass(this%weights, this%velocity_modes, mass_local)
    this%viscosity = real(simulation%fluid%mu%x(1, 1, 1, 1), rp)
    this%density = real(simulation%fluid%rho%x(1, 1, 1, 1), rp)
    nu = this%viscosity / this%density
    call rom_assemble_lift_convection(this%weights, this%velocity_modes, &
         this%velocity_gradients, this%lift_velocity, this%lift_gradients, &
         linear_local, constant_local)
    call rom_assemble_diffusion(this%weights, nu, this%velocity_gradients, &
         implicit_local, this%lift_gradients, constant_local)
    call rom_assemble_convection(this%weights, this%velocity_modes, &
         this%velocity_gradients, quad_local)
    do i = 1, r
       quad_local(i)%x = -quad_local(i)%x
    end do
    ! The ROM assemblers use vector_glsc2/vector_glsc3.  Those Neko
    ! contractions already reduce over NEKO_COMM, so their reduced operators
    ! must be copied directly rather than all-reduced a second time.
    call this%mass%init(r,r); call this%fixed_linear%init(r,r)
    call this%fixed_implicit%init(r,r); call this%fixed_constant%init(r)
    this%mass = mass_local
    this%fixed_linear = linear_local
    this%fixed_implicit = implicit_local
    this%fixed_constant = constant_local
    allocate(this%quadratic(r))
    do i = 1, r
       call this%quadratic(i)%init(r,r)
       this%quadratic(i) = quad_local(i)
    end do
    call this%viscous_constant%init(r); call this%viscous_linear%init(r)
    call this%viscous_quadratic%init(r,r); this%viscous_constant = 0.0_rp
    this%viscous_linear = 0.0_rp; this%viscous_quadratic = 0.0_rp
    ! The Neko forward initial condition is projected after subtracting the
    ! same spatial lift used to homogenise the streamed POD snapshots.
    call neko_scratch_registry%request_vector(rhs_local, rhs_local_index, r, .true.)
    call neko_scratch_registry%request_vector(initial_difference, &
         initial_difference_index, n, .true.)
    rhs_local = 0.0_rp
    do c = 1, 3
       call field_copy(field, this%lift%reference(c))
       call field_add2s2(field, this%lift%velocity(c), -1.0_rp, n)
       if (NEKO_BCKND_DEVICE .eq. 1) call field%copy_from(DEVICE_TO_HOST, sync = .true.)
       initial_difference%x = reshape(field%x, [n])
       do i = 1, r
          rhs_local%x(i) = rhs_local%x(i) + vector_glsc3( &
               this%velocity_modes(c,i), this%weights, initial_difference)
       end do
    end do
    call this%initial_state%init(r)
    call this%rom%init(this%mass, this%fixed_constant, this%fixed_linear, &
         this%quadratic, this%fixed_implicit)
    call this%rom%solve_mass(rhs_local, this%initial_state)
    this%rom_ready = .true.
    call neko_scratch_registry%relinquish_field([field_index, gx_index, gy_index, gz_index])
    nullify(field); nullify(gx); nullify(gy); nullify(gz)
    call neko_scratch_registry%relinquish_matrix([mass_local_index, linear_local_index, &
         implicit_local_index])
    nullify(mass_local); nullify(linear_local); nullify(implicit_local)
    do i = 1, r; call quad_local(i)%free(); end do
    deallocate(quad_local)
    call neko_scratch_registry%relinquish_vector(rhs_local_index)
    call neko_scratch_registry%relinquish_vector([constant_local_index, &
         initial_difference_index])
    nullify(rhs_local); nullify(constant_local); nullify(initial_difference)
  end subroutine rom_problem_assemble_fixed

  !> Resolve objective and constraint point-zone masks on the SEM mesh.
  subroutine rom_problem_setup_masks(this, simulation)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    integer :: n
    n = simulation%fluid%u%dof%size()
    call get_zone_mask(this%viscous_mask_name, n, this%viscous_mask)
    call get_zone_mask(this%brinkman_mask_name, n, this%brinkman_mask)
    call get_zone_mask(this%volume_mask_name, n, this%volume_mask)
    this%viscous_domain = vector_glsc2(this%weights, this%viscous_mask)
    this%brinkman_domain = vector_glsc2(this%weights, this%brinkman_mask)
    this%volume_domain = vector_glsc2(this%weights, this%volume_mask)
    ! vector_glsc2 performs the global MPI reduction internally.
    if (this%volume_domain <= 0.0_rp) call neko_error('ROM volume mask has zero measure')
    if (this%has_viscous .and. this%viscous_domain <= 0.0_rp) &
         call neko_error('ROM viscous objective mask has zero measure')
    if (this%has_brinkman .and. this%brinkman_domain <= 0.0_rp) &
         call neko_error('ROM Brinkman objective mask has zero measure')
  end subroutine rom_problem_setup_masks

  !> Convert a registered point zone into a binary local vector mask.
  subroutine get_zone_mask(name, n, mask)
    character(len=*), intent(in) :: name
    integer, intent(in) :: n
    type(vector_t), intent(inout) :: mask
    class(point_zone_t), pointer :: zone
    call mask%init(n); mask = 1.0_rp
    if (len_trim(name) == 0) return
    zone => neko_point_zone_registry%get_point_zone(trim(name))
    mask = 0.0_rp
    mask%x(zone%mask%get()) = 1.0_rp
  end subroutine get_zone_mask

  !> Assemble reduced viscous and design-dependent objective operators.
  subroutine rom_problem_assemble_objective_operators(this)
    class(rom_problem_t), intent(inout) :: this
    integer :: r
    r = this%mass%get_nrows()
    call this%viscous_constant%init(1)
    call this%viscous_linear%init(r)
    call this%viscous_quadratic%init(r,r)
    this%viscous_constant = 0.0_rp
    this%viscous_linear = 0.0_rp
    this%viscous_quadratic = 0.0_rp
    if (.not. this%has_viscous) return
    call rom_assemble_viscous_terms(this%weights, this%viscous_mask, &
         this%viscosity, this%velocity_gradients, this%lift_gradients, &
         this%viscous_constant, this%viscous_linear, this%viscous_quadratic)
    ! rom_assemble_viscous_terms is globally assembled through vector_glsc3.
  end subroutine rom_problem_assemble_objective_operators

  !> Check the prescribed velocity trace of an arbitrary ROM reconstruction.
  !!
  !! A copy of the reconstructed field is passed through Neko's actual
  !! velocity BC list and compared with the unmodified reconstruction.  Since
  !! the homogeneous modes are restricted by `rom_project_basis`, a zero
  !! difference is a direct boundary-trace check rather than a check against a
  !! separately recreated mask.
  !! @param[in] this ROM problem containing the spatial lift and modes.
  !! @param[inout] simulation Neko simulation carrying the native BC list.
  !! @param[in] coefficients Arbitrary reduced coefficients.
  !! @param[out] residual Global maximum boundary correction.
  subroutine rom_problem_verify_dirichlet_reconstruction(this, simulation, &
       coefficients, residual)
    class(rom_problem_t), intent(in) :: this
    type(simulation_t), intent(inout) :: simulation
    real(kind=rp), intent(in) :: coefficients(:)
    real(kind=rp), intent(out) :: residual
    type(field_t), pointer :: candidate_u, candidate_v, candidate_w
    type(field_t), pointer :: corrected_u, corrected_v, corrected_w, mode_field
    integer :: indices(7), i, n
    real(kind=rp) :: local_residual

    if (.not. this%lift%enabled) call neko_error( &
         'Cannot verify ROM BCs before constructing the spatial lift')
    if (size(coefficients) /= this%basis%n_modes) call neko_error( &
         'ROM BC verification coefficient count does not match the basis')
    n = this%weights%size()
    call neko_scratch_registry%request_field(candidate_u, indices(1), &
         simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(candidate_v, indices(2), &
         simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(candidate_w, indices(3), &
         simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(corrected_u, indices(4), &
         simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(corrected_v, indices(5), &
         simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(corrected_w, indices(6), &
         simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(mode_field, indices(7), &
         simulation%fluid%u%dof, .true.)
    call field_copy(candidate_u, this%lift%velocity(1))
    call field_copy(candidate_v, this%lift%velocity(2))
    call field_copy(candidate_w, this%lift%velocity(3))
    do i = 1, size(coefficients)
       mode_field%x = reshape(this%velocity_modes(1,i)%x, shape(mode_field%x))
       call field_add2s2(candidate_u, mode_field, coefficients(i), n)
       mode_field%x = reshape(this%velocity_modes(2,i)%x, shape(mode_field%x))
       call field_add2s2(candidate_v, mode_field, coefficients(i), n)
       mode_field%x = reshape(this%velocity_modes(3,i)%x, shape(mode_field%x))
       call field_add2s2(candidate_w, mode_field, coefficients(i), n)
    end do
    call field_copy(corrected_u, candidate_u)
    call field_copy(corrected_v, candidate_v)
    call field_copy(corrected_w, candidate_w)
    call simulation%fluid%bcs_vel%apply(corrected_u, corrected_v, corrected_w, &
         time = simulation%neko_case%time, strong = .true.)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call candidate_u%copy_from(DEVICE_TO_HOST, sync = .true.)
       call candidate_v%copy_from(DEVICE_TO_HOST, sync = .true.)
       call candidate_w%copy_from(DEVICE_TO_HOST, sync = .true.)
       call corrected_u%copy_from(DEVICE_TO_HOST, sync = .true.)
       call corrected_v%copy_from(DEVICE_TO_HOST, sync = .true.)
       call corrected_w%copy_from(DEVICE_TO_HOST, sync = .true.)
    end if
    local_residual = max(maxval(abs(corrected_u%x - candidate_u%x)), &
         max(maxval(abs(corrected_v%x - candidate_v%x)), &
         maxval(abs(corrected_w%x - candidate_w%x))))
    call MPI_Allreduce(local_residual, residual, 1, MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM)
    call neko_scratch_registry%relinquish_field(indices)
    nullify(candidate_u); nullify(candidate_v); nullify(candidate_w)
    nullify(corrected_u); nullify(corrected_v); nullify(corrected_w); nullify(mode_field)
  end subroutine rom_problem_verify_dirichlet_reconstruction

  !> Map one design and assemble its reduced Brinkman operators.
  subroutine rom_problem_assemble_candidate(this, design, brinkman, brinkman_lift, &
       objective_brinkman, objective_lift, objective_constant, volume)
    class(rom_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design
    type(matrix_t), intent(inout) :: brinkman, objective_brinkman
    type(vector_t), intent(inout) :: brinkman_lift, objective_lift
    real(kind=rp), intent(out) :: objective_constant
    real(kind=rp), intent(out), optional :: volume
    type(vector_t), pointer :: mapped, raw_design, local_lift
    type(matrix_t), pointer :: local
    integer :: mapped_index, raw_design_index, local_index, local_lift_index
    integer :: n, r

    select type (des => design)
    type is (brinkman_design_t)
       n = des%size()
       if (n /= this%weights%size()) call neko_error('ROM/design local size mismatch')
       call neko_scratch_registry%request_vector(mapped, mapped_index, n, .true.)
       call des%get_mapped_values(mapped)
    class default
       call neko_error('rom_mma requires a Brinkman design')
    end select
    r = this%mass%get_nrows()
    call neko_scratch_registry%request_matrix(local, local_index, r, r, .true.)
    call neko_scratch_registry%request_vector(local_lift, local_lift_index, r, .true.)
    call rom_assemble_brinkman(this%weights, mapped, this%velocity_modes, &
         this%lift_velocity, local, local_lift)
    ! Both assembly routines below use Neko's global vector_glsc* contractions.
    call brinkman%init(r,r); call brinkman_lift%init(r)
    brinkman = local
    brinkman_lift = local_lift
    call rom_assemble_brinkman_objective(this%weights, this%brinkman_mask, mapped, &
         this%velocity_modes, this%lift_velocity, local, local_lift, objective_constant)
    call objective_brinkman%init(r,r); call objective_lift%init(r)
    objective_brinkman = local
    objective_lift = local_lift
    if (present(volume)) then
       ! The production volume constraint integrates the unfiltered design
       ! variable, not the mapped Brinkman coefficient used by the PDE.
       call neko_scratch_registry%request_vector(raw_design, raw_design_index, &
            n, .true.)
       call design%get_values(raw_design)
       volume = vector_glsc3(this%weights, this%volume_mask, raw_design)
       call neko_scratch_registry%relinquish_vector(raw_design_index)
       nullify(raw_design)
    end if
    call neko_scratch_registry%relinquish_vector(mapped_index)
    call neko_scratch_registry%relinquish_matrix(local_index)
    call neko_scratch_registry%relinquish_vector(local_lift_index)
    nullify(mapped); nullify(raw_design); nullify(local); nullify(local_lift)
  end subroutine rom_problem_assemble_candidate

  !> Project a reduced density onto the physical design representation.
  subroutine rom_problem_project_density(this, design, density, output)
    class(rom_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design
    type(vector_t), intent(in) :: density
    type(vector_t), intent(inout) :: output
    output%x = density%x
    if (NEKO_BCKND_DEVICE .eq. 1) call output%copy_from(HOST_TO_DEVICE, sync = .true.)
    call design%map_backward(output)
    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(output)
       ! ROM density is a pointwise SEM sensitivity.  Convert it to the
       ! directional derivative used by MMA with the actual diagonal B,
       ! rather than the generic average-weight design smoothing.
       call des%convert_to_directional_derivative(output)
    class default
       call neko_error('rom_mma requires a Brinkman design')
    end select
  end subroutine rom_problem_project_density

  !> Project one sensitivity row through the design mapping.
  subroutine rom_problem_project_density_row(this, design, density, output)
    class(rom_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design
    real(kind=rp), intent(in) :: density(:)
    real(kind=rp), intent(out) :: output(:)
    type(vector_t), pointer :: input, result
    integer :: input_index, result_index
    call neko_scratch_registry%request_vector(input, input_index, size(density), .true.)
    call neko_scratch_registry%request_vector(result, result_index, size(density), .true.)
    input%x = density
    if (NEKO_BCKND_DEVICE .eq. 1) call input%copy_from(HOST_TO_DEVICE, sync = .true.)
    call design%map_backward(input)
    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(result)
       call des%convert_to_directional_derivative(result)
    class default
       call neko_error('rom_mma requires a Brinkman design')
    end select
    if (NEKO_BCKND_DEVICE .eq. 1) call result%copy_from(DEVICE_TO_HOST, sync = .true.)
    output = result%x
    call neko_scratch_registry%relinquish_vector([input_index, result_index])
    nullify(input); nullify(result)
  end subroutine rom_problem_project_density_row

  !> Accumulate a scaled outer product into a reduced matrix.
  subroutine add_outer_product(target, left, right, scale)
    type(matrix_t), intent(inout) :: target
    type(vector_t), intent(in) :: left, right
    real(kind=rp), intent(in) :: scale
    integer :: i, j
    do i = 1, left%size()
       do j = 1, right%size()
          target%x(i,j) = target%x(i,j) + scale * left%x(i) * right%x(j)
       end do
    end do
  end subroutine add_outer_product

  !> Initialise ROM field-output streams for a candidate iteration.
  subroutine rom_problem_initialize_outputs(this, design, simulation, iter)
    class(rom_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design
    type(simulation_t), intent(inout) :: simulation
    integer, intent(in) :: iter
    character(len=1024) :: output_path
    character(len=80) :: suffix
    call rom_problem_free_outputs(this)
    call filename_path(simulation%output_forward%file_%get_base_fname(), output_path)
    call this%rom_velocity(1)%init(simulation%fluid%u%dof, 'ROM_u')
    call this%rom_velocity(2)%init(simulation%fluid%v%dof, 'ROM_v')
    call this%rom_velocity(3)%init(simulation%fluid%w%dof, 'ROM_w')
    if (this%write_trajectory_fields) then
       call this%rom_adjoint_velocity(1)%init(simulation%fluid%u%dof, 'ROM_adj_u')
       call this%rom_adjoint_velocity(2)%init(simulation%fluid%v%dof, 'ROM_adj_v')
       call this%rom_adjoint_velocity(3)%init(simulation%fluid%w%dof, 'ROM_adj_w')
    end if
    write(suffix, '(i0)') iter
    call this%rom_forward_output%init('ROM_forward_fields_iter_'//trim(suffix), 3, &
         precision=sp, path=trim(output_path), format='fld')
    call this%rom_forward_output%fields%assign(1, this%rom_velocity(1)%field)
    call this%rom_forward_output%fields%assign(2, this%rom_velocity(2)%field)
    call this%rom_forward_output%fields%assign(3, this%rom_velocity(3)%field)
    call this%rom_design_output%init('ROM_design_iter_'//trim(suffix), 1, &
         precision=sp, path=trim(output_path), format='fld')
    call this%rom_sensitivity_output%init('ROM_sensitivity_iter_'//trim(suffix), 1, &
         precision=sp, path=trim(output_path), format='fld')
    if (this%write_trajectory_fields) then
       call this%rom_forward_trajectory_output%init( &
            'ROM_forward_trajectory_iter_'//trim(suffix), 3, precision=sp, &
            path=trim(output_path), format='fld')
       call this%rom_forward_trajectory_output%fields%assign(1, &
            this%rom_velocity(1)%field)
       call this%rom_forward_trajectory_output%fields%assign(2, &
            this%rom_velocity(2)%field)
       call this%rom_forward_trajectory_output%fields%assign(3, &
            this%rom_velocity(3)%field)
       call this%rom_adjoint_trajectory_output%init( &
            'ROM_adjoint_trajectory_iter_'//trim(suffix), 3, precision=sp, &
            path=trim(output_path), format='fld')
       call this%rom_adjoint_trajectory_output%fields%assign(1, &
            this%rom_adjoint_velocity(1)%field)
       call this%rom_adjoint_trajectory_output%fields%assign(2, &
            this%rom_adjoint_velocity(2)%field)
       call this%rom_adjoint_trajectory_output%fields%assign(3, &
            this%rom_adjoint_velocity(3)%field)
    end if
    select type (des => design)
    type is (brinkman_design_t)
       call des%configure_rom_outputs(this%rom_design_output, this%rom_sensitivity_output)
    class default
       call neko_error('rom_mma requires a Brinkman design')
    end select
    this%rom_outputs_initialized = .true.
  end subroutine rom_problem_initialize_outputs

  !> Write reconstructed ROM velocity, design and sensitivity fields.
  subroutine rom_problem_write_outputs(this, design, simulation, inner, state)
    class(rom_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design
    type(simulation_t), intent(inout) :: simulation
    integer, intent(in) :: inner
    real(kind=rp), intent(in) :: state(:)
    type(vector_t), pointer :: state_vector, reconstruction
    integer :: state_index, reconstruction_index
    integer :: n, r, c, i
    if (.not. this%rom_outputs_initialized) return
    n = simulation%fluid%u%dof%size(); r = size(state)
    call neko_scratch_registry%request_vector(state_vector, state_index, r, .true.)
    state_vector%x = state
    call neko_scratch_registry%request_vector(reconstruction, reconstruction_index, n, .true.)
    do c = 1, 3
       reconstruction = this%lift_velocity(c)
       do i = 1, r
          call vector_add2s2(reconstruction, this%velocity_modes(c,i), state_vector%x(i))
       end do
       this%rom_velocity(c)%field%x = reshape(reconstruction%x, &
            shape(this%rom_velocity(c)%field%x))
       if (NEKO_BCKND_DEVICE .eq. 1) call this%rom_velocity(c)%field%copy_from(HOST_TO_DEVICE, sync=.true.)
    end do
    call this%rom_forward_output%sample(real(inner, kind=dp))
    call this%rom_design_output%sample(real(inner, kind=dp))
    call this%rom_sensitivity_output%sample(real(inner, kind=dp))
    call neko_scratch_registry%relinquish_vector([state_index, reconstruction_index])
    nullify(state_vector); nullify(reconstruction)
  end subroutine rom_problem_write_outputs

  !> Reconstruct a three-component physical velocity from ROM coefficients.
  !! @param this ROM problem containing the common basis and lift.
  !! @param coefficients Reduced coefficients to reconstruct.
  !! @param fields Destination field wrappers for the three velocity components.
  !! @param include_lift Whether to add the inhomogeneous forward lift.
  subroutine rom_problem_reconstruct_velocity(this, coefficients, fields, &
       include_lift)
    class(rom_problem_t), intent(in) :: this
    real(kind=rp), intent(in) :: coefficients(:)
    type(field_wrapper_t), intent(inout) :: fields(3)
    logical, intent(in) :: include_lift
    type(vector_t), pointer :: reconstruction
    integer :: reconstruction_index
    integer :: n, r, c, i

    r = size(coefficients)
    n = this%weights%size()
    if (r /= this%basis%n_modes) call neko_error( &
         'ROM trajectory coefficient count does not match the basis')
    call neko_scratch_registry%request_vector(reconstruction, &
         reconstruction_index, n, .true.)
    do c = 1, 3
       if (include_lift) then
          reconstruction = this%lift_velocity(c)
       else
          reconstruction = 0.0_rp
       end if
       do i = 1, r
          call vector_add2s2(reconstruction, this%velocity_modes(c,i), &
               coefficients(i))
       end do
       fields(c)%field%x = reshape(reconstruction%x, &
            shape(fields(c)%field%x))
       if (NEKO_BCKND_DEVICE .eq. 1) &
            call fields(c)%field%copy_from(HOST_TO_DEVICE, sync=.true.)
    end do
    call neko_scratch_registry%relinquish_vector(reconstruction_index)
    nullify(reconstruction)
  end subroutine rom_problem_reconstruct_velocity

  !> Write the reconstructed ROM forward trajectory at physical-time samples.
  subroutine rom_problem_write_forward_trajectory(this)
    class(rom_problem_t), intent(inout) :: this
    integer :: stride, step
    real(kind=rp) :: sample_time, tolerance

    stride = nint(this%trajectory_output_interval / this%dt)
    tolerance = 1.0e-8_rp * max(1.0_rp, this%trajectory_output_interval)
    if (stride < 1 .or. abs(real(stride, rp) * this%dt - &
         this%trajectory_output_interval) > tolerance) call neko_error( &
         'ROM trajectory output interval is not an integer multiple of dt')
    do step = 0, this%n_steps
       if (mod(step, stride) == 0 .or. step == this%n_steps) then
          call rom_problem_reconstruct_velocity(this, &
               this%states%x(:,step+1), this%rom_velocity, .true.)
          sample_time = this%time_start + real(step, rp) * this%dt
          call this%rom_forward_trajectory_output%sample( &
               real(sample_time, kind=dp))
       end if
    end do
  end subroutine rom_problem_write_forward_trajectory

  !> Write the reconstructed ROM adjoint trajectory at physical-time samples.
  !! The reduced adjoint matrix is stored in ascending physical-time index;
  !! the integration that produced it runs backward through that matrix.
  subroutine rom_problem_write_adjoint_trajectory(this)
    class(rom_problem_t), intent(inout) :: this
    integer :: stride, step
    real(kind=rp) :: sample_time, tolerance

    stride = nint(this%trajectory_output_interval / this%dt)
    tolerance = 1.0e-8_rp * max(1.0_rp, this%trajectory_output_interval)
    if (stride < 1 .or. abs(real(stride, rp) * this%dt - &
         this%trajectory_output_interval) > tolerance) call neko_error( &
         'ROM trajectory output interval is not an integer multiple of dt')
    do step = 0, this%n_steps
       if (mod(step, stride) == 0 .or. step == this%n_steps) then
          call rom_problem_reconstruct_velocity(this, &
               this%adjoint%x(:,step+1), this%rom_adjoint_velocity, .false.)
          sample_time = this%time_start + real(step, rp) * this%dt
          call this%rom_adjoint_trajectory_output%sample( &
               real(sample_time, kind=dp))
       end if
    end do
  end subroutine rom_problem_write_adjoint_trajectory

  !> Close and release ROM field-output streams.
  subroutine rom_problem_free_outputs(this)
    class(rom_problem_t), intent(inout) :: this
    integer :: i
    call this%rom_forward_output%free(); call this%rom_design_output%free()
    call this%rom_sensitivity_output%free()
    call this%rom_forward_trajectory_output%free()
    call this%rom_adjoint_trajectory_output%free()
    do i = 1, 3
       call this%rom_velocity(i)%free()
       call this%rom_adjoint_velocity(i)%free()
    end do
    this%rom_outputs_initialized = .false.
  end subroutine rom_problem_free_outputs

  !> Release vector and matrix containers owned by the ROM problem.
  subroutine rom_problem_free_arrays(this)
    class(rom_problem_t), intent(inout) :: this
    integer :: i, j, k
    call this%weights%free(); call this%viscous_mask%free()
    call this%brinkman_mask%free(); call this%volume_mask%free()
    if (allocated(this%velocity_modes)) then
       do i = 1, size(this%velocity_modes,1)
          do j = 1, size(this%velocity_modes,2)
             call this%velocity_modes(i,j)%free()
          end do
       end do
       deallocate(this%velocity_modes)
    end if
    if (allocated(this%velocity_gradients)) then
       do i = 1, size(this%velocity_gradients,1)
          do j = 1, size(this%velocity_gradients,2)
             do k = 1, size(this%velocity_gradients,3)
                call this%velocity_gradients(i,j,k)%free()
             end do
          end do
       end do
       deallocate(this%velocity_gradients)
    end if
    do i = 1, 3; call this%lift_velocity(i)%free(); end do
    do i = 1, 3
       do j = 1, 3
          call this%lift_gradients(i,j)%free()
       end do
    end do
    call this%mass%free(); call this%fixed_constant%free()
    call this%fixed_linear%free(); call this%fixed_implicit%free()
    if (allocated(this%quadratic)) then
       do i = 1, size(this%quadratic); call this%quadratic(i)%free(); end do
       deallocate(this%quadratic)
    end if
    call this%viscous_constant%free(); call this%viscous_linear%free()
    call this%viscous_quadratic%free(); call this%initial_state%free()
  end subroutine rom_problem_free_arrays

end module rom_problem
