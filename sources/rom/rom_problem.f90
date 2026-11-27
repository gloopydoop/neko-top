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
  use fluid_pnpn, only: fluid_pnpn_t
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use field_output, only: field_output_t
  use field, only: field_t, field_wrapper_t
  use field_math, only: field_add2s2, field_copy, field_col3, field_invcol2, &
       field_sub3, field_rzero
  use matrix, only: matrix_t
  use rom_model, only: rom_model_t, rom_bdf1_export_context_t
  use rom_config, only: rom_config_t
  use rom_lift, only: rom_lift_t
  use rom_vector_potential, only: rom_vector_potential_t
       use rom_assembly, only: rom_assemble_mass, rom_assemble_brinkman, &
       rom_assemble_mass_coupled, rom_assemble_scalar_transport, &
       rom_assemble_lift_convection, rom_assemble_diffusion, &
       rom_assemble_convection, rom_assemble_cross_mass, &
       rom_assemble_brinkman_petrov, rom_assemble_lift_convection_petrov, &
       rom_assemble_diffusion_petrov, rom_assemble_convection_petrov, &
       rom_assemble_scalar_transport_petrov, rom_assemble_linearized_convection_tensor, &
       rom_assemble_scalar_mass, rom_assemble_scalar_diffusion_petrov, &
       rom_glsc3_local, rom_reduce_buffer
  use rom_objectives, only: rom_assemble_viscous_terms, &
       rom_assemble_viscous_cross_terms, rom_assemble_brinkman_objective, &
       rom_assemble_brinkman_objective_cross, rom_assemble_scalar_mixing_terms, &
       rom_assemble_scalar_mixing_cross_terms, &
       rom_viscous_value, rom_scalar_mixing_value, rom_brinkman_value
  use rom_basis_handoff, only: rom_receive_basis, rom_receive_petrov_bases, &
       rom_project_basis, rom_project_test_basis, rom_restrict_scalar_basis, &
       rom_orthonormalize_projected_device, &
       rom_orthonormalize_projected_coupled_device, &
       rom_project_velocity_column, rom_project_velocity_cotangent, &
       rom_divergence_norm, rom_bc_residual_norm, rom_weak_continuity_norm
  use operators, only: grad, cdtp
  use num_types, only: rp, sp, dp
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default, json_get_or_lookup, &
       json_extract_item
  use vector, only: vector_t, vector_ptr_t
  use scratch_registry, only: neko_scratch_registry
  use vector_math, only: vector_glsc2, vector_glsc3, vector_col3, vector_add2, &
       vector_add2s2, vector_cmult, vector_cmult2, vector_vdot3, vector_copy, &
       vector_sub3, vector_absval, vector_glmax, vector_rzero, vector_cadd, &
       vector_glmin, vector_invcol2, vector_col2
  use device, only: DEVICE_TO_HOST, HOST_TO_DEVICE
  use device_math, only: device_copy
  use neko_ext, only: field_to_vector, vector_to_field
  use registry, only: neko_registry
  use neko_config, only: NEKO_BCKND_DEVICE, NEKO_BCKND_XSMM
  use comm, only: pe_rank
  use mpi_f08, only: MPI_WTIME
  use point_zone, only: point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  use utils, only: neko_error, filename_path, mkdir
  use logger, only: neko_log
  use lapack_interfaces, only: dgesvd
  use mxm_wrapper, only: mxm
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
  implicit none
  private

  public :: rom_problem_t

  type :: rom_basis_storage_t
     integer :: n_rows = 0
     integer :: n_modes = 0
     integer :: n_fields = 0
     type(matrix_t) :: modes
  end type rom_basis_storage_t

  type :: rom_problem_t
     type(rom_model_t) :: rom
     type(rom_basis_storage_t) :: basis
     type(rom_basis_storage_t) :: test_basis

     character(len=32) :: basis_forward = 'forward'
     character(len=32) :: basis_adjoint = 'forward'
     character(len=64) :: formulation = 'forward_galerkin'
     character(len=32) :: sensitivity_source = 'exact_galerkin'
     logical :: petrov_galerkin = .false.
     logical :: forward_constrained_enriched = .false.
     logical :: forward_galerkin = .true.
     logical :: separate_pod = .false.
     logical :: forward_enriched_adjoint = .false.
     logical :: enriched_galerkin = .false.
     logical :: uses_dual_basis = .false.
     ! Opt-in design-dependent physical velocity modes.  Reference modes are
     ! retained separately so each candidate is re-anchored to the accepted
     ! outer basis rather than compounded from the previous candidate.
     logical :: basis_transform_adaptive = .false.
     character(len=32) :: basis_transform = 'none'
     logical :: differentiate_basis = .false.
     logical :: experimental_algebraic = .false.
     ! Diagnostic mode: evaluate candidate design terms with spatial modes held
     ! fixed at the current candidate, then resume adaptation when cleared.
     logical :: adaptive_modes_frozen = .false.
     logical :: diagnostic_audit = .false.
     real(kind=rp) :: basis_transform_tau = 1.0_rp
     type(vector_t), allocatable :: reference_velocity_modes(:, :)
     type(vector_t), allocatable :: reference_test_velocity_modes(:, :)
     type(vector_t), allocatable :: reference_velocity_gradients(:, :, :)
     ! Reference vector potentials are constructed once at the outer basis
     ! handoff and reused for every candidate in that inner cycle.
     type(rom_vector_potential_t) :: reference_forward_potential
     integer, allocatable :: reference_potential_forward_map(:)
     integer, allocatable :: reference_potential_test_map(:)
     type(vector_t) :: reference_chi, reference_D
     logical :: reference_modes_saved = .false.
     type(simulation_t), pointer :: simulation => null()

     integer :: inner_iterations = 1
     ! Actual ranks after the received POD blocks have been projected and
     ! orthonormalised.
     integer :: rank_forward = 0
     integer :: rank_adjoint = 0
     integer :: rank_enriched = 0
     integer :: history_outer_iterations = 3
     real(kind=rp) :: scalar_metric_weight = 1.0_rp
     real(kind=rp) :: admissibility_tolerance = 1.0e-8_rp
     logical :: enable_output = .true.
     logical :: rom_ready = .false.
     logical :: write_trajectory_fields = .false.
     real(kind=rp) :: trajectory_output_interval = 0.5_rp
     logical :: export_reduced_system = .false.
     character(len=:), allocatable :: reduced_system_export_directory
     logical :: reduced_system_exported = .false.
     integer :: reduced_system_export_count = 0
     logical :: reference_trajectory_exported = .false.

     integer :: substeps = 1
     integer :: n_steps = 0
     real(kind=rp) :: dt = 0.0_rp
     real(kind=rp) :: time_start = 0.0_rp
     real(kind=rp) :: time_end = 0.0_rp
     real(kind=rp) :: viscosity = 0.0_rp
     real(kind=rp) :: density = 1.0_rp

     logical :: has_viscous = .false.
     logical :: has_brinkman = .false.
     logical :: has_scalar = .false.
     logical :: has_scalar_mixing = .false.
     real(kind=rp) :: viscous_weight = 0.0_rp
     real(kind=rp) :: brinkman_weight = 0.0_rp
     real(kind=rp) :: scalar_mixing_weight = 0.0_rp
     real(kind=rp) :: viscous_start = 0.0_rp
     real(kind=rp) :: viscous_end = huge(0.0_rp)
     real(kind=rp) :: brinkman_start = 0.0_rp
     real(kind=rp) :: brinkman_end = huge(0.0_rp)
     real(kind=rp) :: scalar_mixing_start = 0.0_rp
     real(kind=rp) :: scalar_mixing_end = huge(0.0_rp)
     character(len=80) :: viscous_mask_name = ''
     character(len=80) :: brinkman_mask_name = ''
     character(len=80) :: scalar_mixing_mask_name = ''
     character(len=80) :: volume_mask_name = ''
     logical :: has_volume_constraint = .false.
     logical :: volume_is_max = .false.
     real(kind=rp) :: volume_limit = 0.0_rp
     real(kind=rp) :: volume_domain = 0.0_rp
     real(kind=rp) :: viscous_domain = 0.0_rp
     real(kind=rp) :: brinkman_domain = 0.0_rp
     real(kind=rp) :: scalar_mixing_domain = 0.0_rp
     real(kind=rp) :: scalar_reference_value = 0.5_rp

     type(vector_t) :: weights, viscous_mask, brinkman_mask, scalar_mixing_mask, volume_mask
     type(vector_t), allocatable :: velocity_modes(:, :)
     type(vector_t), allocatable :: velocity_gradients(:, :, :)
     type(vector_t), allocatable :: test_velocity_modes(:, :)
     type(vector_t), allocatable :: test_velocity_gradients(:, :, :)
     type(vector_t), allocatable :: scalar_modes(:)
     type(vector_t), allocatable :: scalar_gradients(:, :)
     type(vector_t), allocatable :: test_scalar_modes(:)
     type(vector_t), allocatable :: test_scalar_gradients(:, :)
     type(vector_t) :: scalar_diffusivity
     type(vector_t) :: scalar_lift, scalar_reference
     type(vector_t) :: scalar_lift_gradients(3)
     type(rom_lift_t) :: lift
     type(vector_t) :: lift_gradients(3,3)
     type(vector_t) :: lift_velocity(3)
     type(matrix_t) :: mass, trial_mass, fixed_linear, fixed_implicit
     type(matrix_t) :: enriched_mass, enriched_fixed_linear, enriched_fixed_implicit
     type(matrix_t) :: candidate_enriched_brinkman
     type(rom_model_t) :: enriched_adjoint_model
     integer :: enriched_adjoint_rank = 0
     real(kind=rp) :: cross_mass_sigma_min = 0.0_rp
     real(kind=rp) :: cross_mass_sigma_max = 0.0_rp
     real(kind=rp) :: cross_mass_condition = huge(0.0_rp)
     real(kind=rp) :: petrov_svd_relative_tolerance = 1.0e-8_rp
     integer :: cross_mass_retained_rank = 0
     real(kind=rp), allocatable :: cross_mass_singular_values(:)
     type(vector_t) :: fixed_constant
     type(matrix_t), allocatable :: quadratic(:)
     type(vector_t) :: enriched_viscous_linear
     type(matrix_t) :: enriched_viscous_cross
     type(vector_t) :: enriched_scalar_mixing_linear
     type(matrix_t) :: enriched_scalar_mixing_cross
     type(vector_t) :: candidate_objective_enriched_lift
     type(matrix_t) :: candidate_objective_enriched_cross
     type(vector_t) :: viscous_constant, viscous_linear
     type(matrix_t) :: viscous_quadratic
     type(vector_t) :: scalar_mixing_constant, scalar_mixing_linear
     type(matrix_t) :: scalar_mixing_quadratic
     type(vector_t) :: initial_state

     type(matrix_t) :: candidate_brinkman, candidate_objective_brinkman
     type(vector_t) :: candidate_brinkman_lift, candidate_objective_lift
     real(kind=rp) :: candidate_objective_constant = 0.0_rp
     real(kind=rp) :: candidate_volume = 0.0_rp
     ! Diagnostics for the explicitly selected algebraic experimental ROM.
     real(kind=rp) :: adaptive_strong_divergence = 0.0_rp
     real(kind=rp) :: adaptive_weak_continuity = 0.0_rp
     real(kind=rp) :: adaptive_boundary_residual = 0.0_rp
     real(kind=rp) :: adaptive_mass_condition = 1.0_rp
     ! Rank-local wall times for enriched/projected adjoint phases. The
     ! optimizer packs these with other phase timings before one max reduction.
     real(kind=rp) :: exact_adjoint_seconds = 0.0_rp
     real(kind=rp) :: projected_tensor_assembly_seconds = 0.0_rp
     real(kind=rp) :: projected_adjoint_integration_seconds = 0.0_rp
     real(kind=rp) :: temporal_correlation_seconds = 0.0_rp
     real(kind=rp) :: basis_adaptation_seconds = 0.0_rp
     real(kind=rp) :: spatial_mode_derivatives_seconds = 0.0_rp
     real(kind=rp) :: forward_local_operator_assembly_seconds = 0.0_rp
     real(kind=rp) :: adjoint_local_operator_assembly_seconds = 0.0_rp
     real(kind=rp) :: objective_operator_assembly_seconds = 0.0_rp
     real(kind=rp) :: initial_state_projection_seconds = 0.0_rp
     real(kind=rp) :: forward_ode_seconds = 0.0_rp
     real(kind=rp) :: spatial_sensitivity_seconds = 0.0_rp
     real(kind=rp) :: velocity_mode_cotangent_seconds = 0.0_rp
     real(kind=rp) :: scalar_mode_cotangent_seconds = 0.0_rp
     real(kind=rp) :: derivative_transpose_seconds = 0.0_rp
     real(kind=rp) :: basis_transform_pullback_seconds = 0.0_rp
     ! Vector-potential work is split by lifecycle: construction belongs to
     ! the accepted outer basis, while transformation and pullback accumulate
     ! over candidate evaluations in the corresponding inner cycle.
     real(kind=rp) :: potential_construction_seconds = 0.0_rp
     real(kind=rp) :: potential_poisson_solves = 0.0_rp
     real(kind=rp) :: potential_transform_seconds = 0.0_rp
     real(kind=rp) :: potential_transpose_seconds = 0.0_rp
     real(kind=rp) :: projected_tensor_cache_hit = 0.0_rp
     logical :: projected_tensor_cache_valid = .false.
     integer :: projected_tensor_cache_forward_rank = 0
     integer :: projected_tensor_cache_adjoint_rank = 0
     type(matrix_t), allocatable :: projected_convection_cache(:)
     logical :: scalar_mixing_cache_valid = .false.
     integer :: scalar_mixing_cache_forward_rank = 0
     integer :: scalar_mixing_cache_adjoint_rank = 0
     real(kind=rp) :: scalar_mixing_cache_hit = 0.0_rp
     logical :: scalar_transport_cache_valid = .false.
     integer :: scalar_transport_cache_forward_rank = 0
     integer :: scalar_transport_cache_adjoint_rank = 0
     type(matrix_t) :: scalar_mass_cache, scalar_adjoint_mass_cache
     type(matrix_t) :: scalar_diffusion_cache, scalar_adjoint_diffusion_cache
     type(vector_t) :: scalar_diffusion_lift_cache, scalar_adjoint_diffusion_lift_cache
     real(kind=rp) :: scalar_transport_cache_hit = 0.0_rp
     real(kind=rp) :: scalar_transport_cache_build_seconds = 0.0_rp
     type(matrix_t) :: states, adjoint, lambda_history, running_gradient, correlation, cost_sum_aa
     type(vector_t) :: terminal_gradient, lambda_sum, cost_sum_a
     type(matrix_t) :: enriched_running_gradient, exact_adjoint, enriched_correlation
     type(vector_t) :: enriched_terminal_gradient, enriched_lambda_sum
     real(kind=rp) :: cost_total = 0.0_rp
     real(kind=rp) :: objective_value = 0.0_rp
     real(kind=rp) :: viscous_value = 0.0_rp
     real(kind=rp) :: brinkman_value = 0.0_rp
     real(kind=rp) :: scalar_mixing_value = 0.0_rp
     logical :: forward_ready = .false.
     logical :: objectives_ready = .false.
     logical :: backward_ready = .false.

     type(field_output_t) :: rom_forward_output
     type(field_output_t) :: rom_design_output
     type(field_output_t) :: rom_sensitivity_output
     type(field_output_t) :: rom_forward_trajectory_output
     type(field_output_t) :: rom_adjoint_trajectory_output
     type(field_wrapper_t) :: rom_velocity(3)
     type(field_wrapper_t) :: rom_scalar
     type(field_wrapper_t) :: rom_adjoint_velocity(3)
     type(field_wrapper_t) :: rom_adjoint_scalar
     character(len=:), allocatable :: rom_output_root
     character(len=:), allocatable :: rom_outer_directory
     character(len=:), allocatable :: rom_inner_directory
     integer :: output_outer_iteration = -1
     integer :: output_inner_iteration = -1
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
     !> Write the initial reduced forward field, corresponding to
     !! simulation%write_forward(0) in the full-order optimizer.
     procedure, pass(this), public :: write_forward_output => &
          rom_problem_write_forward_output
     procedure, pass(this), public :: write_outputs => rom_problem_write_outputs
     procedure, pass(this), public :: write_inner_diagnostics => &
          rom_problem_write_inner_diagnostics
     procedure, pass(this), public :: compute_objectives => rom_problem_compute_objectives
     procedure, pass(this), public :: compute_density => rom_problem_compute_density
     procedure, pass(this), public :: compute_volume_constraint => &
          rom_problem_compute_volume_constraint
     procedure, pass(this), public :: verify_dirichlet_reconstruction => &
          rom_problem_verify_dirichlet_reconstruction
  end type rom_problem_t

contains

  !> Read a spatially uniform physical coefficient without using a stale host
  !! mirror on an accelerator backend.
  !!
  !! The reduced equations currently store viscosity and density as scalars.
  !! This helper therefore validates that the corresponding Neko field is
  !! actually uniform, while doing the validation with the same
  !! device-resident weighted contractions used by ROM assembly.  A
  !! nonuniform field is rejected explicitly rather than silently sampling
  !! its first host entry.
  subroutine rom_uniform_coefficient(field, weights, name, value)
    type(field_t), intent(in) :: field
    type(vector_t), intent(in) :: weights
    character(len=*), intent(in) :: name
    real(kind=rp), intent(out) :: value
    type(vector_t), pointer :: values, ones, residual
    integer :: values_index, ones_index, residual_index, n
    real(kind=rp) :: measure, mean_square, relative_variation

    n = weights%size()
    call neko_scratch_registry%request_vector(values, values_index, n, .true.)
    call neko_scratch_registry%request_vector(ones, ones_index, n, .true.)
    call neko_scratch_registry%request_vector(residual, residual_index, n, .true.)
    call field_to_vector(values, field)
    ones = 1.0_rp
    measure = vector_glsc2(weights, ones)
    if (measure <= 0.0_rp) then
       call neko_scratch_registry%relinquish_vector([values_index, ones_index, residual_index])
       nullify(values); nullify(ones); nullify(residual)
       call neko_error('ROM coefficient check encountered a zero SEM measure')
    end if
    value = vector_glsc3(values, weights, ones) / measure
    call vector_copy(residual, values)
    call vector_add2s2(residual, ones, -value)
    mean_square = vector_glsc3(residual, weights, residual)
    relative_variation = sqrt(max(0.0_rp, mean_square)) / &
         max(abs(value) * sqrt(measure), sqrt(tiny(1.0_rp)))
    call neko_scratch_registry%relinquish_vector([values_index, ones_index, residual_index])
    nullify(values); nullify(ones); nullify(residual)
    if (relative_variation > 1.0e-12_rp) then
       call neko_error('ROM requires spatially uniform '//trim(name)// &
            '; use a coefficient-aware ROM formulation for nonuniform fields')
    end if
  end subroutine rom_uniform_coefficient

  !> Read ROM time-integration, rank, objective and admissibility settings.
  subroutine rom_problem_init_from_json(this, parameters, problem)
    class(rom_problem_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    class(problem_t), intent(in) :: problem
    type(rom_config_t) :: config
    logical :: dealias, variable_dt
    integer :: time_order

    call this%free()
    call config%init(parameters)
    this%basis_forward = trim(config%basis_forward)
    this%basis_adjoint = trim(config%basis_adjoint)
    this%formulation = trim(config%formulation)
    this%sensitivity_source = trim(config%sensitivity_source)
    this%inner_iterations = config%inner_iterations
    this%history_outer_iterations = config%history_outer_iterations
    this%rank_forward = 0
    this%rank_adjoint = 0
    this%rank_enriched = 0
    this%substeps = config%substeps
    this%scalar_metric_weight = config%scalar_metric_weight
    this%admissibility_tolerance = config%admissibility_tolerance
    this%petrov_svd_relative_tolerance = config%petrov_svd_relative_tolerance
    this%enable_output = config%enable_output
    this%write_trajectory_fields = config%write_trajectory_fields
    this%trajectory_output_interval = config%trajectory_output_interval
    this%export_reduced_system = config%export_reduced_system
    this%reduced_system_export_directory = config%reduced_system_export_directory
    this%forward_galerkin = this%formulation == 'forward_galerkin'
    this%separate_pod = this%formulation == 'separate_pod'
    this%forward_enriched_adjoint = this%formulation == 'forward_enriched_adjoint'
    this%enriched_galerkin = this%formulation == 'enriched_galerkin'
    this%forward_constrained_enriched = this%separate_pod .or. &
         this%forward_enriched_adjoint
    this%uses_dual_basis = this%separate_pod .or. &
         this%forward_enriched_adjoint .or. this%enriched_galerkin
    this%basis_transform = trim(config%basis_transform)
    this%basis_transform_adaptive = trim(config%basis_transform) /= 'none'
    this%differentiate_basis = config%differentiate_basis
    this%experimental_algebraic = config%experimental_algebraic
    this%basis_transform_tau = config%basis_transform_tau
    this%diagnostic_audit = config%diagnostic_audit
    call neko_log%message('ROM effective sensitivity source: '// &
         trim(this%sensitivity_source))
    ! Canonical formulations use Galerkin forward dynamics.  The old flag is
    ! retained only as a compatibility field for diagnostics; no new case
    ! silently selects a Petrov forward operator.
    this%petrov_galerkin = .false.
    call json_get(parameters, 'case.numerics.dealias', dealias)
    if (dealias) call neko_error('rom_mma requires case.numerics.dealias=false')
    call json_get_or_lookup(parameters, 'case.numerics.time_order', time_order)
    if (time_order /= 1) call neko_error('rom_mma currently requires case.numerics.time_order=1')
    call json_get_or_default(parameters, 'case.time.variable_timestep', &
         variable_dt, .false.)
    if (variable_dt) call neko_error('rom_mma does not support variable timestep yet')
    call read_rom_objectives(this, parameters)
    call read_rom_constraint(this, parameters)
  end subroutine rom_problem_init_from_json

  !> Initialise the ROM time grid and problem-dependent physical metadata.
  subroutine rom_problem_initialize(this, problem, design, simulation)
    class(rom_problem_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), target, intent(inout) :: simulation
    integer :: dns_steps
    real(kind=rp) :: dns_dt
    this%simulation => simulation
    ! Use the pointer association, since simulation_t's legacy convenience
    ! flag is not populated by every Neko scalar initialisation path.
    this%has_scalar = associated(simulation%scalars)
    if (this%has_scalar) then
       if (.not. allocated(simulation%scalars%scalar_fields) .or. &
            size(simulation%scalars%scalar_fields) /= 1) then
          call neko_error('rom_mma supports exactly one passive scalar')
       end if
    end if
    if (this%has_scalar_mixing .and. .not. this%has_scalar) then
       call neko_error('scalar_mixing objective requires one passive scalar')
    end if
    ! Reset Neko before constructing the lift so its retained reference fields
    ! are taken from the same canonical initial-condition path used by the
    ! forward solver, including user-defined fluid and scalar ICs.
    call simulation%reset()
    ! Construct the lift before the full-order run starts.  State recovery has
    ! already used the same ROM-local construction for its initial stream.
    call this%lift%init(simulation%neko_case)
    if (.not. this%lift%enabled) call neko_error( &
         'ROM requires a spatial lift constructed from the case velocity BCs')
    call problem%compute(design, simulation)
    ! The ROM optimizer uses this existing full-order initialization solve to
    ! establish the time grid and objective state.  Mirror mma_optimizer's
    ! iteration-zero forward output while those physical fields are current.
    if (this%enable_output) call simulation%write_forward(0)
    call problem%compute_sensitivity(design, simulation)
    if (this%enable_output) then
       call simulation%write_adjoint(0)
    end if
    dns_steps = simulation%n_timesteps
    dns_dt = real(simulation%neko_case%time%dt, rp)
    this%time_start = real(simulation%neko_case%time%start_time, rp)
    this%time_end = real(simulation%neko_case%time%end_time, rp)
    if (dns_steps < 1 .or. dns_dt <= 0.0_rp) call neko_error('invalid Neko ROM time grid')
    if (abs(this%time_start + dns_dt * real(dns_steps, rp) - this%time_end) &
         > 100.0_rp * epsilon(this%time_end) * max(1.0_rp, abs(this%time_end))) then
       call neko_error('rom_mma requires a fixed timestep ending exactly at case.time.end_time')
    end if
    ! Refine only the reduced time grid.  The full-order/DNS timestep and
    ! simulation timestep count remain unchanged for the surrounding Neko
    ! lifecycle; every ROM step uses the same BDF1 map with dt/substeps.
    this%n_steps = dns_steps * this%substeps
    this%dt = dns_dt / real(this%substeps, rp)
    if (abs(this%time_start + this%dt * real(this%n_steps, rp) - this%time_end) &
         > 1000.0_rp * epsilon(this%time_end) * max(1.0_rp, abs(this%time_end))) then
       call neko_error('ROM substepping does not preserve the configured physical time interval')
    end if
    call simulation%reset()
  end subroutine rom_problem_initialize

  !> Release ROM operators, trajectories, outputs and basis storage.
  subroutine rom_problem_free(this)
    class(rom_problem_t), intent(inout) :: this
    call this%rom%free()
    call this%basis%modes%free()
    call this%test_basis%modes%free()
    call rom_problem_free_reference_basis(this)
    call rom_problem_free_arrays(this)
    call rom_problem_free_outputs(this)
    call rom_problem_free_trajectories(this)
    call this%lift%free()
    nullify(this%simulation)
    this%rom_ready = .false.
    this%reduced_system_exported = .false.
    this%reference_trajectory_exported = .false.
  end subroutine rom_problem_free

  !> Release only forward and adjoint trajectory storage.
  subroutine rom_problem_free_trajectories(this)
    class(rom_problem_t), intent(inout) :: this
    call this%candidate_brinkman%free()
    call this%candidate_enriched_brinkman%free()
    call this%candidate_brinkman_lift%free()
    call this%candidate_objective_brinkman%free()
    call this%candidate_objective_lift%free()
    call this%candidate_objective_enriched_lift%free()
    call this%candidate_objective_enriched_cross%free()
    call this%states%free(); call this%adjoint%free(); call this%lambda_history%free()
    call this%exact_adjoint%free()
    call this%running_gradient%free(); call this%terminal_gradient%free()
    call this%enriched_running_gradient%free()
    call this%enriched_terminal_gradient%free()
    call this%lambda_sum%free(); call this%correlation%free()
    call this%enriched_lambda_sum%free(); call this%enriched_correlation%free()
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
    type(simulation_t), target, intent(inout) :: simulation
    integer :: expected_rows
    this%projected_tensor_cache_valid = .false.
    this%scalar_mixing_cache_valid = .false.
    this%scalar_transport_cache_valid = .false.
    this%simulation => simulation
    if (this%uses_dual_basis) call neko_error( &
         'This ROM formulation requires separate streamed forward and adjoint bases')
    if (n_modes < 1 .or. size(basis, 2) /= n_modes .or. &
         (n_fields /= 3 .and. n_fields /= 4)) then
       call neko_error('ROM diagnostic basis must contain velocity or velocity/scalar fields')
    end if
    expected_rows = 3 * simulation%fluid%u%dof%size()
    if (n_fields == 4) expected_rows = 4 * simulation%fluid%u%dof%size()
    if (size(basis, 1) /= expected_rows) then
       call neko_error('ROM diagnostic basis rows do not match the velocity fields')
    end if
    ! The supplied-basis path can be entered without rom_problem_initialize.
    ! Reset the live Neko case explicitly so the initial-state projection below
    ! uses the case-defined IC rather than whatever fields a previous caller
    ! left in the simulation object.
    call simulation%reset()
    if (.not. this%lift%enabled) call this%lift%init(simulation%neko_case)
    if (.not. this%lift%enabled) call neko_error( &
         'ROM diagnostic basis requires a spatial velocity lift')
    call this%basis%modes%free()
    call this%test_basis%modes%free()
    this%basis%n_rows = size(basis, 1)
    this%basis%n_modes = n_modes
    this%basis%n_fields = n_fields
    this%rank_forward = n_modes
    this%rank_adjoint = n_modes
    this%rank_enriched = n_modes
    call this%basis%modes%init(size(basis,1), size(basis,2))
    this%basis%modes%x = basis
    call rom_problem_restrict_scalar_basis(this, simulation)
    ! This matrix remains host staging.  assemble_fixed performs the one
    ! basis-refresh upload for each retained mesh mode vector; a second device
    ! copy of the complete basis would duplicate transfer and storage.
    this%has_scalar = n_fields == 4
    if (this%has_scalar .and. .not. associated(simulation%scalars)) then
       call neko_error('scalar ROM basis supplied to a scalar-free simulation')
    end if
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
    this%basis_adaptation_seconds = 0.0_rp
    this%spatial_mode_derivatives_seconds = 0.0_rp
    this%forward_local_operator_assembly_seconds = 0.0_rp
    this%adjoint_local_operator_assembly_seconds = 0.0_rp
    this%objective_operator_assembly_seconds = 0.0_rp
    this%initial_state_projection_seconds = 0.0_rp
    this%forward_ode_seconds = 0.0_rp
    this%exact_adjoint_seconds = 0.0_rp
    this%projected_tensor_assembly_seconds = 0.0_rp
    this%projected_adjoint_integration_seconds = 0.0_rp
    this%temporal_correlation_seconds = 0.0_rp
    this%spatial_sensitivity_seconds = 0.0_rp
    this%velocity_mode_cotangent_seconds = 0.0_rp
    this%scalar_mode_cotangent_seconds = 0.0_rp
    this%derivative_transpose_seconds = 0.0_rp
    this%basis_transform_pullback_seconds = 0.0_rp
    this%potential_transform_seconds = 0.0_rp
    this%potential_transpose_seconds = 0.0_rp
    this%projected_tensor_cache_hit = 0.0_rp
    this%scalar_mixing_cache_hit = 0.0_rp
    this%scalar_transport_cache_hit = 0.0_rp
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
    type(vector_t), intent(inout) :: density
    real(kind=rp), intent(out), optional :: constraint_density(:, :)
    type(vector_t), pointer :: density_vector
    integer :: density_index
    real(kind=rp) :: timer_start
    if (.not. this%backward_ready) call neko_error( &
         'ROM sensitivity requested before run_backward_unsteady')
    call neko_scratch_registry%request_vector(density_vector, density_index, &
         this%weights%size(), .true.)
    timer_start = real(MPI_WTIME(),rp)
    call this%compute_density(density_vector)
    this%spatial_sensitivity_seconds = real(MPI_WTIME(),rp) - timer_start
    call vector_copy(density, density_vector)
    if (present(constraint_density)) then
       if (.not. this%has_volume_constraint) then
          constraint_density = 0.0_rp
       else
          ! This legacy array output is a standalone diagnostic adapter, not
          ! part of the production MMA path. Its host-array contract requires
          ! one explicit readback; production callers use the device matrix
          ! assembled by compute_volume_constraint below.
          if (NEKO_BCKND_DEVICE .eq. 1) call this%volume_mask%copy_from( &
               DEVICE_TO_HOST, sync = .true.)
          if (this%volume_is_max) then
             constraint_density(1,:) = this%volume_mask%x / this%volume_domain
          else
             constraint_density(1,:) = -this%volume_mask%x / this%volume_domain
          end if
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
  !! the production ROM case.  The returned derivative is the raw constant
  !! sensitivity (with the optimization-domain mask).  The ROM optimizer then
  !! applies brinkman_design_t%project_sensitivity, matching the production
  !! MMA convention that uses the constant average SEM weight avg_B.
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

    if (.not. this%has_volume_constraint) then
       value = -1.0_rp
       ! This is a host-small dummy return.  Keep its host mirror explicit;
       ! the unconstrained MMA path never consumes a device derivative row.
       sensitivity%x = 0.0_rp
       return
    end if

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
    ! Match the production volume_constraint_t representation.  Its
    ! sensitivity is a constant function over the optimization domain; the
    ! standard MMA path subsequently applies brinkman_design_t%project_sensitivity
    ! and therefore multiplies it by the constant average SEM weight avg_B.
    ! Do not insert the local SEM diagonal here: that would make the ROM use
    ! the mathematically exact quadrature derivative while the full-order MMA
    ! path intentionally uses the smoothed constant-weight convention.
    call vector_copy(derivative, this%volume_mask)
    if (this%volume_is_max) then
       call vector_cmult(derivative, 1.0_rp / this%volume_domain)
    else
       call vector_cmult(derivative, -1.0_rp / this%volume_domain)
    end if
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! The ROM currently exposes one design-only volume constraint.  The
       ! one-row matrix is contiguous, so this is a direct device copy.
       call device_copy(sensitivity%x_d, derivative%x_d, n)
    else
       sensitivity%x(1,:) = derivative%x
    end if

    call neko_scratch_registry%relinquish_vector([raw_index, derivative_index])
    nullify(raw_design); nullify(derivative)
  end subroutine rom_problem_compute_volume_constraint

  !> Integrate the reduced forward equations through the configured time grid.
  subroutine rom_problem_run_forward_unsteady(this, design, compute_constraint)
    class(rom_problem_t), intent(inout) :: this
    type(vector_t), pointer :: candidate_constant
    type(matrix_t), pointer :: candidate_linear
    type(matrix_t), pointer :: candidate_implicit
    type(matrix_t) :: enriched_implicit
    type(rom_bdf1_export_context_t) :: export_context
    integer :: candidate_constant_index, candidate_linear_index
    integer :: candidate_implicit_index
    integer :: export_inner
    character(len=32) :: export_inner_name
    character(len=32) :: export_eval_name
    class(design_t), intent(inout) :: design
    logical, intent(in), optional :: compute_constraint
    integer :: r, e
    real(kind=rp) :: timer_start
    logical :: need_constraint
    need_constraint = .false.
    if (present(compute_constraint)) need_constraint = compute_constraint
    call rom_problem_free_trajectories(this)
    r = this%mass%get_nrows()
    e = r
    if (this%forward_constrained_enriched) e = this%enriched_adjoint_rank
    call this%candidate_brinkman%init(r,r)
    call this%candidate_brinkman_lift%init(r)
    call this%candidate_objective_brinkman%init(r,r)
    call this%candidate_objective_lift%init(r)
    if (this%forward_constrained_enriched) then
       call this%candidate_enriched_brinkman%init(e,e)
       call this%candidate_objective_enriched_cross%init(e,r)
       call this%candidate_objective_enriched_lift%init(e)
    end if
    if (need_constraint) then
       call rom_problem_assemble_candidate(this, design, this%candidate_brinkman, &
            this%candidate_brinkman_lift, this%candidate_objective_brinkman, &
            this%candidate_objective_lift, this%candidate_objective_constant, &
            this%candidate_enriched_brinkman, this%candidate_objective_enriched_cross, &
            this%candidate_objective_enriched_lift, this%candidate_volume)
    else
       call rom_problem_assemble_candidate(this, design, this%candidate_brinkman, &
            this%candidate_brinkman_lift, this%candidate_objective_brinkman, &
            this%candidate_objective_lift, this%candidate_objective_constant, &
            this%candidate_enriched_brinkman, this%candidate_objective_enriched_cross, &
            this%candidate_objective_enriched_lift)
    end if
    call neko_scratch_registry%request_vector(candidate_constant, candidate_constant_index, &
         r, .true.)
    call neko_scratch_registry%request_matrix(candidate_linear, candidate_linear_index, &
         r, r, .true.)
    call neko_scratch_registry%request_matrix(candidate_implicit, candidate_implicit_index, &
         r, r, .true.)
    ! Reduced operators are host-owned dense algebra.  Do not use matrix/vector
    ! container assignment here: on an accelerator it copies the device mirror,
    ! while assembly intentionally populated the host mirror with scalar
    ! contractions.
    candidate_constant%x = this%fixed_constant%x - this%candidate_brinkman_lift%x
    candidate_linear%x = this%fixed_linear%x
    ! Brinkman damping is part of the implicit BDF1 operator.  The positive
    ! assembled damping matrix therefore enters with a minus sign in the
    ! residual operator, so M - dt*(I_fixed - B) adds +dt*B to the solved
    ! step matrix.
    candidate_implicit%x = this%fixed_implicit%x - this%candidate_brinkman%x
    if (this%export_reduced_system .and. .not. this%reduced_system_exported) then
       call rom_problem_export_reduced_system(this, candidate_constant, &
            candidate_linear, candidate_implicit)
       this%reduced_system_exported = .true.
    end if
    call this%rom%init(this%mass, candidate_constant, candidate_linear, &
         this%quadratic, candidate_implicit)
    if (this%forward_constrained_enriched) then
       call enriched_implicit%init(e,e)
       enriched_implicit%x = this%enriched_fixed_implicit%x - &
            this%candidate_enriched_brinkman%x
       call this%enriched_adjoint_model%init_linearized_adjoint( &
            this%enriched_mass, enriched_implicit)
       call enriched_implicit%free()
    end if
    call this%states%init(r, this%n_steps + 1)
    call this%running_gradient%init(r, this%n_steps)
    call this%terminal_gradient%init(r)
    if (this%forward_constrained_enriched) then
       call this%enriched_running_gradient%init(e, this%n_steps)
       call this%enriched_terminal_gradient%init(e)
    end if
    call this%cost_sum_a%init(r)
    call this%cost_sum_aa%init(r,r)
    timer_start = real(MPI_WTIME(),rp)
    if (this%export_reduced_system) then
       export_context%directory = trim(this%reduced_system_export_directory)
       export_context%outer_iteration = this%output_outer_iteration
       export_context%basis_forward = trim(this%basis_forward)
       export_context%basis_adjoint = trim(this%basis_adjoint)
       export_context%trial_basis_rank = r
       export_context%test_basis_rank = r
       if (this%uses_dual_basis) &
            export_context%test_basis_rank = this%test_basis%n_modes
       export_context%forward_basis_rank = this%basis%n_modes
       export_context%adjoint_basis_rank = this%test_basis%n_modes
       export_context%enriched_basis_rank = this%enriched_adjoint_rank
       export_context%time_start = this%time_start
       export_context%evaluation_index = this%reduced_system_export_count
       if (this%rom_outputs_initialized .and. this%output_inner_iteration < 0) then
          export_context%directory = trim(export_context%directory)//'/reference'
          export_context%design_id = 'reference'
          export_context%inner_iteration = -1
       else if (this%rom_outputs_initialized) then
          ! The previous output sample is the current optimizer iterate. The
          ! forward call evaluates the next candidate, recorded after this
          ! integration by write_outputs.
          export_inner = this%output_inner_iteration + 1
          write(export_inner_name, '("inner_",I3.3)') export_inner
          export_context%directory = trim(export_context%directory)//'/'// &
               trim(export_inner_name)
          export_context%design_id = trim(export_inner_name)
          export_context%inner_iteration = export_inner
       else
          export_inner = this%reduced_system_export_count
          write(export_inner_name, '("evaluation_",I3.3)') export_inner
          export_context%directory = trim(export_context%directory)//'/'// &
               trim(export_inner_name)
          export_context%design_id = trim(export_inner_name)
          export_context%outer_iteration = -1
          export_context%inner_iteration = -1
       end if
       if (this%rom_outputs_initialized .and. &
            this%reduced_system_export_count > 0) then
          write(export_eval_name, '("evaluation_",I3.3)') &
               this%reduced_system_export_count
          export_context%directory = trim(export_context%directory)//'/'// &
               trim(export_eval_name)
       end if
       call this%rom%forward_bdf1(this%initial_state, this%dt, this%states, &
            export_context)
       this%reduced_system_export_count = this%reduced_system_export_count + 1
    else
       call this%rom%forward_bdf1(this%initial_state, this%dt, this%states)
    end if
    this%forward_ode_seconds = this%forward_ode_seconds + &
         real(MPI_WTIME(),rp) - timer_start
    call neko_scratch_registry%relinquish_vector(candidate_constant_index)
    call neko_scratch_registry%relinquish_matrix(candidate_linear_index)
    call neko_scratch_registry%relinquish_matrix(candidate_implicit_index)
    nullify(candidate_constant); nullify(candidate_linear); nullify(candidate_implicit)
    this%forward_ready = .true.
    if (this%export_reduced_system .and. .not. this%reference_trajectory_exported) then
       call rom_problem_export_reference_trajectory(this)
       this%reference_trajectory_exported = .true.
    end if
    if (this%write_trajectory_fields .and. this%rom_outputs_initialized) &
         call rom_problem_write_forward_trajectory(this)
  end subroutine rom_problem_run_forward_unsteady

  !> Export the exact reduced operators used by the current BDF1 solve.
  !!
  !! This is a reference-design inspection boundary.  It writes host-owned
  !! reduced arrays only; it does not alter any operator or trajectory.  The
  !! exported equation is
  !!
  !!   (H-dt*L_i) a_(n+1) = H a_n + dt*(b+L_e a_n+q(a_n)),
  !!
  !! with q_i(a) = a^T Q_i a.  `b`, `L_e` and `L_i` are the actual candidate
  !! objects passed to rom_model_t, including the reference-design Brinkman
  !! lift and damping contributions.
  subroutine rom_problem_export_reduced_system(this, constant, linear, implicit_linear)
    class(rom_problem_t), intent(inout) :: this
    type(vector_t), intent(inout) :: constant
    type(matrix_t), intent(in) :: linear, implicit_linear
    integer :: unit, i
    integer :: r, trial_rank, test_rank
    character(len=512) :: filename
    character(len=32) :: logical_string

    if (pe_rank /= 0) return
    call mkdir(trim(this%reduced_system_export_directory))
    r = this%mass%get_nrows()
    trial_rank = this%basis%n_modes
    test_rank = trial_rank
    if (this%uses_dual_basis) &
         test_rank = this%test_basis%n_modes

    filename = trim(this%reduced_system_export_directory)//'/H.csv'
    call rom_problem_write_matrix_csv(trim(filename), this%mass%x)
    filename = trim(this%reduced_system_export_directory)//'/trial_mass.csv'
    call rom_problem_write_matrix_csv(trim(filename), this%trial_mass%x)
    filename = trim(this%reduced_system_export_directory)//'/b.csv'
    call rom_problem_write_vector_csv(trim(filename), constant%x)
    filename = trim(this%reduced_system_export_directory)//'/linear_explicit.csv'
    call rom_problem_write_matrix_csv(trim(filename), linear%x)
    filename = trim(this%reduced_system_export_directory)//'/linear_implicit.csv'
    call rom_problem_write_matrix_csv(trim(filename), implicit_linear%x)
    filename = trim(this%reduced_system_export_directory)//'/Q.csv'
    call rom_problem_write_quadratic_csv(trim(filename), this%quadratic)
    filename = trim(this%reduced_system_export_directory)//'/a0.csv'
    call rom_problem_write_vector_csv(trim(filename), this%initial_state%x)

    ! These split files make the lift and design-dependent pieces inspectable
    ! without changing the compact runtime files above.
    filename = trim(this%reduced_system_export_directory)//'/fixed_constant.csv'
    call rom_problem_write_vector_csv(trim(filename), this%fixed_constant%x)
    filename = trim(this%reduced_system_export_directory)//'/brinkman_lift.csv'
    call rom_problem_write_vector_csv(trim(filename), this%candidate_brinkman_lift%x)
    filename = trim(this%reduced_system_export_directory)//'/fixed_linear.csv'
    call rom_problem_write_matrix_csv(trim(filename), this%fixed_linear%x)
    filename = trim(this%reduced_system_export_directory)//'/fixed_implicit.csv'
    call rom_problem_write_matrix_csv(trim(filename), this%fixed_implicit%x)
    filename = trim(this%reduced_system_export_directory)//'/brinkman_operator.csv'
    call rom_problem_write_matrix_csv(trim(filename), this%candidate_brinkman%x)

    ! The final paired bases are written in Fortran column-major order.  The
    ! metadata records the shape and layout for numpy.reshape(order='F').
    filename = trim(this%reduced_system_export_directory)//'/trial_basis.bin'
    call rom_problem_write_matrix_binary(trim(filename), this%basis%modes%x)
    if (this%uses_dual_basis) then
       filename = trim(this%reduced_system_export_directory)//'/test_basis.bin'
       call rom_problem_write_matrix_binary(trim(filename), this%test_basis%modes%x)
    end if
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call this%weights%copy_from(DEVICE_TO_HOST, sync=.true.)
       do i = 1, 3
          call this%lift_velocity(i)%copy_from(DEVICE_TO_HOST, sync=.true.)
       end do
       if (this%has_scalar) call this%scalar_lift%copy_from(DEVICE_TO_HOST, sync=.true.)
    end if
    filename = trim(this%reduced_system_export_directory)//'/weights.bin'
    call rom_problem_write_vector_binary(trim(filename), this%weights%x)
    filename = trim(this%reduced_system_export_directory)//'/lift_velocity.bin'
    call rom_problem_write_vector_matrix_binary(trim(filename), this%lift_velocity)
    if (this%has_scalar) then
       filename = trim(this%reduced_system_export_directory)//'/scalar_lift.bin'
       call rom_problem_write_vector_binary(trim(filename), this%scalar_lift%x)
    end if

    filename = trim(this%reduced_system_export_directory)//'/metadata.json'
    open(newunit=unit, file=trim(filename), status='replace', action='write')
    write(unit, '(A)') '{'
    write(unit, '(A,I0,A)') '  "n_state": ', r, ','
    write(unit, '(A,I0,A)') '  "retained_rank": ', r, ','
    write(unit, '(A,I0,A)') '  "trial_basis_rank": ', trial_rank, ','
    write(unit, '(A,I0,A)') '  "test_basis_rank": ', test_rank, ','
    write(unit, '(A,I0,A)') '  "enriched_adjoint_rank": ', &
         this%enriched_adjoint_rank, ','
    write(unit, '(A,A,A)') '  "basis_forward": "', &
         trim(this%basis_forward), '",'
    write(unit, '(A,A,A)') '  "basis_adjoint": "', &
         trim(this%basis_adjoint), '",'
    write(unit, '(A,I0,A)') '  "forward_retained_rank": ', &
         this%rank_forward, ','
    write(unit, '(A,I0,A)') '  "adjoint_retained_rank": ', &
         this%rank_adjoint, ','
    write(unit, '(A,I0,A)') '  "enriched_retained_rank": ', &
         this%rank_enriched, ','
    write(unit, '(A,I0,A)') '  "trial_basis_rows": ', this%basis%n_rows, ','
    write(unit, '(A,I0,A)') '  "test_basis_rows": ', this%test_basis%n_rows, ','
    write(unit, '(A,A,A)') '  "formulation": "', trim(this%formulation), '",'
    write(unit, '(A,A,A)') '  "sensitivity_source": "', &
         trim(this%sensitivity_source), '",'
    write(unit, '(A,A,A)') '  "basis_transform": "', &
         trim(this%basis_transform), '",'
    logical_string = rom_problem_json_logical(this%differentiate_basis)
    write(unit, '(A,A,A)') '  "differentiate_basis": ', trim(logical_string), ','
    write(unit, '(A,ES24.16,A)') '  "basis_transform_tau": ', &
         this%basis_transform_tau, ','
    logical_string = rom_problem_json_logical(this%petrov_galerkin)
    write(unit, '(A,A,A)') '  "petrov_galerkin": ', trim(logical_string), ','
    logical_string = rom_problem_json_logical(this%forward_constrained_enriched)
    write(unit, '(A,A,A)') '  "forward_constrained_enriched": ', trim(logical_string), ','
    logical_string = rom_problem_json_logical(this%forward_galerkin)
    write(unit, '(A,A,A)') '  "forward_galerkin": ', trim(logical_string), ','
    logical_string = rom_problem_json_logical(this%separate_pod)
    write(unit, '(A,A,A)') '  "separate_pod": ', trim(logical_string), ','
    logical_string = rom_problem_json_logical(this%forward_enriched_adjoint)
    write(unit, '(A,A,A)') '  "forward_enriched_adjoint": ', trim(logical_string), ','
    logical_string = rom_problem_json_logical(this%enriched_galerkin)
    write(unit, '(A,A,A)') '  "enriched_galerkin": ', trim(logical_string), ','
    logical_string = rom_problem_json_logical(this%has_scalar)
    write(unit, '(A,A,A)') '  "has_scalar": ', trim(logical_string), ','
    write(unit, '(A,ES24.16,A)') '  "scalar_metric_weight": ', &
         this%scalar_metric_weight, ','
    write(unit, '(A,ES24.16,A)') '  "time_start": ', this%time_start, ','
    write(unit, '(A,ES24.16,A)') '  "time_end": ', this%time_end, ','
    write(unit, '(A,ES24.16,A)') '  "dt": ', this%dt, ','
    write(unit, '(A,I0,A)') '  "n_steps": ', this%n_steps, ','
    write(unit, '(A,I0,A)') '  "substeps": ', this%substeps, ','
    write(unit, '(A,ES24.16,A)') '  "viscosity": ', this%viscosity, ','
    write(unit, '(A,ES24.16,A)') '  "density": ', this%density, ','
    write(unit, '(A,ES24.16,A)') '  "cross_mass_sigma_min": ', &
         this%cross_mass_sigma_min, ','
    write(unit, '(A,ES24.16,A)') '  "cross_mass_sigma_max": ', &
         this%cross_mass_sigma_max, ','
    write(unit, '(A,ES26.16E3,A)') '  "cross_mass_condition": ', &
         this%cross_mass_condition, ','
    write(unit, '(A,I0,A)') '  "cross_mass_retained_rank": ', &
         this%cross_mass_retained_rank, ','
    write(unit, '(A)', advance='no') '  "cross_mass_singular_values": ['
    if (allocated(this%cross_mass_singular_values)) then
       do i = 1, size(this%cross_mass_singular_values)
          if (i > 1) write(unit, '(A)', advance='no') ', '
          write(unit, '(ES24.16)', advance='no') this%cross_mass_singular_values(i)
       end do
    end if
    write(unit, '(A)') '],'
    write(unit, '(A)') '  "basis_layout": "field-major, column-major float64 binary",'
    write(unit, '(A)') '  "q_indexing": "Q[i,j,k] = sum_jk coefficient for residual i",'
    write(unit, '(A)') '  "b_definition": "fixed_constant - candidate_brinkman_lift",'
    write(unit, '(A)') '  "linear_explicit_definition": "fixed_linear",'
    write(unit, '(A)') '  "linear_implicit_definition": "fixed_implicit - candidate_brinkman",'
    write(unit, '(A)') '  "files": {'
    write(unit, '(A)') '    "H": "H.csv", "trial_mass": "trial_mass.csv",'
    write(unit, '(A)') '    "b": "b.csv", "linear_explicit": "linear_explicit.csv",'
    write(unit, '(A)') '    "linear_implicit": "linear_implicit.csv", "Q": "Q.csv",'
    write(unit, '(A)') '    "a0": "a0.csv", "trial_basis": "trial_basis.bin",'
    if (this%uses_dual_basis) then
       write(unit, '(A)') '    "test_basis": "test_basis.bin",'
    end if
    write(unit, '(A)') '    "weights": "weights.bin", "lift_velocity": "lift_velocity.bin",'
    if (this%has_scalar) then
       write(unit, '(A)') '    "scalar_lift": "scalar_lift.bin",'
    end if
    write(unit, '(A)') '    "reference_trajectory": "reference_trajectory.csv"'
    write(unit, '(A)') '  }'
    write(unit, '(A)') '}'
    close(unit)
  end subroutine rom_problem_export_reduced_system

  !> Export the trajectory produced by the same reduced BDF1 call.
  subroutine rom_problem_export_reference_trajectory(this)
    class(rom_problem_t), intent(in) :: this
    integer :: unit, step, mode, r
    character(len=512) :: filename
    if (pe_rank /= 0) return
    r = this%states%get_nrows()
    filename = trim(this%reduced_system_export_directory)//'/reference_trajectory.csv'
    open(newunit=unit, file=trim(filename), status='replace', action='write')
    write(unit, '(A)', advance='no') 'step,time'
    do mode = 1, r
       write(unit, '(A,I0)', advance='no') ',a_', mode
    end do
    write(unit, '(A)') ''
    do step = 0, this%n_steps
       write(unit, '(I0,A,ES24.16)', advance='no') step, ',', &
            this%time_start + real(step, rp) * this%dt
       do mode = 1, r
          write(unit, '(A,ES24.16)', advance='no') ',', this%states%x(mode,step+1)
       end do
       write(unit, '(A)') ''
    end do
    close(unit)
  end subroutine rom_problem_export_reference_trajectory

  !> Write a dense reduced matrix as comma-separated rows.
  subroutine rom_problem_write_matrix_csv(filename, values)
    character(len=*), intent(in) :: filename
    real(kind=rp), intent(in) :: values(:, :)
    integer :: unit, i, j
    open(newunit=unit, file=trim(filename), status='replace', action='write')
    do i = 1, size(values,1)
       do j = 1, size(values,2)
          if (j > 1) write(unit, '(A)', advance='no') ','
          write(unit, '(ES24.16)', advance='no') values(i,j)
       end do
       write(unit, '(A)') ''
    end do
    close(unit)
  end subroutine rom_problem_write_matrix_csv

  !> Write a reduced vector as one value per CSV row.
  subroutine rom_problem_write_vector_csv(filename, values)
    character(len=*), intent(in) :: filename
    real(kind=rp), intent(in) :: values(:)
    integer :: unit, i
    open(newunit=unit, file=trim(filename), status='replace', action='write')
    do i = 1, size(values)
       write(unit, '(ES24.16)') values(i)
    end do
    close(unit)
  end subroutine rom_problem_write_vector_csv

  !> Write the residual-row quadratic tensor as i,j,k,value CSV records.
  subroutine rom_problem_write_quadratic_csv(filename, values)
    character(len=*), intent(in) :: filename
    type(matrix_t), intent(in) :: values(:)
    integer :: unit, i, j, k
    open(newunit=unit, file=trim(filename), status='replace', action='write')
    write(unit, '(A)') 'i,j,k,value'
    do i = 1, size(values)
       do j = 1, values(i)%get_nrows()
          do k = 1, values(i)%get_ncols()
             write(unit, '(I0,A,I0,A,I0,A,ES24.16)') i, ',', j, ',', k, ',', &
                  values(i)%x(j,k)
          end do
       end do
    end do
    close(unit)
  end subroutine rom_problem_write_quadratic_csv

  !> Write a dense matrix in raw Fortran column-major float64 form.
  subroutine rom_problem_write_matrix_binary(filename, values)
    character(len=*), intent(in) :: filename
    real(kind=rp), intent(in) :: values(:, :)
    integer :: unit
    open(newunit=unit, file=trim(filename), status='replace', &
         access='stream', form='unformatted', action='write')
    write(unit) values
    close(unit)
  end subroutine rom_problem_write_matrix_binary

  !> Write a reduced vector in raw float64 form.
  subroutine rom_problem_write_vector_binary(filename, values)
    character(len=*), intent(in) :: filename
    real(kind=rp), intent(in) :: values(:)
    integer :: unit
    open(newunit=unit, file=trim(filename), status='replace', &
         access='stream', form='unformatted', action='write')
    write(unit) values
    close(unit)
  end subroutine rom_problem_write_vector_binary

  !> Write three mesh vectors consecutively in component-major raw form.
  subroutine rom_problem_write_vector_matrix_binary(filename, values)
    character(len=*), intent(in) :: filename
    type(vector_t), intent(in) :: values(:)
    integer :: unit, i
    open(newunit=unit, file=trim(filename), status='replace', &
         access='stream', form='unformatted', action='write')
    do i = 1, size(values)
       write(unit) values(i)%x
    end do
    close(unit)
  end subroutine rom_problem_write_vector_matrix_binary

  !> Return a JSON-compatible lowercase logical literal.
  function rom_problem_json_logical(value) result(string)
    logical, intent(in) :: value
    character(len=5) :: string
    if (value) then
       string = 'true '
    else
       string = 'false'
    end if
  end function rom_problem_json_logical

  !> Evaluate time-integrated ROM objectives and their state derivatives.
  subroutine rom_problem_compute_objectives(this, objective, constraints)
    class(rom_problem_t), intent(inout) :: this
    real(kind=rp), intent(out) :: objective
    real(kind=rp), intent(out), optional :: constraints(:)
    type(vector_t), pointer :: state, state_gradient, enriched_state_gradient
    integer :: state_index, state_gradient_index, enriched_state_gradient_index
    integer :: step
    real(kind=rp) :: state_cost, wv, wb, ws, viscous_time_weight
    real(kind=rp) :: brinkman_time_weight, scalar_mixing_time_weight
    real(kind=rp) :: viscous_volume, brinkman_volume, scalar_mixing_volume, time_value
    if (.not. this%forward_ready) call neko_error( &
         'ROM objectives requested before run_forward_unsteady')
    this%running_gradient%x = 0.0_rp
    this%terminal_gradient%x = 0.0_rp
    if (this%forward_constrained_enriched) then
       this%enriched_running_gradient%x = 0.0_rp
       this%enriched_terminal_gradient%x = 0.0_rp
    end if
    this%cost_sum_a%x = 0.0_rp
    this%cost_sum_aa%x = 0.0_rp
    this%cost_total = 0.0_rp
    this%objective_value = 0.0_rp
    this%viscous_value = 0.0_rp
    this%brinkman_value = 0.0_rp
    this%scalar_mixing_value = 0.0_rp
    viscous_volume = window_measure(this, this%viscous_start, this%viscous_end)
    brinkman_volume = window_measure(this, this%brinkman_start, this%brinkman_end)
    scalar_mixing_volume = window_measure(this, this%scalar_mixing_start, &
         this%scalar_mixing_end)
    if (this%has_viscous .and. viscous_volume <= 0.0_rp) call neko_error( &
         'ROM viscous objective time window contains no timestep')
    if (this%has_brinkman .and. brinkman_volume <= 0.0_rp) call neko_error( &
         'ROM Brinkman objective time window contains no timestep')
    if (this%has_scalar_mixing .and. scalar_mixing_volume <= 0.0_rp) &
         call neko_error('ROM scalar mixing objective time window contains no timestep')
    call neko_scratch_registry%request_vector(state, state_index, &
         this%states%get_nrows(), .true.)
    call neko_scratch_registry%request_vector(state_gradient, state_gradient_index, &
         this%states%get_nrows(), .true.)
    if (this%forward_constrained_enriched) then
       call neko_scratch_registry%request_vector(enriched_state_gradient, &
            enriched_state_gradient_index, this%enriched_adjoint_rank, .true.)
    end if
    do step = 1, this%n_steps
       state%x = this%states%x(:,step+1)
       time_value = this%time_start + real(step, rp) * this%dt
       viscous_time_weight = 0.0_rp
       brinkman_time_weight = 0.0_rp
       scalar_mixing_time_weight = 0.0_rp
       if (this%has_viscous) viscous_time_weight = objective_time_weight( &
            this%viscous_start, this%viscous_end, time_value, viscous_volume, &
            this%dt)
       if (this%has_brinkman) brinkman_time_weight = objective_time_weight( &
            this%brinkman_start, this%brinkman_end, time_value, brinkman_volume, &
            this%dt)
       if (this%has_scalar_mixing) scalar_mixing_time_weight = objective_time_weight( &
            this%scalar_mixing_start, this%scalar_mixing_end, time_value, &
            scalar_mixing_volume, this%dt)
       wv = viscous_time_weight * this%viscous_weight
       wb = brinkman_time_weight * this%brinkman_weight
       ws = scalar_mixing_time_weight * this%scalar_mixing_weight
       if (this%has_viscous) then
          state_cost = rom_viscous_value(state, &
               this%viscous_constant%x(1), this%viscous_linear, &
               this%viscous_quadratic, this%viscous_domain)
          this%viscous_value = this%viscous_value + this%dt * &
               viscous_time_weight * state_cost
          this%objective_value = this%objective_value + this%dt * wv * state_cost
          state_gradient%x = (this%viscous_linear%x + &
               matmul(this%viscous_quadratic%x, state%x)) / &
               max(this%viscous_domain, tiny(1.0_rp))
          if (this%forward_constrained_enriched) then
             enriched_state_gradient%x = (this%enriched_viscous_linear%x + &
                  matmul(this%enriched_viscous_cross%x, state%x)) / &
                  max(this%viscous_domain, tiny(1.0_rp))
             if (step < this%n_steps) this%enriched_running_gradient%x(:,step+1) = &
                  this%enriched_running_gradient%x(:,step+1) + wv * &
                  enriched_state_gradient%x
             if (step == this%n_steps) this%enriched_terminal_gradient%x = &
                  this%enriched_terminal_gradient%x + this%dt * wv * &
                  enriched_state_gradient%x
          end if
          if (step < this%n_steps) this%running_gradient%x(:,step+1) = &
               this%running_gradient%x(:,step+1) + wv * state_gradient%x
          if (step == this%n_steps) this%terminal_gradient%x = &
               this%terminal_gradient%x + this%dt * wv * state_gradient%x
       end if
       if (this%has_brinkman) then
          state_cost = rom_brinkman_value(state, &
               this%candidate_objective_constant, this%candidate_objective_lift, &
               this%candidate_objective_brinkman, this%brinkman_domain)
          this%brinkman_value = this%brinkman_value + this%dt * &
               brinkman_time_weight * state_cost
          this%objective_value = this%objective_value + this%dt * wb * state_cost
          state_gradient%x = (this%candidate_objective_lift%x + &
               matmul(this%candidate_objective_brinkman%x, state%x)) / &
               this%brinkman_domain
          if (this%forward_constrained_enriched) then
             enriched_state_gradient%x = (this%candidate_objective_enriched_lift%x + &
                  matmul(this%candidate_objective_enriched_cross%x, state%x)) / &
                  this%brinkman_domain
             if (step < this%n_steps) this%enriched_running_gradient%x(:,step+1) = &
                  this%enriched_running_gradient%x(:,step+1) + wb * &
                  enriched_state_gradient%x
             if (step == this%n_steps) this%enriched_terminal_gradient%x = &
                  this%enriched_terminal_gradient%x + this%dt * wb * &
                  enriched_state_gradient%x
          end if
          if (step < this%n_steps) this%running_gradient%x(:,step+1) = &
               this%running_gradient%x(:,step+1) + wb * state_gradient%x
          if (step == this%n_steps) this%terminal_gradient%x = &
               this%terminal_gradient%x + this%dt * wb * state_gradient%x
          this%cost_sum_a%x = this%cost_sum_a%x + this%dt * wb * state%x
          call add_outer_product(this%cost_sum_aa, state, state, this%dt * wb)
          this%cost_total = this%cost_total + this%dt * wb
       end if
       if (this%has_scalar_mixing) then
          state_cost = rom_scalar_mixing_value(state, this%scalar_mixing_constant%x(1), &
               this%scalar_mixing_linear, this%scalar_mixing_quadratic, &
               this%scalar_mixing_domain)
          this%scalar_mixing_value = this%scalar_mixing_value + this%dt * &
               scalar_mixing_time_weight * state_cost
          this%objective_value = this%objective_value + this%dt * ws * state_cost
          state_gradient%x = (this%scalar_mixing_linear%x + &
               matmul(this%scalar_mixing_quadratic%x, state%x)) / &
               max(this%scalar_mixing_domain, tiny(1.0_rp))
          if (this%forward_constrained_enriched) then
             enriched_state_gradient%x = (this%enriched_scalar_mixing_linear%x + &
                  matmul(this%enriched_scalar_mixing_cross%x, state%x)) / &
                  max(this%scalar_mixing_domain, tiny(1.0_rp))
             if (step < this%n_steps) this%enriched_running_gradient%x(:,step+1) = &
                  this%enriched_running_gradient%x(:,step+1) + ws * &
                  enriched_state_gradient%x
             if (step == this%n_steps) this%enriched_terminal_gradient%x = &
                  this%enriched_terminal_gradient%x + this%dt * ws * &
                  enriched_state_gradient%x
          end if
          if (step < this%n_steps) this%running_gradient%x(:,step+1) = &
               this%running_gradient%x(:,step+1) + ws * state_gradient%x
          if (step == this%n_steps) this%terminal_gradient%x = &
               this%terminal_gradient%x + this%dt * ws * state_gradient%x
       end if
    end do
    if (present(constraints)) then
       if (.not. this%has_volume_constraint) then
          constraints(1) = -1.0_rp
       else if (this%volume_is_max) then
          constraints(1) = this%candidate_volume / this%volume_domain - this%volume_limit
       else
          constraints(1) = this%volume_limit - this%candidate_volume / this%volume_domain
       end if
    end if
    objective = this%objective_value
    this%objectives_ready = .true.
    call neko_scratch_registry%relinquish_vector([state_index, state_gradient_index])
    if (this%forward_constrained_enriched) then
       call neko_scratch_registry%relinquish_vector(enriched_state_gradient_index)
    end if
    nullify(state); nullify(state_gradient); nullify(enriched_state_gradient)
  end subroutine rom_problem_compute_objectives

  !> Integrate the reduced discrete adjoint backward through stored states.
  subroutine rom_problem_run_backward_unsteady(this)
    class(rom_problem_t), intent(inout) :: this
    type(vector_t), pointer :: lambda, state
    integer :: lambda_index, state_index
    integer :: step
    real(kind=rp) :: timer_start
    if (.not. this%objectives_ready) call neko_error( &
         'ROM adjoint requested before compute_objectives')
    if (this%forward_constrained_enriched) then
       call rom_problem_run_enriched_backward(this)
       return
    end if
    call this%adjoint%free(); call this%lambda_sum%free(); call this%correlation%free()
    call this%adjoint%init(this%states%get_nrows(), this%n_steps + 1)
    call this%lambda_sum%init(this%states%get_nrows())
    call this%correlation%init(this%states%get_nrows(), this%states%get_nrows())
    timer_start = real(MPI_WTIME(),rp)
    call this%rom%discrete_adjoint_bdf1(this%states, this%dt, &
         this%running_gradient, this%terminal_gradient, this%adjoint, this%lambda_history)
    this%exact_adjoint_seconds = this%exact_adjoint_seconds + &
         real(MPI_WTIME(),rp) - timer_start
    timer_start = real(MPI_WTIME(),rp)
    this%lambda_sum%x = 0.0_rp
    this%correlation%x = 0.0_rp
    do step = 1, this%n_steps
       if (step == 1) then
          call neko_scratch_registry%request_vector(lambda, lambda_index, &
               this%states%get_nrows(), .true.)
          call neko_scratch_registry%request_vector(state, state_index, &
               this%states%get_nrows(), .true.)
       end if
       ! Reuse the transpose-solve vector accumulated during the backward
       ! recursion; this is the coefficient in the implicit design derivative.
       state%x = this%states%x(:,step+1)
       lambda%x = this%lambda_history%x(:,step)
       this%lambda_sum%x = this%lambda_sum%x + lambda%x
       call add_outer_product(this%correlation, lambda, state, 1.0_rp)
    end do
    this%temporal_correlation_seconds = this%temporal_correlation_seconds + &
         real(MPI_WTIME(),rp) - timer_start
    this%backward_ready = .true.
    if (this%write_trajectory_fields .and. this%rom_outputs_initialized) &
         call rom_problem_write_adjoint_trajectory(this)
    call neko_scratch_registry%relinquish_vector([lambda_index, state_index])
    nullify(lambda); nullify(state)
  end subroutine rom_problem_run_backward_unsteady

  !> Run both enriched and exact reference adjoints on the stored trajectory.
  !! `adjoint` is the independently restricted adjoint-POD solve used when
  !! `sensitivity_source=projected_adjoint`; `exact_adjoint` remains the
  !! exact discrete adjoint of the forward Galerkin ROM for diagnostics and
  !! the selectable `exact_galerkin` sensitivity source.
  subroutine rom_problem_run_enriched_backward(this)
    class(rom_problem_t), intent(inout) :: this
    type(matrix_t), allocatable :: convection_jacobian(:)
    type(matrix_t) :: projected_lambda_history
    type(vector_t), pointer :: lambda, state
    type(vector_t), pointer :: exact_lambda
    integer :: lambda_index, state_index, exact_lambda_index
    integer :: step, e, r, k
    logical :: need_exact_adjoint, need_projected_adjoint
    logical :: need_projected_correlations
    logical :: can_reuse_tensor
    real(kind=rp) :: timer_start

    e = this%enriched_adjoint_rank
    r = this%states%get_nrows()
    call this%adjoint%free(); call this%exact_adjoint%free()
    call this%enriched_lambda_sum%free(); call this%enriched_correlation%free()
    call this%lambda_sum%free(); call this%correlation%free()
    need_exact_adjoint = this%diagnostic_audit .or. &
         this%sensitivity_source == 'exact_algebraic' .or. &
         this%sensitivity_source == 'exact_galerkin'
    need_projected_adjoint = this%diagnostic_audit .or. &
         this%sensitivity_source == 'projected_adjoint' .or. &
         this%write_trajectory_fields .or. this%enable_output
    need_projected_correlations = this%diagnostic_audit .or. &
         this%sensitivity_source == 'projected_adjoint'
    this%exact_adjoint_seconds = 0.0_rp
    this%projected_tensor_assembly_seconds = 0.0_rp
    this%projected_tensor_cache_hit = 0.0_rp
    this%projected_adjoint_integration_seconds = 0.0_rp
    this%temporal_correlation_seconds = 0.0_rp
    if (need_exact_adjoint) then
       call this%exact_adjoint%init(r, this%n_steps + 1)
       timer_start = real(MPI_WTIME(),rp)
       call this%rom%discrete_adjoint_bdf1(this%states, this%dt, &
            this%running_gradient, this%terminal_gradient, this%exact_adjoint, &
            this%lambda_history)
       this%exact_adjoint_seconds = real(MPI_WTIME(),rp) - timer_start
       call this%lambda_sum%free(); call this%correlation%free()
       call this%lambda_sum%init(r)
       call this%correlation%init(r, r)
    end if

    if (need_projected_adjoint) then
       call this%adjoint%init(e, this%n_steps + 1)
       if (need_projected_correlations) then
          call this%enriched_lambda_sum%init(e)
          call this%enriched_correlation%init(e, r)
       end if

       ! Assemble each dynamic Jacobian slice C_k directly. The helper forms
       ! every (forward mode k, adjoint perturbation j) field once and reduces
       ! the complete packed tensor. Lift-linearized terms remain in the
       ! fixed E-space linear operator.
       allocate(convection_jacobian(r))
       do k = 1, r
          call convection_jacobian(k)%init(e,e)
       end do
       can_reuse_tensor = this%projected_tensor_cache_valid .and. &
            .not. this%basis_transform_adaptive .and. &
            NEKO_BCKND_DEVICE .eq. 0 .and. &
            this%projected_tensor_cache_forward_rank == r .and. &
            this%projected_tensor_cache_adjoint_rank == e
       if (can_reuse_tensor) then
          if (allocated(this%projected_convection_cache)) then
             can_reuse_tensor = size(this%projected_convection_cache) == r
          else
             can_reuse_tensor = .false.
          end if
       end if
       if (can_reuse_tensor) then
          do k = 1, r
             convection_jacobian(k)%x = this%projected_convection_cache(k)%x
          end do
          this%projected_tensor_cache_hit = 1.0_rp
       else
          this%projected_tensor_cache_valid = .false.
          timer_start = real(MPI_WTIME(),rp)
          if (this%has_scalar) then
             call rom_assemble_linearized_convection_tensor( &
                  this%weights, this%velocity_modes, this%velocity_gradients, &
                  this%scalar_modes, this%scalar_gradients, this%test_velocity_modes, &
                  this%test_velocity_gradients, this%test_scalar_modes, &
                  this%test_scalar_gradients, this%scalar_metric_weight, &
                  convection_jacobian)
          else
             call rom_assemble_linearized_convection_tensor( &
                  this%weights, this%velocity_modes, this%velocity_gradients, &
                  adjoint_velocity_modes=this%test_velocity_modes, &
                  adjoint_velocity_gradients=this%test_velocity_gradients, &
                  scalar_weight=this%scalar_metric_weight, tensor=convection_jacobian)
          end if
          this%projected_tensor_assembly_seconds = real(MPI_WTIME(),rp) - timer_start
          if (.not. this%basis_transform_adaptive .and. NEKO_BCKND_DEVICE .eq. 0) then
             if (allocated(this%projected_convection_cache)) then
                do k = 1, size(this%projected_convection_cache)
                   call this%projected_convection_cache(k)%free()
                end do
                deallocate(this%projected_convection_cache)
             end if
             allocate(this%projected_convection_cache(r))
             do k = 1, r
                call this%projected_convection_cache(k)%init(e,e)
                this%projected_convection_cache(k)%x = convection_jacobian(k)%x
             end do
             this%projected_tensor_cache_forward_rank = r
             this%projected_tensor_cache_adjoint_rank = e
             this%projected_tensor_cache_valid = .true.
          end if
       end if

       ! Integrate without retaining Nt dense Jacobians. Keep residual
       ! multipliers only when projected sensitivity or audit needs them.
       timer_start = real(MPI_WTIME(),rp)
       if (need_projected_correlations) then
          call this%enriched_adjoint_model%discrete_adjoint_bdf1_tensor( &
               this%states, this%dt, this%enriched_running_gradient, &
               this%enriched_terminal_gradient, this%enriched_fixed_linear, &
               convection_jacobian, this%adjoint, projected_lambda_history)
       else
          call this%enriched_adjoint_model%discrete_adjoint_bdf1_tensor( &
               this%states, this%dt, this%enriched_running_gradient, &
               this%enriched_terminal_gradient, this%enriched_fixed_linear, &
               convection_jacobian, this%adjoint)
       end if
       this%projected_adjoint_integration_seconds = real(MPI_WTIME(),rp) - timer_start
    end if

    if (need_projected_correlations .or. need_exact_adjoint) then
       timer_start = real(MPI_WTIME(),rp)
       if (need_projected_correlations) then
          this%enriched_lambda_sum%x = 0.0_rp
          this%enriched_correlation%x = 0.0_rp
       end if
       if (this%n_steps > 0) then
          if (need_projected_correlations) &
               call neko_scratch_registry%request_vector(lambda, lambda_index, e, .true.)
          call neko_scratch_registry%request_vector(state, state_index, r, .true.)
          if (need_exact_adjoint) call neko_scratch_registry%request_vector( &
               exact_lambda, exact_lambda_index, r, .true.)
       end if
       do step = 1, this%n_steps
          state%x = this%states%x(:,step+1)
          if (need_projected_correlations) then
             lambda%x = projected_lambda_history%x(:,step)
             this%enriched_lambda_sum%x = this%enriched_lambda_sum%x + lambda%x
             call add_outer_product(this%enriched_correlation, lambda, state, 1.0_rp)
          end if
          if (need_exact_adjoint) then
             exact_lambda%x = this%lambda_history%x(:,step)
             this%lambda_sum%x = this%lambda_sum%x + exact_lambda%x
             call add_outer_product(this%correlation, exact_lambda, state, 1.0_rp)
          end if
       end do
       if (this%n_steps > 0) then
          if (need_projected_correlations) &
               call neko_scratch_registry%relinquish_vector(lambda_index)
          if (need_exact_adjoint) then
             call neko_scratch_registry%relinquish_vector([state_index, exact_lambda_index])
          else
             call neko_scratch_registry%relinquish_vector(state_index)
          end if
       end if
       nullify(lambda); nullify(state); nullify(exact_lambda)
       this%temporal_correlation_seconds = real(MPI_WTIME(),rp) - timer_start
    end if
    this%backward_ready = .true.
    if (this%write_trajectory_fields .and. this%rom_outputs_initialized) &
         call rom_problem_write_adjoint_trajectory(this)
    if (allocated(convection_jacobian)) then
       do k = 1, size(convection_jacobian)
          call convection_jacobian(k)%free()
       end do
       deallocate(convection_jacobian)
    end if
    call projected_lambda_history%free()
  end subroutine rom_problem_run_enriched_backward

  !> Assemble the mapped full-dimensional ROM design sensitivity density.
  subroutine rom_problem_compute_density(this, density, source)
    class(rom_problem_t), intent(inout) :: this
    type(vector_t), intent(inout) :: density
    character(len=*), intent(in), optional :: source
    type(vector_t), pointer :: lift_square, direct, masked_direct, product, component_sum
    type(vector_ptr_t) :: mode_sum(3)
    integer :: lift_square_index, direct_index, masked_direct_index, product_index
    integer :: component_sum_index, mode_sum_index(3)
    integer :: i, j, c, r, e, n
    logical :: use_projected
    character(len=32) :: source_name

    source_name = this%sensitivity_source
    if (present(source)) source_name = trim(source)
    if (source_name /= 'exact_galerkin' .and. source_name /= 'projected_adjoint' .and. &
         source_name /= 'lagged_basis' .and. source_name /= 'exact_algebraic') &
         call neko_error('unknown ROM sensitivity source')
    if (source_name == 'exact_algebraic' .and. &
         (.not. this%experimental_algebraic .or. .not. this%basis_transform_adaptive)) &
         call neko_error('exact_algebraic sensitivity requires algebraic adaptive modes')
    use_projected = this%forward_constrained_enriched .and. &
         source_name == 'projected_adjoint'
    if (source_name == 'projected_adjoint' .and. .not. this%forward_constrained_enriched) &
         call neko_error('projected_adjoint requires independent forward and adjoint spaces')
    if (source_name == 'lagged_basis' .and. this%forward_constrained_enriched) &
         call neko_error('lagged_basis is unavailable for independently projected adjoints')
    if (this%experimental_algebraic .and. source_name == 'exact_galerkin') &
         call neko_error('adaptive Galerkin sensitivities must select exact_algebraic or lagged_basis')
    r = this%mass%get_nrows(); n = this%weights%size()
    e = this%enriched_adjoint_rank
    if (use_projected) then
       if (this%enriched_lambda_sum%size() /= e .or. &
            this%enriched_correlation%get_nrows() /= e .or. &
            this%enriched_correlation%get_ncols() /= r) call neko_error( &
            'projected-adjoint ROM sensitivity has inconsistent reduced shapes')
    end if
    if (density%size() /= n) call density%init(n)
       call vector_cmult(density, 0.0_rp)
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
    if (this%has_brinkman) then
       call vector_cmult2(direct, lift_square, this%cost_total)
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
    ! The state-dependent design term is differentiated through the selected
    ! reduced adjoint.  The exact branch uses the adjoint of the restricted
    ! forward Galerkin map.  The projected branch contracts the independently
    ! integrated V_a adjoint against the V_a--V_f Brinkman cross terms.
    if (use_projected) then
       do i = 1, e
          component_sum = 0.0_rp
          do c = 1, 3
             call vector_col3(product, this%test_velocity_modes(c,i), &
                  this%lift_velocity(c))
             call vector_add2(component_sum, product)
          end do
          call vector_add2s2(density, component_sum, -this%dt * &
               this%enriched_lambda_sum%x(i))
          do j = 1, r
             component_sum = 0.0_rp
             do c = 1, 3
                call vector_col3(product, this%test_velocity_modes(c,i), &
                     this%velocity_modes(c,j))
                call vector_add2(component_sum, product)
             end do
             call vector_add2s2(density, component_sum, -this%dt * &
                  this%enriched_correlation%x(i,j))
          end do
       end do
    else
       do j = 1, r
          component_sum = 0.0_rp
          do c = 1, 3
             if (this%forward_constrained_enriched) then
                call vector_col3(product, this%velocity_modes(c,j), &
                     this%lift_velocity(c))
             else
                call vector_col3(product, this%test_velocity_modes(c,j), &
                     this%lift_velocity(c))
             end if
             call vector_add2(component_sum, product)
          end do
          call vector_add2s2(density, component_sum, -this%dt * this%lambda_sum%x(j))
          do c = 1, r
             component_sum = 0.0_rp
             do n = 1, 3
                if (this%forward_constrained_enriched) then
                   call vector_col3(product, this%velocity_modes(n,j), &
                        this%velocity_modes(n,c))
                else
                   call vector_col3(product, this%test_velocity_modes(n,j), &
                        this%velocity_modes(n,c))
                end if
                call vector_add2(component_sum, product)
             end do
             call vector_add2s2(density, component_sum, -this%dt * &
                  this%correlation%x(j,c))
          end do
       end do
    end if
    if (this%basis_transform_adaptive .and. &
         (source_name == 'exact_galerkin' .or. source_name == 'exact_algebraic') .and. &
         abs(this%basis_transform_tau) > &
         10.0_rp * epsilon(1.0_rp)) then
       call rom_problem_add_adaptive_mode_density(this, density)
    end if
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

  !> Add the exact frozen-basis mode-motion contribution to the mapped design
  !! density for the design-dependent Galerkin path.
  !!
  !! The reduced adjoint supplies the cotangents of the BDF1 residual.  This
  !! routine reverse-assembles those cotangents first with respect to the
  !! physical velocity modes and their SEM gradients, then applies the selected
  !! transformation's exact discrete transpose after the native weighted
  !! derivative transpose.
  !!
  !! The transformed POD factors are frozen.  No SVD or reorthogonalisation is
  !! differentiated here. In mixed basis selections this uses the exact
  !! forward-ROM adjoint, independently of the projected-adjoint search path.
  subroutine rom_problem_add_adaptive_mode_density(this, density)
    class(rom_problem_t), intent(inout) :: this
    type(vector_t), intent(inout) :: density
    type(vector_t), allocatable :: mode_value(:, :), mode_gradient(:, :, :)
    type(vector_t), allocatable :: initial_hom(:)
    real(kind=rp), allocatable :: mass_moment(:, :), lambda_moment(:), old_moment(:, :)
    real(kind=rp), allocatable :: quadratic_moment(:, :, :), new_moment(:, :)
    real(kind=rp), allocatable :: viscous_first(:), viscous_second(:, :)
    real(kind=rp), allocatable :: brinkman_first(:), brinkman_second(:, :)
    real(kind=rp), allocatable :: velocity_panel(:, :, :), scalar_panel(:, :)
    real(kind=rp), allocatable :: test_panel(:, :, :), advect_panel(:, :, :)
    real(kind=rp), allocatable :: scalar_contract_panel(:, :)
    real(kind=rp), allocatable :: k_panel(:, :), kt_panel(:, :)
    type(vector_t), pointer :: lambda, state_old, state_new, delta_state
    type(vector_t), pointer :: adjoint_initial, initial_cotangent
    type(vector_t), pointer :: product, product_two, current_chi, candidate_D, inv_D2
    type(vector_t), pointer :: total, projected_cotangent, transform_factor
    type(vector_t), pointer :: factor_cotangent
    type(vector_t) :: correction_cotangent(3)
    type(field_t), pointer :: field_work, transpose_work, cot_x, cot_y, cot_z
    type(field_t), pointer :: chi_field
    type(fluid_pnpn_t), pointer :: fluid
    integer :: n, r, c, d, i, j, k, p, pp, nb, panel_start, step
    integer :: lambda_index, state_old_index, state_new_index
    integer :: delta_state_index, adjoint_initial_index, initial_cotangent_index
    integer :: product_index, product_two_index, current_chi_index, candidate_D_index
    integer :: inv_D2_index, total_index, projected_cotangent_index
    integer :: transform_factor_index, field_work_index, transpose_work_index
    integer :: factor_cotangent_index
    integer :: cot_x_index, cot_y_index, cot_z_index
    real(kind=rp) :: nu, dt_step, timer_start, scalar_contract
    real(kind=rp) :: viscous_volume, brinkman_volume, viscous_time_weight
    real(kind=rp) :: brinkman_time_weight, wv, wb
    real(kind=rp) :: projection_residual, projection_linear_residual
    integer :: projection_iterations
    logical :: projection_converged
    logical :: cpu_contractions
    integer, parameter :: contraction_block = 256

    if (.not. this%basis_transform_adaptive) return
    if (this%petrov_galerkin) &
         call neko_error('exact adaptive pullback requires an equal-space Galerkin forward model')
    if (.not. this%reference_modes_saved) call neko_error( &
         'adaptive mode pullback has no immutable reference modes')
    if (.not. associated(this%simulation)) call neko_error( &
         'adaptive mode pullback has no bound simulation')

    select type (fluid_in => this%simulation%fluid)
    type is (fluid_pnpn_t)
       fluid => fluid_in
    class default
       call neko_error('adaptive mode pullback requires a Pn/Pn fluid')
    end select
    n = this%weights%size()
    r = this%mass%get_nrows()
    ! mxm_wrapper dispatches through LIBXSMM when Neko was configured with the
    ! XSMM backend. Use the packed BLAS path only when that backend is absent;
    ! the native vector path below remains the safe fallback otherwise.
    cpu_contractions = (NEKO_BCKND_DEVICE .eq. 0) .and. &
         (NEKO_BCKND_XSMM .eq. 0) .and. (rp == sp .or. rp == dp)
    if (size(this%velocity_modes, 2) /= r .or. &
         size(this%reference_velocity_modes, 2) /= r .or. &
         (size(this%test_velocity_modes,2) /= r .and. &
         .not. this%forward_constrained_enriched)) &
         call neko_error('adaptive mode pullback rank mismatch')

    nu = this%viscosity / this%density
    allocate(mode_value(3,r), mode_gradient(3,r,3), initial_hom(3))
    allocate(mass_moment(r,r), lambda_moment(r), old_moment(r,r), &
         quadratic_moment(r,r,r), new_moment(r,r), viscous_first(r), &
         viscous_second(r,r), brinkman_first(r), brinkman_second(r,r))
    mass_moment = 0.0_rp
    lambda_moment = 0.0_rp
    old_moment = 0.0_rp
    quadratic_moment = 0.0_rp
    new_moment = 0.0_rp
    viscous_first = 0.0_rp
    viscous_second = 0.0_rp
    brinkman_first = 0.0_rp
    brinkman_second = 0.0_rp
    do i = 1, r
       do c = 1, 3
          call mode_value(c,i)%init(n); mode_value(c,i) = 0.0_rp
          do d = 1, 3
             call mode_gradient(c,i,d)%init(n)
             mode_gradient(c,i,d) = 0.0_rp
          end do
       end do
    end do
    do c = 1, 3; call initial_hom(c)%init(n); end do

    call neko_scratch_registry%request_vector(lambda, lambda_index, r, .true.)
    call neko_scratch_registry%request_vector(state_old, state_old_index, r, .true.)
    call neko_scratch_registry%request_vector(state_new, state_new_index, r, .true.)
    call neko_scratch_registry%request_vector(delta_state, delta_state_index, r, .true.)
    call neko_scratch_registry%request_vector(adjoint_initial, adjoint_initial_index, r, .true.)
    call neko_scratch_registry%request_vector(initial_cotangent, initial_cotangent_index, r, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    call neko_scratch_registry%request_vector(product_two, product_two_index, n, .true.)
    call neko_scratch_registry%request_vector(current_chi, current_chi_index, n, .true.)
    call neko_scratch_registry%request_vector(candidate_D, candidate_D_index, n, .true.)
    call neko_scratch_registry%request_vector(inv_D2, inv_D2_index, n, .true.)
    call neko_scratch_registry%request_vector(total, total_index, n, .true.)
    call neko_scratch_registry%request_vector(projected_cotangent, projected_cotangent_index, n, .true.)
    call neko_scratch_registry%request_vector(transform_factor, transform_factor_index, n, .true.)
    if (this%basis_transform == 'vector_potential') &
         call neko_scratch_registry%request_vector(factor_cotangent, &
         factor_cotangent_index, n, .true.)
    if (this%basis_transform == 'vector_potential') then
       do c = 1, 3
          call correction_cotangent(c)%init(n)
       end do
    end if
    call neko_scratch_registry%request_field(field_work, field_work_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(transpose_work, transpose_work_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(cot_x, cot_x_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(cot_y, cot_y_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(cot_z, cot_z_index, fluid%dm_Xh, .true.)

    chi_field => neko_registry%get_field('brinkman_amplitude')
    call field_to_vector(current_chi, chi_field)
    call vector_cmult2(candidate_D, current_chi, this%basis_transform_tau)
    call vector_cadd(candidate_D, 1.0_rp)
    inv_D2 = 1.0_rp
    product = 1.0_rp
    call vector_col2(product, candidate_D)
    call vector_col2(product, candidate_D)
    call vector_invcol2(inv_D2, product)
    call vector_cmult(inv_D2, this%basis_transform_tau)
    if (this%basis_transform == 'vector_potential') then
       call vector_copy(transform_factor, this%reference_D)
       call vector_invcol2(transform_factor, candidate_D)
       call vector_cadd(transform_factor, -1.0_rp)
    end if
    nullify(chi_field)

    ! The canonical initial condition is the same Neko field used for the
    ! reference reduced projection.  Only its homogeneous velocity part enters
    ! the velocity-mode pullback.
    do c = 1, 3
       select case (c)
       case (1); call field_copy(field_work, fluid%u)
       case (2); call field_copy(field_work, fluid%v)
       case (3); call field_copy(field_work, fluid%w)
       end select
       call field_add2s2(field_work, this%lift%velocity(c), -1.0_rp, n)
       call field_to_vector(initial_hom(c), field_work)
    end do
    ! With an independently selected adjoint basis, use the separately stored
    ! exact forward-ROM adjoint. Equal-space FF and E:E run that exact adjoint
    ! in `adjoint` directly.
    if (this%forward_constrained_enriched) then
       adjoint_initial%x = this%exact_adjoint%x(:,1)
    else
       adjoint_initial%x = this%adjoint%x(:,1)
    end if
    call this%rom%solve_mass_transpose(adjoint_initial, initial_cotangent)
    do i = 1, r
       do c = 1, 3
          call vector_add2s2(mode_value(c,i), initial_hom(c), initial_cotangent%x(i))
          do j = 1, r
             call vector_add2s2(mode_value(c,i), this%velocity_modes(c,j), &
                  -initial_cotangent%x(i) * this%initial_state%x(j))
             call vector_add2s2(mode_value(c,i), this%velocity_modes(c,j), &
                  -this%initial_state%x(i) * initial_cotangent%x(j))
          end do
       end do
    end do

    viscous_volume = window_measure(this, this%viscous_start, this%viscous_end)
    brinkman_volume = window_measure(this, this%brinkman_start, this%brinkman_end)
    ! Accumulate only reduced temporal moments in this sweep.  The adjoint and
    ! ROM states are reduced vectors, so no mesh-sized vector operation belongs
    ! in this loop.  The spatial cotangents below apply these moments once.
    do step = 1, this%n_steps
       state_old%x = this%states%x(:,step)
       state_new%x = this%states%x(:,step+1)
       delta_state%x = state_new%x - state_old%x
       lambda%x = this%lambda_history%x(:,step)

       do i = 1, r
          lambda_moment(i) = lambda_moment(i) + lambda%x(i)
          do j = 1, r
             mass_moment(i,j) = mass_moment(i,j) + lambda%x(i) * delta_state%x(j)
             old_moment(i,j) = old_moment(i,j) + lambda%x(i) * state_old%x(j)
             new_moment(i,j) = new_moment(i,j) + lambda%x(i) * state_new%x(j)
             do k = 1, r
                quadratic_moment(i,j,k) = quadratic_moment(i,j,k) + &
                     lambda%x(i) * state_old%x(j) * state_old%x(k)
             end do
          end do
       end do

       dt_step = this%time_start + real(step, rp) * this%dt
       viscous_time_weight = 0.0_rp
       brinkman_time_weight = 0.0_rp
       if (this%has_viscous) viscous_time_weight = objective_time_weight( &
            this%viscous_start, this%viscous_end, dt_step, viscous_volume, this%dt)
       if (this%has_brinkman) brinkman_time_weight = objective_time_weight( &
            this%brinkman_start, this%brinkman_end, dt_step, brinkman_volume, this%dt)
       wv = this%dt * viscous_time_weight * this%viscous_weight * this%viscosity / &
            max(this%viscous_domain, tiny(1.0_rp))
       wb = this%dt * brinkman_time_weight * this%brinkman_weight / &
            max(this%brinkman_domain, tiny(1.0_rp))
       if (this%has_viscous) then
          do i = 1, r
             viscous_first(i) = viscous_first(i) + wv * state_new%x(i)
             do j = 1, r
                viscous_second(i,j) = viscous_second(i,j) + &
                     wv * state_new%x(i) * state_new%x(j)
             end do
          end do
       end if
       if (this%has_brinkman) then
          do i = 1, r
             brinkman_first(i) = brinkman_first(i) + wb * state_new%x(i)
             do j = 1, r
                brinkman_second(i,j) = brinkman_second(i,j) + &
                     wb * state_new%x(i) * state_new%x(j)
             end do
          end do
       end if
    end do

    ! BDF mass matrix dependence: both test and trial mode occurrences.
    do i = 1, r
       do j = 1, r
          do c = 1, 3
             call vector_add2s2(mode_value(c,i), this%velocity_modes(c,j), -mass_moment(i,j))
             call vector_add2s2(mode_value(c,j), this%velocity_modes(c,i), -mass_moment(i,j))
          end do
       end do
    end do

    ! Explicit velocity lift convection and diffusion.
    do i = 1, r
       do c = 1, 3
          do d = 1, 3
             call rom_problem_add_product(mode_value(c,i), this%lift_velocity(d), &
                  this%lift_gradients(c,d), -this%dt * lambda_moment(i), product)
             call vector_add2s2(mode_gradient(c,i,d), this%lift_gradients(c,d), &
                  -this%dt * lambda_moment(i) * nu)
          end do
       end do
    end do
    do i = 1, r
       do j = 1, r
          do c = 1, 3
             do d = 1, 3
                call rom_problem_add_product(mode_value(c,i), this%lift_velocity(d), &
                     this%velocity_gradients(c,j,d), -this%dt * old_moment(i,j), product)
                call rom_problem_add_product(mode_value(c,i), this%velocity_modes(d,j), &
                     this%lift_gradients(c,d), -this%dt * old_moment(i,j), product)
                call rom_problem_add_product(mode_value(d,j), this%velocity_modes(c,i), &
                     this%lift_gradients(c,d), -this%dt * old_moment(i,j), product)
                call rom_problem_add_product(mode_gradient(c,j,d), this%velocity_modes(c,i), &
                     this%lift_velocity(d), -this%dt * old_moment(i,j), product)
             end do
          end do
       end do
    end do

    ! Quadratic velocity convection, contracted over time before spatial work.
    ! On CPU, fuse the three exact mode cotangents into one node-local reduced
    ! contraction. This avoids issuing three mesh-vector product/update pairs
    ! for every (i,j,k,c,d), while retaining all test, advecting and transported
    ! mode terms. The device path keeps the existing backend-aware vector ops.
    timer_start = real(MPI_WTIME(),rp)
    if (cpu_contractions) then
       allocate(velocity_panel(r,contraction_block,3), &
            test_panel(r,contraction_block,3), advect_panel(r,contraction_block,3), &
            k_panel(r,r), kt_panel(r,r))
       do panel_start = 1, n, contraction_block
          nb = min(contraction_block, n-panel_start+1)
          do d = 1, 3
             do j = 1, r
                do pp = 1, nb
                   p = panel_start + pp - 1
                   velocity_panel(j,pp,d) = this%velocity_modes(d,j)%x(p)
                end do
             end do
          end do
          do k = 1, r
             k_panel = quadratic_moment(:,:,k)
             kt_panel = transpose(k_panel)
             do d = 1, 3
                call mxm(k_panel, r, velocity_panel(:,:,d), r, test_panel(:,:,d), nb)
             end do
             do c = 1, 3
                call mxm(kt_panel, r, velocity_panel(:,:,c), r, advect_panel(:,:,c), nb)
             end do
             do c = 1, 3
                do d = 1, 3
                   ! K_k V_d and K_k^T V_c are shared across the other
                   ! component index, so each bounded panel GEMM is formed once.
                   do pp = 1, nb
                      p = panel_start + pp - 1
                      do i = 1, r
                         mode_value(c,i)%x(p) = mode_value(c,i)%x(p) - this%dt * &
                              this%velocity_gradients(c,k,d)%x(p) * test_panel(i,pp,d)
                      end do
                      do j = 1, r
                         mode_value(d,j)%x(p) = mode_value(d,j)%x(p) - this%dt * &
                              this%velocity_gradients(c,k,d)%x(p) * advect_panel(j,pp,c)
                      end do
                      scalar_contract = 0.0_rp
                      do j = 1, r
                         scalar_contract = scalar_contract + &
                              velocity_panel(j,pp,d) * advect_panel(j,pp,c)
                      end do
                      mode_gradient(c,k,d)%x(p) = mode_gradient(c,k,d)%x(p) - &
                           this%dt * scalar_contract
                   end do
                end do
             end do
          end do
       end do
       deallocate(velocity_panel, test_panel, advect_panel, k_panel, kt_panel)
    else
       do i = 1, r
          do j = 1, r
             do k = 1, r
                do c = 1, 3
                   do d = 1, 3
                      call rom_problem_add_product(mode_value(c,i), this%velocity_modes(d,j), &
                           this%velocity_gradients(c,k,d), -this%dt * quadratic_moment(i,j,k), product)
                      call rom_problem_add_product(mode_value(d,j), this%velocity_modes(c,i), &
                           this%velocity_gradients(c,k,d), -this%dt * quadratic_moment(i,j,k), product)
                      call rom_problem_add_product(mode_gradient(c,k,d), this%velocity_modes(c,i), &
                           this%velocity_modes(d,j), -this%dt * quadratic_moment(i,j,k), product)
                   end do
                end do
             end do
          end do
       end do
    end if
    this%velocity_mode_cotangent_seconds = this%velocity_mode_cotangent_seconds + &
         real(MPI_WTIME(),rp) - timer_start

    ! Implicit viscous and Brinkman terms at the new BDF1 state.
    do i = 1, r
       do j = 1, r
          do c = 1, 3
             do d = 1, 3
                call vector_add2s2(mode_gradient(c,i,d), this%velocity_gradients(c,j,d), &
                     -this%dt * nu * new_moment(i,j))
                call vector_add2s2(mode_gradient(c,j,d), this%velocity_gradients(c,i,d), &
                     -this%dt * nu * new_moment(i,j))
             end do
             call vector_col3(product, current_chi, this%velocity_modes(c,j))
             call vector_add2s2(mode_value(c,i), product, -this%dt * new_moment(i,j))
             call vector_col3(product, current_chi, this%velocity_modes(c,i))
             call vector_add2s2(mode_value(c,j), product, -this%dt * new_moment(i,j))
          end do
       end do
    end do
    do i = 1, r
       do c = 1, 3
          call vector_col3(product, current_chi, this%lift_velocity(c))
          call vector_add2s2(mode_value(c,i), product, -this%dt * lambda_moment(i))
       end do
    end do

    ! Windowed objective moments, applied once to each spatial mode.
    do i = 1, r
       do c = 1, 3
          do d = 1, 3
             call vector_col3(product, this%viscous_mask, this%lift_gradients(c,d))
             call vector_add2s2(mode_gradient(c,i,d), product, viscous_first(i))
             do j = 1, r
                call vector_col3(product, this%viscous_mask, this%velocity_gradients(c,j,d))
                call vector_add2s2(mode_gradient(c,i,d), product, viscous_second(i,j))
             end do
          end do
          call vector_col3(product, current_chi, this%lift_velocity(c))
          call vector_col3(product_two, this%brinkman_mask, product)
          call vector_add2s2(mode_value(c,i), product_two, brinkman_first(i))
          do j = 1, r
             call vector_col3(product, current_chi, this%velocity_modes(c,j))
             call vector_col3(product_two, this%brinkman_mask, product)
             call vector_add2s2(mode_value(c,i), product_two, brinkman_second(i,j))
          end do
       end do
    end do

    ! Scalar advection depends on velocity modes but not their gradients.
    if (this%has_scalar) then
       timer_start = real(MPI_WTIME(),rp)
       do i = 1, r
          do j = 1, r
             do d = 1, 3
                call rom_problem_add_product(mode_value(d,j), this%scalar_modes(i), &
                     this%scalar_lift_gradients(d), -this%dt * this%scalar_metric_weight * &
                     old_moment(i,j), product)
             end do
          end do
       end do
       if (cpu_contractions) then
          allocate(scalar_panel(r,contraction_block), &
               scalar_contract_panel(r,contraction_block), k_panel(r,r), kt_panel(r,r))
          do panel_start = 1, n, contraction_block
             nb = min(contraction_block, n-panel_start+1)
             do i = 1, r
                do pp = 1, nb
                   scalar_panel(i,pp) = this%scalar_modes(i)%x(panel_start+pp-1)
                end do
             end do
             do j = 1, r
                do i = 1, r
                   do k = 1, r
                      k_panel(i,k) = quadratic_moment(i,j,k)
                   end do
                end do
                kt_panel = transpose(k_panel)
                call mxm(kt_panel, r, scalar_panel, r, scalar_contract_panel, nb)
                do d = 1, 3
                   do pp = 1, nb
                      p = panel_start + pp - 1
                      scalar_contract = 0.0_rp
                      do k = 1, r
                         scalar_contract = scalar_contract + scalar_contract_panel(k,pp) * &
                              this%scalar_gradients(k,d)%x(p)
                      end do
                      mode_value(d,j)%x(p) = mode_value(d,j)%x(p) - this%dt * &
                           this%scalar_metric_weight * scalar_contract
                   end do
                end do
             end do
          end do
          deallocate(scalar_panel, scalar_contract_panel, k_panel, kt_panel)
       else
          do i = 1, r
             do j = 1, r
                do d = 1, 3
                   do k = 1, r
                      call rom_problem_add_product(mode_value(d,j), this%scalar_modes(i), &
                           this%scalar_gradients(k,d), -this%dt * this%scalar_metric_weight * &
                           quadratic_moment(i,j,k), product)
                   end do
                end do
             end do
          end do
       end if
       this%scalar_mode_cotangent_seconds = this%scalar_mode_cotangent_seconds + &
            real(MPI_WTIME(),rp) - timer_start
    end if

    ! Reverse the SEM gradient.  cdtp is Neko's device-aware transpose of the
    ! elemental derivative and includes the SEM quadrature weight used by the
    ! reduced contractions.
    do i = 1, r
       timer_start = real(MPI_WTIME(),rp)
       do c = 1, 3
          total = 0.0_rp
          call vector_col3(product, this%weights, mode_value(c,i))
          call vector_copy(total, product)
          do d = 1, 3
             call vector_to_field(field_work, mode_gradient(c,i,d))
             if (d == 1) then
                call cdtp(transpose_work%x, field_work%x, fluid%c_Xh%drdx, &
                     fluid%c_Xh%dsdx, fluid%c_Xh%dtdx, fluid%c_Xh)
             else if (d == 2) then
                call cdtp(transpose_work%x, field_work%x, fluid%c_Xh%drdy, &
                     fluid%c_Xh%dsdy, fluid%c_Xh%dtdy, fluid%c_Xh)
             else
                call cdtp(transpose_work%x, field_work%x, fluid%c_Xh%drdz, &
                     fluid%c_Xh%dsdz, fluid%c_Xh%dtdz, fluid%c_Xh)
             end if
             call field_to_vector(product, transpose_work)
             call vector_add2(total, product)
          end do
          select case (c)
          case (1); call vector_to_field(cot_x, total)
          case (2); call vector_to_field(cot_y, total)
          case (3); call vector_to_field(cot_z, total)
          end select
       end do
       this%derivative_transpose_seconds = this%derivative_transpose_seconds + &
            real(MPI_WTIME(),rp) - timer_start
       timer_start = real(MPI_WTIME(),rp)
       if (this%basis_transform == 'projected') then
          call rom_project_velocity_cotangent(fluid, cot_x, cot_y, cot_z, &
               projection_residual, linear_residual = projection_linear_residual, &
               iterations = projection_iterations, converged = projection_converged)
          if (projection_residual > this%admissibility_tolerance .or. &
               .not. projection_converged) then
             call neko_error('adaptive mode pullback projection failed its continuity gate')
          end if
       end if
       if (this%basis_transform == 'vector_potential') then
          call field_to_vector(correction_cotangent(1), cot_x)
          call field_to_vector(correction_cotangent(2), cot_y)
          call field_to_vector(correction_cotangent(3), cot_z)
          factor_cotangent = 0.0_rp
          call this%reference_forward_potential%pullback( &
               this%reference_potential_forward_map(i), transform_factor, &
               correction_cotangent, factor_cotangent)
          call vector_col3(total, factor_cotangent, inv_D2)
          call vector_col3(total, total, this%reference_D)
          call vector_invcol2(total, this%weights)
          call vector_add2s2(density, total, -1.0_rp)
       else
          call vector_cmult2(transform_factor, this%reference_D, 1.0_rp)
          do c = 1, 3
             select case (c)
             case (1); call field_to_vector(projected_cotangent, cot_x)
             case (2); call field_to_vector(projected_cotangent, cot_y)
             case (3); call field_to_vector(projected_cotangent, cot_z)
             end select
             call vector_col3(product, transform_factor, this%reference_velocity_modes(c,i))
             call vector_col3(product_two, product, inv_D2)
             call vector_col3(total, product_two, projected_cotangent)
             if (this%experimental_algebraic) then
                ! The design mapping consumes a pointwise sensitivity density.
                ! The accumulated spatial cotangent is the assembled Euclidean
                ! covector, so remove one SEM weight here; map_backward then
                ! applies its metric once, matching the direct design derivative.
                call vector_invcol2(total, this%weights)
             end if
             call vector_add2s2(density, total, -1.0_rp)
          end do
       end if
       if (this%basis_transform == 'vector_potential') then
          this%potential_transpose_seconds = this%potential_transpose_seconds + &
               real(MPI_WTIME(),rp) - timer_start
       end if
       this%basis_transform_pullback_seconds = &
            this%basis_transform_pullback_seconds + real(MPI_WTIME(),rp) - timer_start
    end do

    call neko_scratch_registry%relinquish_field([field_work_index, transpose_work_index, &
         cot_x_index, cot_y_index, cot_z_index])
    call neko_scratch_registry%relinquish_vector([lambda_index, &
         state_old_index, state_new_index, delta_state_index, adjoint_initial_index, &
         initial_cotangent_index, product_index, product_two_index, current_chi_index, &
         candidate_D_index, inv_D2_index, total_index, projected_cotangent_index, &
         transform_factor_index])
    if (this%basis_transform == 'vector_potential') &
         call neko_scratch_registry%relinquish_vector(factor_cotangent_index)
    if (this%basis_transform == 'vector_potential') then
       do c = 1, 3
          call correction_cotangent(c)%free()
       end do
    end if
    do i = 1, r
       do c = 1, 3
          call mode_value(c,i)%free()
          do d = 1, 3; call mode_gradient(c,i,d)%free(); end do
       end do
    end do
    do c = 1, 3
       call initial_hom(c)%free()
    end do
    deallocate(mode_value, mode_gradient, initial_hom)
    nullify(lambda); nullify(state_old); nullify(state_new)
    nullify(delta_state); nullify(adjoint_initial); nullify(initial_cotangent)
    nullify(product); nullify(product_two); nullify(current_chi); nullify(candidate_D)
    nullify(inv_D2); nullify(total); nullify(projected_cotangent); nullify(transform_factor)
    nullify(factor_cotangent)
    nullify(field_work); nullify(transpose_work); nullify(cot_x); nullify(cot_y); nullify(cot_z)
    nullify(fluid)
  end subroutine rom_problem_add_adaptive_mode_density

  !> Add a scaled pointwise product to a mesh-sized mode cotangent.
  subroutine rom_problem_add_product(target, left, right, scale, work)
    type(vector_t), intent(inout) :: target
    type(vector_t), intent(in) :: left, right
    real(kind=rp), intent(in) :: scale
    type(vector_t), intent(inout) :: work
    if (abs(scale) <= tiny(1.0_rp)) return
    call vector_col3(work, left, right)
    call vector_add2s2(target, work, scale)
  end subroutine rom_problem_add_product

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
    this%has_scalar_mixing = .false.
    if (.not. parameters%valid_path('optimization.objectives')) then
       call neko_error('rom_mma requires a supported objective')
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
       case ('scalar_mixing')
          if (this%has_scalar_mixing) call neko_error('duplicate scalar mixing objective')
          this%has_scalar_mixing = .true.
          this%scalar_mixing_weight = weight
          this%scalar_mixing_start = start_time
          this%scalar_mixing_end = end_time
          this%scalar_mixing_mask_name = trim(mask_name)
          call json_get_or_default(item, 'target_concentration', &
               this%scalar_reference_value, 0.5_rp)
       case default
          call neko_error('unsupported objective in rom_mma: ' // trim(typ))
       end select
    end do
    if (.not. this%has_viscous .and. .not. this%has_brinkman .and. &
         .not. this%has_scalar_mixing) then
       call neko_error('rom_mma found no supported objective')
    end if
  end subroutine read_rom_objectives

  !> Parse the optional supported volume constraint configuration.
  subroutine read_rom_constraint(this, parameters)
    class(rom_problem_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(json_file) :: item
    character(len=:), allocatable :: typ, mask_name
    logical :: is_max
    real(kind=rp) :: limit
    integer :: n_items

    this%has_volume_constraint = .false.
    this%volume_mask_name = ''
    this%volume_limit = 0.0_rp
    this%volume_is_max = .false.
    if (.not. parameters%valid_path('optimization.constraints')) return
    call parameters%info('optimization.constraints', n_children = n_items)
    if (n_items == 0) return
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
    this%has_volume_constraint = .true.
  end subroutine read_rom_constraint

  !> Rebuild basis and fixed operators after a full-order refresh.
  subroutine rom_problem_refresh(this, simulation)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), target, intent(inout) :: simulation
    this%simulation => simulation
    call rom_problem_load_local_basis(this, simulation)
    call rom_problem_assemble_fixed(this, simulation)
    call rom_problem_setup_masks(this, simulation)
    call rom_problem_assemble_objective_operators(this)
  end subroutine rom_problem_refresh

  !> Receive and construct the selected forward/adjoint basis spaces.
  subroutine rom_problem_load_local_basis(this, simulation)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    integer :: n_modes, n_fields
    real(kind=rp), allocatable :: local_basis(:, :)
    character(len=256) :: log_message
    this%projected_tensor_cache_valid = .false.
    this%scalar_mixing_cache_valid = .false.
    if (this%uses_dual_basis) then
       call rom_problem_load_local_forward_constrained_basis(this, simulation)
       return
    end if
    if (this%petrov_galerkin) then
       call rom_problem_load_local_petrov_basis(this, simulation)
       return
    end if
    call rom_receive_basis(simulation%state_recover, simulation%fluid%u%x, &
         local_basis, n_modes, n_fields)
    if ((n_fields /= 3 .and. n_fields /= 4) .or. &
         (n_fields == 4 .and. .not. associated(simulation%scalars))) then
       call neko_error('ROM basis fields do not match the configured simulation')
    end if
    call rom_project_basis(simulation, local_basis, n_modes, n_fields)
    this%basis%n_fields = n_fields
    if (n_modes < 1) call neko_error('Python ROM peer returned an empty basis')
    this%basis%n_rows = size(local_basis, 1)
    this%basis%n_modes = n_modes
    this%rank_forward = n_modes
    this%rank_adjoint = n_modes
    this%rank_enriched = n_modes
    write(log_message, '("ROM basis selection: forward=",A,", adjoint=",A,", retained forward rank=",I0)') &
         trim(this%basis_forward), trim(this%basis_adjoint), this%rank_forward
    call neko_log%message(trim(log_message))
    if (this%basis%n_rows /= n_fields * simulation%fluid%u%dof%size()) then
       call neko_error('ROM basis rows do not match local velocity fields')
    end if
    call this%basis%modes%free()
    call this%basis%modes%init(size(local_basis,1), size(local_basis,2))
    this%basis%modes%x = local_basis
    ! The received matrix is host staging.  Retained mode vectors are uploaded
    ! once by assemble_fixed at the basis-refresh boundary.
    deallocate(local_basis)
  end subroutine rom_problem_load_local_basis

  !> Receive forward and adjoint POD spaces and construct
  !! the forward POD space V_f and the independently streamed adjoint POD
  !! space V_a.  The selected formulation then chooses between V_f, V_a and
  !! V_E = orth_M([V_f,V_a]) without allowing adjoint modes into a restricted
  !! forward trajectory.
  !!
  !! The forward and adjoint POD blocks are mass-orthonormalised
  !! independently.  `separate_pod` retains those two spaces.  The enriched
  !! formulations additionally build V_E from their concatenation.  The
  !! original V_a overlap is retained in the span; it is not replaced by an
  !! orthogonal-complement-only construction.
  subroutine rom_problem_load_local_forward_constrained_basis(this, simulation)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    integer :: trial_modes, trial_fields, adjoint_modes, adjoint_fields
    integer :: n, trial_retained, adjoint_retained, enriched_retained, fields
    integer :: enriched_input_modes, trial_use, adjoint_use, enriched_use
    integer, allocatable :: trial_mode_index(:), adjoint_mode_index(:)
    real(kind=rp), allocatable :: trial_basis(:, :), adjoint_basis(:, :)
    real(kind=rp), allocatable :: trial_velocity(:, :), adjoint_velocity(:, :)
    real(kind=rp), allocatable :: trial_scalar(:,:), adjoint_scalar(:,:), mass(:)
    real(kind=rp), allocatable :: enriched_velocity(:, :), enriched_scalar(:,:)
    real(kind=rp), allocatable :: trial_norms(:), trial_crosses(:)
    real(kind=rp), allocatable :: adjoint_norms(:), adjoint_crosses(:)
    real(kind=rp) :: scalar_weight, tolerance
    logical :: trial_ok, adjoint_ok, enriched_ok
    character(len=256) :: log_message
    this%projected_tensor_cache_valid = .false.
    this%scalar_mixing_cache_valid = .false.
    this%scalar_transport_cache_valid = .false.

    call rom_receive_petrov_bases(simulation%state_recover, simulation%fluid%u%x, &
         trial_basis, trial_modes, trial_fields, adjoint_basis, adjoint_modes, &
         adjoint_fields)
    if (trial_fields /= adjoint_fields .or. &
         (trial_fields /= 3 .and. trial_fields /= 4)) call neko_error( &
         'forward_constrained_enriched basis field counts do not match')
    if (trial_modes < 1 .or. adjoint_modes < 1) call neko_error( &
         'forward_constrained_enriched requires nonempty forward and adjoint POD spaces')
    if (trial_fields == 4 .and. .not. associated(simulation%scalars)) call neko_error( &
         'forward_constrained_enriched scalar basis supplied to a scalar-free simulation')
    n = simulation%fluid%u%dof%size(); fields = trial_fields
    if (size(trial_basis, 1) /= fields*n .or. size(adjoint_basis, 1) /= fields*n .or. &
         size(trial_basis, 2) /= trial_modes .or. size(adjoint_basis, 2) /= adjoint_modes) &
         call neko_error('forward_constrained_enriched received inconsistent basis shapes')
    call rom_project_basis(simulation, trial_basis, trial_modes, trial_fields)
    call rom_project_test_basis(simulation, adjoint_basis, adjoint_modes, adjoint_fields)
    if (trial_modes < 1 .or. adjoint_modes < 1) call neko_error( &
         'forward_constrained_enriched admissibility processing removed every basis mode')
    if (size(trial_basis, 1) /= fields*n .or. size(adjoint_basis, 1) /= fields*n .or. &
         size(trial_basis, 2) /= trial_modes .or. size(adjoint_basis, 2) /= adjoint_modes) &
         call neko_error('forward_constrained_enriched projection changed basis shapes')

    scalar_weight = this%scalar_metric_weight
    tolerance = this%admissibility_tolerance
    if (scalar_weight <= 0.0_rp) call neko_error( &
         'ROM scalar_metric_weight must be positive')
    allocate(mass(n), trial_mode_index(trial_modes), &
         adjoint_mode_index(trial_modes + adjoint_modes), &
         trial_norms(trial_modes), trial_crosses(trial_modes), &
         adjoint_norms(trial_modes + adjoint_modes), &
         adjoint_crosses(trial_modes + adjoint_modes))
    mass = reshape(simulation%fluid%c_Xh%B, [n])
    ! The streamed forward and adjoint POD blocks define two subspaces of the
    ! same enriched physical space V_E = span([V_f,V_a]).  Orthonormalise each
    ! block independently; do not replace V_a by its complement to V_f.
    if (fields == 4) then
       allocate(trial_velocity(3*n,trial_modes), adjoint_velocity(3*n,adjoint_modes), &
            trial_scalar(n,trial_modes), adjoint_scalar(n,adjoint_modes))
       if (size(trial_scalar,1) /= n .or. size(adjoint_scalar,1) /= n) &
            call neko_error('forward_constrained_enriched scalar workspace shape mismatch')
       trial_velocity = trial_basis(1:3*n,:)
       adjoint_velocity = adjoint_basis(1:3*n,:)
       trial_scalar = trial_basis(3*n+1:4*n,:)
       adjoint_scalar = adjoint_basis(3*n+1:4*n,:)
       call rom_orthonormalize_projected_coupled_device(trial_velocity, trial_scalar, n, &
            trial_modes, mass, scalar_weight, tolerance, trial_retained, &
            trial_mode_index, trial_norms, trial_crosses, trial_ok)
       call rom_orthonormalize_projected_coupled_device(adjoint_velocity, adjoint_scalar, n, &
            adjoint_modes, mass, scalar_weight, tolerance, adjoint_retained, &
            adjoint_mode_index(1:adjoint_modes), adjoint_norms(1:adjoint_modes), &
            adjoint_crosses(1:adjoint_modes), adjoint_ok)
       trial_basis(1:3*n,1:trial_retained) = trial_velocity
       adjoint_basis(1:3*n,1:adjoint_retained) = adjoint_velocity
       trial_basis(3*n+1:4*n,1:trial_retained) = trial_scalar
       adjoint_basis(3*n+1:4*n,1:adjoint_retained) = adjoint_scalar
    else
       call rom_orthonormalize_projected_device(trial_basis, n, trial_modes, mass, &
            tolerance, trial_retained, trial_mode_index, trial_norms, trial_crosses, trial_ok)
       call rom_orthonormalize_projected_device(adjoint_basis, n, adjoint_modes, mass, &
            tolerance, adjoint_retained, adjoint_mode_index(1:adjoint_modes), &
            adjoint_norms(1:adjoint_modes), adjoint_crosses(1:adjoint_modes), &
            adjoint_ok)
    end if
    if (.not. trial_ok .or. .not. adjoint_ok .or. trial_retained < 1 .or. &
         adjoint_retained < 1) call neko_error( &
         'forward_constrained_enriched failed independent coupled mass orthonormalisation')
    if (size(trial_basis,1) /= fields*n .or. size(adjoint_basis,1) /= fields*n .or. &
         size(trial_basis,2) < trial_retained .or. size(adjoint_basis,2) < adjoint_retained) &
         call neko_error('forward_constrained_enriched independent basis shape mismatch')
    trial_use = trial_retained
    adjoint_use = adjoint_retained
    if (trial_use < 1 .or. adjoint_use < 1) call neko_error( &
         'ROM numerical rank removed every forward or adjoint mode')
    trial_modes = trial_use
    adjoint_modes = adjoint_use
    trial_basis = trial_basis(:,1:trial_modes)
    adjoint_basis = adjoint_basis(:,1:adjoint_modes)
    if (fields == 4) then
       trial_scalar = trial_scalar(:,1:trial_modes)
       adjoint_scalar = adjoint_scalar(:,1:adjoint_modes)
    end if

    ! Build the full enriched space only for formulations that request it.
    ! The concatenation is mass-orthonormalised as one coupled basis, so the
    ! scalar component follows exactly the same mode rotations and rejection
    ! decisions as the velocity component.
    if (this%forward_enriched_adjoint .or. this%enriched_galerkin) then
       enriched_input_modes = trial_modes + adjoint_modes
       allocate(enriched_velocity(3*n,enriched_input_modes))
       enriched_velocity(:,1:trial_modes) = trial_basis(1:3*n,:)
       enriched_velocity(:,trial_modes+1:enriched_input_modes) = &
            adjoint_basis(1:3*n,:)
       if (fields == 4) then
          allocate(enriched_scalar(n,enriched_input_modes))
          enriched_scalar(:,1:trial_modes) = trial_scalar
          enriched_scalar(:,trial_modes+1:enriched_input_modes) = adjoint_scalar
          call rom_orthonormalize_projected_coupled_device(enriched_velocity, &
               enriched_scalar, n, enriched_input_modes, mass, scalar_weight, &
               tolerance, enriched_retained, adjoint_mode_index, adjoint_norms, &
               adjoint_crosses, enriched_ok)
       else
          call rom_orthonormalize_projected_device(enriched_velocity, n, &
               enriched_input_modes, mass, tolerance, enriched_retained, &
               adjoint_mode_index, adjoint_norms, adjoint_crosses, enriched_ok)
       end if
       if (.not. enriched_ok .or. enriched_retained < 1) call neko_error( &
            'ROM enriched basis lost every mode during mass orthonormalisation')
       enriched_use = enriched_retained
       if (this%forward_enriched_adjoint .and. enriched_use < trial_modes) &
            call neko_error( &
            'numerical rank loss broke forward:enriched basis containment')
    else
       enriched_use = 0
    end if

    call this%basis%modes%free(); call this%test_basis%modes%free()
    if (this%enriched_galerkin) then
       this%basis%n_rows = fields*n; this%basis%n_modes = enriched_use
       this%test_basis%n_rows = fields*n; this%test_basis%n_modes = enriched_use
    else if (this%forward_enriched_adjoint) then
       this%basis%n_rows = fields*n; this%basis%n_modes = trial_modes
       this%test_basis%n_rows = fields*n; this%test_basis%n_modes = enriched_use
    else
       this%basis%n_rows = fields*n; this%basis%n_modes = trial_modes
       this%test_basis%n_rows = fields*n; this%test_basis%n_modes = adjoint_modes
    end if
    this%basis%n_fields = fields
    this%test_basis%n_fields = fields
    call this%basis%modes%init(fields*n, this%basis%n_modes)
    call this%test_basis%modes%init(fields*n, this%test_basis%n_modes)
    if (size(this%basis%modes%x,1) /= fields*n .or. &
         size(this%test_basis%modes%x,1) /= fields*n) call neko_error( &
         'forward_constrained_enriched output basis shape mismatch')
    if (this%enriched_galerkin) then
       this%basis%modes%x(1:3*n,:) = enriched_velocity(:,1:enriched_use)
       this%test_basis%modes%x(1:3*n,:) = enriched_velocity(:,1:enriched_use)
       if (fields == 4) then
          this%basis%modes%x(3*n+1:4*n,:) = enriched_scalar(:,1:enriched_use)
          this%test_basis%modes%x(3*n+1:4*n,:) = enriched_scalar(:,1:enriched_use)
       end if
    else if (this%forward_enriched_adjoint) then
       this%basis%modes%x = trial_basis
       this%test_basis%modes%x(1:3*n,:) = enriched_velocity(:,1:enriched_use)
       if (fields == 4) this%test_basis%modes%x(3*n+1:4*n,:) = &
            enriched_scalar(:,1:enriched_use)
    else
       this%basis%modes%x = trial_basis
       this%test_basis%modes%x = adjoint_basis
    end if
    this%has_scalar = fields == 4
    this%rank_forward = trial_modes
    this%rank_adjoint = adjoint_modes
    this%rank_enriched = enriched_use
    this%enriched_adjoint_rank = this%test_basis%n_modes
    write(log_message, '("ROM basis selection: forward=",A,", adjoint=",A,", retained Vf/Va/VE=",I0,"/",I0,"/",I0)') &
         trim(this%basis_forward), trim(this%basis_adjoint), this%rank_forward, &
         this%rank_adjoint, this%rank_enriched
    call neko_log%message(trim(log_message))

    deallocate(trial_basis, adjoint_basis, mass, trial_mode_index, adjoint_mode_index, &
         trial_norms, trial_crosses, adjoint_norms, adjoint_crosses)
    if (fields == 4) deallocate(trial_velocity, adjoint_velocity, trial_scalar, adjoint_scalar)
    if (allocated(enriched_velocity)) deallocate(enriched_velocity)
    if (allocated(enriched_scalar)) deallocate(enriched_scalar)
  end subroutine rom_problem_load_local_forward_constrained_basis

  !> Receive and independently prepare the forward trial and adjoint test
  !! bases used by the Petrov--Galerkin formulation.
  subroutine rom_problem_load_local_petrov_basis(this, simulation)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    integer :: trial_modes, trial_fields, test_modes, test_fields, common_modes
    integer :: n
    real(kind=rp), allocatable :: trial_basis(:, :), test_basis(:, :)

    this%projected_tensor_cache_valid = .false.
    this%scalar_mixing_cache_valid = .false.
    this%scalar_transport_cache_valid = .false.

    call rom_receive_petrov_bases(simulation%state_recover, simulation%fluid%u%x, &
         trial_basis, trial_modes, trial_fields, test_basis, test_modes, test_fields)
    if (trial_fields /= test_fields .or. (trial_fields /= 3 .and. trial_fields /= 4)) &
         call neko_error('Petrov trial/test basis field counts do not match')
    if (trial_modes /= test_modes) call neko_error( &
         'Petrov trial/test basis ranks must match before admissibility processing')
    if (trial_fields == 4 .and. .not. associated(simulation%scalars)) then
       call neko_error('Petrov scalar basis supplied to a scalar-free simulation')
    end if
    call rom_project_basis(simulation, trial_basis, trial_modes, trial_fields)
    call rom_project_test_basis(simulation, test_basis, test_modes, test_fields)
    call rom_problem_verify_adjoint_basis_bc(this, simulation, test_basis, test_modes)
    ! The independent admissibility operators can reject different numbers of
    ! nearly dependent columns.  Petrov assembly requires equal trial/test
    ! dimensions, so retain their common leading subspace rather than
    ! aborting after a valid basis handoff.  The spaces remain independent;
    ! this only enforces the rectangular dimensions required by the reduced
    ! cross mass.
    common_modes = min(trial_modes, test_modes)
    if (common_modes < 1) call neko_error( &
         'independent Petrov admissibility processing removed every common mode')
    if (trial_modes /= common_modes) then
       trial_basis = trial_basis(:, 1:common_modes)
       trial_modes = common_modes
    end if
    if (test_modes /= common_modes) then
       test_basis = test_basis(:, 1:common_modes)
       test_modes = common_modes
    end if
    n = simulation%fluid%u%dof%size()
    if (size(trial_basis,1) /= trial_fields*n .or. size(test_basis,1) /= test_fields*n) &
         call neko_error('Petrov basis rows do not match local velocity fields')
    call this%basis%modes%free(); call this%test_basis%modes%free()
    this%basis%n_rows = size(trial_basis,1)
    this%basis%n_modes = trial_modes
    this%basis%n_fields = trial_fields
    this%test_basis%n_rows = size(test_basis,1)
    this%test_basis%n_modes = test_modes
    this%test_basis%n_fields = test_fields
    this%rank_forward = trial_modes
    this%rank_adjoint = test_modes
    this%rank_enriched = max(trial_modes, test_modes)
    call this%basis%modes%init(size(trial_basis,1), trial_modes)
    call this%test_basis%modes%init(size(test_basis,1), test_modes)
    this%basis%modes%x = trial_basis
    this%test_basis%modes%x = test_basis
    this%has_scalar = trial_fields == 4
    deallocate(trial_basis, test_basis)
  end subroutine rom_problem_load_local_petrov_basis

  !> Verify that independently received adjoint test modes satisfy Neko's
  !! actual homogeneous adjoint velocity restriction.  This is a reception
  !! boundary check; it is deliberately not repeated during ROM integration.
  subroutine rom_problem_verify_adjoint_basis_bc(this, simulation, basis, n_modes)
    class(rom_problem_t), intent(in) :: this
    type(simulation_t), intent(inout) :: simulation
    real(kind=rp), intent(in) :: basis(:, :)
    integer, intent(in) :: n_modes
    type(field_t), pointer :: u, v, w, uc, vc, wc
    integer :: indices(6), n, mode
    real(kind=rp) :: residual

    if (.not. allocated(simulation%adjoint_case%fluid_adj)) call neko_error( &
         'Petrov test-basis BC check requires an initialized adjoint fluid')
    n = simulation%fluid%u%dof%size()
    call neko_scratch_registry%request_field(u, indices(1), simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(v, indices(2), simulation%fluid%v%dof, .true.)
    call neko_scratch_registry%request_field(w, indices(3), simulation%fluid%w%dof, .true.)
    call neko_scratch_registry%request_field(uc, indices(4), simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(vc, indices(5), simulation%fluid%v%dof, .true.)
    call neko_scratch_registry%request_field(wc, indices(6), simulation%fluid%w%dof, .true.)
    select type (adjoint_fluid => simulation%adjoint_case%fluid_adj)
    type is (adjoint_fluid_pnpn_t)
       do mode = 1, n_modes
          u%x = reshape(basis(1:n,mode), shape(u%x))
          v%x = reshape(basis(n+1:2*n,mode), shape(v%x))
          w%x = reshape(basis(2*n+1:3*n,mode), shape(w%x))
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call u%copy_from(HOST_TO_DEVICE, sync=.true.)
             call v%copy_from(HOST_TO_DEVICE, sync=.true.)
             call w%copy_from(HOST_TO_DEVICE, sync=.true.)
          end if
          ! Copy only after the source mirror is current.  On an accelerator
          ! field_copy is device-aware; copying before the H2D update would
          ! preserve a stale device source and then overwrite the projected
          ! check field with stale host data.
          call field_copy(uc, u); call field_copy(vc, v); call field_copy(wc, w)
          call adjoint_fluid%bcs_vel_projector%apply(uc%x, vc%x, wc%x, n)
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call uc%copy_from(DEVICE_TO_HOST, sync=.true.)
             call vc%copy_from(DEVICE_TO_HOST, sync=.true.)
             call wc%copy_from(DEVICE_TO_HOST, sync=.true.)
          end if
          residual = max(maxval(abs(uc%x-u%x)), max(maxval(abs(vc%x-v%x)), &
               maxval(abs(wc%x-w%x))))
          if (residual > this%admissibility_tolerance) call neko_error( &
               'Petrov adjoint test basis violates homogeneous adjoint velocity BCs')
       end do
    class default
       call neko_error('Petrov test-basis BC check requires adjoint_fluid_pnpn')
    end select
    call neko_scratch_registry%relinquish_field(indices)
    nullify(u); nullify(v); nullify(w); nullify(uc); nullify(vc); nullify(wc)
  end subroutine rom_problem_verify_adjoint_basis_bc

  !> Apply scalar boundary restrictions to the owned basis matrix.
  subroutine rom_problem_restrict_scalar_basis(this, simulation)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    real(kind=rp), allocatable :: scalar_basis(:, :)
    integer :: n_velocity, n_scalar
    if (this%basis%n_fields /= 4) return
    if (.not. associated(simulation%scalars)) call neko_error( &
         'scalar ROM basis supplied to a scalar-free simulation')
    n_velocity = simulation%fluid%u%dof%size()
    n_scalar = simulation%scalars%scalar_fields(1)%scalar%s%dof%size()
    if (n_scalar /= n_velocity) call neko_error( &
         'coupled ROM basis requires matching velocity and scalar local sizes')
    allocate(scalar_basis(n_scalar, this%basis%n_modes))
    scalar_basis = this%basis%modes%x(3*n_velocity+1:4*n_velocity, :)
    call rom_restrict_scalar_basis(simulation, scalar_basis, this%basis%n_modes)
    this%basis%modes%x(3*n_velocity+1:4*n_velocity, :) = scalar_basis
  end subroutine rom_problem_restrict_scalar_basis

  !> Assemble mass, lift, convection and diffusion operators independent of design.
  subroutine rom_problem_assemble_fixed(this, simulation)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    integer :: n, r, e, i, c, d, offset
    type(matrix_t), pointer :: mass_local, trial_mass_local, linear_local, implicit_local
    type(matrix_t), pointer :: enriched_mass_local, enriched_linear_local, enriched_implicit_local
    type(matrix_t), allocatable :: quad_local(:)
    type(vector_t), pointer :: rhs_local, constant_local, enriched_constant_local, initial_difference
    type(field_t), pointer :: field, gx, gy, gz
    type(field_t), pointer :: scalar_field, scalar_gx, scalar_gy, scalar_gz
    type(field_t), pointer :: scalar_diff_field, scalar_rho_cp
    integer :: mass_local_index, trial_mass_local_index, linear_local_index, implicit_local_index
    integer :: enriched_mass_local_index, enriched_linear_local_index
    integer :: enriched_implicit_local_index
    integer :: rhs_local_index, constant_local_index, enriched_constant_local_index
    integer :: initial_difference_index
    integer :: field_index, gx_index, gy_index, gz_index
    integer :: scalar_field_index, scalar_gx_index, scalar_gy_index, scalar_gz_index
    integer :: scalar_diff_index, scalar_rho_cp_index
    real(kind=rp) :: nu
    logical :: shared_trial_test_modes

    if (this%basis_transform_adaptive) call rom_problem_free_reference_basis(this)
    call rom_problem_free_arrays(this)
    n = simulation%fluid%u%dof%size()
    if (this%enriched_galerkin) then
       r = this%basis%n_modes
    else
       r = this%basis%n_modes
    end if
    e = r
    if (this%separate_pod) e = this%test_basis%n_modes
    if (this%forward_enriched_adjoint .or. this%enriched_galerkin) then
       e = this%test_basis%n_modes
    end if
    if (r < 1) call neko_error('ROM basis contains no retained trial/test modes')
    ! In the equal-space Galerkin paths, test modes are explicitly copied from
    ! the trial modes above. Reuse their derivatives rather than applying the
    ! same native SEM grad operator a second time. Mixed/enriched spaces are
    ! left independent unless identity is guaranteed by that construction.
    shared_trial_test_modes = NEKO_BCKND_DEVICE .eq. 0 .and. &
         .not. this%petrov_galerkin .and. &
         .not. this%forward_constrained_enriched
    call this%weights%init(n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(this%weights%x_d, simulation%fluid%c_Xh%B_d, n)
    else
       this%weights%x = reshape(simulation%fluid%c_Xh%B, [n])
    end if
    allocate(this%velocity_modes(3,r))
    allocate(this%test_velocity_modes(3,e))
    do i = 1, r
       do c = 1, 3
          call this%velocity_modes(c,i)%init(n)
          if (.not. this%forward_constrained_enriched) then
             call this%test_velocity_modes(c,i)%init(n)
          end if
       end do
    end do
    if (this%forward_constrained_enriched) then
       do i = 1, e
          do c = 1, 3
             call this%test_velocity_modes(c,i)%init(n)
          end do
       end do
    end if
    if (this%has_scalar) then
       allocate(this%scalar_modes(r))
       allocate(this%test_scalar_modes(e))
       do i = 1, r
          call this%scalar_modes(i)%init(n)
          if (.not. this%forward_constrained_enriched) call this%test_scalar_modes(i)%init(n)
       end do
       if (this%forward_constrained_enriched) then
          do i = 1, e
             call this%test_scalar_modes(i)%init(n)
          end do
       end if
       call this%scalar_lift%free(); call this%scalar_reference%free()
       call this%lift%copy_scalar(this%scalar_lift)
       call this%scalar_reference%init(n)
       call field_to_vector(this%scalar_reference, this%lift%scalar_reference)
       do d = 1, 3
          call this%scalar_lift_gradients(d)%free()
          call this%scalar_lift_gradients(d)%init(n)
          call field_to_vector(this%scalar_lift_gradients(d), &
               this%lift%scalar_gradient(d))
       end do
    end if
    do c = 1, 3
       call this%lift_velocity(c)%free()
    end do
    call this%lift%copy_velocity(this%lift_velocity)
    do c = 1, 3
       do d = 1, 3
          call this%lift_gradients(c,d)%init(n)
          call field_to_vector(this%lift_gradients(c,d), &
               this%lift%gradient(c,d))
       end do
    end do
    do i = 1, r
       do c = 1, 3
          offset = (c - 1) * n
          this%velocity_modes(c, i)%x = &
               this%basis%modes%x(offset + 1:offset + n, i)
          if (.not. this%forward_constrained_enriched) then
             if (this%petrov_galerkin) then
                this%test_velocity_modes(c, i)%x = &
                     this%test_basis%modes%x(offset + 1:offset + n, i)
             else
                this%test_velocity_modes(c, i)%x = this%velocity_modes(c, i)%x
             end if
          end if
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call this%velocity_modes(c, i)%copy_from(HOST_TO_DEVICE, sync = .true.)
             if (.not. this%forward_constrained_enriched) then
                call this%test_velocity_modes(c, i)%copy_from(HOST_TO_DEVICE, sync = .true.)
             end if
          end if
       end do
    end do
    if (this%forward_constrained_enriched) then
       do i = 1, e
          do c = 1, 3
             offset = (c - 1) * n
             this%test_velocity_modes(c,i)%x = this%test_basis%modes%x( &
                  offset + 1:offset + n, i)
             if (NEKO_BCKND_DEVICE .eq. 1) call this%test_velocity_modes(c,i)%copy_from( &
                  HOST_TO_DEVICE, sync = .true.)
          end do
       end do
    end if
    if (this%has_scalar) then
       do i = 1, r
          this%scalar_modes(i)%x = this%basis%modes%x(3*n+1:4*n, i)
          if (.not. this%forward_constrained_enriched) then
             if (this%petrov_galerkin) then
                this%test_scalar_modes(i)%x = this%test_basis%modes%x(3*n+1:4*n, i)
             else
                this%test_scalar_modes(i)%x = this%scalar_modes(i)%x
             end if
          end if
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call this%scalar_modes(i)%copy_from(HOST_TO_DEVICE, sync = .true.)
             if (.not. this%forward_constrained_enriched) call this%test_scalar_modes(i)%copy_from( &
                  HOST_TO_DEVICE, sync = .true.)
          end if
       end do
       if (this%forward_constrained_enriched) then
          do i = 1, e
             this%test_scalar_modes(i)%x = this%test_basis%modes%x(3*n+1:4*n, i)
             if (NEKO_BCKND_DEVICE .eq. 1) call this%test_scalar_modes(i)%copy_from( &
                  HOST_TO_DEVICE, sync = .true.)
          end do
       end if
    end if

    ! Pair the complete trial and test spaces before any gradients or reduced
    ! operators are formed.  The pairing rotates device-resident mesh vectors
    ! and truncates weakly overlapping directions; all later assembly therefore
    ! sees exactly the final Petrov spaces.
    if (this%petrov_galerkin) call rom_problem_check_cross_mass(this)
    r = size(this%velocity_modes, 2)
    allocate(this%velocity_gradients(3,r,3))
    allocate(this%test_velocity_gradients(3,e,3))
    do i = 1, r
       do c = 1, 3
          do d = 1, 3
             call this%velocity_gradients(c,i,d)%init(n)
          end do
       end do
    end do
    do i = 1, e
       do c = 1, 3
          do d = 1, 3
             call this%test_velocity_gradients(c,i,d)%init(n)
          end do
       end do
    end do
    if (this%has_scalar) then
       allocate(this%scalar_gradients(r,3))
       allocate(this%test_scalar_gradients(e,3))
       do i = 1, r
          do d = 1, 3
             call this%scalar_gradients(i,d)%init(n)
          end do
       end do
       do i = 1, e
          do d = 1, 3
             call this%test_scalar_gradients(i,d)%init(n)
          end do
       end do
    end if

    call neko_scratch_registry%request_field(field, field_index, &
         simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(gx, gx_index, simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(gy, gy_index, simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(gz, gz_index, simulation%fluid%u%dof, .true.)
    if (this%has_scalar) then
       call neko_scratch_registry%request_field(scalar_field, scalar_field_index, &
            simulation%scalars%scalar_fields(1)%scalar%s%dof, .true.)
       call neko_scratch_registry%request_field(scalar_gx, scalar_gx_index, &
            simulation%scalars%scalar_fields(1)%scalar%s%dof, .true.)
       call neko_scratch_registry%request_field(scalar_gy, scalar_gy_index, &
            simulation%scalars%scalar_fields(1)%scalar%s%dof, .true.)
       call neko_scratch_registry%request_field(scalar_gz, scalar_gz_index, &
            simulation%scalars%scalar_fields(1)%scalar%s%dof, .true.)
    end if
    do i = 1, r
       do c = 1, 3
          call vector_to_field(field, this%velocity_modes(c,i))
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call grad(gx%x_d, gy%x_d, gz%x_d, field%x_d, simulation%fluid%c_Xh)
          else
             call grad(gx%x, gy%x, gz%x, field%x, simulation%fluid%c_Xh)
          end if
          call field_to_vector(this%velocity_gradients(c,i,1), gx)
          call field_to_vector(this%velocity_gradients(c,i,2), gy)
          call field_to_vector(this%velocity_gradients(c,i,3), gz)
      end do
   end do
    do i = 1, e
       if (shared_trial_test_modes .and. i <= r) then
          do c = 1, 3
             do d = 1, 3
                call vector_copy(this%test_velocity_gradients(c,i,d), &
                     this%velocity_gradients(c,i,d))
             end do
          end do
          cycle
       end if
       do c = 1, 3
          call vector_to_field(field, this%test_velocity_modes(c,i))
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call grad(gx%x_d, gy%x_d, gz%x_d, field%x_d, simulation%fluid%c_Xh)
          else
             call grad(gx%x, gy%x, gz%x, field%x, simulation%fluid%c_Xh)
          end if
          call field_to_vector(this%test_velocity_gradients(c,i,1), gx)
          call field_to_vector(this%test_velocity_gradients(c,i,2), gy)
          call field_to_vector(this%test_velocity_gradients(c,i,3), gz)
       end do
    end do
    if (this%has_scalar) then
       do i = 1, r
          call vector_to_field(scalar_field, this%scalar_modes(i))
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call grad(scalar_gx%x_d, scalar_gy%x_d, scalar_gz%x_d, &
                  scalar_field%x_d, simulation%scalars%scalar_fields(1)%scalar%c_Xh)
          else
             call grad(scalar_gx%x, scalar_gy%x, scalar_gz%x, scalar_field%x, &
                  simulation%scalars%scalar_fields(1)%scalar%c_Xh)
          end if
          call field_to_vector(this%scalar_gradients(i,1), scalar_gx)
          call field_to_vector(this%scalar_gradients(i,2), scalar_gy)
          call field_to_vector(this%scalar_gradients(i,3), scalar_gz)
       end do
       do i = 1, e
          if (shared_trial_test_modes .and. i <= r) then
             do d = 1, 3
                call vector_copy(this%test_scalar_gradients(i,d), &
                     this%scalar_gradients(i,d))
             end do
             cycle
          end if
          call vector_to_field(scalar_field, this%test_scalar_modes(i))
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call grad(scalar_gx%x_d, scalar_gy%x_d, scalar_gz%x_d, &
                  scalar_field%x_d, simulation%scalars%scalar_fields(1)%scalar%c_Xh)
          else
             call grad(scalar_gx%x, scalar_gy%x, scalar_gz%x, scalar_field%x, &
                  simulation%scalars%scalar_fields(1)%scalar%c_Xh)
          end if
          call field_to_vector(this%test_scalar_gradients(i,1), scalar_gx)
          call field_to_vector(this%test_scalar_gradients(i,2), scalar_gy)
          call field_to_vector(this%test_scalar_gradients(i,3), scalar_gz)
       end do
    end if

    call neko_scratch_registry%request_matrix(mass_local, mass_local_index, r, r, .true.)
    call neko_scratch_registry%request_matrix(trial_mass_local, trial_mass_local_index, r, r, .true.)
    call neko_scratch_registry%request_matrix(linear_local, linear_local_index, r, r, .true.)
    call neko_scratch_registry%request_matrix(implicit_local, implicit_local_index, &
         r, r, .true.)
    call neko_scratch_registry%request_vector(constant_local, constant_local_index, &
         r, .true.)
    allocate(quad_local(r))
    if (this%forward_constrained_enriched) then
       call neko_scratch_registry%request_matrix(enriched_mass_local, enriched_mass_local_index, &
            e, e, .true.)
       call neko_scratch_registry%request_matrix(enriched_linear_local, enriched_linear_local_index, &
            e, e, .true.)
       call neko_scratch_registry%request_matrix(enriched_implicit_local, &
            enriched_implicit_local_index, e, e, .true.)
       call neko_scratch_registry%request_vector(enriched_constant_local, &
            enriched_constant_local_index, e, .true.)
    end if
    if (this%has_scalar) then
       call rom_assemble_mass_coupled(this%weights, this%velocity_modes, &
            this%scalar_modes, this%scalar_metric_weight, trial_mass_local)
       if (this%petrov_galerkin) then
          call rom_assemble_cross_mass(this%weights, this%velocity_modes, &
               this%test_velocity_modes, mass_local, this%scalar_modes, &
               this%test_scalar_modes, this%scalar_metric_weight)
       else
          mass_local%x = trial_mass_local%x
       end if
    else
       call rom_assemble_mass(this%weights, this%velocity_modes, trial_mass_local)
       if (this%petrov_galerkin) then
          call rom_assemble_cross_mass(this%weights, this%velocity_modes, &
               this%test_velocity_modes, mass_local)
       else
          mass_local%x = trial_mass_local%x
       end if
    end if
    call rom_uniform_coefficient(simulation%fluid%mu, this%weights, &
         'viscosity', this%viscosity)
    call rom_uniform_coefficient(simulation%fluid%rho, this%weights, &
         'density', this%density)
    nu = this%viscosity / this%density
    if (this%petrov_galerkin) then
       call rom_assemble_lift_convection_petrov(this%weights, this%velocity_modes, &
            this%velocity_gradients, this%test_velocity_modes, this%lift_velocity, &
            this%lift_gradients, linear_local, constant_local)
       call rom_assemble_diffusion_petrov(this%weights, nu, this%test_velocity_gradients, &
            this%velocity_gradients, implicit_local, this%lift_gradients, constant_local)
       call rom_assemble_convection_petrov(this%weights, this%test_velocity_modes, &
            this%velocity_modes, this%velocity_gradients, quad_local)
    else
       call rom_assemble_lift_convection(this%weights, this%velocity_modes, &
            this%velocity_gradients, this%lift_velocity, this%lift_gradients, &
            linear_local, constant_local)
       call rom_assemble_diffusion(this%weights, nu, this%velocity_gradients, &
            implicit_local, this%lift_gradients, constant_local)
       call rom_assemble_convection(this%weights, this%velocity_modes, &
            this%velocity_gradients, quad_local)
    end if
    if (this%forward_constrained_enriched) then
       if (this%has_scalar) then
          call rom_assemble_mass_coupled(this%weights, this%test_velocity_modes, &
               this%test_scalar_modes, this%scalar_metric_weight, enriched_mass_local)
       else
          call rom_assemble_mass(this%weights, this%test_velocity_modes, enriched_mass_local)
       end if
       call rom_assemble_lift_convection_petrov(this%weights, this%test_velocity_modes, &
            this%test_velocity_gradients, this%test_velocity_modes, this%lift_velocity, &
            this%lift_gradients, enriched_linear_local, enriched_constant_local)
       call rom_assemble_diffusion_petrov(this%weights, nu, this%test_velocity_gradients, &
            this%test_velocity_gradients, enriched_implicit_local, this%lift_gradients, &
            enriched_constant_local)
    end if
    do i = 1, r
       quad_local(i)%x = -quad_local(i)%x
    end do
    if (this%has_scalar) then
       call neko_scratch_registry%request_field(scalar_diff_field, scalar_diff_index, &
            simulation%scalars%scalar_fields(1)%scalar%s%dof, .true.)
       call neko_scratch_registry%request_field(scalar_rho_cp, scalar_rho_cp_index, &
            simulation%scalars%scalar_fields(1)%scalar%s%dof, .true.)
       call field_col3(scalar_rho_cp, &
            simulation%scalars%scalar_fields(1)%scalar%rho, &
            simulation%scalars%scalar_fields(1)%scalar%cp)
       call field_copy(scalar_diff_field, &
            simulation%scalars%scalar_fields(1)%scalar%lambda_tot)
       call field_invcol2(scalar_diff_field, scalar_rho_cp)
       call this%scalar_diffusivity%init(n)
       call field_to_vector(this%scalar_diffusivity, scalar_diff_field)
       call rom_problem_build_scalar_transport_cache(this)
       if (this%petrov_galerkin) then
          call rom_assemble_scalar_transport_petrov(this%weights, this%scalar_metric_weight, &
               this%velocity_modes, this%scalar_modes, this%velocity_gradients, &
               this%scalar_gradients, this%test_scalar_modes, this%test_scalar_gradients, &
               this%lift_velocity, this%scalar_lift, this%scalar_lift_gradients, &
               this%scalar_diffusivity, linear_local, constant_local, quad_local, implicit_local)
       else
          call rom_assemble_scalar_transport(this%weights, this%scalar_metric_weight, &
               this%velocity_modes, this%scalar_modes, this%velocity_gradients, &
               this%scalar_gradients, this%lift_velocity, &
               this%scalar_lift, this%scalar_lift_gradients, this%scalar_diffusivity, &
               linear_local, constant_local, quad_local, implicit_local)
       end if
       if (this%forward_constrained_enriched) then
          call rom_assemble_scalar_transport_petrov(this%weights, this%scalar_metric_weight, &
               this%test_velocity_modes, this%test_scalar_modes, this%test_velocity_gradients, &
               this%test_scalar_gradients, this%test_scalar_modes, this%test_scalar_gradients, &
               this%lift_velocity, this%scalar_lift, this%scalar_lift_gradients, &
               this%scalar_diffusivity, enriched_linear_local, enriched_constant_local, &
               implicit_linear=enriched_implicit_local)
       end if
    end if
    call this%mass%init(r,r); call this%trial_mass%init(r,r)
    call this%fixed_linear%init(r,r)
    call this%fixed_implicit%init(r,r); call this%fixed_constant%init(r)
    ! The reduced operators are host-owned dense data.  Keep the host mirror
    ! authoritative instead of assigning the (unused) device mirrors.
    this%mass%x = mass_local%x
    this%trial_mass%x = trial_mass_local%x
    this%fixed_linear%x = linear_local%x
    this%fixed_implicit%x = implicit_local%x
    this%fixed_constant%x = constant_local%x
    if (this%forward_constrained_enriched) then
       call this%enriched_mass%init(e,e)
       call this%enriched_fixed_linear%init(e,e)
       call this%enriched_fixed_implicit%init(e,e)
       this%enriched_mass%x = enriched_mass_local%x
       this%enriched_fixed_linear%x = enriched_linear_local%x
       this%enriched_fixed_implicit%x = enriched_implicit_local%x
       this%enriched_adjoint_rank = e
    end if
    allocate(this%quadratic(r))
    do i = 1, r
       call this%quadratic(i)%init(r,r)
       this%quadratic(i)%x = quad_local(i)%x
    end do
    call this%viscous_constant%init(r); call this%viscous_linear%init(r)
    call this%viscous_quadratic%init(r,r); this%viscous_constant%x = 0.0_rp
    this%viscous_linear%x = 0.0_rp; this%viscous_quadratic%x = 0.0_rp
    ! Project Neko's canonical initial condition after subtracting the same
    ! spatial lift used to homogenise the streamed POD snapshots.  The fields
    ! are read from the reset Neko case itself rather than from a cached lift
    ! reference, so user-defined and scalar initial conditions are identical
    ! to those used by the forward solver.
    call neko_scratch_registry%request_vector(rhs_local, rhs_local_index, r, .true.)
    call neko_scratch_registry%request_vector(initial_difference, &
         initial_difference_index, n, .true.)
    rhs_local%x = 0.0_rp
    do c = 1, 3
       select case (c)
       case (1)
          call field_copy(field, simulation%fluid%u)
       case (2)
          call field_copy(field, simulation%fluid%v)
       case (3)
          call field_copy(field, simulation%fluid%w)
       end select
       call field_add2s2(field, this%lift%velocity(c), -1.0_rp, n)
       call field_to_vector(initial_difference, field)
       do i = 1, r
          rhs_local%x(i) = rhs_local%x(i) + rom_glsc3_local( &
               this%velocity_modes(c,i), this%weights, initial_difference)
       end do
    end do
    if (this%has_scalar) then
       call field_copy(scalar_field, simulation%scalars%scalar_fields(1)%scalar%s)
       call field_add2s2(scalar_field, this%lift%scalar, -1.0_rp, n)
       call field_to_vector(initial_difference, scalar_field)
       do i = 1, r
          rhs_local%x(i) = rhs_local%x(i) + this%scalar_metric_weight * &
               rom_glsc3_local(this%scalar_modes(i), this%weights, initial_difference)
       end do
    end if
    call rom_reduce_buffer(rhs_local%x)
    call this%initial_state%init(r)
    call this%rom%init(this%mass, this%fixed_constant, this%fixed_linear, &
         this%quadratic, this%fixed_implicit)
    this%rom%mass%x = this%trial_mass%x
    call this%rom%solve_mass(rhs_local, this%initial_state)
    this%rom%mass%x = this%mass%x
    this%rom_ready = .true.
    call neko_scratch_registry%relinquish_field([field_index, gx_index, gy_index, gz_index])
    if (this%has_scalar) call neko_scratch_registry%relinquish_field([ &
         scalar_field_index, scalar_gx_index, scalar_gy_index, scalar_gz_index, &
         scalar_diff_index, scalar_rho_cp_index])
    nullify(field); nullify(gx); nullify(gy); nullify(gz)
    nullify(scalar_field); nullify(scalar_gx); nullify(scalar_gy); nullify(scalar_gz)
    nullify(scalar_diff_field); nullify(scalar_rho_cp)
    call neko_scratch_registry%relinquish_matrix([mass_local_index, trial_mass_local_index, &
         linear_local_index, implicit_local_index])
    nullify(mass_local); nullify(linear_local); nullify(implicit_local)
    do i = 1, r; call quad_local(i)%free(); end do
    deallocate(quad_local)
    call neko_scratch_registry%relinquish_vector(rhs_local_index)
    call neko_scratch_registry%relinquish_vector([constant_local_index, &
         initial_difference_index])
    if (this%forward_constrained_enriched) then
       call neko_scratch_registry%relinquish_matrix([enriched_mass_local_index, &
            enriched_linear_local_index, enriched_implicit_local_index])
       call neko_scratch_registry%relinquish_vector(enriched_constant_local_index)
       nullify(enriched_mass_local); nullify(enriched_linear_local)
       nullify(enriched_implicit_local); nullify(enriched_constant_local)
    end if
    nullify(rhs_local); nullify(constant_local); nullify(initial_difference)
    if (this%basis_transform_adaptive) then
       call rom_problem_save_reference_basis(this)
       if (this%experimental_algebraic) then
          if (this%sensitivity_source == 'exact_algebraic') then
             call neko_log%warning('EXPERIMENTAL ALGEBRAIC ADAPTIVE ROM selected: '&
                  //'candidate modes are not pressure projected; exact_algebraic '&
                  //'sensitivity is UNVERIFIED and must not be treated as accepted')
          else
             call neko_log%warning('EXPERIMENTAL ALGEBRAIC ADAPTIVE ROM selected: '&
                  //'candidate velocity modes are not pressure projected; '&
                  //'sensitivity_source=lagged_basis is an approximate search direction')
          end if
       else if (this%basis_transform == 'vector_potential') then
          call neko_log%message('vector-potential ROM basis saved; candidates use the '// &
               'native gathered curl and homogeneous velocity boundary projector; '// &
               'discrete continuity residuals are measured and not presumed zero')
       else
          call neko_log%warning('design-dependent ROM reference basis saved; '&
               //'candidate admissibility projection and exact mode pullback remain '&
               //'a guarded scientific gate')
       end if
    end if
  end subroutine rom_problem_assemble_fixed

  !> Resolve objective and constraint point-zone masks on the SEM mesh.
  subroutine rom_problem_setup_masks(this, simulation)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    integer :: n
    this%scalar_mixing_cache_valid = .false.
    n = simulation%fluid%u%dof%size()
    call get_zone_mask(this%viscous_mask_name, n, this%viscous_mask)
    call get_zone_mask(this%brinkman_mask_name, n, this%brinkman_mask)
    call get_zone_mask(this%scalar_mixing_mask_name, n, this%scalar_mixing_mask)
    this%viscous_domain = vector_glsc2(this%weights, this%viscous_mask)
    this%brinkman_domain = vector_glsc2(this%weights, this%brinkman_mask)
    this%scalar_mixing_domain = vector_glsc2(this%weights, this%scalar_mixing_mask)
    ! vector_glsc2 performs the global MPI reduction internally.
    if (this%has_volume_constraint) then
       call get_zone_mask(this%volume_mask_name, n, this%volume_mask)
       this%volume_domain = vector_glsc2(this%weights, this%volume_mask)
       if (this%volume_domain <= 0.0_rp) call neko_error('ROM volume mask has zero measure')
    else
       call this%volume_mask%free()
       this%volume_domain = 0.0_rp
    end if
    if (this%has_viscous .and. this%viscous_domain <= 0.0_rp) &
         call neko_error('ROM viscous objective mask has zero measure')
    if (this%has_brinkman .and. this%brinkman_domain <= 0.0_rp) &
         call neko_error('ROM Brinkman objective mask has zero measure')
    if (this%has_scalar_mixing .and. this%scalar_mixing_domain <= 0.0_rp) &
         call neko_error('ROM scalar mixing objective mask has zero measure')
  end subroutine rom_problem_setup_masks

  !> Cache scalar metric and diffusion operators for adaptive CPU candidates.
  !! These terms depend on the accepted scalar basis, scalar metric, diffusivity
  !! and scalar lift, all of which are fixed throughout the inner candidate
  !! loop. Candidate scalar transport assembly then computes advection only.
  subroutine rom_problem_build_scalar_transport_cache(this)
    class(rom_problem_t), intent(inout) :: this
    real(kind=rp) :: timer_start

    this%scalar_transport_cache_valid = .false.
    this%scalar_transport_cache_forward_rank = 0
    this%scalar_transport_cache_adjoint_rank = 0
    this%scalar_transport_cache_build_seconds = 0.0_rp
    call this%scalar_mass_cache%free()
    call this%scalar_adjoint_mass_cache%free()
    call this%scalar_diffusion_cache%free()
    call this%scalar_adjoint_diffusion_cache%free()
    call this%scalar_diffusion_lift_cache%free()
    call this%scalar_adjoint_diffusion_lift_cache%free()
    if (.not. this%has_scalar .or. .not. this%basis_transform_adaptive .or. &
         NEKO_BCKND_DEVICE /= 0) return

    timer_start = real(MPI_WTIME(), rp)
    call rom_assemble_scalar_mass(this%weights, this%scalar_metric_weight, &
         this%scalar_modes, this%scalar_modes, this%scalar_mass_cache)
    call rom_assemble_scalar_diffusion_petrov(this%weights, &
         this%scalar_metric_weight, this%scalar_diffusivity, &
         this%scalar_gradients, this%scalar_gradients, &
         this%scalar_lift_gradients, this%scalar_diffusion_cache, &
         this%scalar_diffusion_lift_cache)
    if (this%forward_constrained_enriched) then
       call rom_assemble_scalar_mass(this%weights, this%scalar_metric_weight, &
            this%test_scalar_modes, this%test_scalar_modes, &
            this%scalar_adjoint_mass_cache)
       call rom_assemble_scalar_diffusion_petrov(this%weights, &
            this%scalar_metric_weight, this%scalar_diffusivity, &
            this%test_scalar_gradients, this%test_scalar_gradients, &
            this%scalar_lift_gradients, this%scalar_adjoint_diffusion_cache, &
            this%scalar_adjoint_diffusion_lift_cache)
       this%scalar_transport_cache_adjoint_rank = size(this%test_scalar_modes)
    end if
    this%scalar_transport_cache_forward_rank = size(this%scalar_modes)
    this%scalar_transport_cache_valid = .true.
    this%scalar_transport_cache_build_seconds = &
         real(MPI_WTIME(), rp) - timer_start
  end subroutine rom_problem_build_scalar_transport_cache

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
    if (NEKO_BCKND_DEVICE .eq. 1) then
       ! Point-zone masks are assembled from host index lists once per basis
       ! refresh, not in the ROM inner loop.
       call mask%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if
  end subroutine get_zone_mask

  !> Assemble reduced viscous and design-dependent objective operators.
  subroutine rom_problem_assemble_objective_operators(this)
    class(rom_problem_t), intent(inout) :: this
    integer :: r
    real(kind=rp) :: timer_start
    logical :: reuse_scalar
    timer_start = real(MPI_WTIME(),rp)
    r = this%mass%get_nrows()
    reuse_scalar = (NEKO_BCKND_DEVICE .eq. 0) .and. &
         this%scalar_mixing_cache_valid .and. &
         this%scalar_mixing_cache_forward_rank == r .and. &
         this%scalar_mixing_cache_adjoint_rank == merge(this%test_basis%n_modes, 0, &
         this%forward_constrained_enriched)
    this%scalar_mixing_cache_hit = merge(1.0_rp, 0.0_rp, reuse_scalar)
    call this%viscous_constant%init(1)
    call this%viscous_linear%init(r)
    call this%viscous_quadratic%init(r,r)
    this%viscous_constant%x = 0.0_rp
    this%viscous_linear%x = 0.0_rp
    this%viscous_quadratic%x = 0.0_rp
    if (.not. reuse_scalar) then
       call this%scalar_mixing_constant%free()
       call this%scalar_mixing_linear%free()
       call this%scalar_mixing_quadratic%free()
       call this%scalar_mixing_constant%init(1)
       call this%scalar_mixing_linear%init(r)
       call this%scalar_mixing_quadratic%init(r,r)
       this%scalar_mixing_constant%x = 0.0_rp
       this%scalar_mixing_linear%x = 0.0_rp
       this%scalar_mixing_quadratic%x = 0.0_rp
    end if
    if (this%has_viscous) then
       call rom_assemble_viscous_terms(this%weights, this%viscous_mask, &
            this%viscosity, this%velocity_gradients, this%lift_gradients, &
            this%viscous_constant, this%viscous_linear, this%viscous_quadratic)
    end if
    if (this%has_scalar_mixing .and. .not. reuse_scalar) then
       call rom_assemble_scalar_mixing_terms(this%weights, this%scalar_mixing_mask, &
            this%scalar_lift, this%scalar_modes, this%scalar_reference_value, &
            this%scalar_mixing_constant, this%scalar_mixing_linear, &
            this%scalar_mixing_quadratic)
    end if
    if (this%forward_constrained_enriched) then
       call this%enriched_viscous_linear%free()
       call this%enriched_viscous_cross%free()
       call this%enriched_viscous_linear%init(this%test_basis%n_modes)
       call this%enriched_viscous_cross%init(this%test_basis%n_modes, &
            this%mass%get_nrows())
       if (.not. reuse_scalar) then
          call this%enriched_scalar_mixing_linear%free()
          call this%enriched_scalar_mixing_cross%free()
          call this%enriched_scalar_mixing_linear%init(this%test_basis%n_modes)
          call this%enriched_scalar_mixing_cross%init(this%test_basis%n_modes, &
               this%mass%get_nrows())
       end if
       if (this%has_viscous) then
          call rom_assemble_viscous_cross_terms(this%weights, this%viscous_mask, &
               this%viscosity, this%test_velocity_gradients, this%velocity_gradients, &
               this%lift_gradients, this%enriched_viscous_linear, &
               this%enriched_viscous_cross)
       else
          this%enriched_viscous_linear%x = 0.0_rp
          this%enriched_viscous_cross%x = 0.0_rp
       end if
       if (this%has_scalar_mixing .and. .not. reuse_scalar) then
          call rom_assemble_scalar_mixing_cross_terms(this%weights, &
               this%scalar_mixing_mask, this%scalar_lift, this%scalar_modes, &
               this%test_scalar_modes, this%scalar_reference_value, &
               this%enriched_scalar_mixing_linear, this%enriched_scalar_mixing_cross)
       else if (.not. reuse_scalar) then
          this%enriched_scalar_mixing_linear%x = 0.0_rp
          this%enriched_scalar_mixing_cross%x = 0.0_rp
       end if
    end if
    this%scalar_mixing_cache_valid = NEKO_BCKND_DEVICE .eq. 0
    this%scalar_mixing_cache_forward_rank = r
    this%scalar_mixing_cache_adjoint_rank = merge(this%test_basis%n_modes, 0, &
         this%forward_constrained_enriched)
    ! rom_assemble_viscous_terms is globally assembled through vector_glsc3.
    this%objective_operator_assembly_seconds = &
         this%objective_operator_assembly_seconds + real(MPI_WTIME(),rp) - timer_start
  end subroutine rom_problem_assemble_objective_operators

  !> Pair and truncate the Petrov trial and test spaces by cross-mass SVD.
  !! For the currently selected equal-rank input spaces this forms
  !! H = W_a^T M V_f, computes H = U Sigma Z^T, and applies V_f <- V_f Z_k
  !! and W_a <- W_a U_k to the complete device-resident mode fields before
  !! any gradients or reduced operators are assembled.  No regularisation is
  !! applied to a weakly overlapping direction: it is removed instead.
  subroutine rom_problem_check_cross_mass(this)
    class(rom_problem_t), intent(inout) :: this
    type(matrix_t) :: cross_mass
    real(kind=rp), allocatable :: h(:, :), singular(:), u(:, :), vt(:, :), work(:)
    type(vector_t), allocatable :: velocity_new(:, :), test_velocity_new(:, :)
    type(vector_t), allocatable :: scalar_new(:), test_scalar_new(:)
    integer :: r, k, info, lwork, i, j, c, n, offset
    real(kind=rp) :: threshold

    if (.not. allocated(this%velocity_modes) .or. &
         .not. allocated(this%test_velocity_modes)) call neko_error( &
         'Petrov cross mass requested before trial/test modes were loaded')
    r = size(this%velocity_modes, 2)
    if (r < 1 .or. size(this%test_velocity_modes, 2) /= r) call neko_error( &
         'Petrov cross mass requires equal trial/test ranks')
    n = this%weights%size()
    allocate(h(r,r), singular(r), u(r,r), vt(r,r))
    if (this%has_scalar) then
       call rom_assemble_cross_mass(this%weights, this%velocity_modes, &
            this%test_velocity_modes, cross_mass, this%scalar_modes, &
            this%test_scalar_modes, this%scalar_metric_weight)
    else
       call rom_assemble_cross_mass(this%weights, this%velocity_modes, &
            this%test_velocity_modes, cross_mass)
    end if
    h = cross_mass%x
    call cross_mass%free()
    lwork = max(1, 5*r)
    allocate(work(lwork))
    call dgesvd('A', 'A', r, r, h, r, singular, u, r, vt, r, work, lwork, info)
    if (info /= 0) call neko_error('LAPACK failed while pairing Petrov cross mass')
    if (singular(1) <= 0.0_rp) call neko_error( &
         'Petrov cross mass has no nonzero singular direction')
    threshold = this%petrov_svd_relative_tolerance * singular(1)
    k = count(singular >= threshold)
    if (k < 1) call neko_error( &
         'Petrov cross mass has no singular direction above its relative threshold')

    if (allocated(this%cross_mass_singular_values)) &
         deallocate(this%cross_mass_singular_values)
    allocate(this%cross_mass_singular_values(r))
    this%cross_mass_singular_values = singular
    this%cross_mass_retained_rank = k
    this%cross_mass_sigma_max = singular(1)
    this%cross_mass_sigma_min = singular(k)
    this%cross_mass_condition = singular(1) / max(singular(k), tiny(1.0_rp))
    call neko_log%message('ROM Petrov cross-mass SVD: relative threshold='// &
         trim(real_to_string(threshold))//', retained rank='// &
         trim(int_to_string(k))//'/'//trim(int_to_string(r))//', condition='// &
         trim(real_to_string(this%cross_mass_condition)))
    do i = 1, r
       call neko_log%message('ROM Petrov cross-mass singular value '// &
            trim(int_to_string(i))//'='//trim(real_to_string(singular(i))))
    end do

    allocate(velocity_new(3,k), test_velocity_new(3,k))
    do j = 1, k
       do c = 1, 3
          call velocity_new(c,j)%init(n)
          call test_velocity_new(c,j)%init(n)
          call vector_rzero(velocity_new(c,j))
          call vector_rzero(test_velocity_new(c,j))
          do i = 1, r
             call vector_add2s2(velocity_new(c,j), this%velocity_modes(c,i), &
                  vt(j,i))
             call vector_add2s2(test_velocity_new(c,j), this%test_velocity_modes(c,i), &
                  u(i,j))
          end do
       end do
    end do
    if (this%has_scalar) then
       allocate(scalar_new(k), test_scalar_new(k))
       do j = 1, k
          call scalar_new(j)%init(n)
          call test_scalar_new(j)%init(n)
          call vector_rzero(scalar_new(j))
          call vector_rzero(test_scalar_new(j))
          do i = 1, r
             call vector_add2s2(scalar_new(j), this%scalar_modes(i), vt(j,i))
             call vector_add2s2(test_scalar_new(j), this%test_scalar_modes(i), u(i,j))
          end do
       end do
    end if

    do i = 1, size(this%velocity_modes,1)
       do j = 1, size(this%velocity_modes,2)
          call this%velocity_modes(i,j)%free()
          call this%test_velocity_modes(i,j)%free()
       end do
    end do
    deallocate(this%velocity_modes, this%test_velocity_modes)
    call move_alloc(velocity_new, this%velocity_modes)
    call move_alloc(test_velocity_new, this%test_velocity_modes)
    if (this%has_scalar) then
       do i = 1, size(this%scalar_modes)
          call this%scalar_modes(i)%free()
          call this%test_scalar_modes(i)%free()
       end do
       deallocate(this%scalar_modes, this%test_scalar_modes)
       call move_alloc(scalar_new, this%scalar_modes)
       call move_alloc(test_scalar_new, this%test_scalar_modes)
    end if
    this%basis%n_modes = k
    this%test_basis%n_modes = k

    ! Keep the owned host basis metadata consistent with the vectors used for
    ! assembly.  This is a basis-refresh boundary, so the readback is not part
    ! of the reduced inner optimisation loop.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       do j = 1, k
          do c = 1, 3
             call this%velocity_modes(c,j)%copy_from(DEVICE_TO_HOST, sync = .true.)
             call this%test_velocity_modes(c,j)%copy_from(DEVICE_TO_HOST, sync = .true.)
          end do
          if (this%has_scalar) then
             call this%scalar_modes(j)%copy_from(DEVICE_TO_HOST, sync = .true.)
             call this%test_scalar_modes(j)%copy_from(DEVICE_TO_HOST, sync = .true.)
          end if
       end do
    end if
    call this%basis%modes%free()
    call this%basis%modes%init(this%basis%n_fields*n, k)
    call this%test_basis%modes%free()
    call this%test_basis%modes%init(this%test_basis%n_fields*n, k)
    do j = 1, k
       do c = 1, 3
          offset = (c-1)*n
          this%basis%modes%x(offset+1:offset+n,j) = this%velocity_modes(c,j)%x
          this%test_basis%modes%x(offset+1:offset+n,j) = &
               this%test_velocity_modes(c,j)%x
       end do
       if (this%has_scalar) then
          this%basis%modes%x(3*n+1:4*n,j) = this%scalar_modes(j)%x
          this%test_basis%modes%x(3*n+1:4*n,j) = this%test_scalar_modes(j)%x
       end if
    end do
    deallocate(h, singular, u, vt, work)
  end subroutine rom_problem_check_cross_mass

  !> Format a reduced diagnostic scalar without introducing a mesh-sized I/O
  !! dependency into the Petrov assembly path.
  function real_to_string(value) result(string)
    real(kind=rp), intent(in) :: value
    character(len=32) :: string
    write(string, '(ES12.5)') value
  end function real_to_string

  !> Format a reduced diagnostic integer without introducing an I/O dependency.
  function int_to_string(value) result(string)
    integer, intent(in) :: value
    character(len=32) :: string
    write(string, '(I0)') value
  end function int_to_string

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
    type(vector_t), pointer :: difference_u, difference_v, difference_w
    integer :: indices(10), i, n
    integer :: difference_indices(3)
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
    call neko_scratch_registry%request_vector(difference_u, difference_indices(1), n, .true.)
    call neko_scratch_registry%request_vector(difference_v, difference_indices(2), n, .true.)
    call neko_scratch_registry%request_vector(difference_w, difference_indices(3), n, .true.)
    call field_copy(candidate_u, this%lift%velocity(1))
    call field_copy(candidate_v, this%lift%velocity(2))
    call field_copy(candidate_w, this%lift%velocity(3))
       do i = 1, size(coefficients)
          call vector_to_field(mode_field, this%velocity_modes(1,i))
          call field_add2s2(candidate_u, mode_field, coefficients(i), n)
          call vector_to_field(mode_field, this%velocity_modes(2,i))
          call field_add2s2(candidate_v, mode_field, coefficients(i), n)
          call vector_to_field(mode_field, this%velocity_modes(3,i))
          call field_add2s2(candidate_w, mode_field, coefficients(i), n)
    end do
    call field_copy(corrected_u, candidate_u)
    call field_copy(corrected_v, candidate_v)
    call field_copy(corrected_w, candidate_w)
    call simulation%fluid%bcs_vel%apply(corrected_u, corrected_v, corrected_w, &
         time = simulation%neko_case%time, strong = .true.)
    call field_sub3(corrected_u, corrected_u, candidate_u, n)
    call field_sub3(corrected_v, corrected_v, candidate_v, n)
    call field_sub3(corrected_w, corrected_w, candidate_w, n)
    call field_to_vector(difference_u, corrected_u)
    call field_to_vector(difference_v, corrected_v)
    call field_to_vector(difference_w, corrected_w)
    call vector_absval(difference_u)
    call vector_absval(difference_v)
    call vector_absval(difference_w)
    local_residual = max(vector_glmax(difference_u), max(vector_glmax(difference_v), &
         vector_glmax(difference_w)))
    ! vector_glmax already performs the global Neko reduction.
    residual = local_residual
    call neko_scratch_registry%relinquish_field(indices)
    call neko_scratch_registry%relinquish_vector(difference_indices)
    nullify(candidate_u); nullify(candidate_v); nullify(candidate_w)
    nullify(corrected_u); nullify(corrected_v); nullify(corrected_w); nullify(mode_field)
    nullify(difference_u); nullify(difference_v); nullify(difference_w)
  end subroutine rom_problem_verify_dirichlet_reconstruction

  !> Map one design and assemble its reduced Brinkman operators.
  subroutine rom_problem_assemble_candidate(this, design, brinkman, brinkman_lift, &
       objective_brinkman, objective_lift, objective_constant, enriched_brinkman, &
       enriched_objective_cross, enriched_objective_lift, volume)
    class(rom_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design
    type(matrix_t), intent(inout) :: brinkman, objective_brinkman
    type(vector_t), intent(inout) :: brinkman_lift, objective_lift
    type(matrix_t), intent(inout) :: enriched_brinkman, enriched_objective_cross
    type(vector_t), intent(inout) :: enriched_objective_lift
    real(kind=rp), intent(out) :: objective_constant
    real(kind=rp), intent(out), optional :: volume
    type(vector_t), pointer :: mapped, raw_design, local_lift
    type(matrix_t), pointer :: local
    integer :: mapped_index, raw_design_index, local_index, local_lift_index
    integer :: n, r
    real(kind=rp) :: timer_start

    select type (des => design)
    type is (brinkman_design_t)
       n = des%size()
       if (n /= this%weights%size()) call neko_error('ROM/design local size mismatch')
       call neko_scratch_registry%request_vector(mapped, mapped_index, n, .true.)
       call des%get_mapped_values(mapped)
    class default
       call neko_error('rom_mma requires a Brinkman design')
    end select
    if (this%basis_transform_adaptive .and. .not. this%adaptive_modes_frozen) then
       call rom_problem_update_adaptive_modes(this, mapped)
       call rom_problem_reassemble_adaptive_operators(this)
       call rom_problem_assemble_objective_operators(this)
    end if
    r = this%mass%get_nrows()
    call neko_scratch_registry%request_matrix(local, local_index, r, r, .true.)
    call neko_scratch_registry%request_vector(local_lift, local_lift_index, r, .true.)
    timer_start = real(MPI_WTIME(),rp)
    if (this%petrov_galerkin) then
       call rom_assemble_brinkman_petrov(this%weights, mapped, this%velocity_modes, &
            this%test_velocity_modes, this%lift_velocity, local, local_lift)
    else
       call rom_assemble_brinkman(this%weights, mapped, this%velocity_modes, &
            this%lift_velocity, local, local_lift)
    end if
    this%forward_local_operator_assembly_seconds = &
         this%forward_local_operator_assembly_seconds + real(MPI_WTIME(),rp) - timer_start
    call brinkman%init(r,r); call brinkman_lift%init(r)
    brinkman%x = local%x
    brinkman_lift%x = local_lift%x
    timer_start = real(MPI_WTIME(),rp)
    call rom_assemble_brinkman_objective(this%weights, this%brinkman_mask, mapped, &
         this%velocity_modes, this%lift_velocity, local, local_lift, objective_constant)
    this%objective_operator_assembly_seconds = &
         this%objective_operator_assembly_seconds + real(MPI_WTIME(),rp) - timer_start
    call objective_brinkman%init(r,r); call objective_lift%init(r)
    objective_brinkman%x = local%x
    objective_lift%x = local_lift%x
    if (this%forward_constrained_enriched) then
       timer_start = real(MPI_WTIME(),rp)
       call rom_assemble_brinkman_petrov(this%weights, mapped, this%test_velocity_modes, &
            this%test_velocity_modes, this%lift_velocity, enriched_brinkman, &
            enriched_objective_lift)
       this%adjoint_local_operator_assembly_seconds = &
            this%adjoint_local_operator_assembly_seconds + &
            real(MPI_WTIME(),rp) - timer_start
       timer_start = real(MPI_WTIME(),rp)
       call rom_assemble_brinkman_objective_cross(this%weights, this%brinkman_mask, &
            mapped, this%velocity_modes, this%test_velocity_modes, this%lift_velocity, &
            enriched_objective_cross, enriched_objective_lift)
       this%objective_operator_assembly_seconds = &
            this%objective_operator_assembly_seconds + &
            real(MPI_WTIME(),rp) - timer_start
    end if
    if (present(volume)) then
       ! The production volume constraint integrates the unfiltered design
       ! variable, not the mapped Brinkman coefficient used by the PDE.
       if (this%has_volume_constraint) then
          call neko_scratch_registry%request_vector(raw_design, raw_design_index, &
               n, .true.)
          call design%get_values(raw_design)
          volume = vector_glsc3(this%weights, this%volume_mask, raw_design)
          call neko_scratch_registry%relinquish_vector(raw_design_index)
          nullify(raw_design)
       else
          volume = 0.0_rp
       end if
    end if
    call neko_scratch_registry%relinquish_vector(mapped_index)
    call neko_scratch_registry%relinquish_matrix(local_index)
    call neko_scratch_registry%relinquish_vector(local_lift_index)
    nullify(mapped); nullify(raw_design); nullify(local); nullify(local_lift)
  end subroutine rom_problem_assemble_candidate

  !> Reassemble all reduced Galerkin operators after a candidate mode update.
  !! Dense reduced arrays are rebuilt on the host from Neko's backend-aware
  !! contractions; mesh-sized mode and gradient data remain in vector_t
  !! containers. Independent-adjoint configurations refresh their E-space
  !! tangent operators separately from the F-space forward model.
  subroutine rom_problem_reassemble_adaptive_operators(this)
    class(rom_problem_t), intent(inout) :: this
    type(matrix_t), pointer :: mass_local, trial_mass_local, linear_local, implicit_local
    type(matrix_t), pointer :: e_mass_local, e_linear_local, e_implicit_local
    type(matrix_t), allocatable :: quad_local(:)
    type(vector_t), pointer :: constant_local, rhs_local, initial_difference
    type(vector_t), pointer :: e_constant_local
    type(field_t), pointer :: field, scalar_field
    integer :: mass_index, trial_mass_index, linear_index, implicit_index
    integer :: constant_index, rhs_index, initial_difference_index
    integer :: field_index, scalar_field_index
    real(kind=rp), allocatable :: mass_svd(:, :), mass_singular(:), mass_u(:, :), &
         mass_vt(:, :), mass_work(:)
    integer :: n, r, e, i, c, d
    integer :: e_mass_index, e_linear_index, e_implicit_index, e_constant_index
    integer :: svd_info, svd_lwork
    logical :: use_scalar_cache, use_scalar_adjoint_cache
    real(kind=rp) :: nu, timer_start

    if (.not. associated(this%simulation)) call neko_error( &
         'adaptive ROM operator reassembly has no simulation')
    if (this%petrov_galerkin) call neko_error( &
         'adaptive algebraic reassembly does not support Petrov forward dynamics')
    ! All candidate velocity modes and their gradients may have changed.
    this%projected_tensor_cache_valid = .false.
    n = this%weights%size()
    r = size(this%velocity_modes, 2)
    e = size(this%test_velocity_modes, 2)
    nu = this%viscosity / this%density
    use_scalar_cache = .false.
    use_scalar_adjoint_cache = .false.
    if (this%has_scalar .and. NEKO_BCKND_DEVICE == 0 .and. &
         this%scalar_transport_cache_valid) then
       if (allocated(this%scalar_modes)) then
          if (this%scalar_transport_cache_forward_rank == size(this%scalar_modes)) &
               use_scalar_cache = .true.
       end if
    end if
    if (use_scalar_cache .and. this%forward_constrained_enriched) then
       if (allocated(this%test_scalar_modes)) then
          if (this%scalar_transport_cache_adjoint_rank == size(this%test_scalar_modes)) &
               use_scalar_adjoint_cache = .true.
       end if
    end if

    call this%mass%free(); call this%trial_mass%free()
    call this%fixed_linear%free(); call this%fixed_implicit%free()
    call this%fixed_constant%free()
    if (allocated(this%quadratic)) then
       do i = 1, size(this%quadratic); call this%quadratic(i)%free(); end do
       deallocate(this%quadratic)
    end if
    call neko_scratch_registry%request_matrix(mass_local, mass_index, r, r, .true.)
    call neko_scratch_registry%request_matrix(trial_mass_local, trial_mass_index, r, r, .true.)
    call neko_scratch_registry%request_matrix(linear_local, linear_index, r, r, .true.)
    call neko_scratch_registry%request_matrix(implicit_local, implicit_index, r, r, .true.)
    call neko_scratch_registry%request_vector(constant_local, constant_index, r, .true.)
    allocate(quad_local(r))

    timer_start = real(MPI_WTIME(),rp)
    if (this%has_scalar) then
       if (use_scalar_cache) then
          call rom_assemble_mass(this%weights, this%velocity_modes, trial_mass_local)
          trial_mass_local%x = trial_mass_local%x + this%scalar_mass_cache%x
       else
          call rom_assemble_mass_coupled(this%weights, this%velocity_modes, &
               this%scalar_modes, this%scalar_metric_weight, trial_mass_local)
       end if
    else
       call rom_assemble_mass(this%weights, this%velocity_modes, trial_mass_local)
    end if
    mass_local%x = trial_mass_local%x
    ! Retain an explicit condition estimate for each adapted candidate.  This
    ! is a diagnostic and rejects only a numerically singular Gram matrix.
    allocate(mass_svd(r,r), mass_singular(r), mass_u(r,r), mass_vt(r,r))
    mass_svd = trial_mass_local%x
    svd_lwork = max(1,5*r)
    allocate(mass_work(svd_lwork))
    call dgesvd('A', 'A', r, r, mass_svd, r, mass_singular, mass_u, r, &
         mass_vt, r, mass_work, svd_lwork, svd_info)
    if (svd_info /= 0 .or. mass_singular(1) <= 0.0_rp .or. &
         mass_singular(r) <= 100.0_rp * epsilon(1.0_rp) * mass_singular(1)) &
         call neko_error( &
         'adaptive ROM candidate mass matrix is singular')
    this%adaptive_mass_condition = mass_singular(1) / mass_singular(r)
    deallocate(mass_svd, mass_singular, mass_u, mass_vt, mass_work)
    call rom_assemble_lift_convection(this%weights, this%velocity_modes, &
         this%velocity_gradients, this%lift_velocity, this%lift_gradients, &
         linear_local, constant_local)
    call rom_assemble_diffusion(this%weights, nu, this%velocity_gradients, &
         implicit_local, this%lift_gradients, constant_local)
    call rom_assemble_convection(this%weights, this%velocity_modes, &
         this%velocity_gradients, quad_local)
    do i = 1, r; quad_local(i)%x = -quad_local(i)%x; end do

    if (this%has_scalar) then
       if (use_scalar_cache) then
          call rom_assemble_scalar_transport_petrov(this%weights, &
               this%scalar_metric_weight, this%velocity_modes, this%scalar_modes, &
               this%velocity_gradients, this%scalar_gradients, this%scalar_modes, &
               this%scalar_gradients, this%lift_velocity, this%scalar_lift, &
               this%scalar_lift_gradients, this%scalar_diffusivity, linear_local, &
               constant_local, quad_local, implicit_local, include_diffusion=.false.)
          implicit_local%x = implicit_local%x + this%scalar_diffusion_cache%x
          constant_local%x = constant_local%x + this%scalar_diffusion_lift_cache%x
          this%scalar_transport_cache_hit = 1.0_rp
       else
          call rom_assemble_scalar_transport(this%weights, this%scalar_metric_weight, &
               this%velocity_modes, this%scalar_modes, this%velocity_gradients, &
               this%scalar_gradients, this%lift_velocity, this%scalar_lift, &
               this%scalar_lift_gradients, this%scalar_diffusivity, linear_local, &
               constant_local, quad_local, implicit_local)
       end if
    end if
    this%forward_local_operator_assembly_seconds = &
         this%forward_local_operator_assembly_seconds + &
         real(MPI_WTIME(),rp) - timer_start

    if (this%forward_constrained_enriched) then
       call neko_scratch_registry%request_matrix(e_mass_local, e_mass_index, e, e, .true.)
       call neko_scratch_registry%request_matrix(e_linear_local, e_linear_index, e, e, .true.)
       call neko_scratch_registry%request_matrix(e_implicit_local, e_implicit_index, e, e, .true.)
       call neko_scratch_registry%request_vector(e_constant_local, e_constant_index, e, .true.)
       timer_start = real(MPI_WTIME(),rp)
       if (this%has_scalar) then
          if (use_scalar_adjoint_cache) then
             call rom_assemble_mass(this%weights, this%test_velocity_modes, e_mass_local)
             e_mass_local%x = e_mass_local%x + this%scalar_adjoint_mass_cache%x
          else
             call rom_assemble_mass_coupled(this%weights, this%test_velocity_modes, &
                  this%test_scalar_modes, this%scalar_metric_weight, e_mass_local)
          end if
       else
          call rom_assemble_mass(this%weights, this%test_velocity_modes, e_mass_local)
       end if
       call rom_assemble_lift_convection_petrov(this%weights, this%test_velocity_modes, &
            this%test_velocity_gradients, this%test_velocity_modes, this%lift_velocity, &
            this%lift_gradients, e_linear_local, e_constant_local)
       call rom_assemble_diffusion_petrov(this%weights, nu, this%test_velocity_gradients, &
            this%test_velocity_gradients, e_implicit_local, this%lift_gradients, &
            e_constant_local)
       if (this%has_scalar) then
          if (use_scalar_adjoint_cache) then
             call rom_assemble_scalar_transport_petrov(this%weights, &
                  this%scalar_metric_weight, this%test_velocity_modes, &
                  this%test_scalar_modes, this%test_velocity_gradients, &
                  this%test_scalar_gradients, this%test_scalar_modes, &
                  this%test_scalar_gradients, this%lift_velocity, this%scalar_lift, &
                  this%scalar_lift_gradients, this%scalar_diffusivity, e_linear_local, &
                  e_constant_local, implicit_linear=e_implicit_local, &
                  include_diffusion=.false.)
             e_implicit_local%x = e_implicit_local%x + &
                  this%scalar_adjoint_diffusion_cache%x
             e_constant_local%x = e_constant_local%x + &
                  this%scalar_adjoint_diffusion_lift_cache%x
             this%scalar_transport_cache_hit = 1.0_rp
          else
             call rom_assemble_scalar_transport_petrov(this%weights, &
                  this%scalar_metric_weight, this%test_velocity_modes, &
                  this%test_scalar_modes, this%test_velocity_gradients, &
                  this%test_scalar_gradients, this%test_scalar_modes, &
                  this%test_scalar_gradients, this%lift_velocity, this%scalar_lift, &
                  this%scalar_lift_gradients, this%scalar_diffusivity, &
                  e_linear_local, e_constant_local, implicit_linear=e_implicit_local)
          end if
       end if
       this%adjoint_local_operator_assembly_seconds = &
            this%adjoint_local_operator_assembly_seconds + &
            real(MPI_WTIME(),rp) - timer_start
       call this%enriched_mass%free()
       call this%enriched_fixed_linear%free()
       call this%enriched_fixed_implicit%free()
       call this%enriched_mass%init(e,e)
       call this%enriched_fixed_linear%init(e,e)
       call this%enriched_fixed_implicit%init(e,e)
       this%enriched_mass%x = e_mass_local%x
       this%enriched_fixed_linear%x = e_linear_local%x
       this%enriched_fixed_implicit%x = e_implicit_local%x
       this%enriched_adjoint_rank = e
    end if

    call this%mass%init(r,r); call this%trial_mass%init(r,r)
    call this%fixed_linear%init(r,r); call this%fixed_implicit%init(r,r)
    call this%fixed_constant%init(r)
    this%mass%x = mass_local%x; this%trial_mass%x = trial_mass_local%x
    this%fixed_linear%x = linear_local%x
    this%fixed_implicit%x = implicit_local%x
    this%fixed_constant%x = constant_local%x
    allocate(this%quadratic(r))
    do i = 1, r
       call this%quadratic(i)%init(r,r)
       this%quadratic(i)%x = quad_local(i)%x
    end do

    ! Reproject the canonical Neko initial state using the candidate Gram
    ! matrix.  The live simulation is reset by the surrounding lifecycle, so
    ! this reads the same prescribed initial condition as the reference path.
    timer_start = real(MPI_WTIME(),rp)
    call neko_scratch_registry%request_vector(rhs_local, rhs_index, r, .true.)
    call neko_scratch_registry%request_vector(initial_difference, &
         initial_difference_index, n, .true.)
    call neko_scratch_registry%request_field(field, field_index, &
         this%simulation%fluid%u%dof, .true.)
    rhs_local%x = 0.0_rp
    do c = 1, 3
       select case (c)
       case (1); call field_copy(field, this%simulation%fluid%u)
       case (2); call field_copy(field, this%simulation%fluid%v)
       case (3); call field_copy(field, this%simulation%fluid%w)
       end select
       call field_add2s2(field, this%lift%velocity(c), -1.0_rp, n)
       call field_to_vector(initial_difference, field)
       do i = 1, r
          rhs_local%x(i) = rhs_local%x(i) + rom_glsc3_local( &
               this%velocity_modes(c,i), this%weights, initial_difference)
       end do
    end do
    if (this%has_scalar) then
       call neko_scratch_registry%request_field(scalar_field, scalar_field_index, &
            this%simulation%scalars%scalar_fields(1)%scalar%s%dof, .true.)
       call field_copy(scalar_field, this%simulation%scalars%scalar_fields(1)%scalar%s)
       call field_add2s2(scalar_field, this%lift%scalar, -1.0_rp, n)
       call field_to_vector(initial_difference, scalar_field)
       do i = 1, r
          rhs_local%x(i) = rhs_local%x(i) + this%scalar_metric_weight * &
               rom_glsc3_local(this%scalar_modes(i), this%weights, initial_difference)
       end do
    end if
    call rom_reduce_buffer(rhs_local%x)
    call this%initial_state%free(); call this%initial_state%init(r)
    call this%rom%free()
    call this%rom%init(this%mass, this%fixed_constant, this%fixed_linear, &
         this%quadratic, this%fixed_implicit)
    this%rom%mass%x = this%trial_mass%x
    call this%rom%solve_mass(rhs_local, this%initial_state)
    this%rom%mass%x = this%mass%x
    this%initial_state_projection_seconds = this%initial_state_projection_seconds + &
         real(MPI_WTIME(),rp) - timer_start

    call neko_scratch_registry%relinquish_field(field_index)
    if (this%has_scalar) call neko_scratch_registry%relinquish_field(scalar_field_index)
    call neko_scratch_registry%relinquish_vector([constant_index, rhs_index, &
         initial_difference_index])
    call neko_scratch_registry%relinquish_matrix([mass_index, trial_mass_index, &
         linear_index, implicit_index])
    if (this%forward_constrained_enriched) then
       call neko_scratch_registry%relinquish_matrix([e_mass_index, e_linear_index, &
            e_implicit_index])
       call neko_scratch_registry%relinquish_vector(e_constant_index)
    end if
    do i = 1, r; call quad_local(i)%free(); end do
    deallocate(quad_local)
    nullify(mass_local, trial_mass_local, linear_local, implicit_local)
    nullify(e_mass_local, e_linear_local, e_implicit_local, e_constant_local)
    nullify(constant_local, rhs_local, initial_difference, field, scalar_field)
  end subroutine rom_problem_reassemble_adaptive_operators

  !> Update physical trial modes from the immutable transformed reference.
  !!
  !! This is the device-resident algebraic part of the design-dependent map:
  !! W=D(chi)^{-1}D(chi_ref)V0.  It intentionally does not re-SVD or change
  !! the retained rank.  The current implementation exposes the unprojected
  !! map as an engineering diagnostic; the pressure/admissibility projector
  !! and its transpose are a separate fail-closed scientific gate.
  subroutine rom_problem_update_adaptive_modes(this, mapped)
    class(rom_problem_t), intent(inout) :: this
    type(vector_t), intent(in) :: mapped
    type(vector_t), pointer :: candidate_D, ratio
    type(vector_t), pointer :: reference_difference
    type(vector_t), pointer :: transform_factor
    type(vector_t) :: potential_correction(3)
    type(field_t), pointer :: field, gx, gy, gz
    type(field_t), pointer :: p_work, p_rhs, f_x, f_y, f_z, u_zero, v_zero, w_zero
    integer :: candidate_D_index, ratio_index, reference_difference_index
    integer :: transform_factor_index
    integer :: field_index, gx_index, gy_index, gz_index
    integer :: diagnostic_field_indices(8)
    integer :: n, r, e, i, c, d
    real(kind=rp) :: projection_residual, projection_linear_residual
    real(kind=rp) :: reference_difference_norm
    real(kind=rp) :: candidate_D_min, candidate_D_max
    real(kind=rp) :: divergence_norm, velocity_norm
    real(kind=rp) :: weak_norm, boundary_norm
    real(kind=rp) :: timer_start
    type(c_ptr) :: event
    integer :: projection_iterations
    logical :: projection_converged
    logical :: shared_trial_test_modes

    if (.not. this%reference_modes_saved) call neko_error( &
         'design-dependent ROM reference modes are not available')
    if (.not. associated(this%simulation)) call neko_error( &
         'design-dependent ROM has no bound simulation mesh')
    n = mapped%size()
    r = size(this%velocity_modes, 2)
    e = size(this%test_velocity_modes, 2)
    shared_trial_test_modes = NEKO_BCKND_DEVICE .eq. 0 .and. &
         .not. this%petrov_galerkin .and. &
         .not. this%forward_constrained_enriched
    if (n /= this%reference_D%size()) call neko_error( &
         'design-dependent ROM mapped field has the wrong size')
    call neko_scratch_registry%request_vector(candidate_D, candidate_D_index, n, .true.)
    call neko_scratch_registry%request_vector(ratio, ratio_index, n, .true.)
    call neko_scratch_registry%request_vector(reference_difference, &
         reference_difference_index, n, .true.)
    if (this%basis_transform == 'vector_potential') &
         call neko_scratch_registry%request_vector(transform_factor, &
         transform_factor_index, n, .true.)
    timer_start = real(MPI_WTIME(),rp)
    call vector_copy(candidate_D, mapped)
    call vector_cmult(candidate_D, this%basis_transform_tau)
    call vector_cadd(candidate_D, 1.0_rp)
    candidate_D_min = vector_glmin(candidate_D)
    candidate_D_max = vector_glmax(candidate_D)
    if (.not. (candidate_D_min > 0.0_rp) .or. .not. (candidate_D_max < huge(1.0_rp))) then
       call neko_error('design-dependent ROM produced an invalid D=1+tau*chi field')
    end if
    call vector_copy(ratio, this%reference_D)
    call vector_invcol2(ratio, candidate_D)
    if (this%basis_transform == 'vector_potential') then
       call vector_copy(transform_factor, ratio)
       call vector_cadd(transform_factor, -1.0_rp)
       do c = 1, 3
          call potential_correction(c)%init(n)
       end do
    end if
    call vector_copy(reference_difference, candidate_D)
    call vector_add2s2(reference_difference, this%reference_D, -1.0_rp)
    call vector_absval(reference_difference)
    reference_difference_norm = vector_glmax(reference_difference)

    timer_start = real(MPI_WTIME(),rp)
    do i = 1, r
       do c = 1, 3
          call vector_copy(this%velocity_modes(c,i), &
               this%reference_velocity_modes(c,i))
       end do
       if (this%basis_transform == 'vector_potential') then
          call this%reference_forward_potential%apply( &
               this%reference_potential_forward_map(i), transform_factor, &
               potential_correction)
          do c = 1, 3
             call vector_add2(this%velocity_modes(c,i), potential_correction(c))
          end do
       else
          do c = 1, 3
             call vector_col2(this%velocity_modes(c,i), ratio)
          end do
       end if
    end do
    if (this%forward_constrained_enriched) then
       do i = 1, e
          do c = 1, 3
             call vector_copy(this%test_velocity_modes(c,i), &
                  this%reference_test_velocity_modes(c,i))
          end do
          if (this%basis_transform == 'vector_potential') then
             call this%reference_forward_potential%apply( &
                  this%reference_potential_test_map(i), transform_factor, &
                  potential_correction)
             do c = 1, 3
                call vector_add2(this%test_velocity_modes(c,i), &
                     potential_correction(c))
             end do
          else
             do c = 1, 3
                call vector_col2(this%test_velocity_modes(c,i), ratio)
             end do
          end if
       end do
    else
       if (e /= r) call neko_error( &
            'adaptive Galerkin basis requires identical forward and adjoint ranks')
       do i = 1, r
          do c = 1, 3
             call vector_copy(this%test_velocity_modes(c,i), this%velocity_modes(c,i))
          end do
       end do
    end if
    if (this%basis_transform == 'vector_potential') then
       this%potential_transform_seconds = this%potential_transform_seconds + &
            real(MPI_WTIME(),rp) - timer_start
       do c = 1, 3
          call potential_correction(c)%free()
       end do
    end if
    call neko_scratch_registry%request_field(field, field_index, &
         this%simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(gx, gx_index, &
         this%simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(gy, gy_index, &
         this%simulation%fluid%u%dof, .true.)
    call neko_scratch_registry%request_field(gz, gz_index, &
         this%simulation%fluid%u%dof, .true.)
    if (this%basis_transform == 'projected' .and. &
         reference_difference_norm > 100.0_rp * epsilon(1.0_rp) * &
         max(1.0_rp, vector_glmax(candidate_D))) then
       select type (fluid => this%simulation%fluid)
       type is (fluid_pnpn_t)
          do i = 1, r
             call vector_to_field(field, this%velocity_modes(1,i))
             call vector_to_field(gx, this%velocity_modes(2,i))
             call vector_to_field(gy, this%velocity_modes(3,i))
             call rom_project_velocity_column(fluid, field, gx, gy, &
                  projection_residual, linear_residual = projection_linear_residual, &
                  iterations = projection_iterations, converged = projection_converged)
             call field_to_vector(this%velocity_modes(1,i), field)
             call field_to_vector(this%velocity_modes(2,i), gx)
             call field_to_vector(this%velocity_modes(3,i), gy)
             if (.not. this%forward_constrained_enriched) then
                call vector_copy(this%test_velocity_modes(1,i), this%velocity_modes(1,i))
                call vector_copy(this%test_velocity_modes(2,i), this%velocity_modes(2,i))
                call vector_copy(this%test_velocity_modes(3,i), this%velocity_modes(3,i))
             end if
             if (projection_residual > this%admissibility_tolerance .or. &
                  .not. projection_converged) call neko_error( &
                  'adaptive ROM candidate mode failed the weak-continuity gate')
          end do
       class default
          call neko_error('adaptive ROM requires a Pn/Pn fluid')
       end select
    end if
    this%basis_adaptation_seconds = this%basis_adaptation_seconds + &
         real(MPI_WTIME(),rp) - timer_start
    timer_start = real(MPI_WTIME(),rp)
    do i = 1, r
       do c = 1, 3
          call vector_to_field(field, this%velocity_modes(c,i))
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call grad(gx%x_d, gy%x_d, gz%x_d, field%x_d, &
                  this%simulation%fluid%c_Xh)
          else
             call grad(gx%x, gy%x, gz%x, field%x, this%simulation%fluid%c_Xh)
          end if
          call field_to_vector(this%velocity_gradients(c,i,1), gx)
          call field_to_vector(this%velocity_gradients(c,i,2), gy)
          call field_to_vector(this%velocity_gradients(c,i,3), gz)
       end do
    end do
    do i = 1, e
       if (shared_trial_test_modes .and. i <= r) then
          do c = 1, 3
             do d = 1, 3
                call vector_copy(this%test_velocity_gradients(c,i,d), &
                     this%velocity_gradients(c,i,d))
             end do
          end do
          cycle
       end if
       do c = 1, 3
          call vector_to_field(field, this%test_velocity_modes(c,i))
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call grad(gx%x_d, gy%x_d, gz%x_d, field%x_d, &
                  this%simulation%fluid%c_Xh)
          else
             call grad(gx%x, gy%x, gz%x, field%x, this%simulation%fluid%c_Xh)
          end if
          call field_to_vector(this%test_velocity_gradients(c,i,1), gx)
          call field_to_vector(this%test_velocity_gradients(c,i,2), gy)
          call field_to_vector(this%test_velocity_gradients(c,i,3), gz)
       end do
    end do
    this%spatial_mode_derivatives_seconds = &
         this%spatial_mode_derivatives_seconds + real(MPI_WTIME(),rp) - timer_start
    this%adaptive_strong_divergence = 0.0_rp
    this%adaptive_weak_continuity = 0.0_rp
    this%adaptive_boundary_residual = 0.0_rp
    if (this%experimental_algebraic .or. &
         this%basis_transform == 'vector_potential') then
       ! Use the same backend-aware discrete diagnostics as basis handoff.
       ! The weak residual uses Neko's pressure residual operator, but performs
       ! no pressure solve or velocity projection.
       call neko_scratch_registry%request_field(p_work, diagnostic_field_indices(1), &
            this%simulation%fluid%u%dof, .true.)
       call neko_scratch_registry%request_field(p_rhs, diagnostic_field_indices(2), &
            this%simulation%fluid%u%dof, .true.)
       call neko_scratch_registry%request_field(f_x, diagnostic_field_indices(3), &
            this%simulation%fluid%u%dof, .true.)
       call neko_scratch_registry%request_field(f_y, diagnostic_field_indices(4), &
            this%simulation%fluid%u%dof, .true.)
       call neko_scratch_registry%request_field(f_z, diagnostic_field_indices(5), &
            this%simulation%fluid%u%dof, .true.)
       call neko_scratch_registry%request_field(u_zero, diagnostic_field_indices(6), &
            this%simulation%fluid%u%dof, .true.)
       call neko_scratch_registry%request_field(v_zero, diagnostic_field_indices(7), &
            this%simulation%fluid%u%dof, .true.)
       call neko_scratch_registry%request_field(w_zero, diagnostic_field_indices(8), &
            this%simulation%fluid%u%dof, .true.)
       call field_rzero(u_zero); call field_rzero(v_zero); call field_rzero(w_zero)
       event = C_NULL_PTR
       select type (fluid => this%simulation%fluid)
       type is (fluid_pnpn_t)
          do i = 1, r
             call vector_to_field(field, this%velocity_modes(1,i))
             call vector_to_field(gx, this%velocity_modes(2,i))
             call vector_to_field(gy, this%velocity_modes(3,i))
             call rom_divergence_norm(fluid, field, gx, gy, gz, n, divergence_norm)
             velocity_norm = sqrt(max(0.0_rp, &
                  vector_glsc3(this%velocity_modes(1,i), this%weights, &
                  this%velocity_modes(1,i)) + &
                  vector_glsc3(this%velocity_modes(2,i), this%weights, &
                  this%velocity_modes(2,i)) + &
                  vector_glsc3(this%velocity_modes(3,i), this%weights, &
                  this%velocity_modes(3,i))))
             this%adaptive_strong_divergence = max( &
                  this%adaptive_strong_divergence, &
                  divergence_norm / max(velocity_norm, tiny(1.0_rp)))
             call rom_weak_continuity_norm(fluid, field, gx, gy, p_work, p_rhs, &
                  f_x, f_y, f_z, u_zero, v_zero, w_zero, event, n, weak_norm)
             this%adaptive_weak_continuity = max( &
                  this%adaptive_weak_continuity, weak_norm)
             call rom_bc_residual_norm(fluid, field, gx, gy, p_work, p_rhs, &
                  f_x, n, boundary_norm)
             this%adaptive_boundary_residual = max( &
                  this%adaptive_boundary_residual, boundary_norm)
          end do
       class default
          call neko_error('adaptive ROM diagnostics require a Pn/Pn fluid')
       end select
       call neko_scratch_registry%relinquish_field(diagnostic_field_indices)
       nullify(p_work, p_rhs, f_x, f_y, f_z, u_zero, v_zero, w_zero)
    end if
    call neko_scratch_registry%relinquish_field([field_index, gx_index, gy_index, gz_index])
    nullify(field, gx, gy, gz)
    call neko_scratch_registry%relinquish_vector([candidate_D_index, ratio_index, &
         reference_difference_index])
    if (this%basis_transform == 'vector_potential') &
         call neko_scratch_registry%relinquish_vector(transform_factor_index)
    nullify(candidate_D, ratio, reference_difference)
    nullify(transform_factor)
  end subroutine rom_problem_update_adaptive_modes

  !> Project a reduced density onto the physical design representation.
  subroutine rom_problem_project_density(this, design, density, output)
    class(rom_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design
    type(vector_t), intent(in) :: density
    type(vector_t), intent(inout) :: output
    call vector_copy(output, density)
    call design%map_backward(output)
    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(output)
       ! Match the production MMA convention.  The mapping cascade produces
       ! the pointwise sensitivity first; project_sensitivity then applies the
       ! constant average SEM weight used by the full-order optimizer.
       call des%project_sensitivity(output)
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
       call des%project_sensitivity(result)
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
    character(len=1024) :: trajectory_path
    integer :: n_output_fields
    call rom_problem_free_outputs(this)
    ! ROM outputs are deliberately independent of Neko's generic field-output
    ! directory.  Every basis/optimization design gets a permanent outer
    ! directory.  The field_output_t streams below stay open for the whole
    ! outer iteration, so Neko supplies the normal f00000, f00001, ...
    ! counters for successive inner samples.
    this%rom_output_root = 'rom_output'
    write(trajectory_path, '("",A,"/outer_",I3.3)') &
         trim(this%rom_output_root), iter
    this%rom_outer_directory = trim(trajectory_path)
    call mkdir(trim(this%rom_output_root))
    call mkdir(trim(this%rom_outer_directory))
    this%output_outer_iteration = iter
    this%output_inner_iteration = -1
    this%rom_inner_directory = ''
    this%reduced_system_export_count = 0

    ! The reduced system export is tied to the basis and design of this outer
    ! iteration.  The configured root name remains a compatibility setting,
    ! while each permanent copy is stored below the outer directory.
    if (this%export_reduced_system) then
       this%reduced_system_export_directory = trim(this%rom_outer_directory)// &
            '/reduced_system'
    end if

    call this%rom_velocity(1)%init(simulation%fluid%u%dof, 'ROM_u')
    call this%rom_velocity(2)%init(simulation%fluid%v%dof, 'ROM_v')
    call this%rom_velocity(3)%init(simulation%fluid%w%dof, 'ROM_w')
    n_output_fields = 3
    if (this%has_scalar) then
       call this%rom_scalar%init(simulation%scalars%scalar_fields(1)%scalar%s%dof, &
            'ROM_scalar')
       n_output_fields = 4
    end if
    ! Keep one standard Neko fld output stream per field family for the
    ! complete outer iteration.  Each sample advances Neko's normal fld
    ! counter, producing e.g.
    !   outer_000/forward_fields0.f00000
    !   outer_000/forward_fields0.f00001
    ! without inventing an inner-iteration directory or filename scheme.
    call this%rom_forward_output%init('forward_fields', n_output_fields, &
         precision=sp, path=trim(this%rom_outer_directory)//'/', format='fld')
    call this%rom_forward_output%fields%assign(1, this%rom_velocity(1)%field)
    call this%rom_forward_output%fields%assign(2, this%rom_velocity(2)%field)
    call this%rom_forward_output%fields%assign(3, this%rom_velocity(3)%field)
    if (this%has_scalar) call this%rom_forward_output%fields%assign(4, &
         this%rom_scalar%field)
    call this%rom_design_output%init('design', 1, precision=sp, &
         path=trim(this%rom_outer_directory)//'/', format='fld')
    call this%rom_sensitivity_output%init('sensitivity', 1, precision=sp, &
         path=trim(this%rom_outer_directory)//'/', format='fld')
    select type (des => design)
    type is (brinkman_design_t)
       call des%configure_rom_outputs(this%rom_design_output, &
            this%rom_sensitivity_output)
    class default
       call neko_error('rom_mma requires a Brinkman design')
    end select
    if (this%write_trajectory_fields) then
       call this%rom_adjoint_velocity(1)%init(simulation%fluid%u%dof, 'ROM_adj_u')
       call this%rom_adjoint_velocity(2)%init(simulation%fluid%v%dof, 'ROM_adj_v')
       call this%rom_adjoint_velocity(3)%init(simulation%fluid%w%dof, 'ROM_adj_w')
       if (this%has_scalar) call this%rom_adjoint_scalar%init( &
            simulation%scalars%scalar_fields(1)%scalar%s%dof, 'ROM_adj_scalar')
    end if
    if (this%write_trajectory_fields) then
       trajectory_path = trim(this%rom_outer_directory)//'/trajectory/'
       call mkdir(trim(trajectory_path))
       call this%rom_forward_trajectory_output%init( &
            'forward_trajectory', n_output_fields, &
            precision=sp, &
            path=trim(trajectory_path), format='fld')
       call this%rom_forward_trajectory_output%fields%assign(1, &
            this%rom_velocity(1)%field)
       call this%rom_forward_trajectory_output%fields%assign(2, &
            this%rom_velocity(2)%field)
       call this%rom_forward_trajectory_output%fields%assign(3, &
            this%rom_velocity(3)%field)
       if (this%has_scalar) call this%rom_forward_trajectory_output%fields%assign(4, &
            this%rom_scalar%field)
       call this%rom_adjoint_trajectory_output%init( &
            'adjoint_trajectory', n_output_fields, &
            precision=sp, &
            path=trim(trajectory_path), format='fld')
       call this%rom_adjoint_trajectory_output%fields%assign(1, &
            this%rom_adjoint_velocity(1)%field)
       call this%rom_adjoint_trajectory_output%fields%assign(2, &
            this%rom_adjoint_velocity(2)%field)
       call this%rom_adjoint_trajectory_output%fields%assign(3, &
            this%rom_adjoint_velocity(3)%field)
       if (this%has_scalar) call this%rom_adjoint_trajectory_output%fields%assign(4, &
            this%rom_adjoint_scalar%field)
    end if
    this%rom_outputs_initialized = .true.
    call rom_problem_write_outer_metadata(this)
  end subroutine rom_problem_initialize_outputs

  !> Select the current inner sample of the outer field-output streams.
  !! @param[inout] this ROM problem owning the outer output directory.
  !! @param[inout] design Design object used to attach design/sensitivity fields.
  !! @param[in] inner ROM candidate index, with zero reserved for the reference.
  subroutine rom_problem_prepare_inner_outputs(this, design, inner)
    class(rom_problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design
    integer, intent(in) :: inner
    if (.not. this%rom_outputs_initialized) return
    ! The field_output_t streams are initialized once per outer iteration and
    ! use Neko's standard counter to distinguish inner samples.  Keep the
    ! outer directory as the CSV destination as well; no inner directory is
    ! created.
    this%rom_inner_directory = trim(this%rom_outer_directory)//'/'
    this%output_inner_iteration = inner
  end subroutine rom_problem_prepare_inner_outputs

  !> Write immutable Fortran-side provenance for one ROM outer iteration.
  !! Python writes the complementary basis_metadata.json after the streamed
  !! POD blocks have been assembled.
  subroutine rom_problem_write_outer_metadata(this)
    class(rom_problem_t), intent(in) :: this
    integer :: unit, ios
    character(len=1024) :: filename
    character(len=8) :: logical_string

    if (pe_rank /= 0 .or. .not. allocated(this%rom_outer_directory)) return
    filename = trim(this%rom_outer_directory)//'/rom_metadata.json'
    open(newunit=unit, file=trim(filename), status='replace', action='write', &
         iostat=ios)
    if (ios /= 0) call neko_error('could not open ROM outer metadata')
    write(unit, '(A)') '{'
    write(unit, '(A,I0,A)') '  "outer_iteration": ', &
         this%output_outer_iteration, ','
    write(unit, '(A,A,A)') '  "basis_forward": "', trim(this%basis_forward), '",'
    write(unit, '(A,A,A)') '  "basis_adjoint": "', trim(this%basis_adjoint), '",'
    write(unit, '(A,A,A)') '  "formulation": "', trim(this%formulation), '",'
    write(unit, '(A,A,A)') '  "sensitivity_source": "', &
         trim(this%sensitivity_source), '",'
    write(unit, '(A,A,A)') '  "basis_transform": "', &
         trim(this%basis_transform), '",'
    logical_string = rom_problem_json_logical(this%differentiate_basis)
    write(unit, '(A,A,A)') '  "differentiate_basis": ', trim(logical_string), ','
    write(unit, '(A,I0,A)') '  "forward_retained_rank": ', this%rank_forward, ','
    write(unit, '(A,I0,A)') '  "adjoint_retained_rank": ', this%rank_adjoint, ','
    write(unit, '(A,I0,A)') '  "enriched_retained_rank": ', this%rank_enriched, ','
    write(unit, '(A,I0,A)') '  "history_outer_iterations": ', &
         this%history_outer_iterations, ','
    write(unit, '(A,ES24.16,A)') '  "time_start": ', this%time_start, ','
    write(unit, '(A,ES24.16,A)') '  "time_end": ', this%time_end, ','
    write(unit, '(A,ES24.16,A)') '  "dt": ', this%dt, ','
    write(unit, '(A,I0,A)') '  "n_steps": ', this%n_steps, ','
    write(unit, '(A,I0,A)') '  "substeps": ', this%substeps, ','
    write(unit, '(A)') '  "source_revision": "not_available_in_fortran",'
    write(unit, '(A)') '  "output_layout": "outer_NNN with standard fld counters and flat inner-indexed CSV files",'
    write(unit, '(A)') '  "latest_aliases": "root-level aliases are optional and explicit"'
    write(unit, '(A)') '}'
    close(unit)
  end subroutine rom_problem_write_outer_metadata

  !> Write the reduced forward field at the initial MMA output sample.
  !!
  !! The full-order MMA optimizer writes the current forward solution with
  !! iteration zero immediately after its initial forward solve and before
  !! the adjoint solve.  Keep the ROM output stream aligned with that lifecycle
  !! without sampling the design or sensitivity streams before the adjoint
  !! sensitivity has been computed.
  !! @param[inout] this ROM problem with initialized output streams.
  !! @param[inout] simulation Neko simulation providing the field layout.
  !! @param[in] inner Output sample index; the initial call uses zero.
  !! @param[in] state Final reduced forward state to reconstruct.
  subroutine rom_problem_write_forward_output(this, simulation, inner, state)
    class(rom_problem_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    integer, intent(in) :: inner
    real(kind=rp), intent(in) :: state(:)

    if (.not. this%rom_outputs_initialized) return
    ! This entry point is retained for callers that only have the initial
    ! forward state.  The complete candidate lifecycle is written by
    ! rom_problem_write_outputs below.
    call rom_problem_reconstruct_velocity(this, state, this%rom_velocity, .true.)
    if (this%has_scalar) call rom_problem_reconstruct_scalar(this, state, &
         this%rom_scalar, .true.)
    if (this%output_inner_iteration /= inner) then
       ! A design is required only for the complete output path.  The
       ! initial-forward compatibility path therefore uses the already
       ! initialised stream when available.
       if (this%output_inner_iteration < 0) return
    end if
    call this%rom_forward_output%sample(real(inner, kind=dp))
  end subroutine rom_problem_write_forward_output

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
    call rom_problem_prepare_inner_outputs(this, design, inner)
    n = simulation%fluid%u%dof%size(); r = size(state)
    call neko_scratch_registry%request_vector(state_vector, state_index, r, .true.)
    state_vector%x = state
    call neko_scratch_registry%request_vector(reconstruction, reconstruction_index, n, .true.)
    do c = 1, 3
       reconstruction = this%lift_velocity(c)
       do i = 1, r
          call vector_add2s2(reconstruction, this%velocity_modes(c,i), state_vector%x(i))
       end do
       call vector_to_field(this%rom_velocity(c)%field, reconstruction)
    end do
    if (this%has_scalar) then
       reconstruction = this%scalar_lift
       do i = 1, r
          call vector_add2s2(reconstruction, this%scalar_modes(i), &
               state_vector%x(i))
       end do
       call vector_to_field(this%rom_scalar%field, reconstruction)
    end if
    call this%rom_forward_output%sample(real(inner, kind=dp))
    call this%rom_design_output%sample(real(inner, kind=dp))
    call this%rom_sensitivity_output%sample(real(inner, kind=dp))
    call rom_problem_write_coefficients(this, inner)
    call neko_scratch_registry%relinquish_vector([state_index, reconstruction_index])
    nullify(state_vector); nullify(reconstruction)
  end subroutine rom_problem_write_outputs

  !> Write reduced forward and adjoint coefficient histories for one inner step.
  !! The field trajectory is the physical inspection output; these compact CSV
  !! files expose the exact common-basis coordinates used to produce it,
  !! including the scalar contribution when a four-field basis is active.
  subroutine rom_problem_write_coefficients(this, inner)
    class(rom_problem_t), intent(in) :: this
    integer, intent(in) :: inner
    integer :: unit, step, mode, r, ra
    character(len=128) :: mode_name

    ! The reduced trajectory is replicated on all ranks, but this is one
    ! shared CSV per candidate.  Restrict the write to Neko's root rank.
    if (pe_rank /= 0 .or. .not. allocated(this%rom_inner_directory)) return
    r = this%states%get_nrows()
    ra = this%adjoint%get_nrows()
    write(mode_name, '("forward_coefficients_inner_",I3.3,".csv")') inner
    open(newunit=unit, file=trim(this%rom_inner_directory)//trim(mode_name), &
         status='replace', action='write')
    write(unit, '(A)', advance='no') 'step,time'
    do mode = 1, r
       write(mode_name, '(A,I0)') ',a_', mode
       write(unit, '(A)', advance='no') trim(mode_name)
    end do
    write(unit, '(A)') ''
    do step = 0, this%n_steps
       write(unit, '(I0,",",ES24.16)', advance='no') step, &
            this%time_start + real(step, rp) * this%dt
       do mode = 1, r
          write(unit, '(",",ES24.16)', advance='no') this%states%x(mode,step+1)
       end do
       write(unit, '(A)') ''
    end do
    close(unit)

    if (.not. this%backward_ready) return
    write(mode_name, '("adjoint_coefficients_inner_",I3.3,".csv")') inner
    open(newunit=unit, file=trim(this%rom_inner_directory)//trim(mode_name), &
         status='replace', action='write')
    write(unit, '(A)', advance='no') 'step,time'
    do mode = 1, ra
       write(mode_name, '(A,I0)') ',lambda_', mode
       write(unit, '(A)', advance='no') trim(mode_name)
    end do
    write(unit, '(A)') ''
    do step = 0, this%n_steps
       write(unit, '(I0,",",ES24.16)', advance='no') step, &
            this%time_start + real(step, rp) * this%dt
       do mode = 1, ra
          write(unit, '(",",ES24.16)', advance='no') this%adjoint%x(mode,step+1)
       end do
       write(unit, '(A)') ''
    end do
    close(unit)
  end subroutine rom_problem_write_coefficients

  !> Write the scalar ROM diagnostics for one candidate design.
  !! @param[in] outer Outer optimization index.
  !! @param[in] inner Inner candidate index, zero for the reference design.
  !! @param[in] objective Weighted ROM objective.
  !! @param[in] constraint ROM constraint value.
  !! @param[in] rom_forward_seconds Forward wall time.
  !! @param[in] rom_adjoint_seconds Adjoint/sensitivity wall time.
  !! @param[in] mma_seconds MMA wall time.
  subroutine rom_problem_write_inner_diagnostics(this, outer, inner, objective, &
       constraint, rom_forward_seconds, rom_adjoint_seconds, mma_seconds)
    class(rom_problem_t), intent(in) :: this
    integer, intent(in) :: outer, inner
    real(kind=rp), intent(in) :: objective, constraint
    real(kind=rp), intent(in) :: rom_forward_seconds, rom_adjoint_seconds
    real(kind=rp), intent(in) :: mma_seconds
    integer :: unit, ios
    character(len=1024) :: filename
    character(len=128) :: diagnostics_file, design_file
    character(len=128) :: forward_file, adjoint_coeff_file

    if (pe_rank /= 0 .or. .not. allocated(this%rom_inner_directory)) return
    write(diagnostics_file, '("objective_diagnostics_inner_",I3.3,".csv")') inner
    filename = trim(this%rom_inner_directory)//trim(diagnostics_file)
    open(newunit=unit, file=trim(filename), status='replace', action='write', &
         iostat=ios)
    if (ios /= 0) call neko_error('could not open ROM objective diagnostics')
    write(design_file, '("design0.f",I5.5)') inner
    write(forward_file, '("forward_coefficients_inner_",I3.3,".csv")') inner
    write(adjoint_coeff_file, '("adjoint_coefficients_inner_",I3.3,".csv")') inner
    write(unit, '(A)') 'outer_iteration,inner_iteration,objective,'// &
         'scalar_mixing,viscous_dissipation,brinkman_dissipation,'// &
         'volume_constraint,rom_forward_seconds,rom_adjoint_seconds,'// &
         'mma_seconds,adaptive_strong_divergence,adaptive_weak_continuity,'// &
         'adaptive_boundary_residual,adaptive_mass_condition,'// &
         'sensitivity_source,design_file,'// &
         'forward_coefficients,adjoint_coefficients,update_status'
    write(unit, '(I0,",",I0,12(",",ES24.16),",",A,",",A,",",A,",",A,",",A)') &
         outer, inner, &
         objective, this%scalar_mixing_value, this%viscous_value, &
         this%brinkman_value, constraint, rom_forward_seconds, &
         rom_adjoint_seconds, mma_seconds, this%adaptive_strong_divergence, &
         this%adaptive_weak_continuity, this%adaptive_boundary_residual, &
         this%adaptive_mass_condition, trim(this%sensitivity_source), &
         trim(design_file), trim(forward_file), trim(adjoint_coeff_file), &
         'accepted'
    close(unit)
  end subroutine rom_problem_write_inner_diagnostics

  !> Reconstruct a three-component physical velocity from ROM coefficients.
  !! @param this ROM problem containing the common basis and lift.
  !! @param coefficients Reduced coefficients to reconstruct.
  !! @param fields Destination field wrappers for the three velocity components.
  !! @param include_lift Whether to add the inhomogeneous forward lift.
  subroutine rom_problem_reconstruct_velocity(this, coefficients, fields, &
       include_lift, test_space)
    class(rom_problem_t), intent(in) :: this
    real(kind=rp), intent(in) :: coefficients(:)
    type(field_wrapper_t), intent(inout) :: fields(3)
    logical, intent(in) :: include_lift
    logical, intent(in), optional :: test_space
    type(vector_t), pointer :: reconstruction
    integer :: reconstruction_index
    integer :: n, r, c, i
    logical :: use_test

    use_test = .false.
    if (present(test_space)) use_test = test_space

    r = size(coefficients)
    n = this%weights%size()
    if (use_test) then
       if (r > this%test_basis%n_modes) call neko_error( &
            'ROM adjoint coefficient count exceeds the enriched basis')
    else if (r /= this%mass%get_nrows()) then
       call neko_error('ROM trajectory coefficient count does not match the retained basis')
    end if
    call neko_scratch_registry%request_vector(reconstruction, &
         reconstruction_index, n, .true.)
    do c = 1, 3
       if (include_lift) then
          reconstruction = this%lift_velocity(c)
       else
          reconstruction = 0.0_rp
       end if
       do i = 1, r
          if (use_test) then
             call vector_add2s2(reconstruction, this%test_velocity_modes(c,i), &
                  coefficients(i))
          else
             call vector_add2s2(reconstruction, this%velocity_modes(c,i), &
                  coefficients(i))
          end if
       end do
       call vector_to_field(fields(c)%field, reconstruction)
    end do
    call neko_scratch_registry%relinquish_vector(reconstruction_index)
    nullify(reconstruction)
  end subroutine rom_problem_reconstruct_velocity

  !> Reconstruct the scalar component from common ROM coefficients.
  subroutine rom_problem_reconstruct_scalar(this, coefficients, field, include_lift, &
       test_space)
    class(rom_problem_t), intent(in) :: this
    real(kind=rp), intent(in) :: coefficients(:)
    type(field_wrapper_t), intent(inout) :: field
    logical, intent(in) :: include_lift
    logical, intent(in), optional :: test_space
    type(vector_t), pointer :: reconstruction
    integer :: reconstruction_index, n, r, i
    logical :: use_test

    use_test = .false.
    if (present(test_space)) use_test = test_space

    if (.not. this%has_scalar) return
    r = size(coefficients); n = this%weights%size()
    if (use_test) then
       if (r > this%test_basis%n_modes) call neko_error( &
            'ROM adjoint scalar coefficient count exceeds the enriched basis')
    else if (r /= this%mass%get_nrows()) then
       call neko_error('ROM scalar coefficient count does not match the retained basis')
    end if
    call neko_scratch_registry%request_vector(reconstruction, &
         reconstruction_index, n, .true.)
    if (include_lift) then
       reconstruction = this%scalar_lift
    else
       reconstruction = 0.0_rp
    end if
    do i = 1, r
       if (use_test) then
          call vector_add2s2(reconstruction, this%test_scalar_modes(i), coefficients(i))
       else
          call vector_add2s2(reconstruction, this%scalar_modes(i), coefficients(i))
       end if
    end do
    call vector_to_field(field%field, reconstruction)
    call neko_scratch_registry%relinquish_vector(reconstruction_index)
    nullify(reconstruction)
  end subroutine rom_problem_reconstruct_scalar

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
          if (this%has_scalar) call rom_problem_reconstruct_scalar(this, &
               this%states%x(:,step+1), this%rom_scalar, .true.)
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
               this%adjoint%x(:,step+1), this%rom_adjoint_velocity, .false., &
               this%petrov_galerkin .or. this%forward_constrained_enriched)
          if (this%has_scalar) call rom_problem_reconstruct_scalar(this, &
               this%adjoint%x(:,step+1), this%rom_adjoint_scalar, .false., &
               this%petrov_galerkin .or. this%forward_constrained_enriched)
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
    call this%rom_scalar%free(); call this%rom_adjoint_scalar%free()
    if (allocated(this%rom_inner_directory)) deallocate(this%rom_inner_directory)
    if (allocated(this%rom_outer_directory)) deallocate(this%rom_outer_directory)
    if (allocated(this%rom_output_root)) deallocate(this%rom_output_root)
    this%output_outer_iteration = -1
    this%output_inner_iteration = -1
    this%rom_outputs_initialized = .false.
  end subroutine rom_problem_free_outputs

  !> Preserve the accepted physical reference velocity basis for adaptation.
  !! The copy is made only at the outer basis handoff.  Candidate evaluations
  !! always start from this immutable reference, preventing compounded
  !! D(chi)-division across inner iterations.
  subroutine rom_problem_save_reference_basis(this)
    class(rom_problem_t), intent(inout) :: this
    type(field_t), pointer :: chi_field
    type(field_t) :: potential_field, grad_x, grad_y, grad_z
    type(vector_t), allocatable :: unique_modes(:,:)
    type(vector_t) :: diagnostic_factor, diagnostic_correction(3)
    type(vector_t) :: diagnostic_difference, diagnostic_divergence
    integer :: n, r, e, c, i, j, k, n_unique
    integer :: source_mode, source_family
    real(kind=rp) :: timer_start
    real(kind=rp) :: reconstruction_error, gauge_residual, reference_norm
    logical :: found
    character(len=256) :: diagnostic_message

    this%projected_tensor_cache_valid = .false.
    this%scalar_mixing_cache_valid = .false.
    if (.not. this%basis_transform_adaptive) return
    this%potential_construction_seconds = 0.0_rp
    this%potential_poisson_solves = 0.0_rp
    n = this%weights%size()
    r = size(this%velocity_modes, 2)
    e = size(this%test_velocity_modes, 2)
    allocate(this%reference_velocity_modes(3, r))
    allocate(this%reference_test_velocity_modes(3, e))
    do i = 1, r
       do c = 1, 3
          call this%reference_velocity_modes(c, i)%init(n)
          call vector_copy(this%reference_velocity_modes(c, i), &
               this%velocity_modes(c, i))
       end do
    end do
    do i = 1, e
       do c = 1, 3
          call this%reference_test_velocity_modes(c, i)%init(n)
          call vector_copy(this%reference_test_velocity_modes(c, i), &
               this%test_velocity_modes(c, i))
       end do
    end do
    call this%reference_chi%init(n)
    call this%reference_D%init(n)
    chi_field => neko_registry%get_field('brinkman_amplitude')
    call field_to_vector(this%reference_chi, chi_field)
    call vector_copy(this%reference_D, this%reference_chi)
    call vector_cmult(this%reference_D, this%basis_transform_tau)
    call vector_cadd(this%reference_D, 1.0_rp)
    nullify(chi_field)
    if (this%basis_transform == 'vector_potential') then
       ! One shared solve set covers forward and independently selected test
       ! directions. Exact duplicates (including an enriched prefix) share the
       ! same stored potential and index.
       allocate(unique_modes(3, r + e))
       allocate(this%reference_potential_forward_map(r), &
            this%reference_potential_test_map(e))
       n_unique = 0
       do i = 1, r
          found = .false.
          do j = 1, n_unique
             if (rom_problem_velocity_modes_equal( &
                  this%reference_velocity_modes(:,i), unique_modes(:,j))) then
                this%reference_potential_forward_map(i) = j
                found = .true.
                exit
             end if
          end do
          if (.not. found) then
             n_unique = n_unique + 1
             do c = 1, 3
                call unique_modes(c,n_unique)%init(n)
                call vector_copy(unique_modes(c,n_unique), &
                     this%reference_velocity_modes(c,i))
             end do
             this%reference_potential_forward_map(i) = n_unique
          end if
       end do
       do i = 1, e
          found = .false.
          do j = 1, n_unique
             if (rom_problem_velocity_modes_equal( &
                  this%reference_test_velocity_modes(:,i), unique_modes(:,j))) then
                this%reference_potential_test_map(i) = j
                found = .true.
                exit
             end if
          end do
          if (.not. found) then
             n_unique = n_unique + 1
             do c = 1, 3
                call unique_modes(c,n_unique)%init(n)
                call vector_copy(unique_modes(c,n_unique), &
                     this%reference_test_velocity_modes(c,i))
             end do
             this%reference_potential_test_map(i) = n_unique
          end if
       end do
       timer_start = real(MPI_WTIME(),rp)
       select type (fluid => this%simulation%fluid)
       type is (fluid_pnpn_t)
          call this%reference_forward_potential%construct( &
               fluid%c_Xh, unique_modes(:,1:n_unique), &
               fluid%bcs_vel_projector, 1.0e-10_rp, 500)
       class default
          call neko_error('vector-potential ROM requires a Pn/Pn fluid')
       end select
       this%potential_construction_seconds = real(MPI_WTIME(),rp) - timer_start
       this%potential_poisson_solves = &
            real(this%reference_forward_potential%solve_count,rp)
       ! Outer-basis diagnostics use the exact candidate correction operator:
       ! native gathered curl followed by the homogeneous velocity projector.
       ! The scalar Poisson gauge is also measured rather than presumed.
       call diagnostic_factor%init(n)
       diagnostic_factor = 1.0_rp
       call diagnostic_difference%init(n)
       call diagnostic_divergence%init(n)
       do c = 1, 3
          call diagnostic_correction(c)%init(n)
       end do
       call potential_field%init(this%simulation%fluid%u%dof, 'vp_diag_potential')
       call grad_x%init(this%simulation%fluid%u%dof, 'vp_diag_grad_x')
       call grad_y%init(this%simulation%fluid%u%dof, 'vp_diag_grad_y')
       call grad_z%init(this%simulation%fluid%u%dof, 'vp_diag_grad_z')
       do j = 1, n_unique
          source_family = 0
          source_mode = 0
          do i = 1, r
             if (this%reference_potential_forward_map(i) == j) then
                source_family = 1
                source_mode = i
                exit
             end if
          end do
          if (source_family == 0) then
             do i = 1, e
                if (this%reference_potential_test_map(i) == j) then
                   source_family = 2
                   source_mode = i
                   exit
                end if
             end do
          end if
          if (source_family == 0) call neko_error( &
               'Vector-potential direction has no reference velocity mode')
          call this%reference_forward_potential%apply( &
               j, diagnostic_factor, diagnostic_correction)
          reference_norm = 0.0_rp
          do c = 1, 3
             if (source_family == 1) then
                call vector_copy(diagnostic_difference, diagnostic_correction(c))
                call vector_add2s2(diagnostic_difference, &
                     this%reference_velocity_modes(c,source_mode), -1.0_rp)
                reference_norm = reference_norm + vector_glsc3( &
                     this%reference_velocity_modes(c,source_mode), this%weights, &
                     this%reference_velocity_modes(c,source_mode))
             else
                call vector_copy(diagnostic_difference, diagnostic_correction(c))
                call vector_add2s2(diagnostic_difference, &
                     this%reference_test_velocity_modes(c,source_mode), -1.0_rp)
                reference_norm = reference_norm + vector_glsc3( &
                     this%reference_test_velocity_modes(c,source_mode), this%weights, &
                     this%reference_test_velocity_modes(c,source_mode))
             end if
             if (c == 1) then
                call vector_to_field(potential_field, &
                     this%reference_forward_potential%potential(c,j))
                if (NEKO_BCKND_DEVICE .eq. 1) then
                   call grad(grad_x%x_d, grad_y%x_d, grad_z%x_d, &
                        potential_field%x_d, this%simulation%fluid%c_Xh)
                else
                   call grad(grad_x%x, grad_y%x, grad_z%x, potential_field%x, &
                        this%simulation%fluid%c_Xh)
                end if
                call field_to_vector(diagnostic_divergence, grad_x)
             else if (c == 2) then
                call vector_to_field(potential_field, &
                     this%reference_forward_potential%potential(c,j))
                if (NEKO_BCKND_DEVICE .eq. 1) then
                   call grad(grad_x%x_d, grad_y%x_d, grad_z%x_d, &
                        potential_field%x_d, this%simulation%fluid%c_Xh)
                else
                   call grad(grad_x%x, grad_y%x, grad_z%x, potential_field%x, &
                        this%simulation%fluid%c_Xh)
                end if
                call field_to_vector(diagnostic_difference, grad_y)
                call vector_add2(diagnostic_divergence, diagnostic_difference)
             else
                call vector_to_field(potential_field, &
                     this%reference_forward_potential%potential(c,j))
                if (NEKO_BCKND_DEVICE .eq. 1) then
                   call grad(grad_x%x_d, grad_y%x_d, grad_z%x_d, &
                        potential_field%x_d, this%simulation%fluid%c_Xh)
                else
                   call grad(grad_x%x, grad_y%x, grad_z%x, potential_field%x, &
                        this%simulation%fluid%c_Xh)
                end if
                call field_to_vector(diagnostic_difference, grad_z)
                call vector_add2(diagnostic_divergence, diagnostic_difference)
             end if
          end do
          ! Compute the full vector reconstruction after the gauge derivative
          ! loop reused the difference scratch.
          reconstruction_error = 0.0_rp
          do c = 1, 3
             call vector_copy(diagnostic_difference, diagnostic_correction(c))
             if (source_family == 1) then
                call vector_add2s2(diagnostic_difference, &
                     this%reference_velocity_modes(c,source_mode), -1.0_rp)
             else
                call vector_add2s2(diagnostic_difference, &
                     this%reference_test_velocity_modes(c,source_mode), -1.0_rp)
             end if
             reconstruction_error = reconstruction_error + vector_glsc3( &
                  diagnostic_difference, this%weights, diagnostic_difference)
          end do
          reconstruction_error = sqrt(max(0.0_rp,reconstruction_error)) / &
               max(sqrt(max(0.0_rp,reference_norm)),tiny(1.0_rp))
          gauge_residual = sqrt(max(0.0_rp,vector_glsc3( &
               diagnostic_divergence,this%weights,diagnostic_divergence)))
          write(diagnostic_message, '(A,I0,A,ES12.4,A,ES12.4)') &
               'VP reference mode ', j, ' assembled-curl reconstruction error=', &
               reconstruction_error, ' Coulomb-gauge residual=', gauge_residual
          call neko_log%message(trim(diagnostic_message))
       end do
       do c = 1, 3
          call diagnostic_correction(c)%free()
       end do
       call diagnostic_factor%free()
       call diagnostic_difference%free()
       call diagnostic_divergence%free()
       call potential_field%free()
       call grad_x%free(); call grad_y%free(); call grad_z%free()
       do j = 1, n_unique
          do c = 1, 3
             call unique_modes(c,j)%free()
          end do
       end do
       deallocate(unique_modes)
    end if
    this%reference_modes_saved = .true.
  end subroutine rom_problem_save_reference_basis

  !> Exact equality is used only to identify vectors that can share a linear
  !! reference-potential solve. Near-duplicate POD directions remain distinct.
  logical function rom_problem_velocity_modes_equal(left, right)
    type(vector_t), intent(in) :: left(3), right(3)
    type(vector_t) :: difference
    integer :: c, n
    rom_problem_velocity_modes_equal = .true.
    n = left(1)%size()
    do c = 1, 3
       if (left(c)%size() /= n .or. right(c)%size() /= n) then
          rom_problem_velocity_modes_equal = .false.
          return
       end if
    end do
    call difference%init(3*n)
    do c = 1, 3
       difference%x((c-1)*n+1:c*n) = left(c)%x - right(c)%x
    end do
    call vector_absval(difference)
    rom_problem_velocity_modes_equal = vector_glmax(difference) == 0.0_rp
    call difference%free()
  end function rom_problem_velocity_modes_equal

  !> Release the immutable adaptive reference basis and its design field.
  subroutine rom_problem_free_reference_basis(this)
    class(rom_problem_t), intent(inout) :: this
    integer :: i, c
    call this%reference_forward_potential%free()
    if (allocated(this%reference_potential_forward_map)) &
         deallocate(this%reference_potential_forward_map)
    if (allocated(this%reference_potential_test_map)) &
         deallocate(this%reference_potential_test_map)
    if (allocated(this%reference_velocity_modes)) then
       do i = 1, size(this%reference_velocity_modes, 2)
          do c = 1, size(this%reference_velocity_modes, 1)
             call this%reference_velocity_modes(c, i)%free()
          end do
       end do
       deallocate(this%reference_velocity_modes)
    end if
    if (allocated(this%reference_test_velocity_modes)) then
       do i = 1, size(this%reference_test_velocity_modes, 2)
          do c = 1, size(this%reference_test_velocity_modes, 1)
             call this%reference_test_velocity_modes(c, i)%free()
          end do
       end do
       deallocate(this%reference_test_velocity_modes)
    end if
    call this%reference_chi%free()
    call this%reference_D%free()
    this%reference_modes_saved = .false.
    this%potential_construction_seconds = 0.0_rp
    this%potential_poisson_solves = 0.0_rp
  end subroutine rom_problem_free_reference_basis

  !> Release vector and matrix containers owned by the ROM problem.
  subroutine rom_problem_free_arrays(this)
    class(rom_problem_t), intent(inout) :: this
    integer :: i, j, k
    if (allocated(this%projected_convection_cache)) then
       do i = 1, size(this%projected_convection_cache)
          call this%projected_convection_cache(i)%free()
       end do
       deallocate(this%projected_convection_cache)
    end if
    this%projected_tensor_cache_valid = .false.
    this%scalar_mixing_cache_valid = .false.
    this%scalar_transport_cache_valid = .false.
    call this%scalar_mass_cache%free()
    call this%scalar_adjoint_mass_cache%free()
    call this%scalar_diffusion_cache%free()
    call this%scalar_adjoint_diffusion_cache%free()
    call this%scalar_diffusion_lift_cache%free()
    call this%scalar_adjoint_diffusion_lift_cache%free()
    call this%weights%free(); call this%viscous_mask%free()
    call this%brinkman_mask%free(); call this%scalar_mixing_mask%free()
    call this%volume_mask%free()
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
    if (allocated(this%test_velocity_modes)) then
       do i = 1, size(this%test_velocity_modes,1)
          do j = 1, size(this%test_velocity_modes,2)
             call this%test_velocity_modes(i,j)%free()
          end do
       end do
       deallocate(this%test_velocity_modes)
    end if
    if (allocated(this%test_velocity_gradients)) then
       do i = 1, size(this%test_velocity_gradients,1)
          do j = 1, size(this%test_velocity_gradients,2)
             do k = 1, size(this%test_velocity_gradients,3)
                call this%test_velocity_gradients(i,j,k)%free()
             end do
          end do
       end do
       deallocate(this%test_velocity_gradients)
    end if
    do i = 1, 3; call this%lift_velocity(i)%free(); end do
    do i = 1, 3
       do j = 1, 3
          call this%lift_gradients(i,j)%free()
       end do
    end do
    if (allocated(this%scalar_modes)) then
       do i = 1, size(this%scalar_modes)
          call this%scalar_modes(i)%free()
          do j = 1, size(this%scalar_gradients,2)
             call this%scalar_gradients(i,j)%free()
          end do
       end do
       deallocate(this%scalar_modes, this%scalar_gradients)
    end if
    if (allocated(this%test_scalar_modes)) then
       do i = 1, size(this%test_scalar_modes)
          call this%test_scalar_modes(i)%free()
          do j = 1, size(this%test_scalar_gradients,2)
             call this%test_scalar_gradients(i,j)%free()
          end do
       end do
       deallocate(this%test_scalar_modes, this%test_scalar_gradients)
    end if
    call this%scalar_diffusivity%free(); call this%scalar_lift%free()
    call this%scalar_reference%free()
    do i = 1, 3; call this%scalar_lift_gradients(i)%free(); end do
    call this%mass%free(); call this%trial_mass%free(); call this%fixed_constant%free()
    call this%fixed_linear%free(); call this%fixed_implicit%free()
    call this%enriched_mass%free(); call this%enriched_fixed_linear%free()
    call this%enriched_fixed_implicit%free()
    call this%candidate_enriched_brinkman%free()
    call this%enriched_adjoint_model%free()
    if (allocated(this%quadratic)) then
       do i = 1, size(this%quadratic); call this%quadratic(i)%free(); end do
       deallocate(this%quadratic)
    end if
    call this%viscous_constant%free(); call this%viscous_linear%free()
    call this%viscous_quadratic%free(); call this%initial_state%free()
    call this%enriched_viscous_linear%free(); call this%enriched_viscous_cross%free()
    call this%enriched_scalar_mixing_linear%free()
    call this%enriched_scalar_mixing_cross%free()
    call this%candidate_objective_enriched_lift%free()
    call this%candidate_objective_enriched_cross%free()
    call this%scalar_mixing_constant%free(); call this%scalar_mixing_linear%free()
    call this%scalar_mixing_quadratic%free()
    if (allocated(this%cross_mass_singular_values)) &
         deallocate(this%cross_mass_singular_values)
    this%cross_mass_retained_rank = 0
    this%cross_mass_sigma_min = 0.0_rp
    this%cross_mass_sigma_max = 0.0_rp
    this%cross_mass_condition = huge(0.0_rp)
    this%enriched_adjoint_rank = 0
  end subroutine rom_problem_free_arrays

end module rom_problem
