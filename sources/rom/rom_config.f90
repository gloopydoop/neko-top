!> @file rom_config.f90
!! @brief Canonical ROM configuration and migration handling.
module rom_config
  use json_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use num_types, only: rp
  use utils, only: neko_error
  implicit none
  private

  public :: rom_config_t, rom_config_get_json

  !> Canonical settings shared by the ROM optimiser and ROM problem.
  type, public :: rom_config_t
     logical :: enabled = .true.
     ! User-facing basis selection.  The formulation and sensitivity source
     ! below are derived implementation metadata, not configuration choices.
     character(len=32) :: basis_forward = 'forward'
     character(len=32) :: basis_adjoint = 'forward'
     character(len=64) :: formulation = 'forward_galerkin'
     character(len=32) :: sensitivity_source = 'exact_galerkin'
     integer :: inner_iterations = 3
     integer :: history_outer_iterations = 3
     integer :: substeps = 1
     real(kind=rp) :: scalar_metric_weight = 1.0_rp
     real(kind=rp) :: admissibility_tolerance = 1.0e-8_rp
     real(kind=rp) :: petrov_svd_relative_tolerance = 1.0e-8_rp
     logical :: enable_output = .true.
     logical :: write_trajectory_fields = .false.
     real(kind=rp) :: trajectory_output_interval = 0.5_rp
     logical :: export_reduced_system = .false.
     character(len=256) :: reduced_system_export_directory = 'rom_reduced_system'
     character(len=32) :: time_scheme = 'bdf1'
     character(len=64) :: pressure_treatment = 'solenoidal_galerkin'
     character(len=64) :: admissibility_projection = 'pressure'
     logical :: strict_admissibility = .false.
     character(len=256) :: projection_diagnostics_file = &
          'rom_projection_diagnostics.csv'
     character(len=64) :: lift_source = 'initial_velocity'
     logical :: validate_full_order = .true.
     real(kind=rp) :: validation_tolerance = 5.0e-2_rp
     real(kind=rp) :: validation_gradient_tolerance = 5.0e-2_rp
     logical :: diagnostic_audit = .false.
     ! ROM-only design-dependent spatial mode map.  The legacy fixed basis is
     ! the default and remains byte-compatible when this is ``none``.
     character(len=32) :: basis_transform = 'none'
     logical :: differentiate_basis = .false.
     logical :: experimental_algebraic = .false.
     real(kind=rp) :: basis_transform_tau = 1.0_rp
   contains
     procedure, pass(this), public :: init => rom_config_init
  end type rom_config_t

contains

  !> Resolve the one supported ROM configuration location.
  !!
  !! ROM cases use only `optimization.solver.rom`.  Reject the former
  !! top-level object so Python and Fortran cannot silently read different
  !! scientific settings.
  !! @param[in] parameters Complete case JSON object.
  !! @param[out] rom_parameters Resolved ROM object.
  !! @param[out] found Whether either supported location was present.
  subroutine rom_config_get_json(parameters, rom_parameters, found)
    type(json_file), intent(inout) :: parameters
    type(json_file), intent(out) :: rom_parameters
    logical, intent(out) :: found
    logical :: nested, legacy

    nested = parameters%valid_path('optimization.solver.rom')
    legacy = 'rom' .in. parameters
    if (legacy) call neko_error( &
         'top-level rom configuration was removed; use optimization.solver.rom')
    found = nested
    if (nested) call json_get(parameters, 'optimization.solver.rom', rom_parameters)
  end subroutine rom_config_get_json

  !> Parse and validate the canonical ROM settings.
  !! @param[inout] this Settings object.
  !! @param[in] parameters Complete case JSON object.
  subroutine rom_config_init(this, parameters)
    class(rom_config_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(json_file) :: rom_parameters
    logical :: found, has_basis_forward, has_basis_adjoint
    character(len=:), allocatable :: basis_forward, basis_adjoint
    character(len=:), allocatable :: reduced_directory, time_scheme
    character(len=:), allocatable :: pressure_treatment, admissibility_projection
    character(len=:), allocatable :: diagnostics_file, lift_source
    character(len=:), allocatable :: basis_transform

    call rom_config_get_json(parameters, rom_parameters, found)
    if (.not. found) call neko_error( &
         'ROM configuration requires optimization.solver.rom')

    has_basis_forward = rom_parameters%valid_path('basis.forward')
    has_basis_adjoint = rom_parameters%valid_path('basis.adjoint')
    call json_get_or_default(rom_parameters, 'enabled', this%enabled, .true.)

    if (rom_parameters%valid_path('formulation') .or. &
         ('sensitivity_source' .in. rom_parameters) .or. &
         rom_parameters%valid_path('rank_forward') .or. &
         rom_parameters%valid_path('rank_adjoint') .or. &
         rom_parameters%valid_path('rank_enriched') .or. &
         rom_parameters%valid_path('rank_total_max')) call neko_error( &
         'removed ROM selection; use basis.forward, basis.adjoint and '// &
         'state_recovery.n_modes')
    if (.not. has_basis_forward .or. .not. has_basis_adjoint) call neko_error( &
         'optimization.solver.rom requires basis.forward and basis.adjoint; '// &
         'supported pairs are (forward,forward), (forward,adjoint), '// &
         '(forward,enriched) and (enriched,enriched)')
    call json_get(rom_parameters, 'basis.forward', basis_forward)
    call json_get(rom_parameters, 'basis.adjoint', basis_adjoint)
    this%basis_forward = trim(basis_forward)
    this%basis_adjoint = trim(basis_adjoint)
    call rom_config_derive_formulation(this)

    call json_get_or_default(rom_parameters, 'inner_iterations', &
         this%inner_iterations, 3)
    call json_get_or_default(rom_parameters, 'history_outer_iterations', &
         this%history_outer_iterations, 3)
    call json_get_or_default(rom_parameters, 'substeps', this%substeps, 1)
    call json_get_or_default(rom_parameters, 'scalar_metric_weight', &
         this%scalar_metric_weight, 1.0_rp)
    call json_get_or_default(rom_parameters, 'admissibility_tolerance', &
         this%admissibility_tolerance, 1.0e-8_rp)
    call json_get_or_default(rom_parameters, 'petrov_svd_relative_tolerance', &
         this%petrov_svd_relative_tolerance, 1.0e-8_rp)
    call json_get_or_default(rom_parameters, 'enable_output', this%enable_output, .true.)
    call json_get_or_default(rom_parameters, 'write_trajectory_fields', &
         this%write_trajectory_fields, .false.)
    call json_get_or_default(rom_parameters, 'trajectory_output_interval', &
         this%trajectory_output_interval, 0.5_rp)
    call json_get_or_default(rom_parameters, 'export_reduced_system', &
         this%export_reduced_system, .false.)
    call json_get_or_default(rom_parameters, 'reduced_system_export_directory', &
         reduced_directory, 'rom_reduced_system')
    this%reduced_system_export_directory = trim(reduced_directory)
    call json_get_or_default(rom_parameters, 'time_scheme', time_scheme, 'bdf1')
    this%time_scheme = trim(time_scheme)
    call json_get_or_default(rom_parameters, 'pressure_treatment', &
         pressure_treatment, 'solenoidal_galerkin')
    this%pressure_treatment = trim(pressure_treatment)
    call json_get_or_default(rom_parameters, 'admissibility_projection', &
         admissibility_projection, 'pressure')
    this%admissibility_projection = trim(admissibility_projection)
    call json_get_or_default(rom_parameters, 'strict_admissibility', &
         this%strict_admissibility, .false.)
    call json_get_or_default(rom_parameters, 'projection_diagnostics_file', &
         diagnostics_file, 'rom_projection_diagnostics.csv')
    this%projection_diagnostics_file = trim(diagnostics_file)
    call json_get_or_default(rom_parameters, 'lift.source', lift_source, &
         'initial_velocity')
    this%lift_source = trim(lift_source)
    call json_get_or_default(rom_parameters, 'validate_full_order', &
         this%validate_full_order, .true.)
    call json_get_or_default(rom_parameters, 'validation_tolerance', &
         this%validation_tolerance, 5.0e-2_rp)
    call json_get_or_default(rom_parameters, 'validation_gradient_tolerance', &
         this%validation_gradient_tolerance, 5.0e-2_rp)
    call json_get_or_default(rom_parameters, 'diagnostic_audit', &
         this%diagnostic_audit, .false.)
    call json_get_or_default(rom_parameters, 'basis_transform', &
         basis_transform, 'none')
    this%basis_transform = trim(basis_transform)
    call json_get_or_default(rom_parameters, 'differentiate_basis', &
         this%differentiate_basis, .false.)
    if (('basis_transform_model' .in. rom_parameters) .or. &
         ('algebraic_sensitivity' .in. rom_parameters) .or. &
         ('basis_pod_source' .in. rom_parameters)) call neko_error( &
         'legacy basis transformation settings were removed; use '// &
         'basis_transform and differentiate_basis')
    this%experimental_algebraic = trim(this%basis_transform) == 'algebraic'
    if (this%experimental_algebraic) then
       if (this%differentiate_basis) then
          this%sensitivity_source = 'exact_algebraic'
       else if (trim(this%formulation) == 'separate_pod' .or. &
            trim(this%formulation) == 'forward_enriched_adjoint') then
          this%sensitivity_source = 'projected_adjoint'
       else
          this%sensitivity_source = 'lagged_basis'
       end if
    else if (this%differentiate_basis .and. &
         trim(this%basis_transform) == 'projected') then
       ! Differentiate the forward Galerkin ROM using its consistent adjoint,
       ! even when the independently selected adjoint basis is enriched.
       this%sensitivity_source = 'exact_galerkin'
    end if
    call json_get_or_default(rom_parameters, 'basis_transform_tau', &
         this%basis_transform_tau, 1.0_rp)

    if (this%inner_iterations < 1 .or. this%history_outer_iterations < 1 .or. &
         this%substeps < 1) call neko_error( &
         'ROM inner_iterations, history_outer_iterations and substeps must be positive')
    if (this%scalar_metric_weight <= 0.0_rp .or. &
         this%admissibility_tolerance <= 0.0_rp) call neko_error( &
         'ROM scalar_metric_weight and admissibility_tolerance must be positive')
    if (this%petrov_svd_relative_tolerance <= 0.0_rp .or. &
         this%petrov_svd_relative_tolerance > 1.0_rp) call neko_error( &
         'ROM petrov_svd_relative_tolerance must be in (0,1]')
    if (this%write_trajectory_fields .and. &
         this%trajectory_output_interval <= 0.0_rp) call neko_error( &
         'ROM trajectory_output_interval must be positive when output is enabled')
    if (trim(this%time_scheme) /= 'bdf1') call neko_error( &
         'ROM requires time_scheme=bdf1')
    if (trim(this%pressure_treatment) /= 'solenoidal_galerkin') call neko_error( &
         'ROM requires pressure_treatment=solenoidal_galerkin')
    select case (trim(this%basis_transform))
    case ('none')
    case ('algebraic', 'projected')
       if (this%basis_transform_tau < 0.0_rp) call neko_error( &
            'ROM basis_transform_tau must be nonnegative')
    case ('vector_potential')
       call neko_error('basis_transform=vector_potential is not implemented yet')
    case default
       call neko_error('basis_transform must be none, algebraic, projected or vector_potential')
    end select
    if (trim(this%basis_transform) == 'projected' .and. &
         (trim(this%basis_forward) /= 'forward' .or. &
          trim(this%basis_adjoint) /= 'forward')) call neko_error( &
         'projected basis transformation currently supports only forward:forward')
  end subroutine rom_config_init

  !> Map a user-selected basis pair to the existing internal implementation.
  subroutine rom_config_derive_formulation(this)
    class(rom_config_t), intent(inout) :: this
    character(len=:), allocatable :: pair
    pair = '('//trim(this%basis_forward)//','//trim(this%basis_adjoint)//')'
    select case (trim(this%basis_forward)//':'//trim(this%basis_adjoint))
    case ('forward:forward')
       this%formulation = 'forward_galerkin'
       this%sensitivity_source = 'exact_galerkin'
    case ('forward:adjoint')
       this%formulation = 'separate_pod'
       this%sensitivity_source = 'projected_adjoint'
    case ('forward:enriched')
       this%formulation = 'forward_enriched_adjoint'
       this%sensitivity_source = 'projected_adjoint'
    case ('enriched:enriched')
       this%formulation = 'enriched_galerkin'
       this%sensitivity_source = 'exact_galerkin'
    case default
       call neko_error('unsupported ROM basis pair '//trim(pair)//'; supported pairs are '// &
            '(forward,forward), (forward,adjoint), (forward,enriched) and (enriched,enriched)')
    end select
  end subroutine rom_config_derive_formulation

end module rom_config
