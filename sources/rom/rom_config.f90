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
     ! Nonzero only for the narrowly scoped legacy top-level `rom` path.
     ! New configurations retain every independent mode received from the
     ! requested POD history and determine ranks after handoff.
     integer :: legacy_rank_forward = 0
     integer :: legacy_rank_adjoint = 0
     integer :: legacy_rank_enriched = 0
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
   contains
     procedure, pass(this), public :: init => rom_config_init
  end type rom_config_t

contains

  !> Resolve the one supported ROM configuration location.
  !!
  !! New cases use `optimization.solver.rom`.  A top-level `rom` object is
  !! accepted only as an explicit legacy migration path.  Supplying both
  !! locations is rejected because silently selecting one creates different
  !! Python and Fortran ROM configurations.
  !! @param[in] parameters Complete case JSON object.
  !! @param[out] rom_parameters Resolved ROM object.
  !! @param[out] found Whether either supported location was present.
  subroutine rom_config_get_json(parameters, rom_parameters, found)
    type(json_file), intent(inout) :: parameters
    type(json_file), intent(out) :: rom_parameters
    logical, intent(out) :: found
    logical :: nested, legacy
    type(json_file) :: solver_parameters

    nested = parameters%valid_path('optimization.solver.rom')
    legacy = 'rom' .in. parameters
    if (nested .and. legacy) call neko_error( &
         'ROM configuration is present in both optimization.solver.rom and '&
         //'legacy top-level rom; remove the legacy object')
    found = nested .or. legacy
    if (nested) then
       call json_get(parameters, 'optimization.solver.rom', rom_parameters)
    else if (legacy) then
       call json_get(parameters, 'rom', rom_parameters)
    else
       call json_get(parameters, 'optimization.solver', solver_parameters)
       rom_parameters = solver_parameters
    end if
  end subroutine rom_config_get_json

  !> Parse and validate the canonical ROM settings.
  !! @param[inout] this Settings object.
  !! @param[in] parameters Complete case JSON object.
  subroutine rom_config_init(this, parameters)
    class(rom_config_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(json_file) :: rom_parameters
    logical :: found, has_legacy_rank, has_forward_rank, has_adjoint_rank
    logical :: has_enriched_rank, has_basis_forward, has_basis_adjoint
    logical :: legacy_location, has_formulation, has_basis
    character(len=:), allocatable :: formulation, sensitivity_source
    character(len=:), allocatable :: basis_forward, basis_adjoint
    character(len=:), allocatable :: reduced_directory, time_scheme
    character(len=:), allocatable :: pressure_treatment, admissibility_projection
    character(len=:), allocatable :: diagnostics_file, lift_source
    integer :: legacy_rank

    call rom_config_get_json(parameters, rom_parameters, found)
    if (.not. found) call neko_error( &
         'ROM configuration requires optimization.solver.rom')

    legacy_location = .not. parameters%valid_path('optimization.solver.rom')
    has_formulation = 'formulation' .in. rom_parameters
    has_basis = rom_parameters%valid_path('basis')
    has_basis_forward = rom_parameters%valid_path('basis.forward')
    has_basis_adjoint = rom_parameters%valid_path('basis.adjoint')
    call json_get_or_default(rom_parameters, 'enabled', this%enabled, .true.)

    if (.not. legacy_location) then
       if (has_formulation .or. ('sensitivity_source' .in. rom_parameters) .or. &
            rom_parameters%valid_path('rank_forward') .or. &
            rom_parameters%valid_path('rank_adjoint') .or. &
            rom_parameters%valid_path('rank_enriched') .or. &
            rom_parameters%valid_path('rank_total_max')) call neko_error( &
            'obsolete ROM selection found under optimization.solver.rom; use '// &
            'basis.forward, basis.adjoint and state_recovery.n_modes')
       if (.not. has_basis_forward .or. .not. has_basis_adjoint) call neko_error( &
            'optimization.solver.rom requires basis.forward and basis.adjoint; '// &
            'supported pairs are (forward,forward), (forward,adjoint), '// &
            '(forward,enriched) and (enriched,enriched)')
       call json_get(rom_parameters, 'basis.forward', basis_forward)
       call json_get(rom_parameters, 'basis.adjoint', basis_adjoint)
       this%basis_forward = trim(basis_forward)
       this%basis_adjoint = trim(basis_adjoint)
       call rom_config_derive_formulation(this)
    else
       if (has_basis .and. has_formulation) call neko_error( &
            'legacy top-level rom cannot contain both formulation and basis')
       if (has_basis) then
          if (.not. has_basis_forward .or. .not. has_basis_adjoint) call neko_error( &
               'legacy ROM basis selection requires both basis.forward and basis.adjoint')
          call json_get(rom_parameters, 'basis.forward', basis_forward)
          call json_get(rom_parameters, 'basis.adjoint', basis_adjoint)
          this%basis_forward = trim(basis_forward)
          this%basis_adjoint = trim(basis_adjoint)
          call rom_config_derive_formulation(this)
       else
          ! The old galerkin peer merged streamed forward and adjoint POD
          ! blocks before sending one common basis.  Preserve that historical
          ! top-level behavior rather than changing its mathematical meaning.
          if (.not. has_formulation) then
             formulation = 'enriched_galerkin'
          else
             call json_get(rom_parameters, 'formulation', formulation)
          end if
          call rom_config_derive_legacy_formulation(this, formulation)
       end if
    end if

    ! Legacy cases may still state their sensitivity source explicitly.  New
    ! cases derive it from the selected basis pair and cannot override it.
    if (legacy_location .and. ('sensitivity_source' .in. rom_parameters)) then
       call json_get(rom_parameters, 'sensitivity_source', sensitivity_source)
       if (trim(sensitivity_source) == 'exact') sensitivity_source = 'exact_galerkin'
       if (trim(sensitivity_source) == 'projected') sensitivity_source = 'projected_adjoint'
       if (trim(sensitivity_source) /= trim(this%sensitivity_source)) call neko_error( &
            'legacy sensitivity_source is inconsistent with the selected ROM formulation')
    end if

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

    ! Explicit ranks are accepted only in the legacy top-level object.  New
    ! cases derive ranks from all independent modes received from history.
    has_legacy_rank = 'rank_total_max' .in. rom_parameters
    if (has_legacy_rank .and. parameters%valid_path('optimization.solver.rom')) &
         call neko_error('rank_total_max is obsolete in nested ROM configuration; use state_recovery.n_modes')
    legacy_rank = 0
    if (has_legacy_rank) call json_get(rom_parameters, 'rank_total_max', legacy_rank)
    has_forward_rank = 'rank_forward' .in. rom_parameters
    has_adjoint_rank = 'rank_adjoint' .in. rom_parameters
    has_enriched_rank = 'rank_enriched' .in. rom_parameters

    if (legacy_location) then
       if (has_legacy_rank) then
          this%legacy_rank_forward = legacy_rank
          this%legacy_rank_adjoint = legacy_rank
          this%legacy_rank_enriched = legacy_rank
       end if
       if (has_forward_rank) call json_get(rom_parameters, 'rank_forward', &
            this%legacy_rank_forward)
       if (has_adjoint_rank) call json_get(rom_parameters, 'rank_adjoint', &
            this%legacy_rank_adjoint)
       if (has_enriched_rank) call json_get(rom_parameters, 'rank_enriched', &
            this%legacy_rank_enriched)
       if (this%legacy_rank_adjoint == 0) this%legacy_rank_adjoint = &
            this%legacy_rank_forward
       if (this%legacy_rank_enriched == 0) this%legacy_rank_enriched = &
            this%legacy_rank_forward + this%legacy_rank_adjoint
    end if

    if (this%inner_iterations < 1 .or. this%history_outer_iterations < 1 .or. &
         this%substeps < 1) call neko_error( &
         'ROM inner_iterations, history_outer_iterations and substeps must be positive')
    if (legacy_location .and. (this%legacy_rank_forward < 0 .or. &
         this%legacy_rank_adjoint < 0 .or. this%legacy_rank_enriched < 0)) &
         call neko_error('legacy ROM ranks must be nonnegative')
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

  !> Preserve the old top-level formulation names without changing semantics.
  subroutine rom_config_derive_legacy_formulation(this, formulation)
    class(rom_config_t), intent(inout) :: this
    character(len=*), intent(in) :: formulation
    select case (trim(formulation))
    case ('forward_galerkin', 'galerkin')
       this%basis_forward = 'forward'; this%basis_adjoint = 'forward'
    case ('separate_pod', 'petrov_galerkin', 'forward_constrained_enriched')
       this%basis_forward = 'forward'; this%basis_adjoint = 'adjoint'
    case ('forward_enriched_adjoint', 'enriched_adjoint_galerkin')
       this%basis_forward = 'forward'; this%basis_adjoint = 'enriched'
    case ('enriched_galerkin')
       this%basis_forward = 'enriched'; this%basis_adjoint = 'enriched'
    case default
       call neko_error('unsupported legacy ROM formulation: '//trim(formulation))
    end select
    call rom_config_derive_formulation(this)
  end subroutine rom_config_derive_legacy_formulation

end module rom_config
