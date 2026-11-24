module rom_basis_handoff
  use simulation_m, only: simulation_t
  use state_recover, only: state_recover_t
  use fluid_pnpn, only: fluid_pnpn_t
  use adjoint_fluid_pnpn, only: adjoint_fluid_pnpn_t
  use scalar_pnpn, only: scalar_pnpn_t
  use vector_bc_projector, only: vector_bc_projector_t
  use field, only: field_t
  use vector, only: vector_t
  use scratch_registry, only: neko_scratch_registry
  use field_math, only: field_add2s2, field_col3, field_rzero, field_copy, &
       field_sub3
  use operators, only: grad, div, ortho
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2, glsc3
  use device_math, only: device_col2, device_glsc3
  use device, only: DEVICE_TO_HOST, HOST_TO_DEVICE
  use neko_config, only: NEKO_BCKND_DEVICE
  use krylov, only: ksp_monitor_t
  use json_file_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use num_types, only: rp
  use comm, only: pe_rank
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
  use logger, only: neko_log
  use utils, only: neko_error, mkdir
  use rom_basis_orthogonalization, only: rom_mass_dot, &
       rom_mass_orthonormalize_device, rom_mass_orthonormalize_coupled_device
  use rom_config, only: rom_config_get_json
  use rom_projection_operator, only: rom_projection_ax_t
  implicit none
  private

  public :: rom_receive_basis, rom_receive_petrov_bases, rom_project_basis, &
       rom_project_test_basis, rom_direct_pod_basis, rom_restrict_scalar_basis, &
       rom_orthonormalize_projected_device, &
       rom_orthonormalize_projected_coupled_device

contains

  !> Orthonormalise the final projected velocity columns on Neko's backend.
  !!
  !! `projected` is the host staging buffer supplied by the ADIOS/Python
  !! handoff.  It is uploaded once after projection, processed with device
  !! vector contractions and axpy operations, and read back once for the
  !! returned host basis.  No invocation occurs in the reduced inner loop.
  subroutine rom_orthonormalize_projected_device(projected, n, n_modes, mass, &
       tolerance, accepted, final_mode_index, mass_norm_final, mass_cross_final, ok)
    real(kind=rp), allocatable, intent(inout) :: projected(:, :)
    integer, intent(in) :: n
    integer, intent(in) :: n_modes
    real(kind=rp), intent(in) :: mass(:), tolerance
    integer, intent(out) :: accepted
    integer, intent(inout) :: final_mode_index(:)
    real(kind=rp), intent(inout) :: mass_norm_final(:), mass_cross_final(:)
    logical, intent(out) :: ok
    type(vector_t), allocatable :: modes(:, :)
    type(vector_t) :: weights
    integer :: mode, component

    accepted = 0
    ok = .false.
    if (n < 1 .or. size(mass) < n .or. size(projected, 1) /= 3*n .or. &
         size(final_mode_index) < n_modes .or. size(mass_norm_final) < n_modes .or. &
         size(mass_cross_final) < n_modes) return
    allocate(modes(3, n_modes))
    do mode = 1, n_modes
       do component = 1, 3
          call modes(component, mode)%init(n)
          modes(component, mode)%x = projected((component-1)*n+1:component*n, mode)
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call modes(component, mode)%copy_from(HOST_TO_DEVICE, sync = .true.)
          end if
       end do
    end do
    call weights%init(n)
    weights%x = mass(1:n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call weights%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if
    call rom_mass_orthonormalize_device(modes, weights, tolerance, accepted, &
         final_mode_index, mass_norm_final, mass_cross_final, ok)
    if (ok) then
       do mode = 1, accepted
          do component = 1, 3
             if (NEKO_BCKND_DEVICE .eq. 1) then
                call modes(component, mode)%copy_from(DEVICE_TO_HOST, sync = .true.)
             end if
             projected((component-1)*n+1:component*n, mode) = modes(component, mode)%x
          end do
       end do
    end if
    call weights%free()
    do mode = 1, size(modes, 2)
       do component = 1, 3
          call modes(component, mode)%free()
       end do
    end do
    deallocate(modes)
  end subroutine rom_orthonormalize_projected_device

  !> Orthonormalise complete velocity--scalar columns on Neko's backend.
  !! The velocity and scalar arrays are one coupled POD basis: every
  !! projection, rejection, subtraction and normalisation is therefore
  !! applied to both arrays with the configured coupled mass product.
  subroutine rom_orthonormalize_projected_coupled_device(projected, scalar_basis, n, &
       n_modes, mass, scalar_weight, tolerance, accepted, final_mode_index, &
       mass_norm_final, mass_cross_final, ok)
    real(kind=rp), allocatable, intent(inout) :: projected(:, :), scalar_basis(:, :)
    integer, intent(in) :: n, n_modes
    real(kind=rp), intent(in) :: mass(:), scalar_weight, tolerance
    integer, intent(out) :: accepted
    integer, intent(inout) :: final_mode_index(:)
    real(kind=rp), intent(inout) :: mass_norm_final(:), mass_cross_final(:)
    logical, intent(out) :: ok
    type(vector_t), allocatable :: velocity_modes(:, :), scalar_modes(:)
    type(vector_t) :: weights
    integer :: mode, component

    accepted = 0
    ok = .false.
    if (n < 1 .or. n_modes < 1 .or. size(projected, 1) /= 3*n .or. &
         size(projected, 2) < n_modes .or. size(scalar_basis, 1) /= n .or. &
         size(scalar_basis, 2) < n_modes .or. size(mass) < n .or. &
         size(final_mode_index) < n_modes .or. size(mass_norm_final) < n_modes .or. &
         size(mass_cross_final) < n_modes) return

    allocate(velocity_modes(3,n_modes), scalar_modes(n_modes))
    do mode = 1, n_modes
       do component = 1, 3
          call velocity_modes(component, mode)%init(n)
          velocity_modes(component, mode)%x = projected((component-1)*n+1:component*n, mode)
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call velocity_modes(component, mode)%copy_from(HOST_TO_DEVICE, sync = .true.)
          end if
       end do
       call scalar_modes(mode)%init(n)
       scalar_modes(mode)%x = scalar_basis(:, mode)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call scalar_modes(mode)%copy_from(HOST_TO_DEVICE, sync = .true.)
       end if
    end do
    call weights%init(n)
    weights%x = mass(1:n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call weights%copy_from(HOST_TO_DEVICE, sync = .true.)
    end if

    call rom_mass_orthonormalize_coupled_device(velocity_modes, scalar_modes, &
         weights, scalar_weight, tolerance, accepted, final_mode_index, &
         mass_norm_final, mass_cross_final, ok)
    if (ok) then
       do mode = 1, accepted
          do component = 1, 3
             if (NEKO_BCKND_DEVICE .eq. 1) then
                call velocity_modes(component, mode)%copy_from(DEVICE_TO_HOST, sync = .true.)
             end if
             projected((component-1)*n+1:component*n, mode) = &
                  velocity_modes(component, mode)%x
          end do
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call scalar_modes(mode)%copy_from(DEVICE_TO_HOST, sync = .true.)
          end if
          scalar_basis(:, mode) = scalar_modes(mode)%x
       end do
       if (accepted < n_modes) then
          projected = projected(:, 1:accepted)
          scalar_basis = scalar_basis(:, 1:accepted)
       end if
    end if
    call weights%free()
    do mode = 1, size(scalar_modes)
       do component = 1, 3
          call velocity_modes(component, mode)%free()
       end do
       call scalar_modes(mode)%free()
    end do
    deallocate(velocity_modes, scalar_modes)
  end subroutine rom_orthonormalize_projected_coupled_device

  !> Multiply a mesh field by Neko's SEM diagonal on the active backend.
  !! The coefficient mirror is selected explicitly so a host `%x` write can
  !! never silently replace the device arithmetic on HIP/CUDA builds.
  subroutine rom_col2_B(field, fluid, n)
    type(field_t), intent(inout) :: field
    type(fluid_pnpn_t), intent(in) :: fluid
    integer, intent(in) :: n
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(field%x_d, fluid%c_Xh%B_d, n)
    else
       call col2(field%x, fluid%c_Xh%B, n)
    end if
  end subroutine rom_col2_B

  !> Divide a mesh field by Neko's SEM diagonal on the active backend.
  !! This is the final step of the mass-weighted gather/average
  !! u_c = B^{-1} GS(B u) used to recover a continuous nodal field.
  subroutine rom_col2_Binv(field, fluid, n)
    type(field_t), intent(inout) :: field
    type(fluid_pnpn_t), intent(in) :: fluid
    integer, intent(in) :: n
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(field%x_d, fluid%c_Xh%Binv_d, n)
    else
       call col2(field%x, fluid%c_Xh%Binv, n)
    end if
  end subroutine rom_col2_Binv

  !> Return a global SEM-weighted squared norm for a mesh field.
  function rom_weighted_norm(field, weight, weight_d, n) result(value)
    type(field_t), intent(in) :: field
    real(kind=rp), intent(in) :: weight(:,:,:,:)
    type(c_ptr), intent(in) :: weight_d
    integer, intent(in) :: n
    real(kind=rp) :: value
    if (NEKO_BCKND_DEVICE .eq. 1) then
       value = device_glsc3(field%x_d, weight_d, field%x_d, n)
    else
       value = glsc3(field%x, weight, field%x, n)
    end if
  end function rom_weighted_norm

  !> Receive one local velocity basis from the ROM state-recovery stream.
  !! Transport belongs to the ROM add-on; the generic simulation only exposes
  !! its existing recovery object and the adjoint stream callback.
  subroutine rom_receive_basis(recover, prototype, basis, n_modes, n_fields)
    class(state_recover_t), intent(inout) :: recover
    real(kind=rp), intent(in) :: prototype(:,:,:,:)
    real(kind=rp), allocatable, intent(out) :: basis(:, :)
    integer, intent(out) :: n_modes, n_fields
    real(kind=rp), allocatable :: field(:, :, :, :)
    integer :: field_size, mode, component, first_row, last_row

    call recover%receive_rom_basis_header(n_modes, n_fields)
    if (n_modes <= 0 .or. n_fields <= 0) then
       call neko_error('Invalid ROM basis dimensions received from Python')
    end if

    allocate(field, mold = prototype)
    field_size = size(field)
    allocate(basis(field_size * n_fields, n_modes))
    do mode = 1, n_modes
       do component = 1, n_fields
          call recover%receive_rom_basis_field(field)
          first_row = (component - 1) * field_size + 1
          last_row = component * field_size
          basis(first_row:last_row, mode) = reshape(field, [field_size])
       end do
    end do
    deallocate(field)
  end subroutine rom_receive_basis

  !> Receive independent forward-trial and adjoint-test bases.
  !! The second header is received only after the first payload has been
  !! consumed, matching the ordered MPI/ADIOS2 streaming protocol used by the
  !! Python peer.
  subroutine rom_receive_petrov_bases(recover, prototype, trial_basis, &
       trial_n_modes, trial_n_fields, test_basis, test_n_modes, test_n_fields)
    class(state_recover_t), intent(inout) :: recover
    real(kind=rp), intent(in) :: prototype(:,:,:,:)
    real(kind=rp), allocatable, intent(out) :: trial_basis(:, :), test_basis(:, :)
    integer, intent(out) :: trial_n_modes, trial_n_fields, test_n_modes, test_n_fields

    call rom_receive_basis(recover, prototype, trial_basis, trial_n_modes, trial_n_fields)
    call rom_receive_basis(recover, prototype, test_basis, test_n_modes, test_n_fields)
  end subroutine rom_receive_petrov_bases

  !> Resolve the ROM configuration from either supported case-file location.
  !! ROM examples historically used a top-level `rom` object, while the
  !! optimisation driver uses `optimization.solver.rom`.  Basis processing
  !! must read the same object as rom_problem_t.
  subroutine rom_get_parameters(simulation, parameters)
    type(simulation_t), intent(inout) :: simulation
    type(json_file), intent(out) :: parameters
    logical :: found

    call rom_config_get_json(simulation%neko_case%params, parameters, found)
    if (.not. found) call neko_error( &
         'ROM basis handoff requires optimization.solver.rom')
  end subroutine rom_get_parameters

  !> Apply the native scalar homogeneous Dirichlet projector to all columns.
  !! The scalar lift has already been subtracted by the basis streamer.  The
  !! host array is used only at the basis-reception boundary; on an accelerator
  !! each column is uploaded once, projected by Neko and downloaded once before
  !! the coupled velocity--scalar orthogonalisation.
  subroutine rom_restrict_scalar_basis(simulation, scalar_basis, n_modes)
    type(simulation_t), intent(inout) :: simulation
    real(kind=rp), allocatable, intent(inout) :: scalar_basis(:, :)
    integer, intent(in) :: n_modes
    type(field_t), pointer :: work
    integer :: work_index, n, mode

    if (n_modes < 1 .or. .not. associated(simulation%scalars)) return
    n = simulation%scalars%scalar_fields(1)%scalar%s%dof%size()
    if (size(scalar_basis, 1) /= n .or. size(scalar_basis, 2) < n_modes) then
       call neko_error('ROM scalar basis shape is inconsistent with scalar dofs')
    end if
    call neko_scratch_registry%request_field(work, work_index, &
         simulation%scalars%scalar_fields(1)%scalar%s%dof, .true.)
    select type (scalar_scheme => simulation%scalars%scalar_fields(1)%scalar)
    type is (scalar_pnpn_t)
       do mode = 1, n_modes
          work%x = reshape(scalar_basis(:, mode), shape(work%x))
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call work%copy_from(HOST_TO_DEVICE, sync = .true.)
          end if
          call scalar_scheme%bc_projector%apply(work%x, n)
          if (NEKO_BCKND_DEVICE .eq. 1) then
             call work%copy_from(DEVICE_TO_HOST, sync = .true.)
          end if
          scalar_basis(:, mode) = reshape(work%x, [n])
       end do
    class default
       call neko_error('ROM scalar basis restriction requires scalar_pnpn')
    end select
    call neko_scratch_registry%relinquish_field(work_index)
    nullify(work)
  end subroutine rom_restrict_scalar_basis

  !> Apply the selected ROM handoff independently to an adjoint-derived test
  !! basis.  No forward lift is subtracted and no forward/test
  !! orthogonalisation is performed.  The selected direct-POD path reports
  !! continuity and BC residuals without pressure correction; the legacy
  !! pressure path remains available when strict admissibility is required.
  subroutine rom_project_test_basis(simulation, basis, n_modes, n_fields)
    type(simulation_t), intent(inout) :: simulation
    real(kind=rp), allocatable, intent(inout) :: basis(:, :)
    integer, intent(inout) :: n_modes
    integer, intent(in), optional :: n_fields
    if (.not. allocated(simulation%adjoint_case%fluid_adj)) call neko_error( &
         'Petrov test-basis projection requires an initialized adjoint fluid')
    select type (adjoint_fluid => simulation%adjoint_case%fluid_adj)
    type is (adjoint_fluid_pnpn_t)
       call rom_project_basis(simulation, basis, n_modes, n_fields, &
            adjoint_fluid%bcs_vel_projector)
    class default
       call neko_error('Petrov test-basis projection requires adjoint_fluid_pnpn')
    end select
  end subroutine rom_project_test_basis

  !> Select the basis handoff used by the ROM.
  !!
  !! `direct_pod` preserves DNS-derived velocity columns apart from the
  !! necessary homogeneous velocity-boundary restriction and the required
  !! coupled mass orthonormalisation.  It does not solve a pressure
  !! correction problem.  The legacy pressure projection remains available
  !! as `pressure` (and as the historical `velocity_inlet_xmin` value).
  !!
  !! The direct path is intentionally not an admissibility waiver.  It
  !! measures Neko's strong and assembled weak continuity residuals before
  !! and after the boundary treatment, and fails when
  !! `strict_admissibility` is enabled and those residuals exceed the
  !! configured tolerance.  With strict mode disabled it emits an explicit
  !! engineering-approximation warning.
  subroutine rom_project_basis(simulation, basis, n_modes, n_fields, velocity_projector)
    type(simulation_t), intent(inout) :: simulation
    real(kind=rp), allocatable, intent(inout) :: basis(:, :)
    integer, intent(inout) :: n_modes
    integer, intent(in), optional :: n_fields
    class(vector_bc_projector_t), target, intent(inout), optional :: velocity_projector
    character(len=:), allocatable :: projection
    type(json_file) :: rom_parameters

    call rom_get_parameters(simulation, rom_parameters)
    call json_get_or_default(rom_parameters, 'admissibility_projection', &
         projection, 'pressure')
    select case (trim(projection))
    case ('direct_pod', 'direct', 'none')
       call rom_direct_pod_basis(simulation, basis, n_modes, n_fields, &
            velocity_projector)
    case ('pressure', 'legacy', 'velocity_inlet_xmin')
       call rom_project_basis_pressure(simulation, basis, n_modes, n_fields, &
            velocity_projector)
    case default
       call neko_error('ROM admissibility_projection must be direct_pod or pressure')
    end select
  end subroutine rom_project_basis

  !> Handoff DNS-derived POD columns without a pressure correction.
  !!
  !! The velocity columns are inspected in their received form first.  The
  !! strong SEM divergence, Neko's assembled weak continuity residual and the
  !! native homogeneous velocity-BC residual are all evaluated before any
  !! boundary restriction.  Only columns that fail the homogeneous BC check
  !! are restricted; no gather-scatter averaging and no pressure solve is
  !! performed on this path.  Coupled velocity--scalar columns are kept in
  !! their original ordering and are orthonormalised together afterwards.
  subroutine rom_direct_pod_basis(simulation, basis, n_modes, n_fields, &
       velocity_projector)
    type(simulation_t), intent(inout) :: simulation
    real(kind=rp), allocatable, intent(inout) :: basis(:, :)
    integer, intent(inout) :: n_modes
    integer, intent(in), optional :: n_fields
    class(vector_bc_projector_t), target, intent(inout), optional :: velocity_projector
    type(field_t), pointer :: u_work, v_work, w_work, divergence
    type(field_t), pointer :: bc_u, bc_v, bc_w
    type(field_t), pointer :: p_work, p_rhs, f_x, f_y, f_z
    type(field_t), pointer :: u_zero, v_zero, w_zero
    integer :: u_index, v_index, w_index, divergence_index
    integer :: bc_u_index, bc_v_index, bc_w_index
    integer :: p_work_index, p_rhs_index, f_x_index, f_y_index, f_z_index
    integer :: u_zero_index, v_zero_index, w_zero_index
    integer :: fields, n, mode, accepted, diagnostic_unit
    real(kind=rp), allocatable :: projected(:, :), scalar_basis(:, :), mass(:)
    real(kind=rp), allocatable :: strong_before(:), weak_before(:), bc_before(:)
    real(kind=rp), allocatable :: strong_after(:), weak_after(:), bc_after(:)
    real(kind=rp), allocatable :: mass_norm_final(:), mass_cross_final(:)
    integer, allocatable :: final_mode_index(:)
    logical, allocatable :: boundary_modified(:)
    real(kind=rp) :: admissibility_tolerance, scalar_weight
    logical :: strict_admissibility, failed, orthonormal, diagnostic_open
    logical :: modified
    character(len=:), allocatable :: diagnostics_file
    type(json_file) :: rom_parameters
    type(c_ptr) :: event
    type(fluid_pnpn_t), pointer :: fluid

    select type (fluid_class => simulation%fluid)
    type is (fluid_pnpn_t)
       fluid => fluid_class
    class default
       call neko_error('direct_pod ROM handoff requires a Pn/Pn fluid')
    end select

    fields = 3
    if (present(n_fields)) fields = n_fields
    if (fields /= 3 .and. fields /= 4) call neko_error( &
         'direct_pod ROM basis must contain three velocity fields or one scalar too')
    if (n_modes < 1) call neko_error('direct_pod ROM basis contains no modes')

    n = fluid%u%dof%size()
    if (size(basis, 1) /= fields*n .or. size(basis, 2) /= n_modes) then
       call neko_error('direct_pod ROM basis shape is inconsistent with the dof map')
    end if
    if (fields == 4 .and. .not. associated(simulation%scalars)) then
       call neko_error('direct_pod scalar basis supplied to a scalar-free simulation')
    end if
    call rom_get_parameters(simulation, rom_parameters)
    call json_get_or_default(rom_parameters, 'admissibility_tolerance', &
         admissibility_tolerance, 1.0e-8_rp)
    call json_get_or_default(rom_parameters, 'strict_admissibility', &
         strict_admissibility, .false.)
    call json_get_or_default(rom_parameters, 'projection_diagnostics_file', &
         diagnostics_file, 'rom_direct_pod_diagnostics.csv')
    call rom_basis_diagnostics_path(simulation, diagnostics_file)
    call json_get_or_default(rom_parameters, 'scalar_metric_weight', &
         scalar_weight, 1.0_rp)
    if (admissibility_tolerance <= 0.0_rp .or. scalar_weight <= 0.0_rp) then
       call neko_error('invalid direct_pod ROM admissibility settings')
    end if

    allocate(projected(3*n,n_modes), mass(n))
    projected = basis(1:3*n,:)
    mass = reshape(fluid%c_Xh%B, [n])
    if (fields == 4) then
       allocate(scalar_basis(n,n_modes))
       scalar_basis = basis(3*n+1:4*n,:)
    end if

    call neko_scratch_registry%request_field(u_work, u_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(v_work, v_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(w_work, w_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(divergence, divergence_index, &
         fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(bc_u, bc_u_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(bc_v, bc_v_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(bc_w, bc_w_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(p_work, p_work_index, &
         fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(p_rhs, p_rhs_index, &
         fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(f_x, f_x_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(f_y, f_y_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(f_z, f_z_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(u_zero, u_zero_index, &
         fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(v_zero, v_zero_index, &
         fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(w_zero, w_zero_index, &
         fluid%dm_Xh, .true.)
    call field_rzero(u_zero); call field_rzero(v_zero); call field_rzero(w_zero)

    allocate(strong_before(n_modes), weak_before(n_modes), bc_before(n_modes), &
         strong_after(n_modes), weak_after(n_modes), bc_after(n_modes), &
         boundary_modified(n_modes))
    strong_before = 0.0_rp; weak_before = 0.0_rp; bc_before = 0.0_rp
    strong_after = 0.0_rp; weak_after = 0.0_rp; bc_after = 0.0_rp
    boundary_modified = .false.
    failed = .false.
    diagnostic_open = .false.
    event = C_NULL_PTR

    if (pe_rank == 0) then
       open(newunit=diagnostic_unit, file=trim(diagnostics_file), &
            status='replace', action='write')
       write(diagnostic_unit, '(A)') &
            'mode,strong_divergence_before,weak_continuity_before,' // &
            'bc_residual_before,strong_divergence_after,weak_continuity_after,' // &
            'bc_residual_after,boundary_modified'
       diagnostic_open = .true.
    end if

    do mode = 1, n_modes
       u_work%x = reshape(projected(1:n,mode), shape(u_work%x))
       v_work%x = reshape(projected(n+1:2*n,mode), shape(v_work%x))
       w_work%x = reshape(projected(2*n+1:3*n,mode), shape(w_work%x))
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call u_work%copy_from(HOST_TO_DEVICE, sync=.true.)
          call v_work%copy_from(HOST_TO_DEVICE, sync=.true.)
          call w_work%copy_from(HOST_TO_DEVICE, sync=.true.)
       end if

       ! These are measurements of the original incoming POD column.  In
       ! particular, no B/Binv averaging is allowed before this evidence is
       ! collected.
       call rom_divergence_norm(fluid, u_work, v_work, w_work, &
            divergence, n, strong_before(mode))
       call rom_bc_residual_norm(fluid, u_work, v_work, w_work, &
            bc_u, bc_v, bc_w, n, bc_before(mode), velocity_projector)
       call rom_weak_continuity_norm(fluid, u_work, v_work, w_work, &
            p_work, p_rhs, f_x, f_y, f_z, u_zero, v_zero, w_zero, event, n, &
            weak_before(mode))

       modified = bc_before(mode) > admissibility_tolerance
       if (modified) then
          call rom_apply_velocity_bc(fluid, velocity_projector, &
               u_work%x, v_work%x, w_work%x, n)
          boundary_modified(mode) = .true.
       end if
       call rom_divergence_norm(fluid, u_work, v_work, w_work, &
            divergence, n, strong_after(mode))
       call rom_bc_residual_norm(fluid, u_work, v_work, w_work, &
            bc_u, bc_v, bc_w, n, bc_after(mode), velocity_projector)
       call rom_weak_continuity_norm(fluid, u_work, v_work, w_work, &
            p_work, p_rhs, f_x, f_y, f_z, u_zero, v_zero, w_zero, event, n, &
            weak_after(mode))
       if (weak_after(mode) > admissibility_tolerance .or. &
            bc_after(mode) > admissibility_tolerance) failed = .true.

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call u_work%copy_from(DEVICE_TO_HOST, sync=.true.)
          call v_work%copy_from(DEVICE_TO_HOST, sync=.true.)
          call w_work%copy_from(DEVICE_TO_HOST, sync=.true.)
       end if
       projected(1:n,mode) = reshape(u_work%x,[n])
       projected(n+1:2*n,mode) = reshape(v_work%x,[n])
       projected(2*n+1:3*n,mode) = reshape(w_work%x,[n])
       if (diagnostic_open) write(diagnostic_unit, &
            '(I0,6(",",ES24.16),",",L1)') mode, strong_before(mode), &
            weak_before(mode), bc_before(mode), strong_after(mode), &
            weak_after(mode), bc_after(mode), boundary_modified(mode)
    end do
    if (diagnostic_open) then
       flush(diagnostic_unit)
       close(diagnostic_unit)
    end if

    call neko_scratch_registry%relinquish_field([u_index, v_index, w_index, &
         divergence_index, bc_u_index, bc_v_index, bc_w_index, p_work_index, &
         p_rhs_index, f_x_index, f_y_index, f_z_index, u_zero_index, &
         v_zero_index, w_zero_index])
    nullify(u_work, v_work, w_work, divergence, bc_u, bc_v, bc_w)
    nullify(p_work, p_rhs, f_x, f_y, f_z, u_zero, v_zero, w_zero)

    if (strict_admissibility .and. failed) then
       if (fields == 4) deallocate(scalar_basis)
       deallocate(projected, mass, strong_before, weak_before, bc_before, &
            strong_after, weak_after, bc_after, boundary_modified)
       call neko_error('direct_pod ROM basis failed the continuity/BC gate; see '&
            //trim(diagnostics_file))
    end if
    if (failed) call neko_log%warning('direct_pod ROM basis has nonzero '&
         //'assembled continuity residual; proceeding as an engineering '&
         //'pressure-free approximation. It is not scientifically admissible; see '&
         //trim(diagnostics_file))

    allocate(mass_norm_final(n_modes), mass_cross_final(n_modes), &
         final_mode_index(n_modes))
    if (fields == 4) then
       ! Scalar restriction is kept separate from the velocity evidence and
       ! preserves the coupled mode ordering used by the objective operators.
       call rom_restrict_scalar_basis(simulation, scalar_basis, n_modes)
       call rom_orthonormalize_projected_coupled_device(projected, scalar_basis, &
            n, n_modes, mass, scalar_weight, admissibility_tolerance, accepted, &
            final_mode_index, mass_norm_final, mass_cross_final, orthonormal)
    else
       call rom_orthonormalize_projected_device(projected, n, n_modes, mass, &
            admissibility_tolerance, accepted, final_mode_index, mass_norm_final, &
            mass_cross_final, orthonormal)
    end if
    if (.not. orthonormal .or. accepted < 1) then
       if (fields == 4) deallocate(scalar_basis)
       deallocate(projected, mass, strong_before, weak_before, bc_before, &
            strong_after, weak_after, bc_after, boundary_modified)
       call neko_error('direct_pod ROM basis is not mass-orthonormal')
    end if

    n_modes = accepted
    if (fields == 4) then
       deallocate(basis)
       allocate(basis(4*n, n_modes))
       basis(1:3*n,:) = projected(:,1:n_modes)
       basis(3*n+1:4*n,:) = scalar_basis(:,1:n_modes)
       deallocate(scalar_basis)
    else
       deallocate(basis)
       allocate(basis(3*n, n_modes))
       basis = projected(:,1:n_modes)
    end if
    deallocate(projected, mass, strong_before, weak_before, bc_before, &
         strong_after, weak_after, bc_after, boundary_modified, mass_norm_final, &
         mass_cross_final, final_mode_index)
  end subroutine rom_direct_pod_basis

  !> Read projection settings and enforce the ROM basis admissibility policy.
  subroutine rom_project_basis_pressure(simulation, basis, n_modes, n_fields, velocity_projector)
    type(simulation_t), intent(inout) :: simulation
    real(kind=rp), allocatable, intent(inout) :: basis(:, :)
    integer, intent(inout) :: n_modes
    integer, intent(in), optional :: n_fields
    class(vector_bc_projector_t), target, intent(inout), optional :: velocity_projector
    real(kind=rp), allocatable :: velocity_basis(:, :), scalar_basis(:, :)
    integer, allocatable :: kept_indices(:)
    real(kind=rp) :: scalar_weight
    integer :: fields, old_modes, n
    type(json_file) :: rom_parameters

    fields = 3
    if (present(n_fields)) fields = n_fields
    if (fields /= 3 .and. fields /= 4) then
       call neko_error('ROM basis must contain three velocity fields or one scalar too')
    end if
    if (fields == 4) then
       n = simulation%fluid%u%dof%size()
       old_modes = n_modes
       if (size(basis,1) /= 4*n) call neko_error('coupled ROM basis has wrong row count')
       allocate(velocity_basis(3*n, old_modes), scalar_basis(n, old_modes))
       velocity_basis = basis(1:3*n,1:old_modes)
       scalar_basis = basis(3*n+1:4*n,1:old_modes)
       call rom_get_parameters(simulation, rom_parameters)
       call json_get_or_default(rom_parameters, 'scalar_metric_weight', &
            scalar_weight, 1.0_rp)
       if (scalar_weight <= 0.0_rp) call neko_error( &
            'ROM scalar_metric_weight must be positive')
       ! Restrict the scalar component before the coupled MGS step.  This is
       ! deliberately done on the isolated scalar columns so the velocity and
       ! scalar columns still have exactly the same mode ordering.
       call rom_restrict_scalar_basis(simulation, scalar_basis, old_modes)
       call rom_project_basis_velocity(simulation, velocity_basis, n_modes, &
            kept_indices, scalar_basis, scalar_weight, velocity_projector)
       if (n_modes < 1) call neko_error('coupled projection removed every ROM mode')
       deallocate(basis)
       allocate(basis(4*n,n_modes))
       basis(1:3*n,:) = velocity_basis(:,1:n_modes)
       basis(3*n+1:4*n,:) = scalar_basis(:,1:n_modes)
       deallocate(velocity_basis, scalar_basis)
       if (allocated(kept_indices)) deallocate(kept_indices)
       return
    end if
    call rom_project_basis_velocity(simulation, basis, n_modes, &
         velocity_projector = velocity_projector)
  end subroutine rom_project_basis_pressure

  !> Apply the velocity-only admissibility handoff used by the common basis.
  subroutine rom_project_basis_velocity(simulation, basis, n_modes, kept_indices, &
       scalar_basis, scalar_weight, velocity_projector)
    type(simulation_t), intent(inout) :: simulation
    real(kind=rp), allocatable, intent(inout) :: basis(:, :)
    integer, intent(inout) :: n_modes
    integer, allocatable, intent(out), optional :: kept_indices(:)
    real(kind=rp), allocatable, intent(inout), optional :: scalar_basis(:, :)
    real(kind=rp), intent(in), optional :: scalar_weight
    class(vector_bc_projector_t), target, intent(inout), optional :: velocity_projector
    real(kind=rp) :: pressure_abs_tol, pressure_rel_tol
    integer :: pressure_max_iter, projection_outer_iterations
    real(kind=rp) :: admissibility_tolerance
    logical :: strict_admissibility
    character(len=:), allocatable :: diagnostics_file
    type(json_file) :: rom_parameters

    call rom_get_parameters(simulation, rom_parameters)
    call json_get_or_default(rom_parameters, &
         'projection_pressure_absolute_tolerance', pressure_abs_tol, &
         1.0e-10_rp)
    call json_get_or_default(rom_parameters, &
         'projection_pressure_relative_tolerance', pressure_rel_tol, &
         0.0_rp)
    call json_get_or_default(rom_parameters, &
         'projection_pressure_max_iterations', pressure_max_iter, &
         simulation%fluid%ksp_prs%max_iter)
    call json_get_or_default(rom_parameters, &
         'projection_outer_iterations', projection_outer_iterations, 1)
    call json_get_or_default(rom_parameters, &
         'admissibility_tolerance', admissibility_tolerance, 1.0e-8_rp)
    call json_get_or_default(rom_parameters, &
         'projection_diagnostics_file', diagnostics_file, &
         'rom_projection_diagnostics.csv')
    call rom_basis_diagnostics_path(simulation, diagnostics_file)
    call json_get_or_default(rom_parameters, 'strict_admissibility', &
         strict_admissibility, .false.)
    if (pressure_abs_tol <= 0.0_rp .or. pressure_rel_tol < 0.0_rp .or. &
         pressure_max_iter < 1 .or. projection_outer_iterations < 1) then
       call neko_error('invalid ROM projection pressure-solver settings')
    end if
    select type (fluid => simulation%fluid)
    type is (fluid_pnpn_t)
       call project_rom_basis_pnpn(fluid, basis, n_modes, pressure_abs_tol, &
            pressure_rel_tol, pressure_max_iter, admissibility_tolerance, &
            projection_outer_iterations, diagnostics_file, strict_admissibility, &
            kept_indices, scalar_basis, scalar_weight, velocity_projector)
    class default
       call neko_error('ROM admissibility projection requires a Pn/Pn fluid')
    end select
  end subroutine rom_project_basis_velocity

  !> Apply the selected primal or adjoint velocity boundary projector.
  !! The optional projector is supplied only for the adjoint test-basis
  !! handoff; the forward trial path continues to use the fluid projector.
  subroutine rom_apply_velocity_bc(fluid, velocity_projector, u, v, w, n)
    type(fluid_pnpn_t), intent(inout) :: fluid
    class(vector_bc_projector_t), target, intent(inout), optional :: velocity_projector
    real(kind=rp), intent(inout) :: u(:,:,:,:), v(:,:,:,:), w(:,:,:,:)
    integer, intent(in) :: n

    if (present(velocity_projector)) then
       call velocity_projector%apply(u, v, w, n)
    else
       call fluid%bcs_vel_projector%apply(u, v, w, n)
    end if
  end subroutine rom_apply_velocity_bc

  !> Project a streamed Pn/Pn basis using Neko's constrained Schur operator.
  subroutine project_rom_basis_pnpn(fluid, basis, n_modes, projection_abs_tol, &
       projection_rel_tol, projection_max_iter, admissibility_tolerance, &
       projection_outer_iterations, diagnostics_file, strict_admissibility, &
       kept_indices, scalar_basis, scalar_weight, velocity_projector)
    type(fluid_pnpn_t), intent(inout), target :: fluid
    real(kind=rp), allocatable, intent(inout) :: basis(:, :)
    integer, intent(inout) :: n_modes
    real(kind=rp), intent(in) :: projection_abs_tol, projection_rel_tol
    integer, intent(in) :: projection_max_iter
    real(kind=rp), intent(in) :: admissibility_tolerance
    integer, intent(in) :: projection_outer_iterations
    character(len=*), intent(in) :: diagnostics_file
    logical, intent(in) :: strict_admissibility
    integer, allocatable, intent(out), optional :: kept_indices(:)
    real(kind=rp), allocatable, intent(inout), optional :: scalar_basis(:, :)
    real(kind=rp), intent(in), optional :: scalar_weight
    class(vector_bc_projector_t), target, intent(inout), optional :: velocity_projector
    type(field_t), pointer :: p_work, p_rhs, p_op, q_op, phi, u_work, v_work, w_work
    type(field_t), pointer :: u_zero, v_zero, w_zero, f_x, f_y, f_z
    type(field_t), pointer :: grad_x, grad_y, grad_z
    type(field_t), pointer :: divergence, bc_u, bc_v, bc_w
    integer :: p_work_index, p_rhs_index, p_op_index, q_op_index, phi_index
    integer :: u_work_index, v_work_index, w_work_index
    integer :: u_zero_index, v_zero_index, w_zero_index
    integer :: f_x_index, f_y_index, f_z_index
    integer :: grad_x_index, grad_y_index, grad_z_index
    integer :: divergence_index, bc_u_index, bc_v_index, bc_w_index
    type(ksp_monitor_t) :: ksp_result
    type(rom_projection_ax_t) :: projection_ax
    real(kind=rp), allocatable :: projected(:, :), mass(:)
    real(kind=rp), allocatable :: divergence_before(:), divergence_bc_input(:), &
         divergence_after(:), bc_before(:), bc_after(:), mass_norm_before(:), &
         mass_norm_raw_after(:), mass_cross_before(:), mass_cross_raw_after(:), &
         weak_divergence_before(:), weak_divergence_after_bc_input(:), &
         pressure_rhs_norm(:), pressure_residual(:), pressure_relative(:), &
         native_pressure_residual(:), native_pressure_relative(:), &
         weak_divergence_after(:), &
         mass_norm_final(:), mass_cross_final(:)
    real(kind=rp), allocatable :: pressure_ksp_start(:), pressure_ksp_final(:)
    integer, allocatable :: pressure_iterations(:), final_mode_index(:)
    logical, allocatable :: pressure_converged(:)
    real(kind=rp) :: mass_dot, divergence_plus
    real(kind=rp) :: pressure_solution_norm, correction_norm
    real(kind=rp) :: divergence_before_final_bc
    real(kind=rp) :: gradient_divergence_norm, operator_consistency_norm
    real(kind=rp) :: pressure_operator_norm
    real(kind=rp) :: saved_abs_tol, saved_rel_tol
    integer :: saved_max_iter
    integer :: n, mode, previous, accepted, kept, diagnostic_unit, projection_pass
    logical :: diagnostic_open, projection_failed, orthonormal
    character(len=512) :: pressure_log
    type(c_ptr) :: event
    real(kind=rp) :: coupled_scalar_weight

    if (n_modes < 1) call neko_error('ROM basis contains no modes')
    n = fluid%u%dof%size()
    if (size(basis, 1) /= 3 * n .or. size(basis, 2) /= n_modes) then
       call neko_error('ROM basis shape is inconsistent with the velocity dof map')
    end if
    if (present(scalar_basis)) then
       if (.not. present(scalar_weight) .or. size(scalar_basis, 1) /= n .or. &
            size(scalar_basis, 2) /= n_modes) then
          call neko_error('coupled ROM scalar basis shape is inconsistent')
       end if
       coupled_scalar_weight = scalar_weight
       if (coupled_scalar_weight <= 0.0_rp) call neko_error( &
            'ROM scalar_metric_weight must be positive')
    end if

    call neko_scratch_registry%request_field(p_work, p_work_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(p_rhs, p_rhs_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(p_op, p_op_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(q_op, q_op_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(phi, phi_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(u_work, u_work_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(v_work, v_work_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(w_work, w_work_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(u_zero, u_zero_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(v_zero, v_zero_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(w_zero, w_zero_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(f_x, f_x_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(f_y, f_y_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(f_z, f_z_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(grad_x, grad_x_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(grad_y, grad_y_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(grad_z, grad_z_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(divergence, divergence_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(bc_u, bc_u_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(bc_v, bc_v_index, fluid%dm_Xh, .true.)
    call neko_scratch_registry%request_field(bc_w, bc_w_index, fluid%dm_Xh, .true.)

    call field_rzero(u_zero)
    call field_rzero(v_zero)
    call field_rzero(w_zero)
    allocate(projected(3 * n, n_modes), mass(n))
    mass = reshape(fluid%c_Xh%B, [n])
    allocate(divergence_before(n_modes), divergence_bc_input(n_modes), &
         divergence_after(n_modes), bc_before(n_modes), bc_after(n_modes), &
         mass_norm_before(n_modes), mass_norm_raw_after(n_modes), &
         mass_cross_before(n_modes), mass_cross_raw_after(n_modes), &
         weak_divergence_before(n_modes), weak_divergence_after_bc_input(n_modes), &
         pressure_rhs_norm(n_modes), pressure_residual(n_modes), &
         pressure_relative(n_modes), native_pressure_residual(n_modes), &
         native_pressure_relative(n_modes), weak_divergence_after(n_modes), &
         mass_norm_final(n_modes), &
         mass_cross_final(n_modes), pressure_ksp_start(n_modes), &
         pressure_ksp_final(n_modes), pressure_iterations(n_modes), &
         pressure_converged(n_modes), final_mode_index(n_modes))
    divergence_before = 0.0_rp; divergence_bc_input = 0.0_rp
    divergence_after = 0.0_rp; bc_before = 0.0_rp; bc_after = 0.0_rp
    mass_norm_before = 0.0_rp; mass_norm_raw_after = 0.0_rp
    mass_cross_before = 0.0_rp; mass_cross_raw_after = 0.0_rp
    weak_divergence_before = 0.0_rp
    weak_divergence_after_bc_input = 0.0_rp
    pressure_rhs_norm = 0.0_rp; pressure_residual = 0.0_rp
    pressure_relative = 0.0_rp; native_pressure_residual = 0.0_rp
    native_pressure_relative = 0.0_rp; weak_divergence_after = 0.0_rp
    mass_norm_final = 0.0_rp
    mass_cross_final = 0.0_rp; pressure_ksp_start = 0.0_rp
    pressure_ksp_final = 0.0_rp; pressure_iterations = 0
    pressure_converged = .false.; final_mode_index = 0
    event = C_NULL_PTR
    diagnostic_open = .false.
    projection_failed = .false.
    if (pe_rank == 0) then
       open(newunit = diagnostic_unit, file = trim(diagnostics_file), &
            status = 'replace', action = 'write')
       write(diagnostic_unit, '(A)') &
            'mode,divergence_before,divergence_after_bc_input,' // &
            'divergence_after,bc_residual_before,bc_residual_after,' // &
            'mass_norm_before,mass_norm_raw_after,mass_cross_before,' // &
            'mass_cross_raw_after,' // &
            'weak_divergence_before,weak_divergence_after_bc_input,' // &
            'pressure_rhs_norm,' // &
            'constrained_schur_residual,constrained_schur_relative,' // &
            'native_pressure_poisson_residual,' // &
            'native_pressure_poisson_relative,pressure_ksp_start,' // &
            'pressure_ksp_final,weak_divergence_after,pressure_iterations,' // &
            'pressure_converged,' // &
            'mass_norm_final,mass_cross_final,final_mode_index'
       diagnostic_open = .true.
    end if

    saved_abs_tol = fluid%ksp_prs%abs_tol
    saved_rel_tol = fluid%ksp_prs%rel_tol
    saved_max_iter = fluid%ksp_prs%max_iter
    ! These are ROM-projection settings, not production time-step settings.
    ! Restore the production values immediately after this basis handoff.
    fluid%ksp_prs%abs_tol = projection_abs_tol
    fluid%ksp_prs%rel_tol = projection_rel_tol
    fluid%ksp_prs%max_iter = projection_max_iter
    call projection_ax%init(fluid, q_op, p_work, p_op, grad_x, grad_y, grad_z, &
         f_x, f_y, f_z, u_zero, v_zero, w_zero, velocity_projector)

    do mode = 1, n_modes
       u_work%x = reshape(basis(1:n, mode), shape(u_work%x))
       v_work%x = reshape(basis(n + 1:2 * n, mode), shape(v_work%x))
       w_work%x = reshape(basis(2 * n + 1:3 * n, mode), shape(w_work%x))
       if (NEKO_BCKND_DEVICE .eq. 1) then
          ! Basis reception is the explicit host/device boundary.  All
          ! subsequent gather, projection and residual work remains on device.
          call u_work%copy_from(HOST_TO_DEVICE, sync = .true.)
          call v_work%copy_from(HOST_TO_DEVICE, sync = .true.)
          call w_work%copy_from(HOST_TO_DEVICE, sync = .true.)
       end if

       ! fld files contain element-local nodal values.  Resolve shared GLL
       ! nodes with Neko's mass-weighted nodal average,
       ! u_c = Binv*GS(B*u), before using any weak operator.  The streamed
       ! values are pointwise, not mass-weighted, so B is applied first.
       ! Without this gather/average step, inter-element jumps enter the
       ! pressure RHS as artificial surface divergence even when the
       ! elementwise strong divergence is small.
       call rom_col2_B(u_work, fluid, n)
       call rom_col2_B(v_work, fluid, n)
       call rom_col2_B(w_work, fluid, n)
       call fluid%gs_Xh%op(u_work%x, v_work%x, w_work%x, n, GS_OP_ADD)
       call rom_col2_Binv(u_work, fluid, n)
       call rom_col2_Binv(v_work, fluid, n)
       call rom_col2_Binv(w_work, fluid, n)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call u_work%copy_from(DEVICE_TO_HOST, sync = .true.)
          call v_work%copy_from(DEVICE_TO_HOST, sync = .true.)
          call w_work%copy_from(DEVICE_TO_HOST, sync = .true.)
       end if
       basis(1:n, mode) = reshape(u_work%x, [n])
       basis(n + 1:2 * n, mode) = reshape(v_work%x, [n])
       basis(2 * n + 1:3 * n, mode) = reshape(w_work%x, [n])

       mass_norm_before(mode) = sqrt(max(0.0_rp, &
            rom_mass_dot(basis(:, mode), basis(:, mode), mass, n)))
       do previous = 1, mode - 1
          mass_dot = abs(rom_mass_dot(basis(:, mode), basis(:, previous), &
               mass, n))
          mass_cross_before(mode) = max(mass_cross_before(mode), mass_dot)
       end do
       call rom_divergence_norm(fluid, u_work, v_work, w_work, divergence, &
            n, divergence_before(mode))
       call rom_bc_residual_norm(fluid, u_work, v_work, w_work, bc_u, bc_v, &
            bc_w, n, bc_before(mode))

       ! Preserve the assembled weak continuity norm of the raw streamed
       ! column before imposing the homogeneous forward velocity restriction.
       ! This is distinct from the pointwise SEM strong divergence below.
       call field_col3(f_x, u_work, fluid%rho)
       call field_col3(f_y, v_work, fluid%rho)
       call field_col3(f_z, w_work, fluid%rho)
       call rom_col2_B(f_x, fluid, n)
       call rom_col2_B(f_y, fluid, n)
       call rom_col2_B(f_z, fluid, n)
       call field_rzero(p_work)
       call fluid%prs_res%compute(p_work, p_rhs, u_work, v_work, w_work, &
            u_zero, v_zero, w_zero, f_x, f_y, f_z, fluid%c_Xh, fluid%gs_Xh, &
            fluid%bc_prs_surface, fluid%bc_sym_surface, fluid%Ax_prs, &
            0.0_rp, 1.0_rp, u_zero, fluid%rho, event)
       if (.not. fluid%prs_dirichlet) then
          call ortho(p_rhs%x, fluid%glb_n_points, n)
       end if
       call fluid%gs_Xh%op(p_rhs%x, n, GS_OP_ADD)
       call fluid%bcs_prs_projector%apply(p_rhs%x, n)
       weak_divergence_before(mode) = sqrt(max(0.0_rp, &
            rom_weighted_norm(p_rhs, fluid%c_Xh%mult, fluid%c_Xh%mult_d, n)))

       ! Restrict the incoming column to the homogeneous forward velocity
       ! space before forming the pressure RHS. This uses Neko's resolved
       ! projector; it is not a boundary-node approximation.
       call rom_apply_velocity_bc(fluid, velocity_projector, u_work%x, v_work%x, &
            w_work%x, n)
       ! Keep the engineering fallback in the same homogeneous velocity space
       ! as the strict projected path.  The basis is deliberately not marked
       ! admissible here: the per-mode diagnostics remain the scientific gate.
       basis(1:n, mode) = reshape(u_work%x, [n])
       basis(n + 1:2 * n, mode) = reshape(v_work%x, [n])
       basis(2 * n + 1:3 * n, mode) = reshape(w_work%x, [n])
       call rom_divergence_norm(fluid, u_work, v_work, w_work, divergence, &
            n, divergence_bc_input(mode))

       ! The pressure residual computes the SEM weak divergence of the
       ! mass-weighted force f/rho.  The outer loop below is a projection
       ! iteration, not an increase of the pressure Krylov iteration limit.
       do projection_pass = 1, projection_outer_iterations
          call rom_divergence_norm(fluid, u_work, v_work, w_work, divergence, &
               n, divergence_after(mode))

       ! The pressure residual computes Neko's assembled weak continuity
       ! vector C(u).  Its production callers provide f with the quadrature
       ! mass already applied.  Therefore use f = rho*B*u here; passing
       ! pointwise rho*u would apply Binv a second time and create an
       ! artificial O(Binv) pressure correction.
       call field_col3(f_x, u_work, fluid%rho)
       call field_col3(f_y, v_work, fluid%rho)
       call field_col3(f_z, w_work, fluid%rho)
       call rom_col2_B(f_x, fluid, n)
       call rom_col2_B(f_y, fluid, n)
       call rom_col2_B(f_z, fluid, n)
       call field_rzero(p_work)
       call field_rzero(phi)
       call fluid%prs_res%compute(p_work, p_rhs, u_work, v_work, w_work, &
            u_zero, v_zero, w_zero, f_x, f_y, f_z, fluid%c_Xh, fluid%gs_Xh, &
            fluid%bc_prs_surface, fluid%bc_sym_surface, fluid%Ax_prs, &
            0.0_rp, 1.0_rp, u_zero, fluid%rho, event)
       ! Keep the native assembled weak continuity residual.  For the Neko
       ! Pn/Pn pair, C(grad_h(phi)) = A_prs(phi) (up to the pressure
       ! boundary restriction), where grad_h is the collocated strong SEM
       ! gradient used by the pressure Schur complement.  Do not replace C
       ! with a separately constructed B*divergence vector.
       if (.not. fluid%prs_dirichlet) then
          call ortho(p_rhs%x, fluid%glb_n_points, n)
       end if
       call fluid%gs_Xh%op(p_rhs%x, n, GS_OP_ADD)
       call fluid%bcs_prs_projector%apply(p_rhs%x, n)
       weak_divergence_after_bc_input(mode) = sqrt(max(0.0_rp, &
            rom_weighted_norm(p_rhs, fluid%c_Xh%mult, fluid%c_Xh%mult_d, n)))
       pressure_rhs_norm(mode) = sqrt(max(0.0_rp, &
            rom_weighted_norm(p_rhs, fluid%c_Xh%mult, fluid%c_Xh%mult_d, n)))
       call fluid%pc_prs%update()
       ksp_result = fluid%ksp_prs%solve(projection_ax, phi, p_rhs%x, n, &
            fluid%c_Xh, fluid%bcs_prs_projector, fluid%gs_Xh, &
            niter = fluid%ksp_prs%max_iter)
       pressure_iterations(mode) = ksp_result%iter
       pressure_converged(mode) = ksp_result%converged
       pressure_ksp_start(mode) = ksp_result%res_start
       pressure_ksp_final(mode) = ksp_result%res_final

       if (.not. fluid%prs_dirichlet) then
          call ortho(phi%x, fluid%glb_n_points, n)
       end if

       ! Measure the constrained Schur-complement residual independently of
       ! Krylov's volume-normalised monitor.  This is the actual pressure
       ! projection equation solved by this ROM path.
       call projection_ax%compute(p_op%x, phi%x, fluid%c_Xh, fluid%msh, &
            fluid%Xh)
       call fluid%gs_Xh%op(p_op%x, n, GS_OP_ADD)
       call fluid%bcs_prs_projector%apply(p_op%x, n)
       call field_add2s2(p_op, p_rhs, -1.0_rp)
       call fluid%bcs_prs_projector%apply(p_op%x, n)
       pressure_residual(mode) = sqrt(max(0.0_rp, &
            rom_weighted_norm(p_op, fluid%c_Xh%mult, fluid%c_Xh%mult_d, n)))
       pressure_relative(mode) = pressure_residual(mode) / &
            max(pressure_rhs_norm(mode), tiny(1.0_rp))

       ! Also report the residual of Neko's unmodified pressure Poisson
       ! operator.  A mismatch here is expected when the essential velocity
       ! restriction changes the Schur complement; it must not be confused
       ! with the residual of the constrained ROM projection equation.
       call fluid%Ax_prs%compute(p_work%x, phi%x, fluid%c_Xh, fluid%msh, &
            fluid%Xh)
       call fluid%gs_Xh%op(p_work%x, n, GS_OP_ADD)
       call fluid%bcs_prs_projector%apply(p_work%x, n)
       call field_add2s2(p_rhs, p_work, -1.0_rp)
       call fluid%bcs_prs_projector%apply(p_rhs%x, n)
       native_pressure_residual(mode) = sqrt(max(0.0_rp, &
            rom_weighted_norm(p_rhs, fluid%c_Xh%mult, fluid%c_Xh%mult_d, n)))
       native_pressure_relative(mode) = native_pressure_residual(mode) / &
            max(pressure_rhs_norm(mode), tiny(1.0_rp))
       pressure_solution_norm = sqrt(max(0.0_rp, &
            rom_weighted_norm(phi, fluid%c_Xh%B, fluid%c_Xh%B_d, n)))
       call fluid%Ax_prs%compute(p_work%x, phi%x, fluid%c_Xh, fluid%msh, &
            fluid%Xh)
       call fluid%gs_Xh%op(p_work%x, n, GS_OP_ADD)
       call fluid%bcs_prs_projector%apply(p_work%x, n)

       ! Form the same restricted strong SEM gradient used by the constrained
       ! Schur complement.  The gradient is pointwise, not mass-weighted, so
       ! use Neko's mass-weighted nodal average, Binv*GS(B*grad(phi)), before
       ! applying the native homogeneous velocity restriction.  This is the
       ! correction operator whose algebra is tested independently, not an
       ! unrelated Laplacian.
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call grad(grad_x%x_d, grad_y%x_d, grad_z%x_d, phi%x_d, fluid%c_Xh)
       else
          call grad(grad_x%x, grad_y%x, grad_z%x, phi%x, fluid%c_Xh)
       end if
       call rom_col2_B(grad_x, fluid, n)
       call rom_col2_B(grad_y, fluid, n)
       call rom_col2_B(grad_z, fluid, n)
       call fluid%gs_Xh%op(grad_x%x, grad_y%x, grad_z%x, n, GS_OP_ADD)
       call rom_col2_Binv(grad_x, fluid, n)
       call rom_col2_Binv(grad_y, fluid, n)
       call rom_col2_Binv(grad_z, fluid, n)
       call rom_apply_velocity_bc(fluid, velocity_projector, grad_x%x, grad_y%x, &
            grad_z%x, n)
       call div(divergence%x, grad_x%x, grad_y%x, grad_z%x, fluid%c_Xh)
       gradient_divergence_norm = sqrt(max(0.0_rp, &
            rom_weighted_norm(divergence, fluid%c_Xh%B, fluid%c_Xh%B_d, n)))
       ! Compare the constrained Schur operator with Neko's unmodified
       ! pressure Poisson operator.  The manufactured diagnostic proves the
       ! latter on compatible fields; this quantity records the boundary /
       ! gather-scatter difference for each adjoint mode.
       pressure_operator_norm = sqrt(max(0.0_rp, &
            rom_weighted_norm(p_work, fluid%c_Xh%mult, fluid%c_Xh%mult_d, n)))
       call projection_ax%compute(p_op%x, phi%x, fluid%c_Xh, fluid%msh, &
            fluid%Xh)
       call fluid%gs_Xh%op(p_op%x, n, GS_OP_ADD)
       call fluid%bcs_prs_projector%apply(p_op%x, n)
       call field_add2s2(p_op, p_work, -1.0_rp)
       call fluid%bcs_prs_projector%apply(p_op%x, n)
       operator_consistency_norm = sqrt(max(0.0_rp, &
            rom_weighted_norm(p_op, fluid%c_Xh%mult, fluid%c_Xh%mult_d, n)))
       correction_norm = sqrt(max(0.0_rp, &
            rom_weighted_norm(grad_x, fluid%c_Xh%B, fluid%c_Xh%B_d, n) + &
            rom_weighted_norm(grad_y, fluid%c_Xh%B, fluid%c_Xh%B_d, n) + &
            rom_weighted_norm(grad_z, fluid%c_Xh%B, fluid%c_Xh%B_d, n)))
       call field_add2s2(u_work, grad_x, -1.0_rp)
       call field_add2s2(v_work, grad_y, -1.0_rp)
       call field_add2s2(w_work, grad_z, -1.0_rp)
       ! A strong elementwise gradient is not single-valued at shared GLL
       ! nodes.  The corrected velocity is pointwise after subtraction, so
       ! resolve it back into Neko's continuous nodal space with the same
       ! mass-weighted average, Binv*GS(B*u), before measuring its divergence
       ! or applying BCs.
       call rom_col2_B(u_work, fluid, n)
       call rom_col2_B(v_work, fluid, n)
       call rom_col2_B(w_work, fluid, n)
       call fluid%gs_Xh%op(u_work%x, v_work%x, w_work%x, n, GS_OP_ADD)
       call rom_col2_Binv(u_work, fluid, n)
       call rom_col2_Binv(v_work, fluid, n)
       call rom_col2_Binv(w_work, fluid, n)
       call rom_divergence_norm(fluid, u_work, v_work, w_work, divergence, &
            n, divergence_before_final_bc)
       call rom_apply_velocity_bc(fluid, velocity_projector, u_work%x, v_work%x, &
            w_work%x, n)

       call rom_divergence_norm(fluid, u_work, v_work, w_work, divergence, &
            n, divergence_after(mode))

       ! Independently evaluate Neko's weak discrete divergence after the
       ! correction.  This is not the Krylov residual from the solved system:
       ! it rebuilds the pressure RHS from the corrected velocity and reports
       ! the operator's actual divergence in its native mass metric.
       call field_col3(f_x, u_work, fluid%rho)
       call field_col3(f_y, v_work, fluid%rho)
       call field_col3(f_z, w_work, fluid%rho)
       call rom_col2_B(f_x, fluid, n)
       call rom_col2_B(f_y, fluid, n)
       call rom_col2_B(f_z, fluid, n)
       call field_rzero(p_work)
       call fluid%prs_res%compute(p_work, p_rhs, u_work, v_work, w_work, &
            u_zero, v_zero, w_zero, f_x, f_y, f_z, fluid%c_Xh, fluid%gs_Xh, &
            fluid%bc_prs_surface, fluid%bc_sym_surface, fluid%Ax_prs, &
            0.0_rp, 1.0_rp, u_zero, fluid%rho, event)
       if (.not. fluid%prs_dirichlet) then
          call ortho(p_rhs%x, fluid%glb_n_points, n)
       end if
       call fluid%gs_Xh%op(p_rhs%x, n, GS_OP_ADD)
       call fluid%bcs_prs_projector%apply(p_rhs%x, n)
       weak_divergence_after(mode) = sqrt(max(0.0_rp, &
            rom_weighted_norm(p_rhs, fluid%c_Xh%mult, fluid%c_Xh%mult_d, n)))

       ! Keep the opposite sign as a diagnostic while the discrete operator
       ! pairing is being verified.  Restore the selected correction before
       ! proceeding; this is not an adaptive tolerance or iteration change.
       call field_add2s2(u_work, grad_x, 2.0_rp)
       call field_add2s2(v_work, grad_y, 2.0_rp)
       call field_add2s2(w_work, grad_z, 2.0_rp)
       call rom_apply_velocity_bc(fluid, velocity_projector, u_work%x, v_work%x, &
            w_work%x, n)
       call rom_divergence_norm(fluid, u_work, v_work, w_work, divergence, &
            n, divergence_plus)
       call field_add2s2(u_work, grad_x, -2.0_rp)
       call field_add2s2(v_work, grad_y, -2.0_rp)
       call field_add2s2(w_work, grad_z, -2.0_rp)
       call rom_apply_velocity_bc(fluid, velocity_projector, u_work%x, v_work%x, &
            w_work%x, n)
       call rom_bc_residual_norm(fluid, u_work, v_work, w_work, bc_u, bc_v, &
            bc_w, n, bc_after(mode))
       write(pressure_log, '(A,I0,A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') &
            'ROM projection mode ', mode, ': pass ', projection_pass, &
            ': divergence,pressure_residual,pressure_relative=', &
            divergence_after(mode), ',', pressure_residual(mode), ',', &
            pressure_relative(mode)
       call neko_log%message(trim(pressure_log))

       end do

       if (NEKO_BCKND_DEVICE .eq. 1) then
          ! Returning the projected basis to the Python/ADIOS-facing host
          ! array is a deliberate basis-reception boundary.
          call u_work%copy_from(DEVICE_TO_HOST, sync = .true.)
          call v_work%copy_from(DEVICE_TO_HOST, sync = .true.)
          call w_work%copy_from(DEVICE_TO_HOST, sync = .true.)
       end if
       projected(1:n, mode) = reshape(u_work%x, [n])
       projected(n + 1:2 * n, mode) = reshape(v_work%x, [n])
       projected(2 * n + 1:3 * n, mode) = reshape(w_work%x, [n])
       mass_norm_raw_after(mode) = sqrt(max(0.0_rp, &
            rom_mass_dot(projected(:, mode), projected(:, mode), mass, n)))
       do previous = 1, mode - 1
          mass_dot = abs(rom_mass_dot(projected(:, mode), projected(:, previous), &
               mass, n))
          mass_cross_raw_after(mode) = max(mass_cross_raw_after(mode), mass_dot)
       end do

       write(pressure_log, '(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') &
            'ROM projection mode ', mode, ': div(raw,bc,after)=', &
            divergence_before(mode), ',', divergence_bc_input(mode), ',', &
            divergence_after(mode)
       call neko_log%message(trim(pressure_log))
       write(pressure_log, '(A,I0,A,ES12.4,A,ES12.4)') &
            'ROM projection mode ', mode, ': pressure_solution/correction=', &
            pressure_solution_norm, ',', correction_norm
       call neko_log%message(trim(pressure_log))
       write(pressure_log, '(A,I0,A,ES12.4,A,ES12.4)') &
            'ROM projection mode ', mode, ': div(grad_phi)/Ax_consistency=', &
            gradient_divergence_norm, ',', operator_consistency_norm
       call neko_log%message(trim(pressure_log))
       write(pressure_log, '(A,I0,A,ES12.4)') &
            'ROM projection mode ', mode, ': pressure_operator_norm=', &
            pressure_operator_norm
       call neko_log%message(trim(pressure_log))
       write(pressure_log, '(A,I0,A,ES12.4,A,ES12.4)') &
            'ROM projection mode ', mode, ': mass_norm(raw/projected)=', &
            mass_norm_before(mode), ',', mass_norm_raw_after(mode)
       call neko_log%message(trim(pressure_log))
       write(pressure_log, '(A,I0,A,ES12.4)') &
            'ROM projection mode ', mode, ': weak_divergence_after=', &
            weak_divergence_after(mode)
       call neko_log%message(trim(pressure_log))
       write(pressure_log, '(A,I0,A,ES12.4,A,ES12.4)') &
            'ROM projection mode ', mode, ': weak_divergence(raw,bc_input)=', &
            weak_divergence_before(mode), ',', weak_divergence_after_bc_input(mode)
       call neko_log%message(trim(pressure_log))
       write(pressure_log, '(A,I0,A,ES12.4)') &
            'ROM projection mode ', mode, ': div(before_final_BC)=', &
            divergence_before_final_bc
       call neko_log%message(trim(pressure_log))
       write(pressure_log, '(A,I0,A,ES12.4)') &
            'ROM projection mode ', mode, ': div(opposite_sign)=', divergence_plus
       call neko_log%message(trim(pressure_log))
       write(pressure_log, '(A,I0,A,ES12.4,A,ES12.4)') &
            'ROM projection mode ', mode, ': BC(before,after)=', &
            bc_before(mode), ',', bc_after(mode)
       call neko_log%message(trim(pressure_log))
       write(pressure_log, '(A,I0,A,ES12.4,A,ES12.4,A,ES12.4)') &
            'ROM projection mode ', mode, ': pressure(rhs,res,rel)=', &
            pressure_rhs_norm(mode), ',', pressure_residual(mode), ',', &
            pressure_relative(mode)
       call neko_log%message(trim(pressure_log))
       write(pressure_log, '(A,I0,A,ES12.4,A,ES12.4)') &
            'ROM projection mode ', mode, ': native_pressure(res,rel)=', &
            native_pressure_residual(mode), ',', native_pressure_relative(mode)
       call neko_log%message(trim(pressure_log))
       ! Flush the per-mode evidence before any fail-closed error.  The final
       ! mass-orthogonalisation fields are intentionally zero on this path;
       ! raw norm/cross/BC/divergence and both pressure residuals remain valid.
       if (diagnostic_open) then
          write(diagnostic_unit, &
               '(I0,19(",",ES24.16),",",I0,",",L1,2(",",ES24.16),",",I0)') &
               mode, divergence_before(mode), divergence_bc_input(mode), &
               divergence_after(mode), bc_before(mode), bc_after(mode), &
               mass_norm_before(mode), mass_norm_raw_after(mode), &
               mass_cross_before(mode), mass_cross_raw_after(mode), &
               weak_divergence_before(mode), weak_divergence_after_bc_input(mode), &
               pressure_rhs_norm(mode), pressure_residual(mode), &
               pressure_relative(mode), native_pressure_residual(mode), &
               native_pressure_relative(mode), pressure_ksp_start(mode), &
               pressure_ksp_final(mode), weak_divergence_after(mode), &
               pressure_iterations(mode), pressure_converged(mode), &
               mass_norm_final(mode), mass_cross_final(mode), &
               final_mode_index(mode)
          flush(diagnostic_unit)
       end if
       if (.not. pressure_converged(mode)) then
          projection_failed = .true.
       end if
       if (weak_divergence_after(mode) > admissibility_tolerance .or. &
            bc_after(mode) > admissibility_tolerance) then
          projection_failed = .true.
       end if
    end do

    if (projection_failed) then
       if (diagnostic_open) then
          flush(diagnostic_unit)
          close(diagnostic_unit)
       end if
       fluid%ksp_prs%abs_tol = saved_abs_tol
       fluid%ksp_prs%rel_tol = saved_rel_tol
       fluid%ksp_prs%max_iter = saved_max_iter
       call neko_scratch_registry%relinquish_field([p_work_index, p_rhs_index, p_op_index, &
            q_op_index, phi_index, u_work_index, v_work_index, w_work_index, &
            u_zero_index, v_zero_index, w_zero_index, f_x_index, f_y_index, f_z_index, &
            grad_x_index, grad_y_index, grad_z_index, divergence_index, bc_u_index, &
            bc_v_index, bc_w_index])
       nullify(p_work, p_rhs, p_op, q_op, phi, u_work, v_work, w_work)
       nullify(u_zero, v_zero, w_zero, f_x, f_y, f_z, grad_x, grad_y, grad_z)
       nullify(divergence, bc_u, bc_v, bc_w)
       if (strict_admissibility) then
          deallocate(projected, mass, divergence_before, divergence_bc_input, &
               divergence_after, bc_before, bc_after, mass_norm_before, &
               mass_norm_raw_after, mass_cross_before, mass_cross_raw_after, &
               weak_divergence_before, weak_divergence_after_bc_input, &
               pressure_rhs_norm, pressure_residual, pressure_relative, &
               native_pressure_residual, native_pressure_relative, &
               weak_divergence_after, mass_norm_final, mass_cross_final, &
               pressure_ksp_start, pressure_ksp_final, pressure_iterations, &
               pressure_converged, final_mode_index)
          call neko_error('ROM admissibility projection failed; see per-mode diagnostics')
       else
          ! The engineering path is deliberately not a scientific admission:
          ! a failed pressure/BC gate remains visible in the warning and CSV.
          ! It must nevertheless remove collapsed columns before the coupled
          ! ROM mass is assembled.  Keeping a zero projected column produces
          ! a singular reduced system and hides the actual gate failure behind
          ! an unrelated linear-algebra error.  Strict mode never reaches this
          ! fallback.
          if (present(scalar_basis)) then
             call rom_orthonormalize_projected_coupled_device(projected, scalar_basis, &
                  n, n_modes, mass, coupled_scalar_weight, admissibility_tolerance, &
                  accepted, final_mode_index, mass_norm_final, mass_cross_final, orthonormal)
          else
             call rom_orthonormalize_projected_device(projected, n, n_modes, mass, &
                  admissibility_tolerance, accepted, final_mode_index, mass_norm_final, &
                  mass_cross_final, orthonormal)
          end if
          if (.not. orthonormal .or. accepted < 1) then
             deallocate(projected, mass, divergence_before, divergence_bc_input, &
                  divergence_after, bc_before, bc_after, mass_norm_before, &
                  mass_norm_raw_after, mass_cross_before, mass_cross_raw_after, &
                  weak_divergence_before, weak_divergence_after_bc_input, &
                  pressure_rhs_norm, pressure_residual, pressure_relative, &
                  native_pressure_residual, native_pressure_relative, &
                  weak_divergence_after, mass_norm_final, mass_cross_final, &
                  pressure_ksp_start, pressure_ksp_final, pressure_iterations, &
                  pressure_converged, final_mode_index)
             call neko_error('ROM engineering fallback removed every projected mode')
          end if
          basis = projected(:, 1:accepted)
          n_modes = accepted
          if (present(kept_indices)) then
             allocate(kept_indices(accepted))
             kept = 0
             do mode = 1, size(final_mode_index)
                if (final_mode_index(mode) > 0) then
                   kept = kept + 1
                   kept_indices(kept) = mode
                end if
             end do
          end if
          call neko_log%warning('ROM admissibility projection did not meet its '&
               //'scientific gate; continuing with the homogeneous-BC basis '&
               //'for engineering experimentation. The ROM is not validated; '&
               //'see ' // trim(diagnostics_file))
       end if
       deallocate(projected, mass, divergence_before, divergence_bc_input, &
            divergence_after, bc_before, bc_after, mass_norm_before, &
            mass_norm_raw_after, mass_cross_before, mass_cross_raw_after, &
            weak_divergence_before, weak_divergence_after_bc_input, &
            pressure_rhs_norm, pressure_residual, pressure_relative, &
            native_pressure_residual, native_pressure_relative, &
            weak_divergence_after, mass_norm_final, mass_cross_final, &
            pressure_ksp_start, pressure_ksp_final, pressure_iterations, &
            pressure_converged, final_mode_index)
       return
    end if

    ! Pressure projection changes the metric, and can make nearly dependent
    ! enrichment columns collapse.  Rebuild the common basis using the same
    ! global SEM mass as the Python POD combination.  The implementation is
    ! in an independent ROM module so this numerical gate is testable without
    ! constructing a simulation object.
    if (present(scalar_basis)) then
       call rom_orthonormalize_projected_coupled_device(projected, scalar_basis, n, &
            n_modes, mass, coupled_scalar_weight, admissibility_tolerance, accepted, &
            final_mode_index, mass_norm_final, mass_cross_final, orthonormal)
    else
       call rom_orthonormalize_projected_device(projected, n, n_modes, mass, &
            admissibility_tolerance, accepted, final_mode_index, mass_norm_final, &
            mass_cross_final, orthonormal)
    end if
    if (.not. orthonormal) then
       call neko_error('ROM projected basis is not mass-orthonormal; '&
            //'gather-scatter normalisation or reorthogonalisation failed')
    end if
    if (accepted /= n_modes) then
       basis = projected(:, 1:accepted)
    else
       basis = projected
    end if
    n_modes = accepted
    if (present(kept_indices)) then
       allocate(kept_indices(accepted))
       kept = 0
       do mode = 1, size(final_mode_index)
          if (final_mode_index(mode) > 0) then
             kept = kept + 1
             kept_indices(kept) = mode
          end if
       end do
    end if

    if (diagnostic_open) then
       do mode = 1, size(divergence_before)
          write(diagnostic_unit, '(I0,19(",",ES24.16),",",I0,",",L1,2(",",ES24.16),",",I0)') &
               mode, divergence_before(mode), divergence_bc_input(mode), &
               divergence_after(mode), bc_before(mode), bc_after(mode), &
               mass_norm_before(mode), mass_norm_raw_after(mode), &
               mass_cross_before(mode), mass_cross_raw_after(mode), &
               weak_divergence_before(mode), weak_divergence_after_bc_input(mode), &
               pressure_rhs_norm(mode), pressure_residual(mode), &
               pressure_relative(mode), native_pressure_residual(mode), &
               native_pressure_relative(mode), pressure_ksp_start(mode), &
               pressure_ksp_final(mode), weak_divergence_after(mode), &
               pressure_iterations(mode), &
               pressure_converged(mode), mass_norm_final(mode), &
               mass_cross_final(mode), final_mode_index(mode)
       end do
       flush(diagnostic_unit)
       close(diagnostic_unit)
    end if

    fluid%ksp_prs%abs_tol = saved_abs_tol
    fluid%ksp_prs%rel_tol = saved_rel_tol
    fluid%ksp_prs%max_iter = saved_max_iter

    call neko_scratch_registry%relinquish_field([p_work_index, p_rhs_index, p_op_index, &
         q_op_index, phi_index, u_work_index, v_work_index, w_work_index, &
         u_zero_index, v_zero_index, w_zero_index, f_x_index, f_y_index, f_z_index, &
         grad_x_index, grad_y_index, grad_z_index, divergence_index, bc_u_index, &
         bc_v_index, bc_w_index])
    nullify(p_work, p_rhs, p_op, q_op, phi, u_work, v_work, w_work)
    nullify(u_zero, v_zero, w_zero, f_x, f_y, f_z, grad_x, grad_y, grad_z)
    nullify(divergence, bc_u, bc_v, bc_w)
    deallocate(projected, mass, divergence_before, divergence_bc_input, &
         divergence_after, bc_before, bc_after, mass_norm_before, &
         mass_norm_raw_after, mass_cross_before, mass_cross_raw_after, &
         weak_divergence_before, weak_divergence_after_bc_input, &
         pressure_rhs_norm, pressure_residual, pressure_relative, &
         native_pressure_residual, native_pressure_relative, &
         weak_divergence_after, mass_norm_final, mass_cross_final, &
         pressure_ksp_start, &
         pressure_ksp_final, pressure_iterations, &
         pressure_converged, final_mode_index)
  end subroutine project_rom_basis_pnpn

  !> Measure the strong SEM divergence in the velocity mass norm.
  subroutine rom_divergence_norm(fluid, u, v, w, divergence, n, value)
    type(fluid_pnpn_t), intent(in) :: fluid
    type(field_t), intent(in) :: u, v, w
    type(field_t), intent(inout) :: divergence
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: value

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call div(divergence%x_d, u%x_d, v%x_d, w%x_d, fluid%c_Xh)
    else
       call div(divergence%x, u%x, v%x, w%x, fluid%c_Xh)
    end if
    value = sqrt(max(0.0_rp, rom_weighted_norm(divergence, fluid%c_Xh%B, &
         fluid%c_Xh%B_d, n)))
  end subroutine rom_divergence_norm

  !> Measure the residual of Neko's homogeneous velocity projector.
  subroutine rom_bc_residual_norm(fluid, u, v, w, bc_u, bc_v, bc_w, n, value, &
       velocity_projector)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(field_t), intent(in) :: u, v, w
    type(field_t), intent(inout) :: bc_u, bc_v, bc_w
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: value
    class(vector_bc_projector_t), target, intent(inout), optional :: velocity_projector
    real(kind=rp) :: norm_squared

    call field_copy(bc_u, u); call field_copy(bc_v, v); call field_copy(bc_w, w)
    call rom_apply_velocity_bc(fluid, velocity_projector, bc_u%x, bc_v%x, bc_w%x, n)
    call field_sub3(bc_u, bc_u, u, n)
    call field_sub3(bc_v, bc_v, v, n)
    call field_sub3(bc_w, bc_w, w, n)
    norm_squared = rom_weighted_norm(bc_u, fluid%c_Xh%B, fluid%c_Xh%B_d, n) + &
         rom_weighted_norm(bc_v, fluid%c_Xh%B, fluid%c_Xh%B_d, n) + &
         rom_weighted_norm(bc_w, fluid%c_Xh%B, fluid%c_Xh%B_d, n)
    value = sqrt(max(0.0_rp, norm_squared))
  end subroutine rom_bc_residual_norm

  !> Measure Neko's assembled weak continuity residual for a velocity field.
  !!
  !! This follows the production Pn/Pn pressure-residual path: the force is
  !! formed as rho*u, weighted by the SEM diagonal B, evaluated by
  !! `prs_res%compute`, gathered and pressure-BC restricted, and finally
  !! measured with Neko's multiplicity metric.  It is deliberately not a
  !! pressure solve and is therefore safe for the direct-POD handoff.
  subroutine rom_weak_continuity_norm(fluid, u, v, w, p_work, p_rhs, f_x, f_y, &
       f_z, u_zero, v_zero, w_zero, event, n, value)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(field_t), intent(inout) :: u, v, w
    type(field_t), intent(inout) :: p_work, p_rhs, f_x, f_y, f_z
    type(field_t), intent(inout) :: u_zero, v_zero, w_zero
    type(c_ptr), intent(inout) :: event
    integer, intent(in) :: n
    real(kind=rp), intent(out) :: value

    call field_col3(f_x, u, fluid%rho)
    call field_col3(f_y, v, fluid%rho)
    call field_col3(f_z, w, fluid%rho)
    call rom_col2_B(f_x, fluid, n)
    call rom_col2_B(f_y, fluid, n)
    call rom_col2_B(f_z, fluid, n)
    call field_rzero(p_work)
    call fluid%prs_res%compute(p_work, p_rhs, u, v, w, u_zero, v_zero, w_zero, &
         f_x, f_y, f_z, fluid%c_Xh, fluid%gs_Xh, fluid%bc_prs_surface, &
         fluid%bc_sym_surface, fluid%Ax_prs, 0.0_rp, 1.0_rp, u_zero, &
         fluid%rho, event)
    if (.not. fluid%prs_dirichlet) then
       call ortho(p_rhs%x, fluid%glb_n_points, n)
    end if
    call fluid%gs_Xh%op(p_rhs%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(p_rhs%x, n)
    value = sqrt(max(0.0_rp, rom_weighted_norm(p_rhs, fluid%c_Xh%mult, &
         fluid%c_Xh%mult_d, n)))
  end subroutine rom_weak_continuity_norm

  !> Place basis-admissibility diagnostics below the matching ROM outer output.
  !! @param[in] simulation Simulation carrying the current outer index.
  !! @param[inout] configured Configured diagnostic filename; its basename is
  !! retained while the ROM output directory is made authoritative.
  subroutine rom_basis_diagnostics_path(simulation, configured)
    type(simulation_t), intent(in) :: simulation
    character(len=:), allocatable, intent(inout) :: configured
    character(len=1024) :: basename, directory
    integer :: slash

    basename = trim(configured)
    slash = scan(trim(basename), '/', back=.true.)
    if (slash > 0) basename = basename(slash+1:)
    write(directory, '("rom_output/outer_",I3.3,"/basis_diagnostics/")') &
         simulation%current_design_iteration
    call mkdir('rom_output')
    call mkdir(trim(directory))
    configured = trim(directory)//trim(basename)
  end subroutine rom_basis_diagnostics_path

  !> Reset the simulation
end module rom_basis_handoff
