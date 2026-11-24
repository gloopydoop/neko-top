!> @file rom_lift.f90
!! @brief Spatially varying velocity lifting for the ROM add-on.
!!
!! The lifting field is constructed from the velocity fields already owned by
!! the Neko case and then passed through the native velocity boundary-condition
!! list.  Keeping the construction here makes the streamed POD state, reduced
!! operators and reconstructed physical state use the same field.
module rom_lift
  use case, only: case_t
  use data_streamer, only: data_streamer_t
  use field, only: field_t
  use rom_config, only: rom_config_t, rom_config_get_json
  use scalar_pnpn, only: scalar_pnpn_t
  use json_module, only: json_file
  use field_math, only: field_add2s2, field_copy
  use operators, only: grad
  use scratch_registry, only: neko_scratch_registry
  use num_types, only: rp
  use vector, only: vector_t
  use neko_ext, only: field_to_vector
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: DEVICE_TO_HOST
  use utils, only: neko_error
  implicit none
  private

  public :: rom_lift_t

  !> Spatial lift and its SEM derivatives.
  type, public :: rom_lift_t
     logical :: enabled = .false.
     logical :: scalar_enabled = .false.
     type(field_t) :: reference(3)
     type(field_t) :: velocity(3)
     type(field_t) :: gradient(3,3)
     type(field_t) :: scalar_reference
     type(field_t) :: scalar
     type(field_t) :: scalar_gradient(3)
   contains
     procedure, pass(this), public :: init => rom_lift_init
     procedure, pass(this), public :: free => rom_lift_free
     procedure, pass(this), public :: sync_host => rom_lift_sync_host
     procedure, pass(this), public :: copy_velocity => rom_lift_copy_velocity
     procedure, pass(this), public :: copy_scalar => rom_lift_copy_scalar
     procedure, pass(this), public :: stream_homogeneous => &
          rom_lift_stream_homogeneous
  end type rom_lift_t

contains

  !> Build the ROM lift from the case reference velocity and native velocity BCs.
  !!
  !! The initial/reference velocity is copied before boundary conditions are
  !! applied.  The resulting lift is therefore spatial, uses any expressions
  !! or field-valued BC objects configured by Neko, and is independent of the
  !! reduced coefficients.  The reference fields are retained for initial-state
  !! projection.
  !! @param[inout] this Lift storage.
  !! @param[inout] neko_case Initialised Neko case.
  subroutine rom_lift_init(this, neko_case)
    class(rom_lift_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(rom_config_t) :: rom_parameters
    type(json_file) :: unused_parameters
    logical :: found
    character(len=64) :: source
    integer :: component, direction

    call this%free()
    call rom_config_get_json(neko_case%params, unused_parameters, found)
    if (.not. found) return
    call rom_parameters%init(neko_case%params)
    source = trim(rom_parameters%lift_source)
    select case (trim(source))
    case ('initial_velocity', 'reference_velocity')
       ! The case velocity at initialisation is the reference solution for the
       ! lift.  The second spelling documents that a precomputed reference
       ! state can be supplied by the case initial-condition machinery.
    case default
       call neko_error('ROM lift source must be initial_velocity or reference_velocity')
    end select

    do component = 1, 3
       call this%reference(component)%init(neko_case%fluid%u%dof, &
            'rom_lift_reference')
       call this%velocity(component)%init(neko_case%fluid%u%dof, &
            'rom_lift_velocity')
    end do
    call field_copy(this%reference(1), neko_case%fluid%u)
    call field_copy(this%reference(2), neko_case%fluid%v)
    call field_copy(this%reference(3), neko_case%fluid%w)
    call field_copy(this%velocity(1), this%reference(1))
    call field_copy(this%velocity(2), this%reference(2))
    call field_copy(this%velocity(3), this%reference(3))

    ! Use the exact BC objects attached to the forward fluid scheme.  This is
    ! deliberately a field-level application so Neko handles host/device data,
    ! masks and time-dependent expression BCs in its normal path.
    call neko_case%fluid%bcs_vel%apply(this%velocity(1), this%velocity(2), &
         this%velocity(3), time = neko_case%time, strong = .true.)

    do component = 1, 3
       do direction = 1, 3
          call this%gradient(component,direction)%init(neko_case%fluid%u%dof, &
               'rom_lift_gradient')
       end do
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call grad(this%gradient(component,1)%x_d, &
               this%gradient(component,2)%x_d, &
               this%gradient(component,3)%x_d, &
               this%velocity(component)%x_d, neko_case%fluid%c_Xh)
       else
          call grad(this%gradient(component,1)%x, &
               this%gradient(component,2)%x, &
               this%gradient(component,3)%x, &
               this%velocity(component)%x, neko_case%fluid%c_Xh)
       end if
    end do

    ! The scalar lift follows the same reference-state convention as the
    ! velocity lift.  Its homogeneous boundary restriction is applied through
    ! scalar_pnpn's resolved Dirichlet projector, not through a ROM-specific
    ! list of boundary nodes.
    this%scalar_enabled = allocated(neko_case%scalars)
    if (this%scalar_enabled) then
       if (size(neko_case%scalars%scalar_fields) /= 1) then
          call neko_error('ROM scalar lift currently supports one scalar')
       end if
       call this%scalar_reference%init( &
            neko_case%scalars%scalar_fields(1)%scalar%s%dof, &
            'rom_scalar_lift_reference')
       call this%scalar%init(neko_case%scalars%scalar_fields(1)%scalar%s%dof, &
            'rom_scalar_lift')
       call field_copy(this%scalar_reference, &
            neko_case%scalars%scalar_fields(1)%scalar%s)
       call field_copy(this%scalar, this%scalar_reference)
       select type (scalar_scheme => neko_case%scalars%scalar_fields(1)%scalar)
       type is (scalar_pnpn_t)
          call scalar_scheme%bcs%apply(this%scalar, time=neko_case%time, &
               strong=.true.)
       class default
          call neko_error('ROM scalar lift requires scalar_pnpn')
       end select
       do direction = 1, 3
          call this%scalar_gradient(direction)%init( &
               neko_case%scalars%scalar_fields(1)%scalar%s%dof, &
               'rom_scalar_lift_gradient')
       end do
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call grad(this%scalar_gradient(1)%x_d, this%scalar_gradient(2)%x_d, &
               this%scalar_gradient(3)%x_d, this%scalar%x_d, &
               neko_case%scalars%scalar_fields(1)%scalar%c_Xh)
       else
          call grad(this%scalar_gradient(1)%x, this%scalar_gradient(2)%x, &
               this%scalar_gradient(3)%x, this%scalar%x, &
               neko_case%scalars%scalar_fields(1)%scalar%c_Xh)
       end if
    end if
    this%enabled = .true.
  end subroutine rom_lift_init

  !> Free all fields owned by the lifting object.
  !! @param[inout] this Lift storage.
  subroutine rom_lift_free(this)
    class(rom_lift_t), intent(inout) :: this
    integer :: component, direction
    do component = 1, 3
       call this%reference(component)%free()
       call this%velocity(component)%free()
       do direction = 1, 3
          call this%gradient(component,direction)%free()
       end do
    end do
    call this%scalar_reference%free()
    call this%scalar%free()
    do direction = 1, 3
       call this%scalar_gradient(direction)%free()
    end do
    this%scalar_enabled = .false.
    this%enabled = .false.
  end subroutine rom_lift_free

  !> Synchronise all lift fields to host memory for vector/reduction consumers.
  !! @param[inout] this Lift storage.
  subroutine rom_lift_sync_host(this)
    class(rom_lift_t), intent(inout) :: this
    integer :: component, direction
    if (.not. this%enabled) return
    if (NEKO_BCKND_DEVICE .eq. 1) then
       do component = 1, 3
          call this%reference(component)%copy_from(DEVICE_TO_HOST, sync = .true.)
          call this%velocity(component)%copy_from(DEVICE_TO_HOST, sync = .true.)
          do direction = 1, 3
             call this%gradient(component,direction)%copy_from( &
                  DEVICE_TO_HOST, sync = .true.)
          end do
       end do
       if (this%scalar_enabled) then
          call this%scalar_reference%copy_from(DEVICE_TO_HOST, sync = .true.)
          call this%scalar%copy_from(DEVICE_TO_HOST, sync = .true.)
          do direction = 1, 3
             call this%scalar_gradient(direction)%copy_from( &
                  DEVICE_TO_HOST, sync = .true.)
          end do
       end if
    end if
  end subroutine rom_lift_sync_host

  !> Copy the three spatial lift components into vector containers.
  !! @param[inout] this Lift storage.
  !! @param[inout] vectors Destination vectors.
  subroutine rom_lift_copy_velocity(this, vectors)
    class(rom_lift_t), intent(inout) :: this
    type(vector_t), intent(inout) :: vectors(3)
    integer :: component, n
    if (.not. this%enabled) call neko_error('ROM lift is not initialised')
    if (size(vectors) /= 3) call neko_error('ROM lift requires three velocity vectors')
    n = this%velocity(1)%dof%size()
    do component = 1, 3
       call vectors(component)%init(n)
       call field_to_vector(vectors(component), this%velocity(component))
    end do
  end subroutine rom_lift_copy_velocity

  !> Copy the scalar lift and its homogeneous reference into a vector.
  subroutine rom_lift_copy_scalar(this, vector)
    class(rom_lift_t), intent(inout) :: this
    type(vector_t), intent(inout) :: vector
    integer :: n
    if (.not. this%scalar_enabled) call neko_error('ROM scalar lift is not enabled')
    n = this%scalar%dof%size()
    call vector%init(n)
    call field_to_vector(vector, this%scalar)
  end subroutine rom_lift_copy_scalar

  !> Stream a physical velocity after subtracting the spatial lift.
  !! @param[inout] this Lift storage.
  !! @param[inout] dstream Neko/Python data streamer.
  !! @param[in] u Forward x-velocity field.
  !! @param[in] v Forward y-velocity field.
  !! @param[in] w Forward z-velocity field.
  subroutine rom_lift_stream_homogeneous(this, dstream, u, v, w, scalar)
    class(rom_lift_t), intent(inout) :: this
    type(data_streamer_t), intent(inout) :: dstream
    type(field_t), intent(inout) :: u, v, w
    type(field_t), intent(inout), optional :: scalar
    type(field_t), pointer :: work_u, work_v, work_w, work_s
    integer :: indices(4), n

    if (.not. this%enabled) then
       ! This is an explicit output boundary. The streamer consumes host
       ! buffers, so synchronise the current fields only when no lift scratch
       ! field is being used to perform the streamed subtraction.
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call u%copy_from(DEVICE_TO_HOST, sync = .true.)
          call v%copy_from(DEVICE_TO_HOST, sync = .true.)
          call w%copy_from(DEVICE_TO_HOST, sync = .true.)
          if (present(scalar)) call scalar%copy_from(DEVICE_TO_HOST, sync = .true.)
       end if
       call dstream%stream(u%x)
       call dstream%stream(v%x)
       call dstream%stream(w%x)
       if (present(scalar)) call dstream%stream(scalar%x)
       return
    end if
    n = u%dof%size()
    call neko_scratch_registry%request_field(work_u, indices(1), u%dof, .true.)
    call neko_scratch_registry%request_field(work_v, indices(2), u%dof, .true.)
    call neko_scratch_registry%request_field(work_w, indices(3), u%dof, .true.)
    if (present(scalar) .and. this%scalar_enabled) then
       call neko_scratch_registry%request_field(work_s, indices(4), scalar%dof, .true.)
    end if
    call field_copy(work_u, u)
    call field_copy(work_v, v)
    call field_copy(work_w, w)
    call field_add2s2(work_u, this%velocity(1), -1.0_rp, n)
    call field_add2s2(work_v, this%velocity(2), -1.0_rp, n)
    call field_add2s2(work_w, this%velocity(3), -1.0_rp, n)
    if (present(scalar) .and. this%scalar_enabled) then
       call field_copy(work_s, scalar)
       call field_add2s2(work_s, this%scalar, -1.0_rp, scalar%dof%size())
    end if
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call work_u%copy_from(DEVICE_TO_HOST, sync = .true.)
       call work_v%copy_from(DEVICE_TO_HOST, sync = .true.)
       call work_w%copy_from(DEVICE_TO_HOST, sync = .true.)
       if (present(scalar) .and. this%scalar_enabled) &
           call work_s%copy_from(DEVICE_TO_HOST, sync = .true.)
    end if
    call dstream%stream(work_u%x)
    call dstream%stream(work_v%x)
    call dstream%stream(work_w%x)
    if (present(scalar) .and. this%scalar_enabled) call dstream%stream(work_s%x)
    call neko_scratch_registry%relinquish_field(indices)
    nullify(work_u); nullify(work_v); nullify(work_w); nullify(work_s)
  end subroutine rom_lift_stream_homogeneous

end module rom_lift
