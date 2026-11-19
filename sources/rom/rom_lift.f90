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
  use field_math, only: field_add2s2, field_copy
  use json_file_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use operators, only: grad
  use scratch_registry, only: neko_scratch_registry
  use num_types, only: rp
  use vector, only: vector_t
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: DEVICE_TO_HOST
  use utils, only: neko_error
  implicit none
  private

  public :: rom_lift_t

  !> Spatial lift and its SEM derivatives.
  type, public :: rom_lift_t
     logical :: enabled = .false.
     type(field_t) :: reference(3)
     type(field_t) :: velocity(3)
     type(field_t) :: gradient(3,3)
   contains
     procedure, pass(this), public :: init => rom_lift_init
     procedure, pass(this), public :: free => rom_lift_free
     procedure, pass(this), public :: sync_host => rom_lift_sync_host
     procedure, pass(this), public :: copy_velocity => rom_lift_copy_velocity
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
    type(json_file) :: lift_parameters
    character(len=:), allocatable :: source
    integer :: component, direction

    call this%free()
    if (.not. neko_case%params%valid_path('rom')) return

    source = 'initial_velocity'
    if ('rom.lift' .in. neko_case%params) then
       call json_get(neko_case%params, 'rom.lift', lift_parameters)
       call json_get_or_default(lift_parameters, 'source', source, &
            'initial_velocity')
    end if
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
       call grad(this%gradient(component,1)%x, this%gradient(component,2)%x, &
            this%gradient(component,3)%x, this%velocity(component)%x, &
            neko_case%fluid%c_Xh)
    end do
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
    call this%sync_host()
    n = this%velocity(1)%dof%size()
    do component = 1, 3
       call vectors(component)%init(n)
       vectors(component)%x = reshape(this%velocity(component)%x, [n])
    end do
  end subroutine rom_lift_copy_velocity

  !> Stream a physical velocity after subtracting the spatial lift.
  !! @param[inout] this Lift storage.
  !! @param[inout] dstream Neko/Python data streamer.
  !! @param[in] u Forward x-velocity field.
  !! @param[in] v Forward y-velocity field.
  !! @param[in] w Forward z-velocity field.
  subroutine rom_lift_stream_homogeneous(this, dstream, u, v, w)
    class(rom_lift_t), intent(inout) :: this
    type(data_streamer_t), intent(inout) :: dstream
    type(field_t), intent(inout) :: u, v, w
    type(field_t), pointer :: work_u, work_v, work_w
    integer :: indices(3), n

    if (.not. this%enabled) then
       call dstream%stream(u%x)
       call dstream%stream(v%x)
       call dstream%stream(w%x)
       return
    end if
    n = u%dof%size()
    call neko_scratch_registry%request_field(work_u, indices(1), u%dof, .true.)
    call neko_scratch_registry%request_field(work_v, indices(2), u%dof, .true.)
    call neko_scratch_registry%request_field(work_w, indices(3), u%dof, .true.)
    call field_copy(work_u, u)
    call field_copy(work_v, v)
    call field_copy(work_w, w)
    call field_add2s2(work_u, this%velocity(1), -1.0_rp, n)
    call field_add2s2(work_v, this%velocity(2), -1.0_rp, n)
    call field_add2s2(work_w, this%velocity(3), -1.0_rp, n)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call work_u%copy_from(DEVICE_TO_HOST, sync = .true.)
       call work_v%copy_from(DEVICE_TO_HOST, sync = .true.)
       call work_w%copy_from(DEVICE_TO_HOST, sync = .true.)
    end if
    call dstream%stream(work_u%x)
    call dstream%stream(work_v%x)
    call dstream%stream(work_w%x)
    call neko_scratch_registry%relinquish_field(indices)
    nullify(work_u); nullify(work_v); nullify(work_w)
  end subroutine rom_lift_stream_homogeneous

end module rom_lift
