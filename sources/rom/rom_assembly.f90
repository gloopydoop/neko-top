!> SEM-weighted reduced Galerkin assembly.
!!
!! Pointwise data are carried by Neko vector_t containers and reduced
!! operators by matrix_t/vector_t.  The routines use vector_math for all
!! pointwise products and weighted contractions, so the same host/device path
!! is used as in the rest of Neko-TOP.
module rom_assembly
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_loc, c_size_t, c_sizeof, c_int
  use num_types, only: rp, xp
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_alloc, device_free, device_memcpy, HOST_TO_DEVICE
  use device_math, only: device_glsc3_many
  use comm, only: NEKO_COMM, MPI_EXTRA_PRECISION
  use mpi_f08, only: MPI_Allreduce, MPI_IN_PLACE, MPI_SUM, MPI_Wtime
  use vector, only: vector_t
  use matrix, only: matrix_t
  use vector_math, only: vector_col3, vector_add2, vector_add2s2, vector_cmult, vector_cmult2, &
       vector_glsc2, vector_glsc3
  use scratch_registry, only: neko_scratch_registry
  implicit none
  private

  integer(kind=8), save :: rom_collective_calls = 0_8
  integer(kind=8), save :: rom_collective_elements = 0_8
  integer(kind=8), save :: rom_device_collective_calls = 0_8
  integer(kind=8), save :: rom_device_collective_elements = 0_8
  real(kind=rp), save :: rom_host_collective_seconds = 0.0_rp
  real(kind=rp), save :: rom_device_collective_seconds = 0.0_rp

  public :: rom_assemble_mass
  public :: rom_assemble_mass_coupled
  public :: rom_assemble_scalar_mass
  public :: rom_assemble_cross_mass
  public :: rom_assemble_brinkman
  public :: rom_assemble_brinkman_petrov
  public :: rom_assemble_lift_convection
  public :: rom_assemble_lift_convection_petrov
  public :: rom_assemble_diffusion
  public :: rom_assemble_diffusion_petrov
  public :: rom_assemble_convection
  public :: rom_assemble_convection_petrov
  public :: rom_add_linearized_convection
  public :: rom_add_linearized_scalar_advection
  public :: rom_assemble_linearized_convection_tensor
  public :: rom_assemble_scalar_transport
  public :: rom_assemble_scalar_transport_petrov
  public :: rom_assemble_scalar_diffusion_petrov
  public :: rom_glsc2_local
  public :: rom_glsc3_local
  public :: rom_reduce_buffer
  public :: rom_reset_collective_counters
  public :: rom_get_collective_counters

contains

  !> Evaluate a batch of weighted contractions against one vector mode family.
  !! The mode vectors remain on the accelerator. Only the short array of their
  !! device addresses is copied to device memory for Neko's native many-kernel.
  !! That kernel performs the backend-global reduction for the complete batch;
  !! do not apply an additional MPI reduction to its results.
  subroutine rom_glsc3_many_modes(a, modes, component, c, values)
    type(vector_t), intent(in) :: a, modes(:, :), c
    integer, intent(in) :: component
    real(kind=rp), intent(out) :: values(:)
    type(c_ptr), allocatable, target :: mode_ptrs(:)
    integer :: mode

    if (component < 1 .or. component > size(modes,1) .or. &
         size(values) /= size(modes,2)) error stop 'ROM batch contraction shape mismatch'
    if (NEKO_BCKND_DEVICE .eq. 1 .and. size(values) > 1 .and. &
         size(values) <= 1024) then
       allocate(mode_ptrs(size(values)))
       do mode = 1, size(values)
          mode_ptrs(mode) = modes(component,mode)%x_d
       end do
       call rom_glsc3_many_ptrs(a, mode_ptrs, c, values)
       deallocate(mode_ptrs)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       do mode = 1, size(values)
          values(mode) = rom_glsc3_local(a, c, modes(component,mode))
       end do
    else
       do mode = 1, size(values)
          values(mode) = rom_glsc3_local(a, c, modes(component,mode))
       end do
    end if
  end subroutine rom_glsc3_many_modes

  !> Batch contractions against one component/direction of a gradient basis.
  subroutine rom_glsc3_many_gradient_modes(a, gradients, component, direction, c, values)
    type(vector_t), intent(in) :: a, gradients(:, :, :), c
    integer, intent(in) :: component, direction
    real(kind=rp), intent(out) :: values(:)
    type(c_ptr), allocatable, target :: mode_ptrs(:)
    integer :: mode

    if (component < 1 .or. component > size(gradients,1) .or. &
         direction < 1 .or. direction > size(gradients,3) .or. &
         size(values) /= size(gradients,2)) &
         error stop 'ROM batch gradient contraction shape mismatch'
    if (NEKO_BCKND_DEVICE .eq. 1 .and. size(values) > 1 .and. &
         size(values) <= 1024) then
       allocate(mode_ptrs(size(values)))
       do mode = 1, size(values)
          mode_ptrs(mode) = gradients(component,mode,direction)%x_d
       end do
       call rom_glsc3_many_ptrs(a, mode_ptrs, c, values)
       deallocate(mode_ptrs)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       do mode = 1, size(values)
          values(mode) = rom_glsc3_local(a, c, gradients(component,mode,direction))
       end do
    else
       do mode = 1, size(values)
          values(mode) = rom_glsc3_local(a, c, gradients(component,mode,direction))
       end do
    end if
  end subroutine rom_glsc3_many_gradient_modes

  !> Scalar-mode counterpart to rom_glsc3_many_modes.
  subroutine rom_glsc3_many_scalar_modes(a, modes, c, values)
    type(vector_t), intent(in) :: a, modes(:), c
    real(kind=rp), intent(out) :: values(:)
    type(c_ptr), allocatable, target :: mode_ptrs(:)
    integer :: mode

    if (size(values) /= size(modes)) error stop 'ROM batch scalar contraction shape mismatch'
    if (NEKO_BCKND_DEVICE .eq. 1 .and. size(values) > 1 .and. &
         size(values) <= 1024) then
       allocate(mode_ptrs(size(values)))
       do mode = 1, size(values)
          mode_ptrs(mode) = modes(mode)%x_d
       end do
       call rom_glsc3_many_ptrs(a, mode_ptrs, c, values)
       deallocate(mode_ptrs)
    else if (NEKO_BCKND_DEVICE .eq. 1) then
       do mode = 1, size(values)
          values(mode) = rom_glsc3_local(a, c, modes(mode))
       end do
    else
       do mode = 1, size(values)
          values(mode) = rom_glsc3_local(a, c, modes(mode))
       end do
    end if
  end subroutine rom_glsc3_many_scalar_modes

  subroutine rom_glsc3_many_ptrs(a, mode_ptrs, c, values)
    type(vector_t), intent(in) :: a, c
    type(c_ptr), intent(in), target :: mode_ptrs(:)
    real(kind=rp), intent(out) :: values(:)
    type(c_ptr) :: mode_ptrs_d, mode_ptrs_h
    integer(c_size_t) :: pointer_bytes
    integer :: n, n_modes
    integer(c_int) :: n_c, n_modes_c
    real(kind=rp) :: timer_start

    n_modes = size(mode_ptrs)
    n = a%size()
    if (n_modes == 0) return
    if (n_modes /= size(values) .or. c%size() /= n) &
         error stop 'ROM batch contraction vector size mismatch'
    if (NEKO_BCKND_DEVICE .eq. 1) then
       timer_start = real(MPI_Wtime(),rp)
       pointer_bytes = int(n_modes,c_size_t) * c_sizeof(mode_ptrs(1))
       mode_ptrs_d = c_null_ptr
       call device_alloc(mode_ptrs_d, pointer_bytes)
       mode_ptrs_h = c_loc(mode_ptrs(1))
       call device_memcpy(mode_ptrs_h, mode_ptrs_d, pointer_bytes, &
            HOST_TO_DEVICE, sync=.true.)
       n_c = int(n,c_int); n_modes_c = int(n_modes,c_int)
       call device_glsc3_many(values, a%x_d, mode_ptrs_d, c%x_d, n_modes_c, n_c)
       call device_free(mode_ptrs_d)
       rom_device_collective_seconds = rom_device_collective_seconds + &
            real(MPI_Wtime(),rp) - timer_start
       rom_device_collective_calls = rom_device_collective_calls + 1_8
       rom_device_collective_elements = rom_device_collective_elements + int(n_modes,8)
    end if
  end subroutine rom_glsc3_many_ptrs

  !> A local contraction for the host backend. Device contractions retain
  !! vector_glsc3, whose backend implementation already performs the global
  !! reduction. Host results are reduced after each complete reduced operator
  !! has been assembled, avoiding one collective per reduced coefficient.
  function rom_glsc3_local(a, b, c) result(value)
    type(vector_t), intent(in) :: a, b, c
    real(kind=rp) :: value
    real(kind=xp) :: local_value
    real(kind=rp) :: timer_start
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then
       timer_start = real(MPI_Wtime(),rp)
       value = vector_glsc3(a, b, c)
       rom_device_collective_seconds = rom_device_collective_seconds + &
            real(MPI_Wtime(),rp) - timer_start
       rom_device_collective_calls = rom_device_collective_calls + 1_8
       rom_device_collective_elements = rom_device_collective_elements + 1_8
       return
    end if
    local_value = 0.0_xp
    do i = 1, a%size()
       local_value = local_value + real(a%x(i),xp) * real(b%x(i),xp) * real(c%x(i),xp)
    end do
    value = real(local_value,rp)
  end function rom_glsc3_local

  !> Local weighted contraction. The device fallback is globally reduced by
  !! Neko's vector backend; host contractions are genuinely rank-local.
  function rom_glsc2_local(a, b) result(value)
    type(vector_t), intent(in) :: a, b
    real(kind=rp) :: value
    real(kind=xp) :: local_value
    real(kind=rp) :: timer_start
    integer :: i

    if (NEKO_BCKND_DEVICE .eq. 1) then
       timer_start = real(MPI_Wtime(),rp)
       value = vector_glsc2(a, b)
       rom_device_collective_seconds = rom_device_collective_seconds + &
            real(MPI_Wtime(),rp) - timer_start
       rom_device_collective_calls = rom_device_collective_calls + 1_8
       rom_device_collective_elements = rom_device_collective_elements + 1_8
       return
    end if
    local_value = 0.0_xp
    do i = 1, a%size()
       local_value = local_value + real(a%x(i),xp) * real(b%x(i),xp)
    end do
    value = real(local_value,rp)
  end function rom_glsc2_local

  !> Reduce a complete local real buffer in one collective.
  !! Device callers use globally reduced vector contractions and must not
  !! call an additional MPI reduction on their host mirrors.
  subroutine rom_reduce_buffer(values)
    real(kind=rp), intent(inout) :: values(:)
    real(kind=xp), allocatable :: buffer(:)
    integer :: i, ierr
    real(kind=rp) :: timer_start

    if (NEKO_BCKND_DEVICE .eq. 1 .or. size(values) == 0) return
    allocate(buffer(size(values)))
    do i = 1, size(values)
       buffer(i) = real(values(i),xp)
    end do
    timer_start = real(MPI_Wtime(),rp)
    call MPI_Allreduce(MPI_IN_PLACE, buffer, size(buffer), MPI_EXTRA_PRECISION, &
         MPI_SUM, NEKO_COMM, ierr)
    rom_host_collective_seconds = rom_host_collective_seconds + &
         real(MPI_Wtime(),rp) - timer_start
    rom_collective_calls = rom_collective_calls + 1_8
    rom_collective_elements = rom_collective_elements + int(size(buffer),8)
    do i = 1, size(values)
       values(i) = real(buffer(i),rp)
    end do
    deallocate(buffer)
  end subroutine rom_reduce_buffer

  !> Reset/read assembly communication counters. Host counts record packed
  !! MPI_Allreduce calls and xp buffer sizes. Device counts separately record
  !! backend-global scalar vector_glsc contractions and xp result sizes; their
  !! elapsed time includes the device contraction/synchronization path.
  subroutine rom_reset_collective_counters()
    rom_collective_calls = 0_8
    rom_collective_elements = 0_8
    rom_device_collective_calls = 0_8
    rom_device_collective_elements = 0_8
    rom_host_collective_seconds = 0.0_rp
    rom_device_collective_seconds = 0.0_rp
  end subroutine rom_reset_collective_counters

  subroutine rom_get_collective_counters(calls, elements, bytes, device_calls, &
       device_elements, device_bytes, host_seconds, device_seconds)
    integer(kind=8), intent(out) :: calls, elements, bytes
    integer(kind=8), intent(out), optional :: device_calls, device_elements, device_bytes
    real(kind=rp), intent(out), optional :: host_seconds, device_seconds
    calls = rom_collective_calls
    elements = rom_collective_elements
    bytes = rom_collective_elements * int(storage_size(0.0_xp)/8,8)
    if (present(device_calls)) device_calls = rom_device_collective_calls
    if (present(device_elements)) device_elements = rom_device_collective_elements
    if (present(device_bytes)) device_bytes = &
         rom_device_collective_elements * int(storage_size(0.0_xp)/8,8)
    if (present(host_seconds)) host_seconds = rom_host_collective_seconds
    if (present(device_seconds)) device_seconds = rom_device_collective_seconds
  end subroutine rom_get_collective_counters

  !> Sum one locally assembled matrix and vector over the Neko communicator.
  !! Call only on the host backend: device vector_glsc3 values are already
  !! globally reduced and the device host mirrors may be stale.
  subroutine rom_reduce_matrix_vector(matrix, vector)
    type(matrix_t), intent(inout) :: matrix
    type(vector_t), intent(inout) :: vector
    call rom_reduce_matrix_array(matrix, vector%x)
  end subroutine rom_reduce_matrix_vector

  !> Sum one locally assembled matrix and real vector buffer.
  subroutine rom_reduce_matrix_array(matrix, values)
    type(matrix_t), intent(inout) :: matrix
    real(kind=rp), intent(inout) :: values(:)
    real(kind=rp), allocatable :: buffer(:)
    integer :: i, n_matrix, nrows

    if (NEKO_BCKND_DEVICE .eq. 1) return
    n_matrix = size(matrix%x)
    nrows = matrix%get_nrows()
    allocate(buffer(n_matrix + size(values)))
    do i = 1, n_matrix
       buffer(i) = matrix%x(mod(i-1,nrows)+1,(i-1)/nrows+1)
    end do
    do i = 1, size(values)
       buffer(n_matrix+i) = values(i)
    end do
    call rom_reduce_buffer(buffer)
    do i = 1, n_matrix
       matrix%x(mod(i-1,nrows)+1,(i-1)/nrows+1) = buffer(i)
    end do
    do i = 1, size(values)
       values(i) = buffer(n_matrix+i)
    end do
    deallocate(buffer)
  end subroutine rom_reduce_matrix_array

  !> Sum a locally assembled matrix over the Neko communicator.
  subroutine rom_reduce_matrix(matrix)
    type(matrix_t), intent(inout) :: matrix
    real(kind=rp), allocatable :: buffer(:)
    integer :: i, nrows
    if (NEKO_BCKND_DEVICE .eq. 1) return
    allocate(buffer(size(matrix%x)))
    nrows = matrix%get_nrows()
    do i = 1, size(matrix%x)
       buffer(i) = matrix%x(mod(i-1,nrows)+1,(i-1)/nrows+1)
    end do
    call rom_reduce_buffer(buffer)
    do i = 1, size(matrix%x)
       matrix%x(mod(i-1,nrows)+1,(i-1)/nrows+1) = buffer(i)
    end do
    deallocate(buffer)
  end subroutine rom_reduce_matrix

  !> Sum a locally assembled tensor represented by reduced matrices.
  subroutine rom_reduce_matrix_tensor(tensor)
    type(matrix_t), intent(inout) :: tensor(:)
    real(kind=rp), allocatable :: buffer(:)
    integer :: i, j, k, n_matrix, nrows

    if (NEKO_BCKND_DEVICE .eq. 1) return
    n_matrix = 0
    do i = 1, size(tensor)
       n_matrix = n_matrix + size(tensor(i)%x)
    end do
    allocate(buffer(n_matrix))
    n_matrix = 0
    do i = 1, size(tensor)
       nrows = tensor(i)%get_nrows()
       do j = 1, size(tensor(i)%x)
          n_matrix = n_matrix + 1
          buffer(n_matrix) = tensor(i)%x(mod(j-1,nrows)+1,(j-1)/nrows+1)
       end do
    end do
    call rom_reduce_buffer(buffer)
    n_matrix = 0
    do i = 1, size(tensor)
       nrows = tensor(i)%get_nrows()
       do j = 1, size(tensor(i)%x)
          n_matrix = n_matrix + 1
          tensor(i)%x(mod(j-1,nrows)+1,(j-1)/nrows+1) = buffer(n_matrix)
       end do
    end do
    deallocate(buffer)
  end subroutine rom_reduce_matrix_tensor

  !> Assemble the SEM mass matrix in the reduced velocity basis.
  subroutine rom_assemble_mass(weights, modes, mass)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: modes(:, :)
    type(matrix_t), intent(inout) :: mass
    integer :: n, n_fields, n_modes, i, j, field

    n = weights%size(); n_fields = size(modes, 1); n_modes = size(modes, 2)
    if (n < 1 .or. n_fields < 1 .or. any([(modes(field,1)%size(), field=1,n_fields)] /= n)) &
         error stop 'ROM mass assembly shape mismatch'
    call rom_assemble_cross_mass(weights, modes, modes, mass)
  end subroutine rom_assemble_mass

  !> Assemble the Petrov--Galerkin cross mass matrix.
  !! The first mode argument is the trial (right) basis and the second is the
  !! test (left) basis.  Optional scalar modes use the same coupled-state
  !! metric as the existing Galerkin assembly.
  subroutine rom_assemble_cross_mass(weights, trial_modes, test_modes, mass, &
       scalar_trial_modes, scalar_test_modes, scalar_weight)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: trial_modes(:, :), test_modes(:, :)
    type(matrix_t), intent(inout) :: mass
    type(vector_t), intent(in), optional :: scalar_trial_modes(:), scalar_test_modes(:)
    real(kind=rp), intent(in), optional :: scalar_weight
    integer :: n, r_trial, r_test, i, j, field
    real(kind=rp) :: weight
    real(kind=rp), allocatable :: values(:)

    n = weights%size(); r_trial = size(trial_modes, 2); r_test = size(test_modes, 2)
    if (r_trial < 1 .or. r_trial /= r_test .or. size(trial_modes, 1) /= size(test_modes, 1)) &
         error stop 'ROM cross mass basis shape mismatch'
    do field = 1, size(trial_modes, 1)
       if (trial_modes(field,1)%size() /= n .or. test_modes(field,1)%size() /= n) &
            error stop 'ROM cross mass vector size mismatch'
    end do
    call mass%init(r_test, r_trial)
    mass%x = 0.0_rp
    allocate(values(r_test))
    do j = 1, r_trial
       do field = 1, size(trial_modes, 1)
          call rom_glsc3_many_modes(trial_modes(field,j), test_modes, field, &
               weights, values)
          mass%x(:,j) = mass%x(:,j) + values
       end do
    end do
    if (present(scalar_trial_modes) .or. present(scalar_test_modes)) then
       if (.not. present(scalar_trial_modes) .or. .not. present(scalar_test_modes) .or. &
            size(scalar_trial_modes) /= r_trial .or. size(scalar_test_modes) /= r_test) &
            error stop 'ROM cross scalar mass shape mismatch'
       weight = 1.0_rp
       if (present(scalar_weight)) weight = scalar_weight
       do j = 1, r_trial
          call rom_glsc3_many_scalar_modes(scalar_trial_modes(j), &
               scalar_test_modes, weights, values)
          mass%x(:,j) = mass%x(:,j) + weight * values
       end do
    end if
    call rom_reduce_matrix(mass)
    deallocate(values)
  end subroutine rom_assemble_cross_mass

  !> Assemble the coupled velocity/scalar SEM mass matrix.
  !! The scalar multiplier is the declared state metric.  It is applied to
  !! both the scalar mass and scalar residual contractions, so the reduced
  !! equations remain the Galerkin projection of the block-diagonal physical
  !! mass operator.
  subroutine rom_assemble_mass_coupled(weights, velocity_modes, scalar_modes, &
       scalar_weight, mass)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: velocity_modes(:, :)
    type(vector_t), intent(in), optional :: scalar_modes(:)
    real(kind=rp), intent(in), optional :: scalar_weight
    type(matrix_t), intent(inout) :: mass
    integer :: r

    r = size(velocity_modes, 2)
    call rom_assemble_cross_mass(weights, velocity_modes, velocity_modes, mass, &
         scalar_modes, scalar_modes, scalar_weight)
    if (.not. present(scalar_modes)) return
    if (size(scalar_modes) /= r) error stop 'ROM scalar mass shape mismatch'
  end subroutine rom_assemble_mass_coupled

  !> Assemble a coupled-metric scalar mass block independently of velocity.
  !! Keeping this block separate lets an adaptive velocity basis refresh its
  !! mass contribution without repeating the invariant scalar contractions.
  subroutine rom_assemble_scalar_mass(weights, scalar_weight, trial_scalar_modes, &
       test_scalar_modes, mass)
    type(vector_t), intent(in) :: weights
    real(kind=rp), intent(in) :: scalar_weight
    type(vector_t), intent(in) :: trial_scalar_modes(:), test_scalar_modes(:)
    type(matrix_t), intent(inout) :: mass
    real(kind=rp), allocatable :: local(:)
    integer :: n, nt, nr, i, j
    real(kind=rp), allocatable :: values(:)

    n = weights%size(); nt = size(test_scalar_modes); nr = size(trial_scalar_modes)
    if (n < 1 .or. nt < 1 .or. nr < 1 .or. scalar_weight <= 0.0_rp) &
         error stop 'ROM scalar mass shape/metric mismatch'
    if (any([(test_scalar_modes(i)%size(), i=1,nt)] /= n) .or. &
         any([(trial_scalar_modes(j)%size(), j=1,nr)] /= n)) &
         error stop 'ROM scalar mass vector size mismatch'

    allocate(local(nt*nr), values(nt)); local = 0.0_rp
    do j = 1, nr
       call rom_glsc3_many_scalar_modes(trial_scalar_modes(j), &
            test_scalar_modes, weights, values)
       local((j-1)*nt+1:j*nt) = scalar_weight * values
    end do
    call rom_reduce_buffer(local)
    call mass%init(nt,nr)
    do j = 1, nr
       do i = 1, nt
          mass%x(i,j) = local((j-1)*nt+i)
       end do
    end do
    deallocate(local, values)
  end subroutine rom_assemble_scalar_mass

  !> Assemble the reduced Brinkman operator and its lift contribution.
  subroutine rom_assemble_brinkman(weights, chi, velocity_modes, lift_velocity, &
       brinkman, lift_vector)
    type(vector_t), intent(in) :: weights, chi
    type(vector_t), intent(in) :: velocity_modes(:, :), lift_velocity(:)
    type(matrix_t), intent(inout) :: brinkman
    type(vector_t), intent(inout) :: lift_vector
    call rom_assemble_brinkman_petrov(weights, chi, velocity_modes, velocity_modes, &
         lift_velocity, brinkman, lift_vector)
  end subroutine rom_assemble_brinkman

  !> Assemble the Petrov--Galerkin Brinkman operator and lift contribution.
  !! Test modes are used only on the left; trial modes are used on the right.
  subroutine rom_assemble_brinkman_petrov(weights, chi, trial_velocity_modes, &
       test_velocity_modes, lift_velocity, brinkman, lift_vector)
    type(vector_t), intent(in) :: weights, chi
    type(vector_t), intent(in) :: trial_velocity_modes(:, :), test_velocity_modes(:, :)
    type(vector_t), intent(in) :: lift_velocity(:)
    type(matrix_t), intent(inout) :: brinkman
    type(vector_t), intent(inout) :: lift_vector
    type(vector_t), pointer :: effective
    integer :: effective_index
    integer :: n, n_trial, n_test, i, j, component
    real(kind=rp), allocatable :: values(:)

    n = weights%size(); n_trial = size(trial_velocity_modes, 2); n_test = size(test_velocity_modes, 2)
    if (chi%size() /= n .or. size(trial_velocity_modes, 1) /= 3 .or. &
         size(test_velocity_modes, 1) /= 3 .or. size(lift_velocity) /= 3 .or. &
         any([(trial_velocity_modes(component,1)%size(), component=1,3)] /= n) .or. &
         any([(test_velocity_modes(component,1)%size(), component=1,3)] /= n)) &
         error stop 'ROM Brinkman assembly shape mismatch'
    call neko_scratch_registry%request_vector(effective, effective_index, n, .true.)
    call vector_col3(effective, weights, chi)
    call brinkman%init(n_test, n_trial); call lift_vector%init(n_test)
    ! The contractions return host scalars and these reduced outputs are
    ! deliberately host-owned.
    brinkman%x = 0.0_rp
    lift_vector%x = 0.0_rp
    allocate(values(n_test))
    do j = 1, n_trial
       do component = 1, 3
          call rom_glsc3_many_modes(trial_velocity_modes(component,j), &
               test_velocity_modes, component, effective, values)
          brinkman%x(:,j) = brinkman%x(:,j) + values
       end do
    end do
    do component = 1, 3
       call rom_glsc3_many_modes(lift_velocity(component), test_velocity_modes, &
            component, effective, values)
       lift_vector%x = lift_vector%x + values
    end do
    call rom_reduce_matrix_vector(brinkman, lift_vector)
    deallocate(values)
    call neko_scratch_registry%relinquish_vector(effective_index); nullify(effective)
  end subroutine rom_assemble_brinkman_petrov

  !> Assemble the convection generated by a spatial physical lift.
  !!
  !! For `u = lift + sum(a_j phi_j)`, this routine returns the negative
  !! weak residual contributions from `lift.grad(phi_j)`,
  !! `phi_j.grad(lift)` and `lift.grad(lift)`.  The latter is returned in the
  !! constant vector, while the two cross terms are returned in the matrix.
  !! @param[in] weights SEM quadrature weights.
  !! @param[in] velocity_modes Homogeneous velocity modes.
  !! @param[in] gradients Mode gradients, indexed by component, mode, direction.
  !! @param[in] lift_velocity Spatial lift components.
  !! @param[in] lift_gradients Spatial lift gradients, indexed by component and direction.
  !! @param[out] linear Reduced cross-convection matrix.
  !! @param[out] constant Reduced lift-convection vector.
  subroutine rom_assemble_lift_convection(weights, velocity_modes, gradients, &
       lift_velocity, lift_gradients, linear, constant)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: velocity_modes(:, :), gradients(:, :, :)
    type(vector_t), intent(in) :: lift_velocity(:), lift_gradients(:, :)
    type(matrix_t), intent(inout) :: linear
    type(vector_t), intent(inout) :: constant
    call rom_assemble_lift_convection_petrov(weights, velocity_modes, gradients, &
         velocity_modes, lift_velocity, lift_gradients, linear, constant)
  end subroutine rom_assemble_lift_convection

  !> Assemble lift convection with independent trial and test spaces.
  subroutine rom_assemble_lift_convection_petrov(weights, trial_velocity_modes, &
       trial_gradients, test_velocity_modes, lift_velocity, lift_gradients, &
       linear, constant)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: trial_velocity_modes(:, :), trial_gradients(:, :, :)
    type(vector_t), intent(in) :: test_velocity_modes(:, :)
    type(vector_t), intent(in) :: lift_velocity(:), lift_gradients(:, :)
    type(matrix_t), intent(inout) :: linear
    type(vector_t), intent(inout) :: constant
    type(vector_t), pointer :: pointwise, product
    integer :: pointwise_index, product_index
    integer :: n, r, j, component, direction
    real(kind=rp), allocatable :: values(:)

    n = weights%size(); r = size(trial_velocity_modes, 2)
    if (size(trial_velocity_modes,1) /= 3 .or. size(test_velocity_modes,1) /= 3 .or. &
         size(test_velocity_modes,2) /= r .or. size(trial_gradients,1) /= 3 .or. &
         size(trial_gradients,2) /= r .or. size(trial_gradients,3) /= 3 .or. &
         size(lift_velocity) /= 3 .or. size(lift_gradients,1) /= 3 .or. &
         size(lift_gradients,2) /= 3) error stop 'ROM lift convection shape mismatch'
    call linear%init(r, r)
    call constant%init(r)
    linear%x = 0.0_rp
    constant%x = 0.0_rp
    call neko_scratch_registry%request_vector(pointwise, pointwise_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    allocate(values(r))
    do j = 1, r
       do component = 1, 3
          pointwise = 0.0_rp
          do direction = 1, 3
             call vector_col3(product, lift_velocity(direction), &
                  trial_gradients(component,j,direction))
             call vector_add2(pointwise, product)
             call vector_col3(product, trial_velocity_modes(direction,j), &
                  lift_gradients(component,direction))
             call vector_add2(pointwise, product)
          end do
          call rom_glsc3_many_modes(weights, test_velocity_modes, component, &
               pointwise, values)
          linear%x(:,j) = linear%x(:,j) - values
       end do
    end do
    do component = 1, 3
       pointwise = 0.0_rp
       do direction = 1, 3
          call vector_col3(product, lift_velocity(direction), &
               lift_gradients(component,direction))
          call vector_add2(pointwise, product)
       end do
       call rom_glsc3_many_modes(weights, test_velocity_modes, component, &
            pointwise, values)
       constant%x = constant%x - values
    end do
    call rom_reduce_matrix_vector(linear, constant)
    deallocate(values)
    call neko_scratch_registry%relinquish_vector([pointwise_index, product_index])
    nullify(pointwise); nullify(product)
  end subroutine rom_assemble_lift_convection_petrov

  !> Assemble the implicit viscous residual matrix and lift diffusion vector.
  !!
  !! The matrix is `-nu <grad(test),grad(state)>`.  When lift gradients and a
  !! constant vector are supplied, the routine additionally accumulates
  !! `-nu <grad(test),grad(lift)>` in that vector.
  subroutine rom_assemble_diffusion(weights, viscosity_over_density, gradients, &
       implicit_linear, lift_gradients, lift_constant)
    type(vector_t), intent(in) :: weights
    real(kind=rp), intent(in) :: viscosity_over_density
    type(vector_t), intent(in) :: gradients(:, :, :)
    type(matrix_t), intent(inout) :: implicit_linear
    type(vector_t), intent(in), optional :: lift_gradients(:, :)
    type(vector_t), intent(inout), optional :: lift_constant
    call rom_assemble_diffusion_petrov(weights, viscosity_over_density, gradients, &
         gradients, implicit_linear, lift_gradients, lift_constant)
  end subroutine rom_assemble_diffusion

  !> Assemble diffusion with independent test and trial gradients.
  subroutine rom_assemble_diffusion_petrov(weights, viscosity_over_density, &
       test_gradients, trial_gradients, implicit_linear, lift_gradients, lift_constant)
    type(vector_t), intent(in) :: weights
    real(kind=rp), intent(in) :: viscosity_over_density
    type(vector_t), intent(in) :: test_gradients(:, :, :), trial_gradients(:, :, :)
    type(matrix_t), intent(inout) :: implicit_linear
    type(vector_t), intent(in), optional :: lift_gradients(:, :)
    type(vector_t), intent(inout), optional :: lift_constant
    real(kind=rp), allocatable :: lift_delta(:)
    integer :: n, r, i, j, component, direction
    real(kind=rp), allocatable :: values(:)

    n = weights%size(); r = size(trial_gradients,2)
    if (size(test_gradients,1) /= 3 .or. size(test_gradients,2) /= r .or. &
         size(test_gradients,3) /= 3 .or. size(trial_gradients,1) /= 3 .or. &
         size(trial_gradients,3) /= 3) &
         error stop 'ROM diffusion assembly shape mismatch'
    if (present(lift_gradients) .neqv. present(lift_constant)) &
         error stop 'ROM diffusion lift arguments must be supplied together'
    if (present(lift_gradients)) then
       if (size(lift_gradients,1) /= 3 .or. size(lift_gradients,2) /= 3) &
            error stop 'ROM diffusion lift shape mismatch'
    end if
    call implicit_linear%init(r, r)
    implicit_linear%x = 0.0_rp
    if (present(lift_gradients)) then
       allocate(lift_delta(r))
       lift_delta = 0.0_rp
    end if
    allocate(values(r))
    if (present(lift_gradients)) then
       do component = 1, 3
          do direction = 1, 3
             call rom_glsc3_many_gradient_modes(lift_gradients(component,direction), &
                  test_gradients, component, direction, weights, values)
             lift_delta = lift_delta - viscosity_over_density * values
          end do
       end do
    end if
    do j = 1, r
       do component = 1, 3
          do direction = 1, 3
             call rom_glsc3_many_gradient_modes(trial_gradients(component,j,direction), &
                  test_gradients, component, direction, weights, values)
             implicit_linear%x(:,j) = implicit_linear%x(:,j) - &
                  viscosity_over_density * values
          end do
       end do
    end do
    if (present(lift_gradients)) then
       call rom_reduce_matrix_array(implicit_linear, lift_delta)
       lift_constant%x = lift_constant%x + lift_delta
       deallocate(lift_delta)
    else
       call rom_reduce_matrix(implicit_linear)
    end if
    deallocate(values)
  end subroutine rom_assemble_diffusion_petrov

  !> Assemble the quadratic convective residual tensor.
  subroutine rom_assemble_convection(weights, modes, gradients, quadratic)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: modes(:, :), gradients(:, :, :)
    type(matrix_t), intent(inout) :: quadratic(:)
    call rom_assemble_convection_petrov(weights, modes, modes, gradients, quadratic)
  end subroutine rom_assemble_convection

  !> Assemble the Petrov--Galerkin quadratic velocity convection tensor.
  !! The advecting and transported fields are trial modes; only the test
  !! function is taken from the independent test basis.
  subroutine rom_assemble_convection_petrov(weights, test_modes, trial_modes, &
       trial_gradients, quadratic)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: test_modes(:, :), trial_modes(:, :)
    type(vector_t), intent(in) :: trial_gradients(:, :, :)
    type(matrix_t), intent(inout) :: quadratic(:)
    type(vector_t), pointer :: pointwise, product
    integer :: pointwise_index, product_index
    integer :: n, r, test_mode, advect_mode, state_mode, field, component
    real(kind=rp), allocatable :: values(:)

    n = weights%size(); r = size(trial_modes,2)
    if (size(test_modes,1) < 3 .or. size(test_modes,2) /= r .or. &
         size(trial_modes,1) < 3 .or. size(trial_gradients,1) < 3 .or. &
         size(trial_gradients,2) /= r .or. size(trial_gradients,3) /= 3 .or. &
         size(quadratic) /= r) error stop 'ROM convection assembly shape mismatch'
    do field = 1, r
       call quadratic(field)%init(r, r)
       quadratic(field)%x = 0.0_rp
    end do
    call neko_scratch_registry%request_vector(pointwise, pointwise_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    allocate(values(r))
    ! The velocity-only ROM has three transported fields.  Additional scalar
    ! fields are intentionally unsupported by rom_problem's configuration gate.
    do advect_mode = 1, r
       do state_mode = 1, r
          do field = 1, 3
             pointwise = 0.0_rp
             do component = 1, 3
                call vector_col3(product, trial_modes(component,advect_mode), &
                     trial_gradients(field,state_mode,component))
                call vector_add2(pointwise, product)
             end do
             call rom_glsc3_many_modes(weights, test_modes, field, pointwise, values)
             do test_mode = 1, r
                quadratic(test_mode)%x(advect_mode,state_mode) = &
                     quadratic(test_mode)%x(advect_mode,state_mode) + values(test_mode)
             end do
          end do
       end do
    end do
    call rom_reduce_matrix_tensor(quadratic)
    deallocate(values)
    call neko_scratch_registry%relinquish_vector([pointwise_index, product_index])
    nullify(pointwise); nullify(product)
  end subroutine rom_assemble_convection_petrov

  !> Add the velocity-convection Jacobian about a forward reduced state.
  !!
  !! The forward state is represented by `coefficients` in the forward trial
  !! basis.  The perturbation and test spaces may be enriched and therefore
  !! have an independent rank.  The two contributions are
  !! `-(delta_u.grad(u_f) + u_f.grad(delta_u))`.
  subroutine rom_add_linearized_convection(weights, forward_modes, forward_gradients, &
       test_modes, perturb_modes, perturb_gradients, coefficients, linearized)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: forward_modes(:, :), forward_gradients(:, :, :)
    type(vector_t), intent(in) :: test_modes(:, :), perturb_modes(:, :)
    type(vector_t), intent(in) :: perturb_gradients(:, :, :)
    real(kind=rp), intent(in) :: coefficients(:)
    type(matrix_t), intent(inout) :: linearized
    type(vector_t), pointer :: pointwise, product
    real(kind=rp), allocatable :: local(:)
    integer :: pointwise_index, product_index
    integer :: n, rf, re, i, j, k, field, direction

    n = weights%size(); rf = size(forward_modes,2); re = size(perturb_modes,2)
    if (size(forward_modes,1) /= 3 .or. size(test_modes,1) /= 3 .or. &
         size(perturb_modes,1) /= 3 .or. size(coefficients) /= rf .or. &
         size(forward_gradients,1) /= 3 .or. size(forward_gradients,2) /= rf .or. &
         size(forward_gradients,3) /= 3 .or. size(perturb_gradients,1) /= 3 .or. &
         size(perturb_gradients,2) /= re .or. size(perturb_gradients,3) /= 3 .or. &
         size(test_modes,2) /= re .or. linearized%get_nrows() /= re .or. &
         linearized%get_ncols() /= re) error stop 'ROM linearized convection shape mismatch'
    do field = 1, 3
       if (forward_modes(field,1)%size() /= n .or. test_modes(field,1)%size() /= n .or. &
            perturb_modes(field,1)%size() /= n) error stop 'ROM linearized convection vector mismatch'
    end do

    call neko_scratch_registry%request_vector(pointwise, pointwise_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    allocate(local(re*re)); local = 0.0_rp
    do j = 1, re
       do field = 1, 3
          pointwise = 0.0_rp
          do k = 1, rf
             do direction = 1, 3
                call vector_col3(product, perturb_modes(direction,j), &
                     forward_gradients(field,k,direction))
                call vector_add2s2(pointwise, product, coefficients(k))
                call vector_col3(product, forward_modes(direction,k), &
                     perturb_gradients(field,j,direction))
                call vector_add2s2(pointwise, product, coefficients(k))
             end do
          end do
          do i = 1, re
             local((j-1)*re+i) = local((j-1)*re+i) - &
                  rom_glsc3_local(test_modes(field,i), weights, pointwise)
          end do
       end do
    end do
    call rom_reduce_buffer(local)
    do j = 1, re
       do i = 1, re
          linearized%x(i,j) = linearized%x(i,j) + local((j-1)*re+i)
       end do
    end do
    deallocate(local)
    call neko_scratch_registry%relinquish_vector([pointwise_index, product_index])
    nullify(pointwise); nullify(product)
  end subroutine rom_add_linearized_convection

  !> Add the scalar-advection Jacobian about a forward coupled state.
  !! The derivative is `-(delta_u.grad(phi_f) + u_f.grad(delta_phi))`, with
  !! the configured scalar metric applied to the scalar residual row.
  subroutine rom_add_linearized_scalar_advection(weights, scalar_weight, &
       forward_velocity_modes, forward_scalar_modes, forward_scalar_gradients, &
       test_scalar_modes, perturb_velocity_modes, perturb_scalar_modes, &
       perturb_scalar_gradients, coefficients, linearized)
    type(vector_t), intent(in) :: weights
    real(kind=rp), intent(in) :: scalar_weight
    type(vector_t), intent(in) :: forward_velocity_modes(:, :)
    type(vector_t), intent(in) :: forward_scalar_modes(:), forward_scalar_gradients(:, :)
    type(vector_t), intent(in) :: test_scalar_modes(:)
    type(vector_t), intent(in) :: perturb_velocity_modes(:, :)
    type(vector_t), intent(in) :: perturb_scalar_modes(:), perturb_scalar_gradients(:, :)
    real(kind=rp), intent(in) :: coefficients(:)
    type(matrix_t), intent(inout) :: linearized
    type(vector_t), pointer :: pointwise, product
    real(kind=rp), allocatable :: local(:)
    integer :: pointwise_index, product_index
    integer :: n, rf, re, i, j, k, direction

    n = weights%size(); rf = size(forward_velocity_modes,2); re = size(perturb_velocity_modes,2)
    if (size(forward_velocity_modes,1) /= 3 .or. size(perturb_velocity_modes,1) /= 3 .or. &
         size(forward_scalar_modes) /= rf .or. size(perturb_scalar_modes) /= re .or. &
         size(forward_scalar_gradients,1) /= rf .or. size(forward_scalar_gradients,2) /= 3 .or. &
         size(perturb_scalar_gradients,1) /= re .or. size(perturb_scalar_gradients,2) /= 3 .or. &
         size(test_scalar_modes) /= re .or. size(coefficients) /= rf .or. &
         linearized%get_nrows() /= re .or. linearized%get_ncols() /= re) &
         error stop 'ROM linearized scalar advection shape mismatch'
    if (scalar_weight <= 0.0_rp) error stop 'ROM scalar metric must be positive'

    call neko_scratch_registry%request_vector(pointwise, pointwise_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    allocate(local(re*re)); local = 0.0_rp
    do j = 1, re
       pointwise = 0.0_rp
       do k = 1, rf
          do direction = 1, 3
             call vector_col3(product, perturb_velocity_modes(direction,j), &
                  forward_scalar_gradients(k,direction))
             call vector_add2s2(pointwise, product, coefficients(k))
             call vector_col3(product, forward_velocity_modes(direction,k), &
                  perturb_scalar_gradients(j,direction))
             call vector_add2s2(pointwise, product, coefficients(k))
          end do
       end do
       do i = 1, re
          local((j-1)*re+i) = -scalar_weight * &
               rom_glsc3_local(test_scalar_modes(i), weights, pointwise)
       end do
    end do
    call rom_reduce_buffer(local)
    do j = 1, re
       do i = 1, re
          linearized%x(i,j) = linearized%x(i,j) + local((j-1)*re+i)
       end do
    end do
    deallocate(local)
    call neko_scratch_registry%relinquish_vector([pointwise_index, product_index])
    nullify(pointwise); nullify(product)
  end subroutine rom_add_linearized_scalar_advection

  !> Assemble the independently projected velocity/scalar convection Jacobian
  !! tensor C(:,:,k) about a forward trajectory mode k. The pair fields are
  !! formed once for each (adjoint perturbation j, forward mode k), then
  !! projected against every adjoint test mode. The tensor is reduced as one
  !! packed local operator on the host backend.
  !!
  !! The returned dynamic contribution is
  !! C_ijk = -<w_i,(v_k.grad)w_j+(w_j.grad)v_k>
  !!         -gamma<zeta_i,v_k.grad(zeta_j)+w_j.grad(psi_k)>.
  !! Lift-linearised terms belong to the fixed adjoint linear operator and
  !! remain there; this tensor contains only terms proportional to a_k.
  subroutine rom_assemble_linearized_convection_tensor(weights, &
       forward_velocity_modes, forward_velocity_gradients, &
       forward_scalar_modes, forward_scalar_gradients, &
       adjoint_velocity_modes, adjoint_velocity_gradients, &
       adjoint_scalar_modes, adjoint_scalar_gradients, scalar_weight, tensor)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: forward_velocity_modes(:, :)
    type(vector_t), intent(in) :: forward_velocity_gradients(:, :, :)
    type(vector_t), intent(in), optional :: forward_scalar_modes(:)
    type(vector_t), intent(in), optional :: forward_scalar_gradients(:, :)
    type(vector_t), intent(in) :: adjoint_velocity_modes(:, :)
    type(vector_t), intent(in) :: adjoint_velocity_gradients(:, :, :)
    type(vector_t), intent(in), optional :: adjoint_scalar_modes(:)
    type(vector_t), intent(in), optional :: adjoint_scalar_gradients(:, :)
    real(kind=rp), intent(in) :: scalar_weight
    type(matrix_t), intent(inout) :: tensor(:)
    type(vector_t), allocatable :: velocity_pair(:)
    type(vector_t), pointer :: scalar_pair, product
    integer :: scalar_pair_index, product_index
    integer :: n, r, e, i, j, k, d, component
    logical :: has_scalar

    n = weights%size()
    r = size(forward_velocity_modes,2)
    e = size(adjoint_velocity_modes,2)
    if (r < 1 .or. e < 1 .or. size(forward_velocity_modes,1) /= 3 .or. &
         size(adjoint_velocity_modes,1) /= 3 .or. &
         size(forward_velocity_gradients,1) /= 3 .or. &
         size(forward_velocity_gradients,2) /= r .or. &
         size(forward_velocity_gradients,3) /= 3 .or. &
         size(adjoint_velocity_gradients,1) /= 3 .or. &
         size(adjoint_velocity_gradients,2) /= e .or. &
         size(adjoint_velocity_gradients,3) /= 3 .or. size(tensor) /= r) &
         error stop 'ROM projected convection tensor shape mismatch'
    has_scalar = present(forward_scalar_modes) .or. present(forward_scalar_gradients) .or. &
         present(adjoint_scalar_modes) .or. present(adjoint_scalar_gradients)
    if (has_scalar) then
       if (.not. present(forward_scalar_modes) .or. &
            .not. present(forward_scalar_gradients) .or. &
            .not. present(adjoint_scalar_modes) .or. &
            .not. present(adjoint_scalar_gradients)) &
            error stop 'ROM projected scalar convection tensor arguments are incomplete'
       if (size(forward_scalar_modes) /= r .or. &
            size(forward_scalar_gradients,1) /= r .or. &
            size(forward_scalar_gradients,2) /= 3 .or. &
            size(adjoint_scalar_modes) /= e .or. &
            size(adjoint_scalar_gradients,1) /= e .or. &
            size(adjoint_scalar_gradients,2) /= 3 .or. scalar_weight <= 0.0_rp) &
            error stop 'ROM projected scalar convection tensor shape mismatch'
    end if
    do k = 1, r
       if (tensor(k)%get_nrows() /= e .or. tensor(k)%get_ncols() /= e) &
            error stop 'ROM projected convection tensor output shape mismatch'
       tensor(k)%x = 0.0_rp
    end do

    allocate(velocity_pair(3))
    do component = 1, 3
       call velocity_pair(component)%init(n)
    end do
    nullify(scalar_pair)
    if (has_scalar) &
         call neko_scratch_registry%request_vector(scalar_pair, scalar_pair_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)

    do k = 1, r
       do j = 1, e
          do component = 1, 3
             velocity_pair(component) = 0.0_rp
          end do
          if (has_scalar) scalar_pair = 0.0_rp
          do d = 1, 3
             do component = 1, 3
                call vector_col3(product, adjoint_velocity_modes(d,j), &
                     forward_velocity_gradients(component,k,d))
                call vector_add2(velocity_pair(component), product)
                call vector_col3(product, forward_velocity_modes(d,k), &
                     adjoint_velocity_gradients(component,j,d))
                call vector_add2(velocity_pair(component), product)
             end do
             if (has_scalar) then
                call vector_col3(product, adjoint_velocity_modes(d,j), &
                     forward_scalar_gradients(k,d))
                call vector_add2(scalar_pair, product)
                call vector_col3(product, forward_velocity_modes(d,k), &
                     adjoint_scalar_gradients(j,d))
                call vector_add2(scalar_pair, product)
             end if
          end do
          do i = 1, e
             do component = 1, 3
                tensor(k)%x(i,j) = tensor(k)%x(i,j) - rom_glsc3_local( &
                     adjoint_velocity_modes(component,i), weights, &
                     velocity_pair(component))
             end do
             if (has_scalar) tensor(k)%x(i,j) = tensor(k)%x(i,j) - scalar_weight * &
                  rom_glsc3_local(adjoint_scalar_modes(i), weights, scalar_pair)
          end do
       end do
    end do
    call rom_reduce_matrix_tensor(tensor)

    if (has_scalar) then
       call neko_scratch_registry%relinquish_vector(scalar_pair_index)
       nullify(scalar_pair)
    end if
    call neko_scratch_registry%relinquish_vector(product_index)
    nullify(product)
    do component = 1, 3
       call velocity_pair(component)%free()
    end do
    deallocate(velocity_pair)
  end subroutine rom_assemble_linearized_convection_tensor

  !> Assemble scalar advection and diffusion in the common ROM basis.
  !!
  !! For
  !!   phi = phi_L + sum_j a_j phi_j and
  !!   u   = u_L   + sum_j a_j u_j,
  !! this adds the negative weak residual of u.grad(phi) and the implicit
  !! diffusion matrix to the same reduced rows used by the velocity equation.
  !! The routine deliberately takes spatial diffusivity so nonuniform material
  !! properties follow Neko's scalar scheme rather than a guessed Péclet
  !! constant.
  subroutine rom_assemble_scalar_transport(weights, scalar_weight, velocity_modes, &
       scalar_modes, velocity_gradients, scalar_gradients, lift_velocity, &
       scalar_lift, lift_scalar_gradient, diffusivity, linear, constant, quadratic, &
       implicit_linear)
    type(vector_t), intent(in) :: weights
    real(kind=rp), intent(in) :: scalar_weight
    type(vector_t), intent(in) :: velocity_modes(:, :), scalar_modes(:)
    type(vector_t), intent(in) :: velocity_gradients(:, :, :)
    type(vector_t), intent(in) :: scalar_gradients(:, :)
    type(vector_t), intent(in) :: lift_velocity(:), scalar_lift, lift_scalar_gradient(:)
    type(vector_t), intent(in) :: diffusivity
    type(matrix_t), intent(inout) :: linear, implicit_linear
    type(vector_t), intent(inout) :: constant
    type(matrix_t), intent(inout) :: quadratic(:)
    call rom_assemble_scalar_transport_petrov(weights, scalar_weight, velocity_modes, &
         scalar_modes, velocity_gradients, scalar_gradients, scalar_modes, scalar_gradients, &
         lift_velocity, &
         scalar_lift, lift_scalar_gradient, diffusivity, linear, constant, quadratic, &
         implicit_linear)
  end subroutine rom_assemble_scalar_transport

  !> Assemble the design-independent scalar diffusion matrix and lift vector.
  !! Diffusivity is a nodal field, so this retains Neko's variable-property
  !! SEM quadrature rather than replacing it with a uniform coefficient.
  !! Both outputs use the negative weak-residual sign of the ROM operators.
  subroutine rom_assemble_scalar_diffusion_petrov(weights, scalar_weight, &
       diffusivity, test_gradients, trial_gradients, lift_gradient, &
       implicit_linear, lift_constant)
    type(vector_t), intent(in) :: weights, diffusivity
    real(kind=rp), intent(in) :: scalar_weight
    type(vector_t), intent(in) :: test_gradients(:, :), trial_gradients(:, :)
    type(vector_t), intent(in) :: lift_gradient(:)
    type(matrix_t), intent(inout) :: implicit_linear
    type(vector_t), intent(inout) :: lift_constant
    type(vector_t), pointer :: effective
    real(kind=rp), allocatable :: local(:)
    integer :: effective_index
    integer :: n, nt, nr, i, j, d, matrix_size

    n = weights%size(); nt = size(test_gradients,1); nr = size(trial_gradients,1)
    if (n < 1 .or. diffusivity%size() /= n .or. nt < 1 .or. nr < 1 .or. &
         size(test_gradients,2) /= 3 .or. size(trial_gradients,2) /= 3 .or. &
         size(lift_gradient) /= 3 .or. scalar_weight <= 0.0_rp) &
         error stop 'ROM scalar diffusion shape/metric mismatch'
    if (any([(test_gradients(i,1)%size(), i=1,nt)] /= n) .or. &
         any([(trial_gradients(j,1)%size(), j=1,nr)] /= n) .or. &
         any([(lift_gradient(d)%size(), d=1,3)] /= n)) &
         error stop 'ROM scalar diffusion vector size mismatch'

    matrix_size = nt*nr
    allocate(local(matrix_size+nt)); local = 0.0_rp
    call neko_scratch_registry%request_vector(effective, effective_index, n, .true.)
    call vector_col3(effective, weights, diffusivity)
    call vector_cmult(effective, scalar_weight)
    do j = 1, nr
       do i = 1, nt
          do d = 1, 3
             local((j-1)*nt+i) = local((j-1)*nt+i) - &
                  rom_glsc3_local(test_gradients(i,d), effective, trial_gradients(j,d))
          end do
       end do
    end do
    do i = 1, nt
       do d = 1, 3
          local(matrix_size+i) = local(matrix_size+i) - &
               rom_glsc3_local(test_gradients(i,d), effective, lift_gradient(d))
       end do
    end do
    call rom_reduce_buffer(local)

    call implicit_linear%init(nt,nr)
    call lift_constant%init(nt)
    do j = 1, nr
       do i = 1, nt
          implicit_linear%x(i,j) = local((j-1)*nt+i)
       end do
    end do
    do i = 1, nt
       lift_constant%x(i) = local(matrix_size+i)
    end do
    call neko_scratch_registry%relinquish_vector(effective_index)
    nullify(effective)
    deallocate(local)
  end subroutine rom_assemble_scalar_diffusion_petrov

  !> Assemble scalar transport with trial velocity/scalar states and test
  !! scalar modes.  The scalar equation shares the Petrov test index with the
  !! coupled reduced state.
  subroutine rom_assemble_scalar_transport_petrov(weights, scalar_weight, &
       trial_velocity_modes, trial_scalar_modes, trial_velocity_gradients, &
       trial_scalar_gradients, test_scalar_modes, test_scalar_gradients, lift_velocity, scalar_lift, &
       lift_scalar_gradient, diffusivity, linear, constant, quadratic, implicit_linear, &
       include_diffusion)
    type(vector_t), intent(in) :: weights
    real(kind=rp), intent(in) :: scalar_weight
    type(vector_t), intent(in) :: trial_velocity_modes(:, :), trial_scalar_modes(:)
    type(vector_t), intent(in) :: trial_velocity_gradients(:, :, :)
    type(vector_t), intent(in) :: trial_scalar_gradients(:, :)
    type(vector_t), intent(in) :: test_scalar_modes(:)
    type(vector_t), intent(in) :: test_scalar_gradients(:, :)
    type(vector_t), intent(in) :: lift_velocity(:), scalar_lift, lift_scalar_gradient(:)
    type(vector_t), intent(in) :: diffusivity
    type(matrix_t), intent(inout) :: linear, implicit_linear
    type(vector_t), intent(inout) :: constant
    type(matrix_t), intent(inout), optional :: quadratic(:)
    logical, intent(in), optional :: include_diffusion
    type(vector_t), pointer :: effective, pointwise, product
    real(kind=rp), allocatable :: local_operators(:)
    integer :: effective_index, pointwise_index, product_index
    integer :: n, r, i, j, k, direction, q_count
    integer :: linear_offset, q_offset, implicit_offset
    logical :: add_diffusion

    n = weights%size(); r = size(trial_scalar_modes)
    if (size(trial_velocity_modes,1) /= 3 .or. size(trial_velocity_modes,2) /= r .or. &
         size(test_scalar_modes) /= r .or. size(trial_velocity_gradients,1) /= 3 .or. &
         size(trial_velocity_gradients,2) /= r .or. size(trial_velocity_gradients,3) /= 3 .or. &
         size(trial_scalar_gradients,1) /= r .or. size(trial_scalar_gradients,2) /= 3 .or. &
         size(test_scalar_gradients,1) /= r .or. size(test_scalar_gradients,2) /= 3 .or. &
         size(lift_velocity) /= 3 .or. &
         scalar_lift%size() /= n .or. &
         size(lift_scalar_gradient) /= 3 .or. diffusivity%size() /= n) then
       error stop 'ROM scalar transport shape mismatch'
    end if
    if (linear%get_nrows() /= r .or. linear%get_ncols() /= r .or. &
         implicit_linear%get_nrows() /= r .or. implicit_linear%get_ncols() /= r) then
       error stop 'ROM scalar transport operator shape mismatch'
    end if
    if (present(quadratic)) then
       if (size(quadratic) /= r) error stop 'ROM scalar transport quadratic shape mismatch'
       if (any([(quadratic(i)%get_nrows(), i=1,r)] /= r) .or. &
            any([(quadratic(i)%get_ncols(), i=1,r)] /= r)) &
            error stop 'ROM scalar transport quadratic shape mismatch'
    end if
    add_diffusion = .true.
    if (present(include_diffusion)) add_diffusion = include_diffusion
    call neko_scratch_registry%request_vector(effective, effective_index, n, .true.)
    call neko_scratch_registry%request_vector(pointwise, pointwise_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    effective = weights
    call vector_cmult(effective, scalar_weight)

    q_count = 0
    if (present(quadratic)) q_count = r*r*r
    linear_offset = r
    q_offset = linear_offset + r*r
    implicit_offset = q_offset + q_count
    allocate(local_operators(implicit_offset + r*r))
    local_operators = 0.0_rp

    ! The outputs are small host-owned reduced coefficients.  This routine
    ! accumulates the scalar block into the velocity operators assembled by
    ! the preceding routines, so their host values are intentionally retained.

    ! Lift--lift and lift--mode advection terms.  The lift--lift term must use
    ! grad(phi_L); using grad(phi_i) here would silently turn a constant
    ! residual into a state-dependent term and is not the Galerkin expansion
    ! of u.grad(phi).
    pointwise = 0.0_rp
    do direction = 1, 3
       call vector_col3(product, lift_velocity(direction), &
            lift_scalar_gradient(direction))
       call vector_add2(pointwise, product)
    end do
    do i = 1, r
       local_operators(i) = -rom_glsc3_local(test_scalar_modes(i), effective, pointwise)
    end do
    do j = 1, r
       pointwise = 0.0_rp
       do direction = 1, 3
          call vector_col3(product, trial_velocity_modes(direction,j), &
               lift_scalar_gradient(direction))
          call vector_add2(pointwise, product)
          call vector_col3(product, lift_velocity(direction), &
               trial_scalar_gradients(j,direction))
          call vector_add2(pointwise, product)
       end do
       do i = 1, r
          local_operators(linear_offset + (j-1)*r + i) = &
               -rom_glsc3_local(test_scalar_modes(i), effective, pointwise)
       end do
    end do

    ! Nonlinear scalar advection u_j.grad(phi_k).
    if (present(quadratic)) then
       do j = 1, r
          do k = 1, r
             pointwise = 0.0_rp
             do direction = 1, 3
                call vector_col3(product, trial_velocity_modes(direction,j), &
                     trial_scalar_gradients(k,direction))
                call vector_add2(pointwise, product)
             end do
             do i = 1, r
                local_operators(q_offset + (k-1)*r*r + (j-1)*r + i) = &
                     -rom_glsc3_local(test_scalar_modes(i), effective, pointwise)
             end do
          end do
       end do
    end if

    ! Scalar diffusion is implicit in the BDF1 map, exactly as the velocity
    ! diffusion assembled by rom_assemble_diffusion.
    if (add_diffusion) then
       call vector_col3(effective, weights, diffusivity)
       call vector_cmult(effective, scalar_weight)
       do j = 1, r
          do i = 1, r
             do direction = 1, 3
                local_operators(implicit_offset + (j-1)*r + i) = &
                     local_operators(implicit_offset + (j-1)*r + i) - &
                     rom_glsc3_local(test_scalar_gradients(i,direction), effective, &
                     trial_scalar_gradients(j,direction))
             end do
          end do
       end do
       do i = 1, r
          do direction = 1, 3
             ! The lift diffusion is a state-independent residual term.  It
             ! belongs in the explicit constant vector; fixed_implicit is
             ! reserved for the homogeneous scalar fluctuation.
             local_operators(i) = local_operators(i) - &
                  rom_glsc3_local(test_scalar_gradients(i,direction), effective, &
                  lift_scalar_gradient(direction))
          end do
       end do
    end if
    call rom_reduce_buffer(local_operators)
    do i = 1, r
       constant%x(i) = constant%x(i) + local_operators(i)
    end do
    do j = 1, r
       do i = 1, r
          linear%x(i,j) = linear%x(i,j) + &
               local_operators(linear_offset + (j-1)*r + i)
          implicit_linear%x(i,j) = implicit_linear%x(i,j) + &
               local_operators(implicit_offset + (j-1)*r + i)
       end do
    end do
    if (present(quadratic)) then
       do k = 1, r
          do j = 1, r
             do i = 1, r
                quadratic(i)%x(j,k) = quadratic(i)%x(j,k) + &
                     local_operators(q_offset + (k-1)*r*r + (j-1)*r + i)
             end do
          end do
       end do
    end if
    call neko_scratch_registry%relinquish_vector([effective_index, pointwise_index, &
         product_index])
    nullify(effective); nullify(pointwise); nullify(product)
    deallocate(local_operators)
  end subroutine rom_assemble_scalar_transport_petrov

end module rom_assembly
