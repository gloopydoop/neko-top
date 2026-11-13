!> @file URANS_POD_state_recover.f90
!! @copyright
!! Copyright (c) 2025-2026, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> @brief Fixed-window URANS-POD state recovery using in-situ statistics.
module simulation_URANS_POD_state_recover
  use num_types, only: rp, sp, dp
  use case, only: case_t
  use json_file_module, only: json_file
  use json_utils, only: json_get, json_get_or_default
  use field, only: field_t
  use registry, only: neko_registry
  use file, only: file_t
  use matrix, only: matrix_t
  use field_output, only: field_output_t
  use coefs, only: coef_t
  use data_streamer, only: data_streamer_t
  use profiler, only: profiler_start_region, profiler_end_region
  use state_recover, only: state_recover_t
  use time_state, only: time_state_t
  use logger, only : neko_log
  use comm, only: neko_comm, mpi_real_precision
  use mpi_f08, only: MPI_Allreduce, MPI_IN_PLACE, MPI_SUM

  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: DEVICE_TO_HOST, device_memcpy, HOST_TO_DEVICE
  use utils, only: neko_error
  use vector, only: vector_t
  use field_math, only: field_add2s2, field_rzero

  use neko_ctrl_mod, only: ctrl_stream_t, &
       MODE_IDLE, MODE_FORWARD, MODE_ADJOINT, MODE_STOP, &
       PHASE_INIT, PHASE_FWD_RUNNING, PHASE_FWD_DONE, &
       PHASE_ADJ_RUNNING, PHASE_ADJ_DONE
  use, intrinsic :: iso_fortran_env, only: int32, real64

  implicit none
  private

  ! Control protocol used by save()/restore() and the Python POD driver:
  ! - MODE_IDLE + PHASE_INIT is a startup tick sent after ctrl init.
  ! - MODE_FORWARD + PHASE_FWD_RUNNING is emitted before each streamed forward
  !   snapshot batch.
  ! - MODE_FORWARD + PHASE_FWD_DONE marks the forward-to-adjoint boundary.
  !   Neko then blocks until Python has updated the POD basis and responds with
  !   MODE_ADJOINT + PHASE_ADJ_RUNNING.
  ! - MODE_ADJOINT + PHASE_ADJ_RUNNING means restore() may reconstruct states
  !   from the received POD basis and time coefficients.
  ! - MODE_ADJOINT + PHASE_ADJ_DONE is emitted on reset after an adjoint pass
  !   so that Python can clear its old POD state before the next forward pass.
  ! - MODE_STOP is emitted during free() to terminate the Python control loop.
  !> URANS-POD state recovery implementation for forward/adjoint runs.
  type, public, extends(state_recover_t) :: URANS_POD_state_recover_t
     private

     logical :: enabled = .true.

     integer :: n_modes
     integer :: n_flds = 3
     integer :: pod_tstep = 0
     integer :: window_steps = 0
     integer :: window_step = 0
     integer :: window_count = 0
     real(kind=rp) :: window_duration = 0.0_rp
     real(kind=rp), allocatable :: sum_u(:,:,:,:), sum_v(:,:,:,:), sum_w(:,:,:,:)
     real(kind=rp), allocatable :: sum_uu(:,:,:,:), sum_vv(:,:,:,:), sum_ww(:,:,:,:)
     real(kind=rp), allocatable :: sum_uv(:,:,:,:), sum_uw(:,:,:,:), sum_vw(:,:,:,:)
     real(kind=rp), allocatable :: sum_s(:,:,:,:), sum_us(:,:,:,:), sum_vs(:,:,:,:), sum_ws(:,:,:,:)
     logical :: include_scalar = .false.
     logical :: include_nut = .false.
     logical :: include_gammat = .false.
     character(len=16) :: dtype = "double"
     type(coef_t), pointer :: coef => null()

     ! SINGLE lifetime streamer (both directions)
     type(data_streamer_t) :: dstream

     ! POD mode storage
     type(field_t), allocatable :: u_modes(:), v_modes(:), w_modes(:)
     type(field_t), allocatable :: s_modes(:)
     type(field_t), allocatable :: nut_modes(:)
     type(field_t), allocatable :: gammat_modes(:)
     type(field_t), pointer :: s => null()
     type(field_t), pointer :: nut => null()
     type(field_t), pointer :: gammat => null()

     ! CSV coeffs
     type(file_t)   :: csv_reader
     type(matrix_t) :: time_coefs
     type(vector_t) :: a_interp

     ! optional mode output configuration
     integer :: mode_output_precision = sp
     character(len=32) :: mode_output_format = 'fld'
     character(len=80) :: mode_file_name = 'POD_modes'
     logical :: write_modes = .false.
     logical :: output_reconstruction = .false.
     type(field_output_t) :: recon_output
     integer :: recon_output_precision = sp
     character(len=32) :: recon_output_format = 'fld'
     character(len=80) :: recon_file_name = 'pod_reconstruction'
     character(len=16) :: recon_output_control = "never"
     real(kind=rp) :: recon_output_value = 0.0_rp
     real(kind=rp) :: recon_time_interval = 0.0_rp
     integer :: recon_nsteps = 0

     ! Control state
     type(ctrl_stream_t) :: ctrl

     logical :: have_received_modes = .false.
     logical :: adjoint_started = .false.
     logical :: adj_running_sent = .false.
   contains
     procedure, public, pass(this) :: init => URANS_POD_state_recover_init_from_json
     procedure, public, pass(this) :: init_from_json => &
          URANS_POD_state_recover_init_from_json
     procedure, public, pass(this) :: init_from_components => &
          URANS_POD_state_recover_init_from_components
     procedure, public, pass(this) :: free => URANS_POD_state_recover_free
     procedure, public, pass(this) :: reset => URANS_POD_state_recover_reset
     procedure, public, pass(this) :: save => URANS_POD_state_recover_save
     procedure, public, pass(this) :: restore => URANS_POD_state_recover_restore
  end type URANS_POD_state_recover_t

contains

  !> Initialize URANS-POD state recovery from JSON parameters.
  !! @param[inout] this URANS-POD state recovery instance.
  !! @param[inout] neko_case Case data structure.
  !! @param[inout] params JSON parameters.
  subroutine URANS_POD_state_recover_init_from_json(this, neko_case, params)
    class(URANS_POD_state_recover_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(json_file), intent(inout) :: params
    integer :: n_modes
    logical :: write_modes
    logical :: debug
    logical :: output_reconstruction
    character(len=:), allocatable :: recon_output_precision_str
    character(len=:), allocatable :: output_control
    real(kind=rp) :: output_value
    character(len=:), allocatable :: dtype
    character(len=:), allocatable :: mode_output_precision_str
    character(len=:), allocatable :: mode_output_format
    character(len=:), allocatable :: mode_file_name
    character(len=:), allocatable :: recon_output_format
    character(len=:), allocatable :: recon_file_name
    integer :: mode_output_precision
    integer :: recon_output_precision
    logical :: fit_eddy_viscosity
    logical :: fit_scalar_diffusivity
    character(len=:), allocatable :: scalar_closure
    real(kind=rp) :: window_duration

    call json_get(params, "n_modes", n_modes)
    call json_get(params, "window_duration", window_duration)
    call json_get_or_default(params, "dtype", dtype, "double")
    call json_get_or_default(params, "write_modes", write_modes, .false.)
    call json_get_or_default(params, "fit_eddy_viscosity", fit_eddy_viscosity, .false.)
    call json_get_or_default(params, "scalar_closure", scalar_closure, "off")
    fit_scalar_diffusivity = trim(lower_string(scalar_closure)) .ne. 'off' .and. &
         trim(lower_string(scalar_closure)) .ne. 'none'
    call json_get_or_default(params, "debug", debug, .false.)
    call json_get_or_default(params, "output_reconstruction", &
         output_reconstruction, .false.)
    call json_get_or_default(params, "output_precision", &
         mode_output_precision_str, "sp")
    call json_get_or_default(params, "output_format", &
         mode_output_format, "fld")
    call json_get_or_default(params, "output_file_name", mode_file_name, &
         "POD_modes")

    mode_output_precision = precision_from_string(mode_output_precision_str, &
         'state_recovery.output_precision')

    recon_output_precision = mode_output_precision
    recon_output_precision_str = trim(mode_output_precision_str)
    recon_output_format = trim(mode_output_format)
    recon_file_name = 'pod_reconstruction'
    output_control = 'never'
    output_value = 0.0_rp

    if (output_reconstruction) then
       if (.not. ('output_precision' .in. params)) then
          call json_get_or_default(neko_case%params, 'case.output_precision', &
               recon_output_precision_str, 'single')
          recon_output_precision = precision_from_string( &
               recon_output_precision_str, &
               'state_recovery.reconstruction_output_precision')
       end if
       if ('reconstruction_output_precision' .in. params) then
          call json_get(params, 'reconstruction_output_precision', &
               recon_output_precision_str)
          recon_output_precision = precision_from_string( &
               recon_output_precision_str, &
               'state_recovery.reconstruction_output_precision')
       end if
       if ('reconstruction_output_format' .in. params) then
          call json_get(params, 'reconstruction_output_format', &
               recon_output_format)
       end if
       call json_get_or_default(params, 'reconstruction_output_file_name', &
            recon_file_name, 'pod_reconstruction')

       call json_get_or_default(neko_case%params, 'case.fluid.output_control', &
            output_control, 'org')
       select case (trim(lower_string(output_control)))
       case ('org')
          call json_get(neko_case%params, 'case.nsamples', output_value)
       case ('nsamples', 'simulationtime', 'tsteps')
          call json_get(neko_case%params, 'case.fluid.output_value', &
               output_value)
       case ('never')
          call json_get_or_default(neko_case%params, &
               'case.fluid.output_value', output_value, 0.0_rp)
       case default
          call json_get_or_default(neko_case%params, &
               'case.fluid.output_value', output_value, 0.0_rp)
       end select
    end if

    call this%init_from_components(neko_case, window_duration, n_modes, dtype, &
         write_modes, output_reconstruction, output_control, output_value, &
         debug, mode_output_precision, &
         mode_output_format, mode_file_name, recon_output_precision, &
         recon_output_format, recon_file_name, fit_eddy_viscosity, &
         fit_scalar_diffusivity)
  end subroutine URANS_POD_state_recover_init_from_json


  !> Initialize URANS-POD state recovery from explicit components.
  !! @param[inout] this URANS-POD state recovery instance.
  !! @param[inout] neko_case Case data structure.
  !! @param[in] window_duration Physical length of each statistics window.
  !! @param[in] n_modes Number of POD modes to keep.
  !! @param[in] dtype POD floating-point precision.
  !! @param[in] write_modes Whether Python should persist POD modes.
  !! @param[in] output_reconstruction Whether to output reconstructions.
  !! @param[in] output_control Reconstruction output cadence control.
  !! @param[in] output_value Reconstruction cadence value.
  !! @param[in] debug Optional debug flag.
  !! @param[in] mode_output_precision Optional POD mode output precision.
  !! @param[in] mode_output_format Optional POD mode output format.
  !! @param[in] mode_file_name Optional POD mode output base name.
  !! @param[in] recon_output_precision Optional reconstruction precision.
  !! @param[in] recon_output_format Optional reconstruction output format.
  !! @param[in] recon_file_name Optional reconstruction output base name.
  subroutine URANS_POD_state_recover_init_from_components(this, neko_case, &
       window_duration, n_modes, dtype, write_modes, output_reconstruction, &
       output_control, output_value, debug, &
       mode_output_precision, mode_output_format, mode_file_name, &
       recon_output_precision, recon_output_format, recon_file_name, &
       fit_eddy_viscosity, fit_scalar_diffusivity)
    class(URANS_POD_state_recover_t), intent(inout), target :: this
    class(case_t), target, intent(inout) :: neko_case
    real(kind=rp), intent(in) :: window_duration
    integer, intent(in) :: n_modes
    character(len=*), intent(in) :: dtype
    logical, intent(in) :: write_modes
    logical, intent(in) :: output_reconstruction
    character(len=*), intent(in) :: output_control
    real(kind=rp), intent(in) :: output_value
    logical, intent(in), optional :: debug
    integer, intent(in), optional :: mode_output_precision
    character(len=*), intent(in), optional :: mode_output_format
    character(len=*), intent(in), optional :: mode_file_name
    integer, intent(in), optional :: recon_output_precision
    character(len=*), intent(in), optional :: recon_output_format
    character(len=*), intent(in), optional :: recon_file_name
    logical, intent(in), optional :: fit_eddy_viscosity
    logical, intent(in), optional :: fit_scalar_diffusivity
    integer :: i
    character(len=80) :: str

    this%neko_case => neko_case
    this%enabled = .true.
    this%n_modes = n_modes
    this%window_duration = window_duration
    this%dtype = adjustl(dtype)
    this%write_modes = write_modes
    this%output_reconstruction = output_reconstruction
    this%recon_output_control = 'never'
    this%recon_output_value = 0.0_rp
    this%recon_time_interval = 0.0_rp
    this%recon_nsteps = 0
    this%recon_output_precision = sp
    this%recon_output_format = 'fld'
    this%recon_file_name = 'pod_reconstruction'
    this%mode_output_precision = sp
    this%mode_output_format = 'fld'
    this%mode_file_name = 'POD_modes'
    this%n_flds = 3
    this%pod_tstep = 0
    this%window_step = 0
    this%window_count = 0
    this%coef => neko_case%fluid%c_Xh

    if (this%window_duration .le. 0.0_rp) then
       call neko_error('URANS-POD requires a positive window_duration.')
    end if
    this%window_steps = nint(this%window_duration / neko_case%time%dt)
    if (this%window_steps .lt. 1 .or. abs(real(this%window_steps, rp) * &
         neko_case%time%dt - this%window_duration) .gt. &
         100.0_rp * epsilon(this%window_duration)) then
       call neko_error('URANS-POD window_duration must be an integer number ' // &
            'of solver timesteps.')
    end if
    if (abs((neko_case%time%end_time - neko_case%time%start_time) / &
         this%window_duration - nint((neko_case%time%end_time - &
         neko_case%time%start_time) / this%window_duration)) .gt. &
         100.0_rp * epsilon(this%window_duration)) then
       call neko_error('URANS-POD final time must be a multiple of ' // &
            'window_duration.')
    end if

    this%include_scalar = allocated(neko_case%scalars)
    this%include_nut = .false.
    if (present(fit_eddy_viscosity)) this%include_nut = fit_eddy_viscosity
    this%include_gammat = .false.
    if (present(fit_scalar_diffusivity)) this%include_gammat = fit_scalar_diffusivity
    this%n_flds = 3 + merge(1, 0, this%include_scalar) + &
         merge(1, 0, this%include_nut) + merge(1, 0, this%include_gammat)

    allocate(this%sum_u, mold=neko_case%fluid%u%x)
    allocate(this%sum_v, mold=neko_case%fluid%v%x)
    allocate(this%sum_w, mold=neko_case%fluid%w%x)
    allocate(this%sum_uu, mold=neko_case%fluid%u%x)
    allocate(this%sum_vv, mold=neko_case%fluid%v%x)
    allocate(this%sum_ww, mold=neko_case%fluid%w%x)
    allocate(this%sum_uv, mold=neko_case%fluid%u%x)
    allocate(this%sum_uw, mold=neko_case%fluid%u%x)
    allocate(this%sum_vw, mold=neko_case%fluid%u%x)
    if (this%include_scalar) then
       allocate(this%sum_s, mold=neko_case%fluid%u%x)
       allocate(this%sum_us, mold=neko_case%fluid%u%x)
       allocate(this%sum_vs, mold=neko_case%fluid%u%x)
       allocate(this%sum_ws, mold=neko_case%fluid%u%x)
    end if
    call URANS_POD_state_recover_reset_window(this)

    if (this%include_nut) then
       call neko_registry%add_field(this%coef%dof, 'adjoint_nut', .true.)
       this%nut => neko_registry%get_field('adjoint_nut')
       call field_rzero(this%nut)
    end if
    if (this%include_gammat) then
       call neko_registry%add_field(this%coef%dof, 'adjoint_gammat', .true.)
       this%gammat => neko_registry%get_field('adjoint_gammat')
       call field_rzero(this%gammat)
    end if

    if (present(mode_output_precision)) then
       if (mode_output_precision .ne. sp .and. &
            mode_output_precision .ne. dp) then
          call neko_error('mode_output_precision must be either sp or dp')
       end if
       this%mode_output_precision = mode_output_precision
    end if
    if (present(mode_output_format)) then
       this%mode_output_format = trim(mode_output_format)
    end if
    if (present(mode_file_name)) then
       this%mode_file_name = trim(mode_file_name)
    end if
    if (present(recon_output_precision)) then
      if (recon_output_precision .ne. sp .and. &
           recon_output_precision .ne. dp) then
         call neko_error('recon_output_precision must be either sp or dp')
      end if
      this%recon_output_precision = recon_output_precision
    end if
    if (present(recon_output_format)) then
       this%recon_output_format = trim(recon_output_format)
    end if
    if (present(recon_file_name)) then
       this%recon_file_name = trim(recon_file_name)
    end if

    select case (trim(lower_string(this%dtype)))
    case ('single')
       if (rp .ne. sp) then
          call neko_error("POD dtype single but code not single precision.")
       end if
    case ('double')
       if (rp .ne. dp) then
          call neko_error("POD dtype double but code not double precision.")
       end if
    case default
       call neko_error("Unsupported POD dtype: " // trim(this%dtype))
    end select

    if (this%include_scalar) then
       if (size(neko_case%scalars%scalar_fields) .ne. 1) then
          call neko_error('POD scalar recovery currently supports ' // &
               'exactly one scalar.')
       end if
       this%s => neko_case%scalars%scalar_fields(1)%scalar%s
    else
       nullify(this%s)
    end if

    allocate(this%u_modes(this%n_modes))
    allocate(this%v_modes(this%n_modes))
    allocate(this%w_modes(this%n_modes))
    if (this%include_scalar) allocate(this%s_modes(this%n_modes))
    if (this%include_nut) allocate(this%nut_modes(this%n_modes))
    if (this%include_gammat) allocate(this%gammat_modes(this%n_modes))

    do i = 1, this%n_modes
       write(str, '(A,I0)') "u_mode_", i
       call this%u_modes(i)%init(this%coef%dof, trim(str))

       write(str, '(A,I0)') "v_mode_", i
       call this%v_modes(i)%init(this%coef%dof, trim(str))

       write(str, '(A,I0)') "w_mode_", i
       call this%w_modes(i)%init(this%coef%dof, trim(str))

       if (this%include_scalar) then
          write(str, '(A,I0)') "s_mode_", i
          call this%s_modes(i)%init(this%s%dof, trim(str))
       end if
       if (this%include_nut) then
          write(str, '(A,I0)') "nut_mode_", i
          call this%nut_modes(i)%init(this%coef%dof, trim(str))
       end if
       if (this%include_gammat) then
          write(str, '(A,I0)') "gammat_mode_", i
          call this%gammat_modes(i)%init(this%coef%dof, trim(str))
       end if
    end do

    if (this%output_reconstruction) then
       call this%recon_output%init(trim(this%recon_file_name), 3 + &
            merge(1, 0, this%include_scalar), &
            precision = this%recon_output_precision, &
            path = trim(neko_case%output_directory), &
            format = trim(this%recon_output_format))
       call this%recon_output%fields%assign_to_field(1, neko_case%fluid%u)
       call this%recon_output%fields%assign_to_field(2, neko_case%fluid%v)
       call this%recon_output%fields%assign_to_field(3, neko_case%fluid%w)
       if (this%include_scalar) then
          call this%recon_output%fields%assign_to_field(4, this%s)
       end if

       this%recon_output_control = trim(lower_string(output_control))
       this%recon_output_value = output_value

       select case (this%recon_output_control)
       case ('org')
          this%recon_output_control = 'nsamples'
          if (output_value .gt. 0.0_rp) then
             this%recon_time_interval = (neko_case%time%end_time - &
                  neko_case%time%start_time) / output_value
          end if
       case ('nsamples')
          if (output_value .gt. 0.0_rp) then
             this%recon_time_interval = (neko_case%time%end_time - &
                  neko_case%time%start_time) / output_value
          end if
       case ('simulationtime')
          this%recon_time_interval = output_value
       case ('tsteps')
          this%recon_nsteps = int(output_value)
       case ('never')
       case default
          call neko_error('Unsupported output_control for reconstruction: ' // &
               trim(output_control))
       end select
    end if

    call this%a_interp%init(this%n_modes)
    call this%csv_reader%init('pod_time_coeffs.csv')

    ! Single lifetime streamer: init once, keep open
    call this%dstream%init(this%coef)

    ! Send mesh once
    call this%dstream%stream(this%coef%dof%x%x)
    call this%dstream%stream(this%coef%dof%y%x)
    call this%dstream%stream(this%coef%dof%z%x)

    ! Stream the initial condition before the forward time loop begins.
    call URANS_POD_state_recover_stream_fields(this, neko_case)

    ! Control init for root-to-root MPI coordination with Python.
    if (present(debug)) then
       this%ctrl%debug = debug
    end if
    call this%ctrl%init()

    ! Fire an init tick (python might miss it; harmless)
    call this%ctrl%send(MODE_IDLE, PHASE_INIT, 0_int32, 0.0_real64)
    call this%set_n_timesteps(0)

  end subroutine URANS_POD_state_recover_init_from_components


  subroutine URANS_POD_state_recover_reset_window(this)
    class(URANS_POD_state_recover_t), intent(inout) :: this

    this%window_step = 0
    this%sum_u = 0.0_rp
    this%sum_v = 0.0_rp
    this%sum_w = 0.0_rp
    this%sum_uu = 0.0_rp
    this%sum_vv = 0.0_rp
    this%sum_ww = 0.0_rp
    this%sum_uv = 0.0_rp
    this%sum_uw = 0.0_rp
    this%sum_vw = 0.0_rp
    if (this%include_scalar) then
       this%sum_s = 0.0_rp
       this%sum_us = 0.0_rp
       this%sum_vs = 0.0_rp
       this%sum_ws = 0.0_rp
    end if
  end subroutine URANS_POD_state_recover_reset_window


  !> Release URANS-POD state recovery resources.
  !! @param[inout] this URANS-POD state recovery instance.
  subroutine URANS_POD_state_recover_free(this)
    class(URANS_POD_state_recover_t), intent(inout) :: this
    integer :: i

    if (this%ctrl%inited) then
       call this%ctrl%send(MODE_STOP, PHASE_ADJ_DONE, 0_int32, &
            0.0_real64)
    end if

    if (allocated(this%u_modes)) then
       do i = 1, size(this%u_modes); call this%u_modes(i)%free(); end do
       deallocate(this%u_modes)
    end if
    if (allocated(this%v_modes)) then
       do i = 1, size(this%v_modes); call this%v_modes(i)%free(); end do
       deallocate(this%v_modes)
    end if
    if (allocated(this%w_modes)) then
       do i = 1, size(this%w_modes); call this%w_modes(i)%free(); end do
       deallocate(this%w_modes)
    end if
    if (allocated(this%s_modes)) then
       do i = 1, size(this%s_modes); call this%s_modes(i)%free(); end do
       deallocate(this%s_modes)
    end if
    if (allocated(this%nut_modes)) then
       do i = 1, size(this%nut_modes); call this%nut_modes(i)%free(); end do
       deallocate(this%nut_modes)
    end if
    if (allocated(this%gammat_modes)) then
       do i = 1, size(this%gammat_modes); call this%gammat_modes(i)%free(); end do
       deallocate(this%gammat_modes)
    end if
    if (allocated(this%sum_u)) deallocate(this%sum_u)
    if (allocated(this%sum_v)) deallocate(this%sum_v)
    if (allocated(this%sum_w)) deallocate(this%sum_w)
    if (allocated(this%sum_uu)) deallocate(this%sum_uu)
    if (allocated(this%sum_vv)) deallocate(this%sum_vv)
    if (allocated(this%sum_ww)) deallocate(this%sum_ww)
    if (allocated(this%sum_uv)) deallocate(this%sum_uv)
    if (allocated(this%sum_uw)) deallocate(this%sum_uw)
    if (allocated(this%sum_vw)) deallocate(this%sum_vw)
    if (allocated(this%sum_s)) deallocate(this%sum_s)
    if (allocated(this%sum_us)) deallocate(this%sum_us)
    if (allocated(this%sum_vs)) deallocate(this%sum_vs)
    if (allocated(this%sum_ws)) deallocate(this%sum_ws)

    call this%recon_output%free()
    call this%dstream%free()
    call this%csv_reader%free()
    call this%a_interp%free()
    nullify(this%coef)
    nullify(this%neko_case)

    call this%ctrl%free()

    this%include_scalar = .false.
    this%include_nut = .false.
    this%include_gammat = .false.
    this%pod_tstep = 0
    this%mode_output_precision = sp
    this%mode_output_format = 'fld'
    this%mode_file_name = 'POD_modes'
    this%recon_output_precision = sp
    this%recon_output_format = 'fld'
    this%recon_file_name = 'pod_reconstruction'
    nullify(this%s)
    nullify(this%nut)
    nullify(this%gammat)
    this%enabled = .false.
  end subroutine URANS_POD_state_recover_free


  !> Reset URANS-POD state recovery to initial control state.
  !! @param[inout] this URANS-POD state recovery instance.
  subroutine URANS_POD_state_recover_reset(this)
    class(URANS_POD_state_recover_t), intent(inout) :: this
    integer :: i

    if (.not. this%enabled) return

    ! If we were in adjoint previously, emit ADJ_DONE once on reset
    if (this%ctrl%inited .and. this%adjoint_started) then
       call this%ctrl%send(MODE_ADJOINT, PHASE_ADJ_DONE, 0_int32, &
            0.0_real64)
    end if

    this%have_received_modes = .false.
    this%adjoint_started = .false.
    this%adj_running_sent = .false.
    this%pod_tstep = 0
    this%window_count = 0
    call URANS_POD_state_recover_reset_window(this)

    call this%set_n_timesteps(0)

    do i = 1, this%n_modes
       call field_rzero(this%u_modes(i))
       call field_rzero(this%v_modes(i))
       call field_rzero(this%w_modes(i))
       if (this%include_scalar) call field_rzero(this%s_modes(i))
       if (this%include_nut) call field_rzero(this%nut_modes(i))
       if (this%include_gammat) call field_rzero(this%gammat_modes(i))
    end do

  end subroutine URANS_POD_state_recover_reset

  subroutine URANS_POD_state_recover_copy_fields(this, neko_case)
    class(URANS_POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    type(field_t), pointer :: u, v, w, s
    integer :: n

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    if (this%include_scalar) s => this%s
    n = u%dof%size()

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u%x, u%x_d, n, DEVICE_TO_HOST, sync = .true.)
       call device_memcpy(v%x, v%x_d, n, DEVICE_TO_HOST, sync = .true.)
       call device_memcpy(w%x, w%x_d, n, DEVICE_TO_HOST, sync = .true.)
       if (this%include_scalar) then
          call device_memcpy(s%x, s%x_d, n, DEVICE_TO_HOST, sync = .true.)
       end if
    end if

  end subroutine URANS_POD_state_recover_copy_fields


  subroutine URANS_POD_state_recover_stream_fields(this, neko_case)
    class(URANS_POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    type(field_t), pointer :: u, v, w, s

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    if (this%include_scalar) s => this%s

    call URANS_POD_state_recover_copy_fields(this, neko_case)
    call this%dstream%stream(u%x)
    call this%dstream%stream(v%x)
    call this%dstream%stream(w%x)
    if (this%include_scalar) call this%dstream%stream(s%x)
  end subroutine URANS_POD_state_recover_stream_fields


  subroutine URANS_POD_state_recover_accumulate_window(this, neko_case)
    class(URANS_POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    type(field_t), pointer :: u, v, w, s

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    if (this%include_scalar) s => this%s

    call URANS_POD_state_recover_copy_fields(this, neko_case)
    this%sum_u = this%sum_u + u%x
    this%sum_v = this%sum_v + v%x
    this%sum_w = this%sum_w + w%x
    this%sum_uu = this%sum_uu + u%x * u%x
    this%sum_vv = this%sum_vv + v%x * v%x
    this%sum_ww = this%sum_ww + w%x * w%x
    this%sum_uv = this%sum_uv + u%x * v%x
    this%sum_uw = this%sum_uw + u%x * w%x
    this%sum_vw = this%sum_vw + v%x * w%x
    if (this%include_scalar) then
       this%sum_s = this%sum_s + s%x
       this%sum_us = this%sum_us + u%x * s%x
       this%sum_vs = this%sum_vs + v%x * s%x
       this%sum_ws = this%sum_ws + w%x * s%x
    end if
    this%window_step = this%window_step + 1
  end subroutine URANS_POD_state_recover_accumulate_window


  subroutine URANS_POD_state_recover_stream_window(this, neko_case)
    class(URANS_POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    real(kind=rp) :: scale

    scale = 1.0_rp / real(this%window_steps, rp)
    this%sum_u = this%sum_u * scale
    this%sum_v = this%sum_v * scale
    this%sum_w = this%sum_w * scale
    this%sum_uu = this%sum_uu * scale - this%sum_u * this%sum_u
    this%sum_vv = this%sum_vv * scale - this%sum_v * this%sum_v
    this%sum_ww = this%sum_ww * scale - this%sum_w * this%sum_w
    this%sum_uv = this%sum_uv * scale - this%sum_u * this%sum_v
    this%sum_uw = this%sum_uw * scale - this%sum_u * this%sum_w
    this%sum_vw = this%sum_vw * scale - this%sum_v * this%sum_w
    if (this%include_scalar) then
       this%sum_s = this%sum_s * scale
       this%sum_us = this%sum_us * scale - this%sum_u * this%sum_s
       this%sum_vs = this%sum_vs * scale - this%sum_v * this%sum_s
       this%sum_ws = this%sum_ws * scale - this%sum_w * this%sum_s
    end if

    this%window_count = this%window_count + 1
    if (this%ctrl%inited) then
       call this%ctrl%send(MODE_FORWARD, PHASE_FWD_RUNNING, &
            int(this%pod_tstep, int32), &
            real(neko_case%time%t - neko_case%time%start_time - &
            0.5_rp * this%window_duration, real64))
    end if

    call this%dstream%stream(this%sum_u)
    call this%dstream%stream(this%sum_v)
    call this%dstream%stream(this%sum_w)
    if (this%include_scalar) call this%dstream%stream(this%sum_s)
    call this%dstream%stream(this%sum_uu)
    call this%dstream%stream(this%sum_vv)
    call this%dstream%stream(this%sum_ww)
    call this%dstream%stream(this%sum_uv)
    call this%dstream%stream(this%sum_uw)
    call this%dstream%stream(this%sum_vw)
    if (this%include_scalar) then
       call this%dstream%stream(this%sum_us)
       call this%dstream%stream(this%sum_vs)
       call this%dstream%stream(this%sum_ws)
    end if
  end subroutine URANS_POD_state_recover_stream_window


  !> Stream forward state for POD updates.
  !! @param[inout] this URANS-POD state recovery instance.
  subroutine URANS_POD_state_recover_save(this)
    class(URANS_POD_state_recover_t), intent(inout) :: this
    class(case_t), pointer :: neko_case
    if (.not. this%enabled) return
    if (.not. associated(this%neko_case)) then
       call neko_error('URANS-POD state recovery is not associated with a case.')
    end if
    neko_case => this%neko_case

    this%pod_tstep = this%pod_tstep + 1
    call this%set_n_timesteps(max(this%get_n_timesteps(), this%pod_tstep))
    call neko_log%message(" ")
    call neko_log%message("----------------")
    call neko_log%message("Accumulating URANS-POD statistics")
    call neko_log%message("----------------")
    call neko_log%message(" ")

    call profiler_start_region("POD save")

    call URANS_POD_state_recover_accumulate_window(this, neko_case)
    if (this%window_step .eq. this%window_steps) then
       call URANS_POD_state_recover_stream_window(this, neko_case)
       call URANS_POD_state_recover_reset_window(this)
    end if

    call profiler_end_region("POD save")
  end subroutine URANS_POD_state_recover_save


  !> Reconstruct and restore state from POD during adjoint.
  !! @param[inout] this URANS-POD state recovery instance.
  !! @param[in] tstep Forward-state index to restore.
  subroutine URANS_POD_state_recover_restore(this, tstep)
    class(URANS_POD_state_recover_t), intent(inout) :: this
    integer, intent(in) :: tstep
    class(case_t), pointer :: neko_case
    type(time_state_t) :: time
    type(time_state_t) :: time_out
    real(kind=rp) :: t_pod
    integer :: n_case_timesteps

    if (.not. this%enabled) return
    if (.not. associated(this%neko_case)) then
       call neko_error('URANS-POD state recovery is not associated with a case.')
    end if
    neko_case => this%neko_case

    n_case_timesteps = nint((neko_case%time%end_time - &
         neko_case%time%start_time) / neko_case%time%dt)
    time = neko_case%time
    time%tstep = tstep
    time%t = neko_case%time%end_time - &
         real(n_case_timesteps - tstep, rp) * neko_case%time%dt

    ! First restore() call is the phase boundary forward->adjoint
    if (.not. this%have_received_modes) then
       call URANS_POD_state_recover_receive_modes(this, time)
    end if

    ! Emit ADJ_RUNNING only once (avoid flooding SST)
    if (this%ctrl%inited .and. .not. this%adj_running_sent) then
       call this%ctrl%send(MODE_ADJOINT, PHASE_ADJ_RUNNING, &
            int(time%tstep, int32), real(time%t, real64))
       this%adj_running_sent = .true.
    end if

    call profiler_start_region("POD restore")
    t_pod = time%t - time%start_time
    call interpolate_time_coeffs_vec(this%a_interp, this%time_coefs, t_pod)
    call reconstruct_from_coeffs(this, neko_case, this%a_interp)
    if (this%output_reconstruction) then
       if (recon_should_output(this, time, time_out)) then
          call this%recon_output%sample(time_out%t)
       end if
    end if
    call profiler_end_region("POD restore")
  end subroutine URANS_POD_state_recover_restore

  !> Receive POD modes at the forward-to-adjoint boundary.
  !! @param[inout] this URANS-POD state recovery instance.
  !! @param[in] time Target time state.
  subroutine URANS_POD_state_recover_receive_modes(this, time)
    class(URANS_POD_state_recover_t), intent(inout) :: this
    type(time_state_t), intent(in) :: time
    integer :: i, ierr, n_lines, nrows, ncols, n
    integer(int32) :: mode_cmd, phase_cmd

    call profiler_start_region("POD receive modes")

    if (this%window_step .ne. 0) then
       call neko_error('URANS-POD ended with an incomplete statistics window.')
    end if

    if (this%ctrl%inited) then
      call this%ctrl%send(MODE_FORWARD, PHASE_FWD_DONE, &
           int(time%tstep, int32), real(time%t, real64))

      mode_cmd  = MODE_FORWARD
      phase_cmd = PHASE_FWD_DONE

      ! BLOCK until Python says "go adjoint"
      call this%ctrl%receive(mode_cmd, phase_cmd)

      if (mode_cmd /= MODE_ADJOINT) then
         call neko_error('Expected MODE_ADJOINT from Python at ' // &
              'forward-to-adjoint boundary.')
      end if
    end if

    this%adjoint_started = .true.

    ! Receive modes (Python streams after sending cmd)
    do i = 1, this%n_modes
       call this%dstream%recieve(this%u_modes(i)%x)
       call this%dstream%recieve(this%v_modes(i)%x)
       call this%dstream%recieve(this%w_modes(i)%x)
       if (this%include_scalar) then
          call this%dstream%recieve(this%s_modes(i)%x)
       end if
       if (this%include_nut) call this%dstream%recieve(this%nut_modes(i)%x)
       if (this%include_gammat) call this%dstream%recieve(this%gammat_modes(i)%x)
    end do

    ! Move modes back to GPU
    n = this%u_modes(1)%dof%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       do i = 1, this%n_modes
          call device_memcpy(this%u_modes(i)%x, this%u_modes(i)%x_d, n, &
               HOST_TO_DEVICE, sync = .true.)
          call device_memcpy(this%v_modes(i)%x, this%v_modes(i)%x_d, n, &
               HOST_TO_DEVICE, sync = .true.)
          call device_memcpy(this%w_modes(i)%x, this%w_modes(i)%x_d, n, &
               HOST_TO_DEVICE, sync = .true.)
          if (this%include_scalar) then
             call device_memcpy(this%s_modes(i)%x, this%s_modes(i)%x_d, &
                  n, HOST_TO_DEVICE, sync = .true.)
          end if
          if (this%include_nut) then
             call device_memcpy(this%nut_modes(i)%x, this%nut_modes(i)%x_d, &
                  n, HOST_TO_DEVICE, sync = .true.)
          end if
          if (this%include_gammat) then
             call device_memcpy(this%gammat_modes(i)%x, this%gammat_modes(i)%x_d, &
                  n, HOST_TO_DEVICE, sync = .true.)
          end if
       end do
    end if

    ! Read CSV once
    n_lines = csv_file_count_lines(this%csv_reader)
    nrows = n_lines - 1
    ncols = 1 + this%n_modes
    call this%time_coefs%init(nrows, ncols)
    call this%csv_reader%read(this%time_coefs)

    call MPI_Allreduce(MPI_IN_PLACE, this%time_coefs%x, &
         this%time_coefs%size(), mpi_real_precision, MPI_SUM, &
         neko_comm, ierr)

    this%have_received_modes = .true.
    call profiler_end_region("POD receive modes")
  end subroutine URANS_POD_state_recover_receive_modes

  logical function recon_should_output(this, time, time_out)
    class(URANS_POD_state_recover_t), intent(in) :: this
    type(time_state_t), intent(in) :: time
    type(time_state_t), intent(out) :: time_out
    real(kind=rp) :: t_rel, tol, interval

    recon_should_output = .false.
    time_out = time
    time_out%t = time%start_time + real(time%tstep, rp) * time%dt

    select case (this%recon_output_control)
    case ('tsteps')
       if (this%recon_nsteps .gt. 0) then
          recon_should_output = mod(time%tstep, this%recon_nsteps) .eq. 0
       end if
    case ('simulationtime', 'nsamples')
       interval = this%recon_time_interval
       if (interval .gt. 0.0_rp) then
          t_rel = time_out%t - time%start_time
          tol = 0.1_rp * abs(time%dt)
          recon_should_output = abs(t_rel - nint(t_rel / interval) * interval) &
               .le. tol
       end if
    case ('never')
       recon_should_output = .false.
    case default
       recon_should_output = .false.
    end select
  end function recon_should_output

  integer function precision_from_string(str, name) result(precision)
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: name

    select case (trim(lower_string(str)))
    case ('sp', 'single')
       precision = sp
    case ('dp', 'double')
       precision = dp
    case default
       call neko_error('Invalid ' // trim(name) // '. Expected ''sp'', ' // &
            '''dp'', ''single'', or ''double''.')
    end select
  end function precision_from_string

  pure function lower_string(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: out
    integer :: i, c

    do i = 1, len(str)
       c = iachar(str(i:i))
       if (c >= iachar('A') .and. c <= iachar('Z')) then
          out(i:i) = achar(c + 32)
       else
          out(i:i) = str(i:i)
       end if
    end do
  end function lower_string


  function csv_file_count_lines(file_in) result(n)
    class(file_t), intent(in) :: file_in
    integer :: n
    integer :: ierr, file_unit

    open(file = trim(file_in%get_fname()), status = 'old', &
         newunit = file_unit, iostat = ierr)
    if (ierr .ne. 0) then
       call neko_error("Error opening " // trim(file_in%get_fname()))
    end if
    rewind(file_unit)

    n = 0
    do
       read(file_unit, *, iostat = ierr)
       if (ierr .ne. 0) exit
       n = n + 1
    end do

    close(unit = file_unit)
  end function csv_file_count_lines


  subroutine find_bracket_time(i0, i1, time_coefs, tq)
    integer, intent(out)       :: i0, i1
    type(matrix_t), intent(in) :: time_coefs
    real(kind=rp), intent(in)  :: tq
    integer :: lo, hi, mid
    real(kind=rp) :: tm

    lo = 1
    hi = time_coefs%get_nrows()

    do while (hi - lo > 1)
       mid = (lo + hi) / 2
       tm = time_coefs%x(mid, 1)
       if (tm <= tq) then
          lo = mid
       else
          hi = mid
       end if
    end do

    i0 = lo
    i1 = lo + 1
  end subroutine find_bracket_time


  subroutine interpolate_time_coeffs_vec(a_out, time_coefs, tq)
    type(vector_t), intent(inout) :: a_out
    type(matrix_t), intent(in)    :: time_coefs
    real(kind=rp), intent(in)     :: tq
    integer :: nrows, ncols, j, i0, i1
    real(kind=rp) :: t0, t1, w

    nrows = time_coefs%get_nrows()
    ncols = time_coefs%get_ncols()

    if (ncols < 2) call neko_error("time_coefs must have (t + coeffs) columns")
    if (size(a_out%x) /= ncols-1) call neko_error("a_out wrong size")

    if (tq <= time_coefs%x(1,1)) then
       do j = 1, ncols-1
          a_out%x(j) = time_coefs%x(1, j+1)
       end do
       return
    end if

    if (tq >= time_coefs%x(nrows,1)) then
       do j = 1, ncols-1
          a_out%x(j) = time_coefs%x(nrows, j+1)
       end do
       return
    end if

    call find_bracket_time(i0, i1, time_coefs, tq)
    t0 = time_coefs%x(i0,1)
    t1 = time_coefs%x(i1,1)

    if (abs(t1 - t0) < tiny(1.0_rp)) then
       w = 0.0_rp
    else
       w = (tq - t0) / (t1 - t0)
    end if

    do j = 1, ncols-1
       a_out%x(j) = (1.0_rp - w) * time_coefs%x(i0, j+1) + &
            w * time_coefs%x(i1, j+1)
    end do
  end subroutine interpolate_time_coeffs_vec


  subroutine reconstruct_from_coeffs(this, neko_case, a)
    class(URANS_POD_state_recover_t), intent(inout) :: this
    class(case_t), intent(inout)              :: neko_case
    type(vector_t), intent(in)                :: a
    integer :: j
    type(field_t), pointer :: u, v, w, s

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    if (this%include_scalar) s => this%s
    if (this%include_nut) then
       this%nut => neko_registry%get_field('adjoint_nut')
    end if
    if (this%include_gammat) then
       this%gammat => neko_registry%get_field('adjoint_gammat')
    end if

    call field_rzero(u)
    call field_rzero(v)
    call field_rzero(w)
    if (this%include_scalar) call field_rzero(s)
    if (this%include_nut) call field_rzero(this%nut)
    if (this%include_gammat) call field_rzero(this%gammat)

    do j = 1, this%n_modes
       call field_add2s2(u, this%u_modes(j), a%x(j))
       call field_add2s2(v, this%v_modes(j), a%x(j))
       call field_add2s2(w, this%w_modes(j), a%x(j))
       if (this%include_scalar) then
          call field_add2s2(s, this%s_modes(j), a%x(j))
       end if
       if (this%include_nut) then
          call field_add2s2(this%nut, this%nut_modes(j), a%x(j))
       end if
       if (this%include_gammat) then
          call field_add2s2(this%gammat, this%gammat_modes(j), a%x(j))
       end if
    end do
  end subroutine reconstruct_from_coeffs

end module simulation_URANS_POD_state_recover
