!> Dense coefficient-space nonlinear Galerkin model.
!!
!! Reduced vectors and matrices deliberately use Neko's vector_t and matrix_t
!! containers.  This keeps the reduced state on the same host/device-aware
!! abstraction used by the rest of Neko-TOP; only the small dense contractions
!! themselves access the container storage.
module rom_model
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t
  use scratch_registry, only: neko_scratch_registry
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  type, public :: rom_model_t
     integer :: n_state = 0
     type(matrix_t) :: mass
     type(vector_t) :: constant
     type(matrix_t) :: linear
     type(matrix_t) :: implicit_linear
     ! The third-order tensor is represented as one matrix per residual row.
     type(matrix_t), allocatable :: quadratic(:)
   contains
     procedure, public :: init => rom_model_init
     procedure, public :: free => rom_model_free
     procedure, public :: rhs => rom_model_rhs
     procedure, public :: jacobian_state => rom_model_jacobian_state
     procedure, public :: forward_euler => rom_model_forward_euler
     procedure, public :: forward_bdf1 => rom_model_forward_bdf1
     procedure, public :: discrete_adjoint => rom_model_discrete_adjoint
     procedure, public :: discrete_adjoint_bdf1 => rom_model_discrete_adjoint_bdf1
     procedure, public :: init_linearized_adjoint => rom_model_init_linearized_adjoint
     procedure, public :: discrete_adjoint_bdf1_timevarying => &
          rom_model_discrete_adjoint_bdf1_timevarying
     procedure, public :: solve_mass => rom_model_solve_mass
     procedure, public :: solve_mass_transpose => rom_model_solve_mass_transpose
     procedure, public :: solve_step_transpose => rom_model_solve_step_transpose
  end type rom_model_t

contains

  !> Initialise the reduced operators and copy their coefficient data.
  subroutine rom_model_init(this, mass, constant, linear, quadratic, implicit_linear)
    class(rom_model_t), intent(inout) :: this
    type(matrix_t), intent(in) :: mass, linear
    type(vector_t), intent(in) :: constant
    type(matrix_t), intent(in) :: quadratic(:)
    type(matrix_t), intent(in), optional :: implicit_linear
    integer :: i

    call this%free()
    this%n_state = constant%size()
    if (mass%get_nrows() /= this%n_state .or. mass%get_ncols() /= this%n_state) &
         error stop 'ROM mass has the wrong shape'
    if (linear%get_nrows() /= this%n_state .or. linear%get_ncols() /= this%n_state) &
         error stop 'ROM linear operator has the wrong shape'
    if (size(quadratic) /= this%n_state) error stop 'ROM quadratic operator has the wrong shape'
    do i = 1, this%n_state
       if (quadratic(i)%get_nrows() /= this%n_state .or. &
            quadratic(i)%get_ncols() /= this%n_state) then
          error stop 'ROM quadratic operator has the wrong shape'
       end if
    end do

    call this%mass%init(this%n_state, this%n_state)
    call this%linear%init(this%n_state, this%n_state)
    call this%implicit_linear%init(this%n_state, this%n_state)
    call this%constant%init(this%n_state)
    ! Reduced dense algebra is intentionally host-owned.  Copy the host
    ! mirrors explicitly; container assignment would copy device mirrors on
    ! accelerator builds even though assembly returns scalar host entries.
    this%mass%x = mass%x
    this%linear%x = linear%x
    this%constant%x = constant%x
    this%implicit_linear%x = 0.0_rp
    if (present(implicit_linear)) then
       if (implicit_linear%get_nrows() /= this%n_state .or. &
            implicit_linear%get_ncols() /= this%n_state) then
          error stop 'ROM implicit linear operator has the wrong shape'
       end if
       this%implicit_linear%x = implicit_linear%x
    end if

    allocate(this%quadratic(this%n_state))
    do i = 1, this%n_state
       call this%quadratic(i)%init(this%n_state, this%n_state)
       this%quadratic(i)%x = quadratic(i)%x
    end do
  end subroutine rom_model_init

  !> Initialise the coefficient containers used by a time-varying linearised
  !! adjoint.  The nonlinear forward operators are not copied into this
  !! model; each timestep supplies its own projected Jacobian.
  subroutine rom_model_init_linearized_adjoint(this, mass, implicit_linear)
    class(rom_model_t), intent(inout) :: this
    type(matrix_t), intent(in) :: mass, implicit_linear
    type(vector_t) :: zero_constant
    type(matrix_t), allocatable :: zero_quadratic(:)
    integer :: n, i

    n = mass%get_nrows()
    if (n < 1 .or. mass%get_ncols() /= n .or. &
         implicit_linear%get_nrows() /= n .or. implicit_linear%get_ncols() /= n) &
         error stop 'ROM linearized adjoint operator has the wrong shape'
    call zero_constant%init(n)
    zero_constant%x = 0.0_rp
    allocate(zero_quadratic(n))
    do i = 1, n
       call zero_quadratic(i)%init(n, n)
       zero_quadratic(i)%x = 0.0_rp
    end do
    call this%init(mass, zero_constant, implicit_linear, zero_quadratic, &
         implicit_linear)
    call zero_constant%free()
    do i = 1, n
       call zero_quadratic(i)%free()
    end do
    deallocate(zero_quadratic)
  end subroutine rom_model_init_linearized_adjoint

  !> Release all reduced operators owned by the model.
  subroutine rom_model_free(this)
    class(rom_model_t), intent(inout) :: this
    integer :: i
    call this%mass%free()
    call this%constant%free()
    call this%linear%free()
    call this%implicit_linear%free()
    if (allocated(this%quadratic)) then
       do i = 1, size(this%quadratic)
          call this%quadratic(i)%free()
       end do
       deallocate(this%quadratic)
    end if
    this%n_state = 0
  end subroutine rom_model_free

  !> Evaluate the explicit reduced nonlinear residual.
  subroutine rom_model_rhs(this, state, rhs)
    class(rom_model_t), intent(in) :: this
    type(vector_t), intent(in) :: state
    type(vector_t), intent(inout) :: rhs
    type(vector_t), pointer :: residual
    integer :: residual_index
    integer :: i

    if (state%size() /= this%n_state) error stop 'ROM rhs state size mismatch'
    if (rhs%size() /= this%n_state) call rhs%init(this%n_state)
    call neko_scratch_registry%request_vector(residual, residual_index, &
         this%n_state, .true.)
    residual%x = this%constant%x + matmul(this%linear%x, state%x)
    do i = 1, this%n_state
       residual%x(i) = residual%x(i) + dot_product(state%x, &
            matmul(this%quadratic(i)%x, state%x))
    end do
    call rom_solve(this%mass, residual, rhs)
    call neko_scratch_registry%relinquish_vector(residual_index)
    nullify(residual)
  end subroutine rom_model_rhs

  !> Assemble the Jacobian of the explicit residual at one reduced state.
  subroutine rom_model_jacobian_state(this, state, jacobian)
    class(rom_model_t), intent(in) :: this
    type(vector_t), intent(in) :: state
    type(matrix_t), intent(inout) :: jacobian
    type(matrix_t), pointer :: residual_jac
    integer :: residual_jac_index
    integer :: i, j, k

    if (state%size() /= this%n_state) error stop 'ROM Jacobian state size mismatch'
    if (any(jacobian%get_dims() /= [this%n_state, this%n_state])) &
         call jacobian%init(this%n_state, this%n_state)
    call neko_scratch_registry%request_matrix(residual_jac, residual_jac_index, &
         this%n_state, this%n_state, .true.)
    residual_jac%x = this%linear%x
    do i = 1, this%n_state
       do j = 1, this%n_state
          do k = 1, this%n_state
             residual_jac%x(i,j) = residual_jac%x(i,j) + &
                  (this%quadratic(i)%x(j,k) + this%quadratic(i)%x(k,j)) * state%x(k)
          end do
       end do
    end do
    call rom_solve_matrix(this%mass, residual_jac, jacobian)
    call neko_scratch_registry%relinquish_matrix(residual_jac_index)
    nullify(residual_jac)
  end subroutine rom_model_jacobian_state

  !> Integrate the reduced model with explicit forward Euler steps.
  subroutine rom_model_forward_euler(this, initial_state, dt, states)
    class(rom_model_t), intent(in) :: this
    type(vector_t), intent(in) :: initial_state
    real(kind=rp), intent(in) :: dt
    type(matrix_t), intent(inout) :: states
    type(vector_t), pointer :: state, rhs
    integer :: state_index, rhs_index
    integer :: n_steps, step

    n_steps = states%get_ncols() - 1
    if (states%get_nrows() /= this%n_state .or. initial_state%size() /= this%n_state) &
         error stop 'ROM trajectory size mismatch'
    call neko_scratch_registry%request_vector(state, state_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(rhs, rhs_index, this%n_state, .true.)
    states%x(:,1) = initial_state%x
    do step = 1, n_steps
       state%x = states%x(:,step)
       call this%rhs(state, rhs)
       states%x(:,step+1) = states%x(:,step) + dt * rhs%x
    end do
    call neko_scratch_registry%relinquish_vector([state_index, rhs_index])
    nullify(state); nullify(rhs)
  end subroutine rom_model_forward_euler

  !> First-order Neko Pn/Pn map for a fixed timestep.
  !! The viscous and design-dependent Brinkman operators are implicit in the
  !! reduced BDF1 Helmholtz solve.  The solved matrix is
  !! `mass - dt*implicit_linear`.
  !> Integrate the reduced model with the supported BDF1 scheme.
  subroutine rom_model_forward_bdf1(this, initial_state, dt, states)
    class(rom_model_t), intent(in) :: this
    type(vector_t), intent(in) :: initial_state
    real(kind=rp), intent(in) :: dt
    type(matrix_t), intent(inout) :: states
    type(vector_t), pointer :: state, residual, rhs
    type(matrix_t), pointer :: step_matrix
    integer :: state_index, residual_index, rhs_index, step_matrix_index
    integer :: n_steps, step

    n_steps = states%get_ncols() - 1
    if (states%get_nrows() /= this%n_state .or. initial_state%size() /= this%n_state) &
         error stop 'ROM BDF1 trajectory size mismatch'
    call neko_scratch_registry%request_vector(state, state_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(residual, residual_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(rhs, rhs_index, this%n_state, .true.)
    call neko_scratch_registry%request_matrix(step_matrix, step_matrix_index, &
         this%n_state, this%n_state, .true.)
    step_matrix%x = this%mass%x - dt * this%implicit_linear%x
    if (any(.not. ieee_is_finite(step_matrix%x))) then
       write(*,'(A)') 'ROM BDF1 failure: non-finite implicit step matrix'
       error stop 'non-finite ROM BDF1 step matrix'
    end if
    states%x(:,1) = initial_state%x
    if (any(.not. ieee_is_finite(states%x(:,1)))) then
       write(*,'(A)') 'ROM BDF1 failure: non-finite initial reduced state'
       error stop 'non-finite ROM BDF1 initial state'
    end if
    do step = 1, n_steps
       state%x = states%x(:,step)
       call rom_model_explicit_residual(this, state, residual)
       if (any(.not. ieee_is_finite(residual%x))) then
          write(*,'(A,I0,A,ES14.6)') 'ROM BDF1 failure: non-finite residual at step ', &
               step, ', max_abs_state=', maxval(abs(state%x))
          error stop 'non-finite ROM BDF1 residual'
       end if
       rhs%x = matmul(this%mass%x, state%x) + dt * residual%x
       if (any(.not. ieee_is_finite(rhs%x))) then
          write(*,'(A,I0)') 'ROM BDF1 failure: non-finite right-hand side at step ', step
          error stop 'non-finite ROM BDF1 right-hand side'
       end if
       call rom_solve(step_matrix, rhs, state)
       if (any(.not. ieee_is_finite(state%x))) then
          write(*,'(A,I0)') 'ROM BDF1 failure: non-finite state at step ', step
          error stop 'non-finite ROM BDF1 state'
       end if
       states%x(:,step+1) = state%x
    end do
    call neko_scratch_registry%relinquish_vector([state_index, residual_index, rhs_index])
    call neko_scratch_registry%relinquish_matrix(step_matrix_index)
    nullify(state); nullify(residual); nullify(rhs); nullify(step_matrix)
  end subroutine rom_model_forward_bdf1

  !> Discrete adjoint of the explicit Euler map used above.
  !> Integrate the discrete adjoint of the explicit reduced model backward.
  subroutine rom_model_discrete_adjoint(this, states, dt, &
       running_state_gradient, terminal_state_gradient, adjoint)
    class(rom_model_t), intent(in) :: this
    type(matrix_t), intent(in) :: states, running_state_gradient
    real(kind=rp), intent(in) :: dt
    type(vector_t), intent(in) :: terminal_state_gradient
    type(matrix_t), intent(inout) :: adjoint
    type(vector_t), pointer :: state, running, terminal, next_adjoint, current_adjoint
    type(matrix_t), pointer :: jacobian, transition
    integer :: state_index, running_index, terminal_index, next_index, current_index
    integer :: jacobian_index, transition_index
    integer :: n_steps, step, i

    n_steps = states%get_ncols() - 1
    if (states%get_nrows() /= this%n_state .or. &
         running_state_gradient%get_nrows() /= this%n_state .or. &
         running_state_gradient%get_ncols() /= n_steps .or. &
         terminal_state_gradient%size() /= this%n_state) then
       error stop 'ROM adjoint state size mismatch'
    end if
    if (any(adjoint%get_dims() /= [this%n_state, n_steps + 1])) &
         call adjoint%init(this%n_state, n_steps + 1)
    call neko_scratch_registry%request_vector(state, state_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(running, running_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(next_adjoint, next_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(current_adjoint, current_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(terminal, terminal_index, this%n_state, .true.)
    terminal%x = terminal_state_gradient%x
    adjoint%x(:,n_steps+1) = terminal%x
    if (any(.not. ieee_is_finite(adjoint%x(:,n_steps+1)))) then
       error stop 'non-finite ROM adjoint terminal state'
    end if
    call neko_scratch_registry%request_matrix(jacobian, jacobian_index, &
         this%n_state, this%n_state, .true.)
    call neko_scratch_registry%request_matrix(transition, transition_index, &
         this%n_state, this%n_state, .true.)
    do step = n_steps, 1, -1
       state%x = states%x(:,step)
       next_adjoint%x = adjoint%x(:,step+1)
       running%x = running_state_gradient%x(:,step)
       call this%jacobian_state(state, jacobian)
       transition%x = dt * jacobian%x
       do i = 1, this%n_state
          transition%x(i,i) = transition%x(i,i) + 1.0_rp
       end do
       current_adjoint%x = matmul(transpose(transition%x), next_adjoint%x) + dt * running%x
       adjoint%x(:,step) = current_adjoint%x
    end do
    call neko_scratch_registry%relinquish_vector([state_index, running_index, &
         terminal_index, next_index, current_index])
    call neko_scratch_registry%relinquish_matrix([jacobian_index, transition_index])
    nullify(state); nullify(running); nullify(terminal); nullify(next_adjoint)
    nullify(current_adjoint); nullify(jacobian); nullify(transition)
  end subroutine rom_model_discrete_adjoint

  !> Discrete adjoint of the first-order Neko Pn/Pn BDF1 map.
  !> Integrate the BDF1 discrete adjoint backward through the trajectory.
  subroutine rom_model_discrete_adjoint_bdf1(this, states, dt, &
       running_state_gradient, terminal_state_gradient, adjoint)
    class(rom_model_t), intent(in) :: this
    type(matrix_t), intent(in) :: states, running_state_gradient
    real(kind=rp), intent(in) :: dt
    type(vector_t), intent(in) :: terminal_state_gradient
    type(matrix_t), intent(inout) :: adjoint
    type(vector_t), pointer :: state, running, terminal, next_adjoint, current_adjoint
    type(matrix_t), pointer :: jacobian, map_jacobian, step_matrix, rhs_matrix
    integer :: state_index, running_index, terminal_index, next_index, current_index
    integer :: jacobian_index, map_index, step_index, rhs_index
    integer :: n_steps, step

    n_steps = states%get_ncols() - 1
    if (states%get_nrows() /= this%n_state .or. &
         running_state_gradient%get_nrows() /= this%n_state .or. &
         running_state_gradient%get_ncols() /= n_steps .or. &
         terminal_state_gradient%size() /= this%n_state) then
       error stop 'ROM BDF1 adjoint state size mismatch'
    end if
    if (any(adjoint%get_dims() /= [this%n_state, n_steps + 1])) &
         call adjoint%init(this%n_state, n_steps + 1)
    call neko_scratch_registry%request_vector(state, state_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(running, running_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(terminal, terminal_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(next_adjoint, next_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(current_adjoint, current_index, this%n_state, .true.)
    terminal%x = terminal_state_gradient%x
    adjoint%x(:,n_steps+1) = terminal%x
    if (any(.not. ieee_is_finite(adjoint%x(:,n_steps+1)))) then
       error stop 'non-finite ROM BDF1 adjoint terminal state'
    end if
    call neko_scratch_registry%request_matrix(jacobian, jacobian_index, &
         this%n_state, this%n_state, .true.)
    call neko_scratch_registry%request_matrix(map_jacobian, map_index, &
         this%n_state, this%n_state, .true.)
    call neko_scratch_registry%request_matrix(step_matrix, step_index, &
         this%n_state, this%n_state, .true.)
    call neko_scratch_registry%request_matrix(rhs_matrix, rhs_index, &
         this%n_state, this%n_state, .true.)
    step_matrix%x = this%mass%x - dt * this%implicit_linear%x
    do step = n_steps, 1, -1
       state%x = states%x(:,step)
       next_adjoint%x = adjoint%x(:,step+1)
       running%x = running_state_gradient%x(:,step)
       call this%jacobian_state(state, jacobian)
       rhs_matrix%x = this%mass%x + dt * matmul(this%mass%x, jacobian%x)
       call rom_solve_matrix(step_matrix, rhs_matrix, map_jacobian)
       current_adjoint%x = matmul(transpose(map_jacobian%x), next_adjoint%x) + dt * running%x
       if (any(.not. ieee_is_finite(current_adjoint%x))) then
          error stop 'non-finite ROM BDF1 adjoint state'
       end if
       adjoint%x(:,step) = current_adjoint%x
    end do
    call neko_scratch_registry%relinquish_vector([state_index, running_index, &
         terminal_index, next_index, current_index])
    call neko_scratch_registry%relinquish_matrix([jacobian_index, map_index, &
         step_index, rhs_index])
    nullify(state); nullify(running); nullify(terminal); nullify(next_adjoint)
    nullify(current_adjoint); nullify(jacobian); nullify(map_jacobian)
    nullify(step_matrix); nullify(rhs_matrix)
  end subroutine rom_model_discrete_adjoint_bdf1

  !> Integrate an enriched adjoint with a time-dependent projected Jacobian.
  !!
  !! For the BDF1/IMEX map represented by this model, the supplied
  !! `time_jacobian(step)` is the explicit residual Jacobian projected onto
  !! the enriched basis at the stored forward state at `step`.  The backward
  !! map is therefore formed as
  !!
  !!   S^T lambda_n = (H + dt J_n)^T lambda_(n+1) + dt g_n,
  !!
  !! where `S = H - dt*implicit_linear`.  This keeps the adjoint timestep
  !! consistent with the forward BDF1 split while allowing the forward state
  !! and the adjoint basis to have different dimensions.
  subroutine rom_model_discrete_adjoint_bdf1_timevarying(this, dt, &
       running_state_gradient, terminal_state_gradient, time_jacobian, adjoint)
    class(rom_model_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt
    type(matrix_t), intent(in) :: running_state_gradient
    type(vector_t), intent(in) :: terminal_state_gradient
    type(matrix_t), intent(in) :: time_jacobian(:)
    type(matrix_t), intent(inout) :: adjoint
    type(vector_t), pointer :: running, next_adjoint, current_adjoint
    type(matrix_t), pointer :: step_matrix, rhs_matrix, map_jacobian
    integer :: running_index, next_index, current_index
    integer :: step_index, rhs_index, map_index
    integer :: n_steps, step

    n_steps = running_state_gradient%get_ncols()
    if (running_state_gradient%get_nrows() /= this%n_state .or. &
         terminal_state_gradient%size() /= this%n_state .or. &
         size(time_jacobian) /= n_steps) then
       error stop 'ROM time-varying adjoint state size mismatch'
    end if
    if (any(adjoint%get_dims() /= [this%n_state, n_steps + 1])) &
         call adjoint%init(this%n_state, n_steps + 1)

    call neko_scratch_registry%request_vector(running, running_index, &
         this%n_state, .true.)
    call neko_scratch_registry%request_vector(next_adjoint, next_index, &
         this%n_state, .true.)
    call neko_scratch_registry%request_vector(current_adjoint, current_index, &
         this%n_state, .true.)
    call neko_scratch_registry%request_matrix(step_matrix, step_index, &
         this%n_state, this%n_state, .true.)
    call neko_scratch_registry%request_matrix(rhs_matrix, rhs_index, &
         this%n_state, this%n_state, .true.)
    call neko_scratch_registry%request_matrix(map_jacobian, map_index, &
         this%n_state, this%n_state, .true.)

    step_matrix%x = this%mass%x - dt * this%implicit_linear%x
    adjoint%x(:,n_steps+1) = terminal_state_gradient%x
    if (any(.not. ieee_is_finite(adjoint%x(:,n_steps+1)))) then
       write(*,'(A)') 'ROM time-varying adjoint failure: non-finite terminal gradient'
       error stop 'non-finite ROM time-varying adjoint terminal state'
    end if
    do step = n_steps, 1, -1
       next_adjoint%x = adjoint%x(:,step+1)
       running%x = running_state_gradient%x(:,step)
       rhs_matrix%x = this%mass%x + dt * time_jacobian(step)%x
       if (any(.not. ieee_is_finite(rhs_matrix%x))) then
          write(*,'(A,I0)') 'ROM time-varying adjoint failure: non-finite map at step ', step
          error stop 'non-finite ROM time-varying adjoint map'
       end if
       call rom_solve_matrix(step_matrix, rhs_matrix, map_jacobian)
       current_adjoint%x = matmul(transpose(map_jacobian%x), next_adjoint%x) + &
            dt * running%x
       if (any(.not. ieee_is_finite(current_adjoint%x))) then
          write(*,'(A,I0)') 'ROM time-varying adjoint failure: non-finite state at step ', step
          error stop 'non-finite ROM time-varying adjoint state'
       end if
       adjoint%x(:,step) = current_adjoint%x
    end do

    call neko_scratch_registry%relinquish_vector([running_index, next_index, &
         current_index])
    call neko_scratch_registry%relinquish_matrix([step_index, rhs_index, map_index])
    nullify(running); nullify(next_adjoint); nullify(current_adjoint)
    nullify(step_matrix); nullify(rhs_matrix); nullify(map_jacobian)
  end subroutine rom_model_discrete_adjoint_bdf1_timevarying

  !> Solve the transpose reduced mass system.
  subroutine rom_model_solve_mass_transpose(this, rhs, solution)
    class(rom_model_t), intent(in) :: this
    type(vector_t), intent(in) :: rhs
    type(vector_t), intent(inout) :: solution
    call rom_solve_transpose(this%mass, rhs, solution)
  end subroutine rom_model_solve_mass_transpose

  !> Solve the reduced mass system.
  subroutine rom_model_solve_mass(this, rhs, solution)
    class(rom_model_t), intent(in) :: this
    type(vector_t), intent(in) :: rhs
    type(vector_t), intent(inout) :: solution
    call rom_solve(this%mass, rhs, solution)
  end subroutine rom_model_solve_mass

  !> Solve one transpose implicit BDF1 step matrix.
  subroutine rom_model_solve_step_transpose(this, dt, rhs, solution)
    class(rom_model_t), intent(in) :: this
    real(kind=rp), intent(in) :: dt
    type(vector_t), intent(in) :: rhs
    type(vector_t), intent(inout) :: solution
    type(matrix_t), pointer :: step_matrix
    integer :: step_matrix_index
    call neko_scratch_registry%request_matrix(step_matrix, step_matrix_index, &
         this%n_state, this%n_state, .true.)
    step_matrix%x = this%mass%x - dt * this%implicit_linear%x
    call rom_solve_transpose(step_matrix, rhs, solution)
    call neko_scratch_registry%relinquish_matrix(step_matrix_index)
    nullify(step_matrix)
  end subroutine rom_model_solve_step_transpose

  !> Evaluate only the explicit portion of the reduced residual.
  subroutine rom_model_explicit_residual(this, state, residual)
    class(rom_model_t), intent(in) :: this
    type(vector_t), intent(in) :: state
    type(vector_t), intent(inout) :: residual
    integer :: i
    if (residual%size() /= this%n_state) call residual%init(this%n_state)
    residual%x = this%constant%x + matmul(this%linear%x, state%x)
    do i = 1, this%n_state
       residual%x(i) = residual%x(i) + dot_product(state%x, &
            matmul(this%quadratic(i)%x, state%x))
    end do
  end subroutine rom_model_explicit_residual

  !> Solve a dense reduced linear system with pivoted elimination.
  subroutine rom_solve(matrix, rhs, solution)
    type(matrix_t), intent(in) :: matrix
    type(vector_t), intent(in) :: rhs
    type(vector_t), intent(inout) :: solution
    type(matrix_t), pointer :: work
    type(vector_t), pointer :: b
    integer :: work_index, b_index
    integer :: n, i, j, pivot
    real(kind=rp) :: factor, largest, tmp

    n = rhs%size()
    if (matrix%get_nrows() /= n .or. matrix%get_ncols() /= n) error stop 'ROM solve shape mismatch'
    if (solution%size() /= n) call solution%init(n)
    call neko_scratch_registry%request_matrix(work, work_index, n, n, .false.)
    call neko_scratch_registry%request_vector(b, b_index, n, .false.)
    work%x = matrix%x; b%x = rhs%x
    do i = 1, n
       pivot = i
       largest = abs(work%x(i, i))
       do j = i + 1, n
          if (abs(work%x(j, i)) > largest) then
             largest = abs(work%x(j, i)); pivot = j
          end if
       end do
       if (largest <= epsilon(1.0_rp)) error stop 'singular ROM mass matrix'
       if (pivot /= i) then
          call swap_rows(work, i, pivot)
          tmp = b%x(i); b%x(i) = b%x(pivot); b%x(pivot) = tmp
       end if
       do j = i + 1, n
          factor = work%x(j, i) / work%x(i, i)
          work%x(j, i:n) = work%x(j, i:n) - factor * work%x(i, i:n)
          b%x(j) = b%x(j) - factor * b%x(i)
       end do
    end do
    solution%x = b%x
    do i = n, 1, -1
       if (i < n) solution%x(i) = solution%x(i) - &
            sum(work%x(i, i + 1:n) * solution%x(i + 1:n))
       solution%x(i) = solution%x(i) / work%x(i, i)
    end do
    call neko_scratch_registry%relinquish_matrix(work_index)
    call neko_scratch_registry%relinquish_vector(b_index)
    nullify(work); nullify(b)
  end subroutine rom_solve

  !> Solve a dense reduced transpose linear system.
  subroutine rom_solve_transpose(matrix, rhs, solution)
    type(matrix_t), intent(in) :: matrix
    type(vector_t), intent(in) :: rhs
    type(vector_t), intent(inout) :: solution
    type(matrix_t), pointer :: transposed
    integer :: transposed_index
    call neko_scratch_registry%request_matrix(transposed, transposed_index, &
         matrix%get_ncols(), matrix%get_nrows(), .false.)
    transposed%x = transpose(matrix%x)
    call rom_solve(transposed, rhs, solution)
    call neko_scratch_registry%relinquish_matrix(transposed_index)
    nullify(transposed)
  end subroutine rom_solve_transpose

  !> Solve a dense reduced matrix system one right-hand side at a time.
  subroutine rom_solve_matrix(matrix, rhs, solution)
    type(matrix_t), intent(in) :: matrix, rhs
    type(matrix_t), intent(inout) :: solution
    integer :: col, b_index, x_index
    type(vector_t), pointer :: b, x
    if (any(solution%get_dims() /= [rhs%get_nrows(), rhs%get_ncols()])) &
         call solution%init(rhs%get_nrows(), rhs%get_ncols())
    call neko_scratch_registry%request_vector(b, b_index, rhs%get_nrows(), .false.)
    call neko_scratch_registry%request_vector(x, x_index, rhs%get_nrows(), .false.)
    do col = 1, rhs%get_ncols()
       b%x = rhs%x(:,col)
       call rom_solve(matrix, b, x)
       solution%x(:,col) = x%x
    end do
    call neko_scratch_registry%relinquish_vector([b_index, x_index])
    nullify(b); nullify(x)
  end subroutine rom_solve_matrix

  !> Swap two rows of a dense reduced matrix during pivoting.
  subroutine swap_rows(matrix, i, j)
    type(matrix_t), intent(inout) :: matrix
    integer, intent(in) :: i, j
    real(kind=rp) :: tmp(matrix%get_ncols())
    tmp = matrix%x(i,:)
    matrix%x(i,:) = matrix%x(j,:)
    matrix%x(j,:) = tmp
  end subroutine swap_rows

end module rom_model
