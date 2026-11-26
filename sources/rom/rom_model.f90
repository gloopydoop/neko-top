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
  use comm, only: pe_rank
  use utils, only: mkdir, neko_error
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  !> Optional pre-integration export metadata for one candidate BDF1 system.
  !! `directory` is the final, candidate-specific output directory.
  type, public :: rom_bdf1_export_context_t
     character(len=1024) :: directory = ''
     character(len=128) :: design_id = ''
     character(len=32) :: basis_forward = ''
     character(len=32) :: basis_adjoint = ''
     integer :: outer_iteration = -1
     integer :: inner_iteration = -1
     integer :: evaluation_index = -1
     integer :: trial_basis_rank = -1
     integer :: test_basis_rank = -1
     integer :: forward_basis_rank = -1
     integer :: adjoint_basis_rank = -1
     integer :: enriched_basis_rank = -1
     real(kind=rp) :: time_start = 0.0_rp
  end type rom_bdf1_export_context_t

  type, public :: rom_model_t
     integer :: n_state = 0
     type(matrix_t) :: mass
     type(vector_t) :: constant
     type(matrix_t) :: linear
     type(matrix_t) :: implicit_linear
     type(matrix_t) :: bdf1_lu
     integer, allocatable :: bdf1_pivots(:)
     real(kind=rp) :: bdf1_factor_dt = -1.0_rp
     logical :: bdf1_factor_valid = .false.
     ! The third-order tensor is represented as one matrix per residual row.
     type(matrix_t), allocatable :: quadratic(:)
   contains
     procedure, public :: init => rom_model_init
     procedure, public :: free => rom_model_free
     procedure, public :: rhs => rom_model_rhs
     procedure, public :: jacobian_state => rom_model_jacobian_state
     procedure, public :: jacobian_residual_state => rom_model_jacobian_residual_state
     procedure, public :: forward_euler => rom_model_forward_euler
     procedure, public :: forward_bdf1 => rom_model_forward_bdf1
     procedure, public :: discrete_adjoint => rom_model_discrete_adjoint
     procedure, public :: discrete_adjoint_bdf1 => rom_model_discrete_adjoint_bdf1
     procedure, public :: init_linearized_adjoint => rom_model_init_linearized_adjoint
     procedure, public :: discrete_adjoint_bdf1_timevarying => &
          rom_model_discrete_adjoint_bdf1_timevarying
     procedure, public :: discrete_adjoint_bdf1_tensor => &
          rom_model_discrete_adjoint_bdf1_tensor
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
    call this%bdf1_lu%free()
    if (allocated(this%bdf1_pivots)) deallocate(this%bdf1_pivots)
    this%bdf1_factor_dt = -1.0_rp
    this%bdf1_factor_valid = .false.
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
    if (state%size() /= this%n_state) error stop 'ROM Jacobian state size mismatch'
    if (any(jacobian%get_dims() /= [this%n_state, this%n_state])) &
         call jacobian%init(this%n_state, this%n_state)
    call neko_scratch_registry%request_matrix(residual_jac, residual_jac_index, &
         this%n_state, this%n_state, .true.)
    call this%jacobian_residual_state(state, residual_jac)
    call rom_solve_matrix(this%mass, residual_jac, jacobian)
    call neko_scratch_registry%relinquish_matrix(residual_jac_index)
    nullify(residual_jac)
  end subroutine rom_model_jacobian_state

  !> Assemble the Jacobian of the weak explicit residual, before mass inversion.
  subroutine rom_model_jacobian_residual_state(this, state, jacobian)
    class(rom_model_t), intent(in) :: this
    type(vector_t), intent(in) :: state
    type(matrix_t), intent(inout) :: jacobian
    integer :: i, j, k
    if (state%size() /= this%n_state) error stop 'ROM residual Jacobian state size mismatch'
    if (any(jacobian%get_dims() /= [this%n_state, this%n_state])) &
         call jacobian%init(this%n_state, this%n_state)
    jacobian%x = this%linear%x
    do i = 1, this%n_state
       do j = 1, this%n_state
          do k = 1, this%n_state
             jacobian%x(i,j) = jacobian%x(i,j) + &
                  (this%quadratic(i)%x(j,k) + this%quadratic(i)%x(k,j)) * state%x(k)
          end do
       end do
    end do
  end subroutine rom_model_jacobian_residual_state

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
  subroutine rom_model_forward_bdf1(this, initial_state, dt, states, export_context)
    class(rom_model_t), intent(inout) :: this
    type(vector_t), intent(in) :: initial_state
    real(kind=rp), intent(in) :: dt
    type(matrix_t), intent(inout) :: states
    type(rom_bdf1_export_context_t), intent(in), optional :: export_context
    type(vector_t), pointer :: state, residual, rhs
    integer :: state_index, residual_index, rhs_index
    integer :: n_steps, step

    n_steps = states%get_ncols() - 1
    if (states%get_nrows() /= this%n_state .or. initial_state%size() /= this%n_state) &
         error stop 'ROM BDF1 trajectory size mismatch'
    ! Persist the exact reduced map before factorisation or time stepping.  A
    ! later singular solve, non-finite residual, or growing state therefore
    ! leaves the candidate operators and identifiers available for diagnosis.
    if (present(export_context)) &
         call rom_model_export_bdf1_system(this, initial_state, dt, n_steps, export_context)
    call neko_scratch_registry%request_vector(state, state_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(residual, residual_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(rhs, rhs_index, this%n_state, .true.)
    call rom_model_prepare_bdf1_factor(this, dt)
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
       call rom_lu_solve(this%bdf1_lu, this%bdf1_pivots, rhs, state)
       if (any(.not. ieee_is_finite(state%x))) then
          write(*,'(A,I0)') 'ROM BDF1 failure: non-finite state at step ', step
          error stop 'non-finite ROM BDF1 state'
       end if
       states%x(:,step+1) = state%x
    end do
    call neko_scratch_registry%relinquish_vector([state_index, residual_index, rhs_index])
    nullify(state); nullify(residual); nullify(rhs)
  end subroutine rom_model_forward_bdf1

  !> Export one candidate reduced system using the standalone ROM reader's
  !! `metadata.json`, `H.csv`, `b.csv`, `linear_explicit.csv`,
  !! `linear_implicit.csv`, `Q.csv`, and `a0.csv` contract.
  subroutine rom_model_export_bdf1_system(this, initial_state, dt, n_steps, context)
    class(rom_model_t), intent(in) :: this
    type(vector_t), intent(in) :: initial_state
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n_steps
    type(rom_bdf1_export_context_t), intent(in) :: context
    integer :: unit, ios, i, j, k
    character(len=1024) :: filename

    if (pe_rank /= 0) return
    if (len_trim(context%directory) == 0) &
         call neko_error('ROM BDF1 export context has an empty directory')
    if (initial_state%size() /= this%n_state .or. n_steps < 0) &
         call neko_error('ROM BDF1 export context has inconsistent state dimensions')
    call mkdir(trim(context%directory))

    filename = trim(context%directory)//'/H.csv'
    call rom_model_write_matrix_csv(trim(filename), this%mass%x)
    filename = trim(context%directory)//'/b.csv'
    call rom_model_write_vector_csv(trim(filename), this%constant%x)
    filename = trim(context%directory)//'/linear_explicit.csv'
    call rom_model_write_matrix_csv(trim(filename), this%linear%x)
    filename = trim(context%directory)//'/linear_implicit.csv'
    call rom_model_write_matrix_csv(trim(filename), this%implicit_linear%x)
    filename = trim(context%directory)//'/a0.csv'
    call rom_model_write_vector_csv(trim(filename), initial_state%x)

    filename = trim(context%directory)//'/Q.csv'
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) call neko_error('could not open ROM BDF1 quadratic export')
    write(unit, '(A)') 'i,j,k,value'
    do i = 1, this%n_state
       do j = 1, this%n_state
          do k = 1, this%n_state
             if (this%quadratic(i)%x(j,k) == 0.0_rp) cycle
             write(unit, '(I0,",",I0,",",I0,",",ES27.16E3)') &
                  i, j, k, this%quadratic(i)%x(j,k)
          end do
       end do
    end do
    close(unit)

    filename = trim(context%directory)//'/metadata.json'
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) call neko_error('could not open ROM BDF1 export metadata')
    write(unit, '(A)') '{'
    write(unit, '(A,I0,A)') '  "n_state": ', this%n_state, ','
    write(unit, '(A,I0,A)') '  "trial_basis_rank": ', context%trial_basis_rank, ','
    write(unit, '(A,I0,A)') '  "test_basis_rank": ', context%test_basis_rank, ','
    write(unit, '(A,I0,A)') '  "forward_basis_rank": ', context%forward_basis_rank, ','
    write(unit, '(A,I0,A)') '  "adjoint_basis_rank": ', context%adjoint_basis_rank, ','
    write(unit, '(A,I0,A)') '  "enriched_basis_rank": ', context%enriched_basis_rank, ','
    write(unit, '(A,I0,A)') '  "outer_iteration": ', context%outer_iteration, ','
    write(unit, '(A,I0,A)') '  "inner_iteration": ', context%inner_iteration, ','
    write(unit, '(A,I0,A)') '  "evaluation_index": ', context%evaluation_index, ','
    write(unit, '(A,A,A)') '  "design_id": "', trim(context%design_id), '",'
    write(unit, '(A,A,A)') '  "basis_forward": "', trim(context%basis_forward), '",'
    write(unit, '(A,A,A)') '  "basis_adjoint": "', trim(context%basis_adjoint), '",'
    write(unit, '(A)') '  "export_stage": "pre_forward_bdf1",'
    write(unit, '(A,ES27.16E3,A)') '  "time_start": ', context%time_start, ','
    write(unit, '(A,ES27.16E3,A)') '  "dt": ', dt, ','
    write(unit, '(A,I0,A)') '  "n_steps": ', n_steps, ','
    write(unit, '(A,ES27.16E3,A)') '  "time_end": ', &
         context%time_start + real(n_steps, rp) * dt, ','
    write(unit, '(A)') '  "basis_layout": "reduced coefficient space",'
    write(unit, '(A)') '  "q_indexing": "Q[i,j,k] = sum_jk coefficient for residual i",'
    write(unit, '(A)') '  "b_definition": "constant residual vector",'
    write(unit, '(A)') '  "linear_explicit_definition": "linear residual matrix",'
    write(unit, '(A)') '  "linear_implicit_definition": "implicit BDF1 residual matrix",'
    write(unit, '(A)') '  "files": {'
    write(unit, '(A)') '    "H": "H.csv", "b": "b.csv",'
    write(unit, '(A)') '    "linear_explicit": "linear_explicit.csv",'
    write(unit, '(A)') '    "linear_implicit": "linear_implicit.csv",'
    write(unit, '(A)') '    "Q": "Q.csv", "a0": "a0.csv"'
    write(unit, '(A)') '  }'
    write(unit, '(A)') '}'
    flush(unit)
    close(unit)
  end subroutine rom_model_export_bdf1_system

  !> Write dense reduced rows using an explicit three-digit exponent.
  subroutine rom_model_write_matrix_csv(filename, values)
    character(len=*), intent(in) :: filename
    real(kind=rp), intent(in) :: values(:, :)
    integer :: unit, ios, i, j
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) call neko_error('could not open ROM BDF1 matrix export')
    do i = 1, size(values,1)
       do j = 1, size(values,2)
          if (j > 1) write(unit, '(A)', advance='no') ','
          write(unit, '(ES27.16E3)', advance='no') values(i,j)
       end do
       write(unit, '(A)') ''
    end do
    close(unit)
  end subroutine rom_model_write_matrix_csv

  !> Write a reduced vector in the reader-compatible one-column CSV form.
  subroutine rom_model_write_vector_csv(filename, values)
    character(len=*), intent(in) :: filename
    real(kind=rp), intent(in) :: values(:)
    integer :: unit, ios, i
    open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
    if (ios /= 0) call neko_error('could not open ROM BDF1 vector export')
    do i = 1, size(values)
       write(unit, '(ES27.16E3)') values(i)
    end do
    close(unit)
  end subroutine rom_model_write_vector_csv

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
       running_state_gradient, terminal_state_gradient, adjoint, lambda_history)
    class(rom_model_t), intent(inout) :: this
    type(matrix_t), intent(in) :: states, running_state_gradient
    real(kind=rp), intent(in) :: dt
    type(vector_t), intent(in) :: terminal_state_gradient
    type(matrix_t), intent(inout) :: adjoint
    type(matrix_t), intent(inout), optional :: lambda_history
    type(vector_t), pointer :: state, running, terminal, next_adjoint, current_adjoint
    type(vector_t), pointer :: lambda, mass_adjoint
    type(matrix_t), pointer :: residual_jacobian
    integer :: state_index, running_index, terminal_index, next_index, current_index
    integer :: lambda_index, mass_adjoint_index, jacobian_index
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
    if (present(lambda_history)) then
       if (any(lambda_history%get_dims() /= [this%n_state, n_steps])) &
            call lambda_history%init(this%n_state, n_steps)
    end if
    call neko_scratch_registry%request_vector(state, state_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(running, running_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(terminal, terminal_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(next_adjoint, next_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(current_adjoint, current_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(lambda, lambda_index, this%n_state, .true.)
    call neko_scratch_registry%request_vector(mass_adjoint, mass_adjoint_index, &
         this%n_state, .true.)
    terminal%x = terminal_state_gradient%x
    adjoint%x(:,n_steps+1) = terminal%x
    if (any(.not. ieee_is_finite(adjoint%x(:,n_steps+1)))) then
       error stop 'non-finite ROM BDF1 adjoint terminal state'
    end if
    call neko_scratch_registry%request_matrix(residual_jacobian, jacobian_index, &
         this%n_state, this%n_state, .true.)
    do step = n_steps, 1, -1
       state%x = states%x(:,step)
       next_adjoint%x = adjoint%x(:,step+1)
       running%x = running_state_gradient%x(:,step)
       call this%jacobian_residual_state(state, residual_jacobian)
       call this%solve_step_transpose(dt, next_adjoint, lambda)
       if (present(lambda_history)) lambda_history%x(:,step) = lambda%x
       mass_adjoint%x = matmul(transpose(this%mass%x), lambda%x)
       current_adjoint%x = mass_adjoint%x + dt * &
            matmul(transpose(residual_jacobian%x), lambda%x) + dt * running%x
       if (any(.not. ieee_is_finite(current_adjoint%x))) then
          error stop 'non-finite ROM BDF1 adjoint state'
       end if
       adjoint%x(:,step) = current_adjoint%x
    end do
    call neko_scratch_registry%relinquish_vector([state_index, running_index, &
         terminal_index, next_index, current_index])
    call neko_scratch_registry%relinquish_vector([lambda_index, mass_adjoint_index])
    call neko_scratch_registry%relinquish_matrix(jacobian_index)
    nullify(state); nullify(running); nullify(terminal); nullify(next_adjoint)
    nullify(current_adjoint); nullify(lambda); nullify(mass_adjoint)
    nullify(residual_jacobian)
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
    class(rom_model_t), intent(inout) :: this
    real(kind=rp), intent(in) :: dt
    type(matrix_t), intent(in) :: running_state_gradient
    type(vector_t), intent(in) :: terminal_state_gradient
    type(matrix_t), intent(in) :: time_jacobian(:)
    type(matrix_t), intent(inout) :: adjoint
    type(vector_t), pointer :: running, next_adjoint, current_adjoint, lambda
    integer :: running_index, next_index, current_index, lambda_index
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
    call neko_scratch_registry%request_vector(lambda, lambda_index, this%n_state, .true.)
    adjoint%x(:,n_steps+1) = terminal_state_gradient%x
    if (any(.not. ieee_is_finite(adjoint%x(:,n_steps+1)))) then
       write(*,'(A)') 'ROM time-varying adjoint failure: non-finite terminal gradient'
       error stop 'non-finite ROM time-varying adjoint terminal state'
    end if
    do step = n_steps, 1, -1
       next_adjoint%x = adjoint%x(:,step+1)
       running%x = running_state_gradient%x(:,step)
       call this%solve_step_transpose(dt, next_adjoint, lambda)
       current_adjoint%x = matmul(transpose(this%mass%x), lambda%x) + &
            dt * matmul(transpose(time_jacobian(step)%x), lambda%x) + dt * running%x
       if (any(.not. ieee_is_finite(current_adjoint%x))) then
          write(*,'(A,I0)') 'ROM time-varying adjoint failure: non-finite state at step ', step
          error stop 'non-finite ROM time-varying adjoint state'
       end if
       adjoint%x(:,step) = current_adjoint%x
    end do

    call neko_scratch_registry%relinquish_vector([running_index, next_index, &
         current_index])
    call neko_scratch_registry%relinquish_vector(lambda_index)
    nullify(running); nullify(next_adjoint); nullify(current_adjoint)
    nullify(lambda)
  end subroutine rom_model_discrete_adjoint_bdf1_timevarying

  !> Integrate a BDF1 discrete adjoint using a preassembled affine Jacobian
  !! tensor without retaining one dense Jacobian per time step.
  !!
  !! `linearized_tensor(k)` is the reduced matrix C_k in
  !! J_n = fixed_linear + sum_k a_(n,k) C_k.  The forward coefficients are
  !! read from `states(:,n)`. If requested, `residual_multipliers(:,n)` stores
  !! lambda_n = (H-dt*L_i)^(-T) p_(n+1) for subsequent spatial correlations.
  subroutine rom_model_discrete_adjoint_bdf1_tensor(this, states, dt, &
       running_state_gradient, terminal_state_gradient, fixed_linear, &
       linearized_tensor, adjoint, residual_multipliers)
    class(rom_model_t), intent(inout) :: this
    type(matrix_t), intent(in) :: states, running_state_gradient, fixed_linear
    real(kind=rp), intent(in) :: dt
    type(vector_t), intent(in) :: terminal_state_gradient
    type(matrix_t), intent(in) :: linearized_tensor(:)
    type(matrix_t), intent(inout) :: adjoint
    type(matrix_t), intent(inout), optional :: residual_multipliers
    type(vector_t), pointer :: running, next_adjoint, current_adjoint, lambda
    integer :: running_index, next_index, current_index, lambda_index
    integer :: n_steps, n_forward, step, k, e

    n_steps = running_state_gradient%get_ncols()
    e = this%n_state
    n_forward = states%get_nrows()
    if (states%get_ncols() /= n_steps + 1 .or. &
         running_state_gradient%get_nrows() /= e .or. &
         terminal_state_gradient%size() /= e .or. &
         fixed_linear%get_nrows() /= e .or. fixed_linear%get_ncols() /= e .or. &
         size(linearized_tensor) /= n_forward) &
         error stop 'ROM tensor adjoint state size mismatch'
    do k = 1, n_forward
       if (linearized_tensor(k)%get_nrows() /= e .or. &
            linearized_tensor(k)%get_ncols() /= e) &
            error stop 'ROM tensor adjoint Jacobian shape mismatch'
    end do
    if (any(adjoint%get_dims() /= [e, n_steps + 1])) &
         call adjoint%init(e, n_steps + 1)
    if (present(residual_multipliers)) then
       if (any(residual_multipliers%get_dims() /= [e, n_steps])) &
            call residual_multipliers%init(e, n_steps)
       residual_multipliers%x = 0.0_rp
    end if

    call neko_scratch_registry%request_vector(running, running_index, e, .true.)
    call neko_scratch_registry%request_vector(next_adjoint, next_index, e, .true.)
    call neko_scratch_registry%request_vector(current_adjoint, current_index, e, .true.)
    call neko_scratch_registry%request_vector(lambda, lambda_index, e, .true.)
    adjoint%x(:,n_steps+1) = terminal_state_gradient%x
    if (any(.not. ieee_is_finite(adjoint%x(:,n_steps+1)))) &
         error stop 'non-finite ROM tensor-adjoint terminal state'

    do step = n_steps, 1, -1
       next_adjoint%x = adjoint%x(:,step+1)
       running%x = running_state_gradient%x(:,step)
       call this%solve_step_transpose(dt, next_adjoint, lambda)
       current_adjoint%x = matmul(transpose(this%mass%x), lambda%x) + &
            dt * matmul(transpose(fixed_linear%x), lambda%x) + dt * running%x
       do k = 1, n_forward
          current_adjoint%x = current_adjoint%x + dt * states%x(k,step) * &
               matmul(transpose(linearized_tensor(k)%x), lambda%x)
       end do
       if (any(.not. ieee_is_finite(current_adjoint%x))) &
            error stop 'non-finite ROM tensor-adjoint state'
       adjoint%x(:,step) = current_adjoint%x
       if (present(residual_multipliers)) &
            residual_multipliers%x(:,step) = lambda%x
    end do

    call neko_scratch_registry%relinquish_vector([running_index, next_index, &
         current_index, lambda_index])
    nullify(running); nullify(next_adjoint); nullify(current_adjoint); nullify(lambda)
  end subroutine rom_model_discrete_adjoint_bdf1_tensor

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
    class(rom_model_t), intent(inout) :: this
    real(kind=rp), intent(in) :: dt
    type(vector_t), intent(in) :: rhs
    type(vector_t), intent(inout) :: solution
    call rom_model_prepare_bdf1_factor(this, dt)
    call rom_lu_solve_transpose(this%bdf1_lu, this%bdf1_pivots, rhs, solution)
  end subroutine rom_model_solve_step_transpose

  !> Build and cache the BDF1 LU factors for the current reduced operators.
  subroutine rom_model_prepare_bdf1_factor(this, dt)
    class(rom_model_t), intent(inout) :: this
    real(kind=rp), intent(in) :: dt
    if (this%bdf1_factor_valid .and. this%bdf1_factor_dt == dt) return
    call this%bdf1_lu%init(this%n_state, this%n_state)
    this%bdf1_lu%x = this%mass%x - dt * this%implicit_linear%x
    if (any(.not. ieee_is_finite(this%bdf1_lu%x))) then
       error stop 'non-finite ROM BDF1 step matrix'
    end if
    call rom_lu_factorize(this%bdf1_lu, this%bdf1_pivots)
    this%bdf1_factor_dt = dt
    this%bdf1_factor_valid = .true.
  end subroutine rom_model_prepare_bdf1_factor

  !> In-place partial-pivot LU factorization, P*B = L*U.
  subroutine rom_lu_factorize(factors, pivots)
    type(matrix_t), intent(inout) :: factors
    integer, allocatable, intent(inout) :: pivots(:)
    integer :: n, i, j, pivot
    real(kind=rp) :: largest, factor
    n = factors%get_nrows()
    if (n < 1 .or. factors%get_ncols() /= n) error stop 'ROM LU shape mismatch'
    if (allocated(pivots)) deallocate(pivots)
    allocate(pivots(n))
    do i = 1, n
       pivot = i
       largest = abs(factors%x(i,i))
       do j = i + 1, n
          if (abs(factors%x(j,i)) > largest) then
             largest = abs(factors%x(j,i))
             pivot = j
          end if
       end do
       if (largest <= epsilon(1.0_rp)) error stop 'singular ROM BDF1 step matrix'
       pivots(i) = pivot
       if (pivot /= i) call swap_rows(factors, i, pivot)
       do j = i + 1, n
          factor = factors%x(j,i) / factors%x(i,i)
          factors%x(j,i) = factor
          factors%x(j,i+1:n) = factors%x(j,i+1:n) - factor * factors%x(i,i+1:n)
       end do
    end do
  end subroutine rom_lu_factorize

  !> Solve B*x=rhs from cached P*B=L*U factors.
  subroutine rom_lu_solve(factors, pivots, rhs, solution)
    type(matrix_t), intent(in) :: factors
    integer, intent(in) :: pivots(:)
    type(vector_t), intent(in) :: rhs
    type(vector_t), intent(inout) :: solution
    real(kind=rp), allocatable :: work(:)
    real(kind=rp) :: tmp
    integer :: n, i, j
    n = rhs%size()
    if (factors%get_nrows() /= n .or. size(pivots) /= n) error stop 'ROM LU solve shape mismatch'
    if (solution%size() /= n) call solution%init(n)
    allocate(work(n)); work = rhs%x
    do i = 1, n
       if (pivots(i) /= i) then
          tmp = work(i); work(i) = work(pivots(i)); work(pivots(i)) = tmp
       end if
    end do
    do i = 1, n
       do j = i + 1, n
          work(j) = work(j) - factors%x(j,i) * work(i)
       end do
    end do
    do i = n, 1, -1
       if (i < n) work(i) = work(i) - sum(factors%x(i,i+1:n) * work(i+1:n))
       work(i) = work(i) / factors%x(i,i)
    end do
    solution%x = work
    deallocate(work)
  end subroutine rom_lu_solve

  !> Solve B^T*x=rhs from cached P*B=L*U factors.
  subroutine rom_lu_solve_transpose(factors, pivots, rhs, solution)
    type(matrix_t), intent(in) :: factors
    integer, intent(in) :: pivots(:)
    type(vector_t), intent(in) :: rhs
    type(vector_t), intent(inout) :: solution
    real(kind=rp), allocatable :: work(:)
    real(kind=rp) :: tmp
    integer :: n, i, j
    n = rhs%size()
    if (factors%get_nrows() /= n .or. size(pivots) /= n) error stop 'ROM transpose LU shape mismatch'
    if (solution%size() /= n) call solution%init(n)
    allocate(work(n)); work = rhs%x
    ! U^T*y=rhs.
    do i = 1, n
       if (i > 1) work(i) = work(i) - sum(factors%x(1:i-1,i) * work(1:i-1))
       work(i) = work(i) / factors%x(i,i)
    end do
    ! L^T*z=y (unit diagonal).
    do i = n, 1, -1
       if (i < n) work(i) = work(i) - sum(factors%x(i+1:n,i) * work(i+1:n))
    end do
    ! Undo the row permutations in reverse order to apply P^T.
    do i = n, 1, -1
       if (pivots(i) /= i) then
          tmp = work(i); work(i) = work(pivots(i)); work(pivots(i)) = tmp
       end if
    end do
    solution%x = work
    deallocate(work)
  end subroutine rom_lu_solve_transpose

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
