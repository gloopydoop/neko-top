!> Actual coupled-POD F:F ROM and vector-potential directional derivative check.
!!
!! Reads a one-rank, column-major [u,v,w,scalar] trial basis exported by the
!! ROM reduced-system writer. At the case's reference design it checks exact
!! mode preservation, checks the local design derivative of one transformed
!! velocity mode against the production Poisson/curl representation, and
!! compares the complete ROM objective directional derivative with a central
!! finite difference. This is a fixed-basis diagnostic; it does not run MMA.
program vector_potential_rom_directional_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use fluid_pnpn, only: fluid_pnpn_t
  use design, only: design_t, design_factory
  use brinkman_design, only: brinkman_design_t
  use problem, only: problem_t
  use rom_optimizer, only: rom_optimizer_t
  use field, only: field_t
  use operators, only: div
  use math, only: glsc3
  use json_module, only: json_file
  use json_utils, only: json_get
  use json_utils_ext, only: json_read_file
  use vector, only: vector_t
  use vector_math, only: vector_glsc2
  use num_types, only: rp
  use comm, only: NEKO_COMM, pe_rank, MPI_REAL_PRECISION
  use mpi_f08, only: MPI_Comm_size, MPI_Allreduce, MPI_SUM
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use utils, only: neko_error
  implicit none

  integer, parameter :: n_modes = 8, n_fields = 4
  real(kind=rp), parameter :: mode_tol = 5.0e-5_rp
  real(kind=rp), parameter :: objective_tol = 2.0e-5_rp
  type(simulation_t) :: simulation
  type(problem_t) :: problem
  type(rom_optimizer_t) :: optimizer
  class(design_t), allocatable :: design
  type(fluid_pnpn_t), pointer :: fluid
  type(json_file) :: parameters, design_parameters
  type(vector_t) :: x0, xplus, xminus, direction, mapped0, mapped_plus, mapped_minus
  type(vector_t) :: ratio_direction, mode_tangent(3), mode_fd(3)
  type(vector_t) :: dot_factor, factor_cotangent, dot_cotangent(3), dot_correction(3)
  type(vector_t) :: objective_gradient, mapped_density, directional_gradient
  type(field_t) :: mode_u, mode_v, mode_w, mode_divergence
  real(kind=rp), allocatable :: basis(:,:)
  real(kind=rp), allocatable :: mode_plus(:,:), mode_minus(:,:)
  real(kind=rp), allocatable :: final_state(:)
  real(kind=rp) :: objective, objective_plus, objective_minus, constraint
  real(kind=rp) :: constraint_plus, constraint_minus, fd, derivative_dot
  real(kind=rp) :: mode_error, mode_scale, mode_relative
  real(kind=rp) :: best_mode_error, best_objective_error, objective_scale
  real(kind=rp) :: poisson_true_residual
  real(kind=rp) :: reference_divergence_ratio, plus_divergence_ratio
  real(kind=rp) :: minus_divergence_ratio
  real(kind=rp) :: epsilons(3), tau
  real(kind=rp) :: configured_tau
  integer :: argc, n, nranks, ierr, unit, file_bytes, expected_bytes
  integer :: i, c, k, mode_map
  logical :: found, configured_logical
  character(len=:), allocatable :: configured_text
  character(len=1024) :: case_file, basis_file

  argc = command_argument_count()
  if (argc /= 2) call neko_error( &
       'Usage: vector_potential_rom_directional_diagnostic CASE TRIAL_BASIS.bin')
  call get_command_argument(1, case_file)
  call get_command_argument(2, basis_file)

  call neko_init()
  call neko_top_register_types()
  parameters = json_read_file(trim(case_file))
  call parameters%update('unsteady', .false., found)
  call parameters%update('optimization.solver.rom.basis.forward', &
       'forward', found)
  call parameters%update('optimization.solver.rom.basis.adjoint', &
       'forward', found)
  call parameters%update('optimization.solver.rom.basis_transform', &
       'vector_potential', found)
  call parameters%update('optimization.solver.rom.basis_transform_tau', &
       1.0_rp, found)
  call parameters%update('optimization.solver.rom.differentiate_basis', &
       .true., found)
  call parameters%update('optimization.solver.rom.enable_output', .false., found)
  call parameters%update('optimization.solver.rom.write_trajectory_fields', &
       .false., found)
  call parameters%update('optimization.solver.rom.export_reduced_system', &
       .false., found)
  call parameters%update('optimization.solver.enable_output', .false., found)
  call json_get(parameters, 'unsteady', configured_logical)
  if (configured_logical) call neko_error( &
       'JSON readback did not disable the unsteady driver')
  call json_get(parameters, 'optimization.solver.rom.basis.forward', &
       configured_text)
  if (trim(configured_text) /= 'forward') call neko_error( &
       'JSON readback did not select a forward trial basis')
  call json_get(parameters, 'optimization.solver.rom.basis.adjoint', &
       configured_text)
  if (trim(configured_text) /= 'forward') call neko_error( &
       'JSON readback did not select a forward test basis')
  call json_get(parameters, 'optimization.solver.rom.basis_transform', &
       configured_text)
  if (trim(configured_text) /= 'vector_potential') call neko_error( &
       'JSON readback did not select vector_potential')
  call json_get(parameters, 'optimization.solver.rom.basis_transform_tau', &
       configured_tau)
  if (abs(configured_tau-1.0_rp) > 10.0_rp*epsilon(1.0_rp)) &
       call neko_error('JSON readback did not set vector-potential tau to one')
  call json_get(parameters, 'optimization.solver.rom.differentiate_basis', &
       configured_logical)
  if (.not. configured_logical) call neko_error( &
       'JSON readback did not enable basis differentiation')
  call json_get(parameters, 'optimization.solver.rom.enable_output', &
       configured_logical)
  if (configured_logical) call neko_error('ROM diagnostic output must be disabled')
  call json_get(parameters, 'optimization.solver.rom.write_trajectory_fields', &
       configured_logical)
  if (configured_logical) call neko_error('Trajectory output must be disabled')
  call json_get(parameters, 'optimization.solver.rom.export_reduced_system', &
       configured_logical)
  if (configured_logical) call neko_error('Reduced-system export must be disabled')
  call json_get(parameters, 'optimization.solver.enable_output', configured_logical)
  if (configured_logical) call neko_error('Optimizer output must be disabled')

  call user_setup(simulation%neko_case%user)
  call simulation%init(parameters)
  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('Vector-potential ROM diagnostic requires Pn/Pn')
  end select
  call json_get(parameters, 'optimization.design', design_parameters)
  call design_factory(design, design_parameters, simulation)
  call problem%init(parameters, design, simulation)
  call optimizer%init_from_json(parameters, problem, design, simulation)
  optimizer%rom_problem%time_start = &
       real(simulation%neko_case%time%start_time, rp)
  optimizer%rom_problem%time_end = real(simulation%neko_case%time%end_time, rp)
  optimizer%rom_problem%dt = real(simulation%neko_case%time%dt, rp)
  optimizer%rom_problem%n_steps = nint((optimizer%rom_problem%time_end - &
       optimizer%rom_problem%time_start) / optimizer%rom_problem%dt)
  if (optimizer%rom_problem%dt <= 0.0_rp .or. &
       optimizer%rom_problem%n_steps < 1) &
       call neko_error('Vector-potential diagnostic has an invalid ROM time grid')
  if (optimizer%rom_problem%basis_transform /= 'vector_potential' .or. &
       .not. optimizer%rom_problem%differentiate_basis .or. &
       optimizer%rom_problem%sensitivity_source /= 'exact_galerkin' .or. &
       optimizer%rom_problem%basis_forward /= 'forward' .or. &
       optimizer%rom_problem%basis_adjoint /= 'forward' .or. &
       abs(optimizer%rom_problem%basis_transform_tau-1.0_rp) > &
       10.0_rp*epsilon(1.0_rp)) &
       call neko_error('Vector-potential diagnostic did not select exact differentiation')
  if (design%size() /= fluid%u%dof%size()) call neko_error( &
       'Vector-potential diagnostic requires an unmasked full local design field')

  call MPI_Comm_size(NEKO_COMM, nranks, ierr)
  if (nranks /= 1) call neko_error( &
       'Archived trial basis is one-rank; use a single MPI rank')
  n = fluid%u%dof%size()
  call mode_u%init(fluid%dm_Xh, 'vp_diag_mode_u')
  call mode_v%init(fluid%dm_Xh, 'vp_diag_mode_v')
  call mode_w%init(fluid%dm_Xh, 'vp_diag_mode_w')
  call mode_divergence%init(fluid%dm_Xh, 'vp_diag_mode_divergence')
  expected_bytes = n_fields * n * n_modes * storage_size(0.0_rp) / 8
  inquire(file=trim(basis_file), size=file_bytes, iostat=ierr)
  if (ierr /= 0 .or. file_bytes /= expected_bytes) &
       call neko_error('Coupled trial basis does not match case mesh/layout/rank')
  allocate(basis(n_fields*n,n_modes))
  open(newunit=unit, file=trim(basis_file), status='old', action='read', &
       access='stream', form='unformatted', iostat=ierr)
  if (ierr /= 0) call neko_error('Could not open coupled trial basis')
  read(unit, iostat=ierr) basis
  close(unit)
  if (ierr /= 0 .or. .not. all(ieee_is_finite(basis))) &
       call neko_error('Could not read finite coupled trial basis')
  call optimizer%initialize_from_basis(basis, n_modes, n_fields, simulation)

  call x0%init(design%size()); call xplus%init(design%size())
  call xminus%init(design%size()); call direction%init(design%size())
  call mapped0%init(design%size()); call mapped_plus%init(design%size())
  call mapped_minus%init(design%size()); call ratio_direction%init(design%size())
  call dot_factor%init(n); call factor_cotangent%init(n)
  do c = 1, 3
     call dot_cotangent(c)%init(n)
     call dot_correction(c)%init(n)
  end do
  call objective_gradient%init(design%size())
  call mapped_density%init(design%size()); call directional_gradient%init(design%size())
  allocate(mode_plus(n,3), mode_minus(n,3))
  do c = 1, 3
     call mode_tangent(c)%init(n)
     call mode_fd(c)%init(n)
  end do
  ! First reproduce the case's captured reference design and verify exact
  ! reference-mode preservation there. The supplied-basis initializer saves
  ! reference_chi while assembling the immutable archived POD basis.
  call design%get_values(x0)
  call design%update_design(x0)
  call optimizer%evaluate_fixed_basis(design, objective, objective_gradient%x, &
       constraint, final_state)
  if (.not. ieee_is_finite(objective) .or. &
       .not. ieee_is_finite(constraint) .or. &
       .not. all(ieee_is_finite(final_state))) &
       call neko_error('Nonfinite reference ROM evaluation')
  deallocate(final_state)
  call assert_reference_modes(optimizer)

  ! The archived starting design is binary. Use an interior test candidate on
  ! the configured optimization mask for central FD; update_design masks all
  ! exterior entries to zero. This leaves the copied case and POD unchanged.
  x0%x = 0.5_rp
  call design%update_design(x0)
  call design%get_values(x0)
  call get_mapped_design(design, mapped0)
  call optimizer%evaluate_fixed_basis(design, objective, objective_gradient%x, &
       constraint, final_state, mapped_density%x)
  call design%map_backward(mapped_density)
  select type (des => design)
  type is (brinkman_design_t)
     call des%get_sensitivity(directional_gradient)
     call des%convert_to_directional_derivative(directional_gradient)
  class default
     call neko_error('Vector-potential diagnostic requires Brinkman design')
  end select
  if (.not. ieee_is_finite(objective) .or. &
       .not. ieee_is_finite(constraint) .or. &
       .not. all(ieee_is_finite(final_state))) &
       call neko_error('Nonfinite interior ROM evaluation')
  deallocate(final_state)
  if (.not. all(ieee_is_finite(objective_gradient%x))) &
       call neko_error('Nonfinite ROM objective gradient at interior design')
  reference_divergence_ratio = mode_divergence_ratio( &
       optimizer%rom_problem%reference_velocity_modes(:,1))
  poisson_true_residual = maxval( &
       optimizer%rom_problem%reference_forward_potential%solve_true_residual)
  if (.not. ieee_is_finite(poisson_true_residual) .or. &
       poisson_true_residual > 1.0e-8_rp) &
       call neko_error('Vector-potential Poisson true residual failed')
  if (pe_rank == 0) write(*,'(A,ES14.6)') &
       'Maximum independent Poisson relative residual = ', poisson_true_residual
  call check_apply_pullback(optimizer, dot_factor, dot_cotangent, &
       dot_correction, factor_cotangent)

  ! Symmetric perturbations stay in [0,1], including coordinates on a bound.
  do i = 1, design%size()
     if (x0%x(i) < 0.0_rp .or. x0%x(i) > 1.0_rp) &
          call neko_error('Reference design is outside feasible [0,1] bounds')
     direction%x(i) = 0.1_rp * sin(real(i,rp)*0.31_rp) * &
          max(0.0_rp, min(x0%x(i), 1.0_rp-x0%x(i)))
  end do
  derivative_dot = vector_glsc2(directional_gradient, direction)
  if (vector_glsc2(direction,direction) <= tiny(1.0_rp)) &
       call neko_error('Reference candidate has no feasible diagnostic direction')
  epsilons = [1.0e-4_rp, 1.0e-5_rp, 1.0e-6_rp]
  tau = optimizer%rom_problem%basis_transform_tau
  best_mode_error = huge(1.0_rp)
  best_objective_error = huge(1.0_rp)

  do k = 1, size(epsilons)
     xplus%x = x0%x + epsilons(k)*direction%x
     xminus%x = x0%x - epsilons(k)*direction%x
     call design%update_design(xplus)
     call optimizer%evaluate_fixed_basis(design, objective_plus, &
          objective_gradient%x, constraint_plus, final_state)
     if (.not. ieee_is_finite(objective_plus) .or. &
          .not. all(ieee_is_finite(final_state))) &
          call neko_error('Nonfinite positive ROM directional evaluation')
     deallocate(final_state)
     call get_mapped_design(design, mapped_plus)
     do c = 1, 3
        mode_plus(:,c) = optimizer%rom_problem%velocity_modes(c,1)%x
     end do
     plus_divergence_ratio = mode_divergence_ratio( &
          optimizer%rom_problem%velocity_modes(:,1))

     call design%update_design(xminus)
     call optimizer%evaluate_fixed_basis(design, objective_minus, &
          objective_gradient%x, constraint_minus, final_state)
     if (.not. ieee_is_finite(objective_minus) .or. &
          .not. all(ieee_is_finite(final_state))) &
          call neko_error('Nonfinite negative ROM directional evaluation')
     deallocate(final_state)
     call get_mapped_design(design, mapped_minus)
     do c = 1, 3
        mode_minus(:,c) = optimizer%rom_problem%velocity_modes(c,1)%x
        mode_fd(c)%x = (mode_plus(:,c)-mode_minus(:,c))/(2.0_rp*epsilons(k))
     end do
     minus_divergence_ratio = mode_divergence_ratio( &
          optimizer%rom_problem%velocity_modes(:,1))
     fd = (objective_plus-objective_minus)/(2.0_rp*epsilons(k))

     ! Local B check: d(R A)/d(direction) is obtained by applying the stored
     ! native potential/curl map to dR, where dchi is measured along the same
     ! feasible design direction and R=D0/D.
     do i = 1, design%size()
        ratio_direction%x(i) = -tau * &
             optimizer%rom_problem%reference_D%x(i) / &
             (1.0_rp + tau*mapped0%x(i))**2 * &
             (mapped_plus%x(i)-mapped_minus%x(i))/(2.0_rp*epsilons(k))
     end do
     mode_map = optimizer%rom_problem%reference_potential_forward_map(1)
     call optimizer%rom_problem%reference_forward_potential%apply( &
          mode_map, ratio_direction, mode_tangent)
     call mode_error_from_fd(mode_fd, mode_tangent, mode_error, mode_scale)
     mode_relative = mode_error / max(mode_scale, tiny(1.0_rp))
     objective_scale = max(abs(fd), abs(derivative_dot))
     objective_scale = max(objective_scale, 1.0e-12_rp)
     if (mode_relative < best_mode_error) best_mode_error = mode_relative
     if (abs(fd-derivative_dot) / (1.0e-12_rp + &
          objective_tol*objective_scale) < best_objective_error) &
          best_objective_error = abs(fd-derivative_dot) / &
          (1.0e-12_rp + objective_tol*objective_scale)
     if (pe_rank == 0) write(*,'(A,ES10.3,A,ES14.6,A,ES14.6,A,ES10.3,A,ES10.3)') &
          'eps=', epsilons(k), ' FD=', fd, ' adjoint=', derivative_dot, &
          ' objective_scaled_error=', abs(fd-derivative_dot) / &
          (1.0e-12_rp + objective_tol*objective_scale), &
          ' mode_relative_error=', mode_relative
     if (pe_rank == 0) write(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
          'mode strong-divergence ratios ref/plus/minus = ', &
          reference_divergence_ratio, ' / ', plus_divergence_ratio, ' / ', &
          minus_divergence_ratio, ' ; plus-minus=', &
          abs(plus_divergence_ratio-minus_divergence_ratio)
  end do

  if (best_mode_error > mode_tol) call neko_error( &
       'Actual-mode vector-potential directional derivative failed')
  if (best_objective_error > 1.0_rp) call neko_error( &
       'Complete ROM objective directional derivative failed')
  if (pe_rank == 0) write(*,'(A,ES10.3,A,ES10.3)') &
       'Best mode relative error / objective scaled error = ', &
       best_mode_error, ' / ', best_objective_error

  do c = 1, 3
     call mode_tangent(c)%free(); call mode_fd(c)%free()
  end do
  call x0%free(); call xplus%free(); call xminus%free(); call direction%free()
  call mapped0%free(); call mapped_plus%free(); call mapped_minus%free()
  call ratio_direction%free()
  call dot_factor%free(); call factor_cotangent%free()
  do c = 1, 3
     call dot_cotangent(c)%free(); call dot_correction(c)%free()
  end do
  call objective_gradient%free(); call mapped_density%free(); call directional_gradient%free()
  call mode_u%free(); call mode_v%free(); call mode_w%free()
  call mode_divergence%free()
  deallocate(basis, mode_plus, mode_minus)
  call optimizer%free(); call problem%free(); call design%free()
  deallocate(design)
  call simulation%free()
  call neko_finalize()

contains

  subroutine get_mapped_design(design_in, values)
    class(design_t), intent(in) :: design_in
    type(vector_t), intent(inout) :: values
    select type (des => design_in)
    type is (brinkman_design_t)
       call des%get_mapped_values(values)
    class default
       call neko_error('Vector-potential ROM diagnostic requires Brinkman design')
    end select
  end subroutine get_mapped_design

  subroutine assert_reference_modes(opt)
    type(rom_optimizer_t), intent(in) :: opt
    real(kind=rp) :: local_err, local_ref, global_err, global_ref
    integer :: comp, mode_index
    local_err = 0.0_rp
    local_ref = 0.0_rp
    do mode_index = 1, n_modes
       do comp = 1, 3
          local_err = local_err + sum(( &
               opt%rom_problem%velocity_modes(comp,mode_index)%x - &
               opt%rom_problem%reference_velocity_modes(comp,mode_index)%x)**2)
          local_ref = local_ref + sum(opt%rom_problem%reference_velocity_modes( &
               comp,mode_index)%x**2)
       end do
    end do
    call MPI_Allreduce(local_err, global_err, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
    call MPI_Allreduce(local_ref, global_ref, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
    if (sqrt(global_err/max(global_ref,tiny(1.0_rp))) > 1.0e-12_rp) &
         call neko_error('Vector-potential map did not preserve reference modes')
  end subroutine assert_reference_modes

  subroutine check_apply_pullback(opt, factor, cotangent, correction, pullback)
    type(rom_optimizer_t), intent(inout) :: opt
    type(vector_t), intent(inout) :: factor, pullback
    type(vector_t), intent(inout) :: cotangent(3), correction(3)
    real(kind=rp) :: lhs, rhs, discrepancy, scale
    integer :: point, comp, potential_index
    do point = 1, factor%size()
       factor%x(point) = sin(0.173_rp*real(point,rp))
       do comp = 1, 3
          cotangent(comp)%x(point) = &
               cos(0.237_rp*real(point+17*comp,rp))
       end do
    end do
    potential_index = opt%rom_problem%reference_potential_forward_map(1)
    call opt%rom_problem%reference_forward_potential%apply( &
         potential_index, factor, correction)
    call opt%rom_problem%reference_forward_potential%pullback( &
         potential_index, factor, cotangent, pullback)
    lhs = 0.0_rp
    do comp = 1, 3
       lhs = lhs + vector_glsc2(cotangent(comp), correction(comp))
    end do
    rhs = vector_glsc2(factor, pullback)
    discrepancy = abs(lhs-rhs)
    scale = max(abs(lhs),abs(rhs),tiny(1.0_rp))
    if (pe_rank == 0) write(*,'(A,ES14.6,A,ES14.6,A,ES10.3)') &
         'Vector-potential dot lhs/rhs/relative-error = ', lhs, ' / ', rhs, &
         ' relative=', discrepancy/scale
    if (.not. ieee_is_finite(lhs) .or. .not. ieee_is_finite(rhs) .or. &
         discrepancy/scale > 5.0e-11_rp) &
         call neko_error('Vector-potential apply/pullback dot test failed')
  end subroutine check_apply_pullback

  function mode_divergence_ratio(mode) result(ratio)
    type(vector_t), intent(in) :: mode(3)
    real(kind=rp) :: ratio, norm_div, norm_velocity
    mode_u%x = reshape(mode(1)%x, shape(mode_u%x))
    mode_v%x = reshape(mode(2)%x, shape(mode_v%x))
    mode_w%x = reshape(mode(3)%x, shape(mode_w%x))
    call div(mode_divergence%x, mode_u%x, mode_v%x, mode_w%x, fluid%c_Xh)
    norm_div = sqrt(max(0.0_rp, glsc3(mode_divergence%x, &
         fluid%c_Xh%B, mode_divergence%x, n)))
    norm_velocity = sqrt(max(0.0_rp, &
         glsc3(mode_u%x,fluid%c_Xh%B,mode_u%x,n) + &
         glsc3(mode_v%x,fluid%c_Xh%B,mode_v%x,n) + &
         glsc3(mode_w%x,fluid%c_Xh%B,mode_w%x,n)))
    ratio = norm_div / max(norm_velocity,tiny(1.0_rp))
  end function mode_divergence_ratio

  subroutine mode_error_from_fd(fd_mode, tangent, error_norm, scale_norm)
    type(vector_t), intent(in) :: fd_mode(3)
    type(vector_t), intent(in) :: tangent(3)
    real(kind=rp), intent(out) :: error_norm, scale_norm
    real(kind=rp) :: local_error, local_scale, global_error, global_scale
    integer :: comp
    local_error = 0.0_rp; local_scale = 0.0_rp
    do comp = 1, 3
       local_error = local_error + sum((fd_mode(comp)%x-tangent(comp)%x)**2)
       local_scale = local_scale + sum(tangent(comp)%x**2)
    end do
    call MPI_Allreduce(local_error, global_error, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
    call MPI_Allreduce(local_scale, global_scale, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
    error_norm = sqrt(global_error)
    scale_norm = sqrt(global_scale)
  end subroutine mode_error_from_fd

end program vector_potential_rom_directional_diagnostic
