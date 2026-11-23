program production_rom_gradient_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use simulation_m, only: simulation_t
  use design, only: design_t, design_factory
  use brinkman_design, only: brinkman_design_t
  use problem, only: problem_t
  use rom_optimizer, only: rom_optimizer_t
  use rom_basis_handoff, only: rom_project_basis
  use fluid_pnpn, only: fluid_pnpn_t
  use file, only: file_t
  use fld_file_data, only: fld_file_data_t
  use field, only: field_t
  use operators, only: grad
  use field_math, only: field_rzero
  use field_math, only: field_rzero
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2
  use json_module, only: json_file
  use json_utils, only: json_get
  use json_utils_ext, only: json_read_file
  use vector, only: vector_t
  use comm, only: NEKO_COMM, pe_rank, MPI_REAL_PRECISION
  use mpi_f08, only: MPI_Allreduce, MPI_SUM
  use num_types, only: rp
  use utils, only: neko_error
  implicit none

  type(simulation_t) :: simulation
  type(problem_t) :: problem
  type(rom_optimizer_t) :: optimizer
  class(design_t), allocatable :: design
  type(fluid_pnpn_t), pointer :: fluid
  type(file_t) :: modes_file
  type(fld_file_data_t) :: modes_data
  type(json_file) :: parameters, design_parameters
  type(field_t) :: u, v, w, psi, grad_x, grad_y, grad_z
  type(vector_t) :: x0, xplus, xminus, direction
  real(kind=rp), allocatable :: basis(:, :), gradient(:), chi_gradient(:), &
       chi_plus(:), chi_minus(:), final_state(:)
  real(kind=rp), allocatable :: epsilons(:)
  real(kind=rp) :: objective, objective_plus, objective_minus
  real(kind=rp) :: constraint, constraint_plus, constraint_minus
  real(kind=rp) :: fd, analytic, chi_analytic, error_value, chi_error
  real(kind=rp) :: local_dot, global_dot, max_error
  real(kind=rp) :: x, y
  integer :: argc, n, n_modes, r, mode, i, k
  logical :: found
  character(len=1024) :: case_file, modes_path, objective_mode, basis_source

  argc = command_argument_count()
  if (argc < 2) call neko_error( &
       'Usage: production_rom_gradient_diagnostic CASE MODES.fld [all|viscous|brinkman]')
  call get_command_argument(1, case_file)
  call get_command_argument(2, modes_path)
  objective_mode = 'all'
  if (argc >= 3) call get_command_argument(3, objective_mode)
  basis_source = 'file'
  if (argc >= 4) call get_command_argument(4, basis_source)

  call neko_init()
  call neko_top_register_types()
  parameters = json_read_file(trim(case_file))
  ! Initialize only the spatial objects.  The normal unsteady simulation
  ! startup would wait for the Python state-recovery peer, so the diagnostic
  ! supplies the already configured time grid directly below.
  call parameters%update('unsteady', .false., found)
  call simulation%init(parameters)
  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('Production ROM gradient diagnostic requires Pn/Pn')
  end select

  call json_get(parameters, 'optimization.design', design_parameters)
  call design_factory(design, design_parameters, simulation)
  call problem%init(parameters, design, simulation)
  call optimizer%init_from_json(parameters, problem, design, simulation)
  select case (trim(objective_mode))
  case ('all')
     continue
  case ('viscous')
     optimizer%rom_problem%has_brinkman = .false.
     optimizer%rom_problem%brinkman_weight = 0.0_rp
  case ('brinkman')
     optimizer%rom_problem%has_viscous = .false.
     optimizer%rom_problem%viscous_weight = 0.0_rp
  case default
     call neko_error('Unknown production ROM objective diagnostic mode')
  end select
  optimizer%rom_problem%time_start = real(simulation%neko_case%time%start_time, rp)
  optimizer%rom_problem%time_end = real(simulation%neko_case%time%end_time, rp)
  optimizer%rom_problem%dt = real(simulation%neko_case%time%dt, rp)
  optimizer%rom_problem%n_steps = nint((optimizer%rom_problem%time_end - &
       optimizer%rom_problem%time_start) / optimizer%rom_problem%dt)
  if (optimizer%rom_problem%n_steps < 1 .or. optimizer%rom_problem%dt <= 0.0_rp) then
     call neko_error('Production ROM diagnostic has an invalid time grid')
  end if

  n = fluid%u%dof%size()
  call u%init(fluid%dm_Xh, 'production_rom_gradient_u')
  call v%init(fluid%dm_Xh, 'production_rom_gradient_v')
  call w%init(fluid%dm_Xh, 'production_rom_gradient_w')
  if (trim(basis_source) == 'manufactured') then
     ! The stream-function columns are the explicitly admissible, forward-only
     ! diagnostic baseline.  The production projection is still called here;
     ! a basis is handed to the ROM only if its native weak/BC/mass gate passes.
     r = 2
     allocate(basis(3*n, r))
     call psi%init(fluid%dm_Xh, 'production_rom_gradient_psi')
     call grad_x%init(fluid%dm_Xh, 'production_rom_gradient_grad_x')
     call grad_y%init(fluid%dm_Xh, 'production_rom_gradient_grad_y')
     call grad_z%init(fluid%dm_Xh, 'production_rom_gradient_grad_z')
     do mode = 1, r
        do i = 1, n
           x = fluid%u%dof%x%x(i, 1, 1, 1)
           y = fluid%u%dof%y%x(i, 1, 1, 1)
           if (mode == 1) then
              psi%x(i,1,1,1) = x*x*(2.0_rp-x)*(2.0_rp-x) * &
                   y*y*(1.0_rp-y)*(1.0_rp-y)
           else
              psi%x(i,1,1,1) = x*x*x*(2.0_rp-x)*(2.0_rp-x) * &
                   y*y*(1.0_rp-y)*(1.0_rp-y)*(1.0_rp+0.25_rp*y)
           end if
        end do
        call fluid%gs_Xh%op(psi%x, n, GS_OP_ADD)
        call col2(psi%x, fluid%c_Xh%mult, n)
        call grad(grad_x%x, grad_y%x, grad_z%x, psi%x, fluid%c_Xh)
        u%x = grad_y%x
        v%x = -grad_x%x
        call field_rzero(w)
        call fluid%gs_Xh%op(u%x, v%x, w%x, n, GS_OP_ADD)
        call col2(u%x, fluid%c_Xh%mult, n)
        call col2(v%x, fluid%c_Xh%mult, n)
        call col2(w%x, fluid%c_Xh%mult, n)
        call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
        basis(1:n, mode) = reshape(u%x, [n])
        basis(n+1:2*n, mode) = reshape(v%x, [n])
        basis(2*n+1:3*n, mode) = reshape(w%x, [n])
     end do
     n_modes = r
     call parameters%update('rom.admissibility_projection', 'pressure', found)
     call rom_project_basis(simulation, basis, n_modes)
     if (n_modes /= r) call neko_error( &
          'Manufactured admissible basis lost a mode during projection')
     call psi%free(); call grad_x%free(); call grad_y%free(); call grad_z%free()
  else
     call modes_data%init()
     call modes_file%init(trim(modes_path))
     call modes_file%set_counter(0)
     call modes_file%read(modes_data)
     if (modes_data%n_scalars < 3 .or. mod(modes_data%n_scalars, 3) /= 0) then
        call neko_error('Production ROM mode file must contain 3*N fields')
     end if
     if (modes_data%s(1)%size() /= n) then
        call neko_error('Production ROM mode partition does not match case')
     end if
     n_modes = modes_data%n_scalars / 3
     r = min(3, n_modes)
     allocate(basis(3*n, r))
     do mode = 1, r
        u%x = reshape(modes_data%s(3*(mode-1)+1)%x, shape(u%x))
        v%x = reshape(modes_data%s(3*(mode-1)+2)%x, shape(v%x))
        w%x = reshape(modes_data%s(3*(mode-1)+3)%x, shape(w%x))
        call fluid%gs_Xh%op(u%x, v%x, w%x, n, GS_OP_ADD)
        call col2(u%x, fluid%c_Xh%mult, n)
        call col2(v%x, fluid%c_Xh%mult, n)
        call col2(w%x, fluid%c_Xh%mult, n)
        call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
        basis(1:n, mode) = reshape(u%x, [n])
        basis(n+1:2*n, mode) = reshape(v%x, [n])
        basis(2*n+1:3*n, mode) = reshape(w%x, [n])
     end do
     call modes_data%free(); call modes_file%free()
  end if
  call u%free(); call v%free(); call w%free()

  ! This hook enters the production assemble_fixed/evaluate_rom/project_density
  ! path with a supplied basis.  It is not a bypass for admissibility: the
  ! forward POD file is used only as a separately diagnosed baseline, while
  ! the adjoint-enriched production handoff remains fail-closed.
  call optimizer%initialize_from_basis(basis, r, 3, simulation)

  call x0%init(design%size()); call xplus%init(design%size())
  call xminus%init(design%size()); call direction%init(design%size())
  x0%x = 0.5_rp
  do i = 1, design%size()
     direction%x(i) = 0.07_rp * sin(real(i, rp) * 0.31_rp)
  end do
  call design%update_design(x0)
  allocate(gradient(design%size()), chi_gradient(design%size()), &
       chi_plus(design%size()), chi_minus(design%size()))
  call optimizer%evaluate_fixed_basis(design, objective, gradient, constraint, &
       final_state, chi_gradient)
  deallocate(final_state)
  call get_mapped_design(design, chi_plus)
  allocate(epsilons(3), source = [1.0e-4_rp, 1.0e-5_rp, 1.0e-6_rp])
  if (pe_rank == 0) write(*,'(A,ES12.4,A,ES12.4)') &
       'Production ROM baseline objective/constraint = ', objective, ' / ', constraint

  max_error = 0.0_rp
  do k = 1, size(epsilons)
     xplus%x = x0%x + epsilons(k) * direction%x
     xminus%x = x0%x - epsilons(k) * direction%x
     call design%update_design(xplus)
     call optimizer%evaluate_fixed_basis(design, objective_plus, gradient, &
          constraint_plus, final_state)
     deallocate(final_state)
     call get_mapped_design(design, chi_plus)
     call design%update_design(xminus)
     call optimizer%evaluate_fixed_basis(design, objective_minus, gradient, &
          constraint_minus, final_state)
     deallocate(final_state)
     call get_mapped_design(design, chi_minus)
     fd = (objective_plus - objective_minus) / (2.0_rp * epsilons(k))

     call design%update_design(x0)
     call optimizer%evaluate_fixed_basis(design, objective, gradient, constraint, &
          final_state)
     deallocate(final_state)
     local_dot = dot_product(gradient, direction%x)
     call MPI_Allreduce(local_dot, global_dot, 1, MPI_REAL_PRECISION, MPI_SUM, &
          NEKO_COMM)
     analytic = global_dot
     error_value = abs(fd - analytic)
     local_dot = sum(optimizer%rom_problem%weights%x * chi_gradient * (chi_plus - chi_minus) / &
          (2.0_rp * epsilons(k)))
     call MPI_Allreduce(local_dot, global_dot, 1, MPI_REAL_PRECISION, MPI_SUM, &
          NEKO_COMM)
     chi_analytic = global_dot
     chi_error = abs(fd - chi_analytic)
     max_error = max(max_error, error_value, chi_error)
     if (pe_rank == 0) write(*,'(A,ES10.3,A,ES14.6,A,ES14.6,A,ES10.3,A,ES14.6,A,ES10.3)') &
          'eps=', epsilons(k), ' FD=', fd, ' adjoint=', analytic, &
          ' mapped_error=', error_value, ' chi_adjoint=', chi_analytic, &
          ' chi_error=', chi_error
  end do

  if (max_error > 2.0e-5_rp * max(1.0_rp, abs(fd))) then
     call neko_error('Production ROM mapped gradient finite-difference gate failed')
  end if

  call x0%free(); call xplus%free(); call xminus%free(); call direction%free()
  deallocate(epsilons, gradient, basis)
  call optimizer%free(); call problem%free(); call design%free()
  if (allocated(design)) deallocate(design)
  call simulation%free()
  call neko_finalize()

contains

  !> Extract the mapped Brinkman field on the diagnostic partition.
  subroutine get_mapped_design(design, values)
    class(design_t), intent(in) :: design
    real(kind=rp), intent(out) :: values(:)
    type(vector_t) :: mapped
    if (size(values) /= design%size()) call neko_error( &
         'Production ROM mapped-design diagnostic size mismatch')
    call mapped%init(size(values))
    select type (des => design)
    type is (brinkman_design_t)
       call des%get_mapped_values(mapped)
    class default
       call neko_error('Production ROM mapped-design diagnostic requires Brinkman')
    end select
    values = mapped%x
    call mapped%free()
  end subroutine get_mapped_design

end program production_rom_gradient_diagnostic
