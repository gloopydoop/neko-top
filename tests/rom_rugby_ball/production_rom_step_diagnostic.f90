program production_rom_step_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use simulation_m, only: simulation_t
  use design, only: design_t, design_factory
  use problem, only: problem_t
  use rom_optimizer, only: rom_optimizer_t
  use rom_basis_handoff, only: rom_project_basis
  use fluid_pnpn, only: fluid_pnpn_t
  use file, only: file_t
  use fld_file_data, only: fld_file_data_t
  use field, only: field_t
  use field_math, only: field_rzero
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2
  use operators, only: grad
  use json_module, only: json_file
  use json_utils, only: json_get
  use json_utils_ext, only: json_read_file
  use vector, only: vector_t
  use num_types, only: rp
  use comm, only: pe_rank
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
  type(vector_t) :: before, after
  real(kind=rp), allocatable :: basis(:, :)
  integer :: argc, n, r
  logical :: found, engineering_mode
  logical :: converged
  real(kind=rp) :: retained_change, engineering_horizon
  character(len=1024) :: case_file, modes_path, basis_source, argument
  integer :: step, n_projected, n_outer, clock_start, clock_end, clock_rate
  integer :: engineering_unit
  real(kind=rp) :: elapsed_seconds, proposed_change, final_change

  argc = command_argument_count()
  if (argc < 2) call neko_error( &
       'Usage: production_rom_step_diagnostic CASE MODES.fld [file|manufactured] '&
       //'[engineering] [horizon] [outer_steps]')
  call get_command_argument(1, case_file)
  call get_command_argument(2, modes_path)
  basis_source = 'file'
  if (argc >= 3) call get_command_argument(3, basis_source)
  engineering_mode = .false.
  if (argc >= 4) then
     call get_command_argument(4, argument)
     engineering_mode = trim(argument) == 'engineering'
  end if
  engineering_horizon = 0.1_rp
  if (argc >= 5) then
     call get_command_argument(5, argument)
     read(argument, *) engineering_horizon
  end if
  n_outer = 3
  if (engineering_mode) n_outer = 1
  if (argc >= 6) then
     call get_command_argument(6, argument)
     read(argument, *) n_outer
  end if
  if (engineering_horizon <= 0.0_rp .or. n_outer < 1) then
     call neko_error('engineering horizon and outer_steps must be positive')
  end if

  call neko_init()
  call neko_top_register_types()
  parameters = json_read_file(trim(case_file))
  if (engineering_mode) then
     call parameters%update('case.time.end_time', engineering_horizon, found)
     if (.not. found) call neko_error('could not set engineering time horizon')
  end if
  call parameters%update('unsteady', .false., found)
  call simulation%init(parameters)
  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('Production ROM step diagnostic requires Pn/Pn')
  end select

  call json_get(parameters, 'optimization.design', design_parameters)
  call design_factory(design, design_parameters, simulation)
  call problem%init(parameters, design, simulation)
  call optimizer%init_from_json(parameters, problem, design, simulation)
  optimizer%enable_output = engineering_mode
  if (.not. engineering_mode) optimizer%inner_iterations = 1
  optimizer%rom_problem%rank_total_max = 1
  optimizer%rom_problem%time_start = real(simulation%neko_case%time%start_time, rp)
  ! Engineering mode exercises the real reduced forward/adjoint, MMA,
  ! full-order refresh, objective and sensitivity path on an explicitly
  ! supplied basis. It intentionally does not call project_rom_basis: this
  ! diagnostic is for an inadmissible basis and is never scientific evidence.
  if (engineering_mode) then
     optimizer%rom_problem%time_end = real(simulation%neko_case%time%end_time, rp)
  else
     optimizer%rom_problem%time_end = min(real(simulation%neko_case%time%end_time, rp), &
          optimizer%rom_problem%time_start + 0.1_rp)
  end if
  simulation%neko_case%time%end_time = optimizer%rom_problem%time_end
  ! The adjoint case is constructed from the original case during
  ! simulation%init.  Keep the bounded diagnostic horizon consistent for
  ! both full-order loops; otherwise validation silently launches the
  ! production three-second adjoint after the shortened forward run.
  simulation%adjoint_case%time%start_time = simulation%neko_case%time%start_time
  simulation%adjoint_case%time%end_time = simulation%neko_case%time%end_time
  optimizer%rom_problem%dt = real(simulation%neko_case%time%dt, rp)
  optimizer%rom_problem%n_steps = nint((optimizer%rom_problem%time_end - &
       optimizer%rom_problem%time_start) / optimizer%rom_problem%dt)

  n = fluid%u%dof%size()
  if (trim(basis_source) == 'manufactured') then
     optimizer%rom_problem%rank_total_max = 2
     r = 2
     allocate(basis(3*n, r))
     call build_manufactured_basis(fluid, basis, n)
     n_projected = r
     call parameters%update('rom.admissibility_projection', 'pressure', found)
     call rom_project_basis(simulation, basis, n_projected)
     if (n_projected /= r) call neko_error( &
          'Manufactured step basis lost a mode during production projection')
  else
     call modes_data%init()
     call modes_file%init(trim(modes_path))
     call modes_file%set_counter(0)
     call modes_file%read(modes_data)
     if (modes_data%n_scalars < 3) call neko_error('Step diagnostic needs one velocity mode')
     if (modes_data%s(1)%size() /= n) call neko_error('Step diagnostic partition mismatch')
     r = 1
     allocate(basis(3*n, r))
     call load_mode(fluid, modes_data, basis(:,1), n)
     call modes_data%free(); call modes_file%free()
  end if

  call before%init(design%size())
  call after%init(design%size())
  if (engineering_mode) then
     call design%get_values(before)
  else
     before%x = 0.5_rp
  end if
  call design%update_design(before)
  retained_change = 0.0_rp
  if (engineering_mode .and. pe_rank == 0) then
     open(newunit=engineering_unit, file='engineering_summary.csv', &
          status='replace', action='write')
     write(engineering_unit, '(A)') &
          'outer_step,proposed_max_design_change,post_validation_design_change,'//&
          'elapsed_seconds,candidate_kept'
  end if
  do step = 1, n_outer
     if (engineering_mode) call system_clock(clock_start, clock_rate)
     converged = optimizer%step_from_basis(step, problem, design, simulation, basis, r, 3)
     if (engineering_mode) call system_clock(clock_end)
     call design%get_values(after)
     final_change = maxval(abs(after%x - before%x))
     retained_change = max(retained_change, final_change)
     proposed_change = optimizer%max_design_change
     if (engineering_mode .and. pe_rank == 0) then
        elapsed_seconds = real(clock_end - clock_start, rp) / &
             real(max(1, clock_rate), rp)
        write(engineering_unit, '(I0,3(",",ES24.16),",",L1)') step, &
             proposed_change, final_change, elapsed_seconds, final_change > 1.0e-12_rp
     end if
     if (.not. engineering_mode .and. &
          maxval(abs(after%x - before%x)) <= 1.0e-12_rp) then
        call neko_error('ROM candidate-retention diagnostic did not keep the update')
     end if
  end do
  if (engineering_mode .and. pe_rank == 0) close(engineering_unit)
  if (engineering_mode) then
     if (pe_rank == 0) then
        open(newunit=step, file='engineering_rom_status.txt', status='replace', action='write')
        write(step, '(A)') 'mode=engineering_supplied_forward_basis'
        write(step, '(A)') 'scientific_validation=NOT_CERTIFIED'
        write(step, '(A)') 'admissibility=NOT_APPLIED_IN_THIS_DIAGNOSTIC'
        write(step, '(A,A)') 'basis_source=', trim(basis_source)
        write(step, '(A,ES24.16)') 'post_validation_max_design_change=', retained_change
        write(step, '(A,ES24.16)') 'proposed_max_design_change=', optimizer%max_design_change
        write(step, '(A)') 'rom_validation.csv records full-order metrics; F means validation failed but the candidate was kept.'
        close(step)
     end if
    write(*,'(A,ES12.4)') &
          'ENGINEERING ROM run completed; not scientifically validated; max design change=', &
          retained_change
  else
     write(*,'(A,ES12.4)') &
          'ROM candidate-retention diagnostic passed for three outer steps; retained max design change=', &
          retained_change
  end if

  call before%free(); call after%free(); deallocate(basis)
  call optimizer%free(); call problem%free(); call design%free()
  if (allocated(design)) deallocate(design)
  call simulation%free()
  call neko_finalize()

contains

  !> Load one velocity mode from the diagnostic field-file data.
  subroutine load_mode(fluid, data, mode, n)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(fld_file_data_t), intent(in) :: data
    real(kind=rp), intent(out) :: mode(:)
    integer, intent(in) :: n
    type(field_t) :: u, v, w
    call u%init(fluid%dm_Xh, 'step_diagnostic_u')
    call v%init(fluid%dm_Xh, 'step_diagnostic_v')
    call w%init(fluid%dm_Xh, 'step_diagnostic_w')
    u%x = reshape(data%s(1)%x, shape(u%x))
    v%x = reshape(data%s(2)%x, shape(v%x))
    w%x = reshape(data%s(3)%x, shape(w%x))
    call fluid%gs_Xh%op(u%x, v%x, w%x, n, GS_OP_ADD)
    call col2(u%x, fluid%c_Xh%mult, n)
    call col2(v%x, fluid%c_Xh%mult, n)
    call col2(w%x, fluid%c_Xh%mult, n)
    call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
    mode(1:n) = reshape(u%x, [n])
    mode(n+1:2*n) = reshape(v%x, [n])
    mode(2*n+1:3*n) = reshape(w%x, [n])
    call u%free(); call v%free(); call w%free()
  end subroutine load_mode

  !> Build the admissible stream-function basis for the step diagnostic.
  subroutine build_manufactured_basis(fluid, basis, n)
    type(fluid_pnpn_t), intent(inout) :: fluid
    real(kind=rp), intent(out) :: basis(:, :)
    integer, intent(in) :: n
    type(field_t) :: psi, gx, gy, gz, u, v, w
    integer :: mode, i
    real(kind=rp) :: x, y

    call psi%init(fluid%dm_Xh, 'step_diagnostic_psi')
    call gx%init(fluid%dm_Xh, 'step_diagnostic_gx')
    call gy%init(fluid%dm_Xh, 'step_diagnostic_gy')
    call gz%init(fluid%dm_Xh, 'step_diagnostic_gz')
    call u%init(fluid%dm_Xh, 'step_diagnostic_u')
    call v%init(fluid%dm_Xh, 'step_diagnostic_v')
    call w%init(fluid%dm_Xh, 'step_diagnostic_w')
    do mode = 1, 2
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
       call grad(gx%x, gy%x, gz%x, psi%x, fluid%c_Xh)
       u%x = gy%x
       v%x = -gx%x
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
    call psi%free(); call gx%free(); call gy%free(); call gz%free()
    call u%free(); call v%free(); call w%free()
  end subroutine build_manufactured_basis

end program production_rom_step_diagnostic
