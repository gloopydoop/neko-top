program rom_projection_preconditioner_diagnostic
  !> Matched bounded PHMG/identity solves for the constrained ROM Schur map.
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use fluid_pnpn, only: fluid_pnpn_t
  use field, only: field_t
  use field_math, only: field_rzero, field_col3
  use operators, only: grad, ortho
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2, glsc3, glmin, glmax
  use num_types, only: rp, xp
  use comm, only: pe_rank, pe_size, NEKO_COMM, MPI_EXTRA_PRECISION
  use mpi_f08, only: MPI_Allreduce, MPI_Barrier, MPI_IN_PLACE, MPI_SUM, MPI_Wtime
  use krylov, only: ksp_monitor_t
  use gmres, only: gmres_t
  use precon, only: pc_t
  use identity, only: ident_t
  use rom_projection_operator, only: rom_projection_ax_t
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
  use utils, only: neko_error
  implicit none

  integer, parameter :: n_caps = 4
  integer, parameter :: iteration_caps(n_caps) = [20, 50, 100, 200]
  integer, parameter :: arnoldi_caps(3) = [160, 199, 200]
  integer, parameter :: arnoldi_restart_caps(3) = [200, 201, 400]
  type(simulation_t) :: simulation
  type(json_file) :: parameters
  type(fluid_pnpn_t), pointer :: fluid
  type(rom_projection_ax_t) :: schur
  type(ident_t), target :: identity_pc
  class(pc_t), pointer :: saved_pc
  type(field_t) :: q, rhs_manufactured, rhs_actual, phi, work
  type(field_t) :: pzero, pres, gx, gy, gz, fx, fy, fz
  type(field_t) :: uzero, vzero, wzero, vu, vv, vw
  type(field_t) :: gq_u, gq_v, gq_w, gp_u, gp_v, gp_w
  type(field_t) :: q_constant, s_constant
  real(kind=rp), allocatable :: rhs(:)
  character(len=1024) :: case_file, csv_file
  character(len=32) :: probe_mode
  integer :: argc, n, unit, rhs_id, pc_id, cap_id, current_cap, i
  integer :: mpi_ierr
  integer :: run_caps(5), n_run_caps, rhs_first, rhs_last, pc_first, pc_last
  logical :: found, restart_probe, restart_compare, arnoldi_audit, &
       arnoldi_reorth_compare, arnoldi_reorth_mfg, arnoldi_reorth_actual, &
       arnoldi_restart_mfg, arnoldi_restart_actual, saved_monitor
  integer :: restart_length
  real(kind=rp) :: saved_abs, saved_rel
  integer :: saved_max
  real(kind=rp) :: xmin, xmax, ymin, ymax, zmin, zmax, x, y, z
  real(kind=rp) :: rhs_norm, true_residual, true_relative
  real(kind=rp) :: compatibility, constant_norm, s_constant_norm
  real(kind=rp) :: correction_norm, correction_difference, correction_relative
  real(kind=rp) :: ksp_residual, ksp_volume_relative, iterations, convergence
  real(kind=rp) :: rho_min, rho_max
  character(len=1200) :: arnoldi_csv
  type(ksp_monitor_t) :: result
  type(c_ptr) :: event
  real(kind=rp), allocatable :: native_x(:, :), native_true(:), native_ksp(:), native_true_abs(:)
  real(kind=rp), allocatable :: audit_x(:, :), audit_true(:), audit_monitor(:)
  real(kind=rp), allocatable :: audit_true_abs(:), audit_grad_diff(:), audit_bc_idempotence(:)
  real(kind=rp), allocatable :: audit_bc_correction(:), audit_continuity_abs(:), &
       audit_continuity_rel(:), audit_grad_relative(:)
  integer, allocatable :: audit_cycle(:)
  real(kind=rp), allocatable :: audit_orth_fro(:), audit_orth_max(:)
  real(kind=rp), allocatable :: audit_relation(:), audit_repeat(:), audit_vnext(:)
  real(kind=rp), allocatable :: reorth_x(:, :), reorth_true(:), reorth_monitor(:)
  real(kind=rp), allocatable :: reorth_true_abs(:), reorth_grad_diff(:), reorth_bc_idempotence(:)
  real(kind=rp), allocatable :: reorth_bc_correction(:), reorth_continuity_abs(:), &
       reorth_continuity_rel(:), reorth_grad_relative(:)
  integer, allocatable :: reorth_cycle(:)
  real(kind=rp), allocatable :: reorth_orth_fro(:), reorth_orth_max(:)
  real(kind=rp), allocatable :: reorth_relation(:), reorth_repeat(:), reorth_vnext(:)
  real(kind=rp), allocatable :: reproduction_solution(:), reproduction_true(:), reproduction_monitor(:)
  real(kind=rp), allocatable :: reorth_solution_difference(:), reorth_true_difference(:)
  real(kind=rp), allocatable :: reorth_monitor_difference(:)
  real(kind=rp) :: clone_solution_error, clone_true_error, clone_monitor_error
  real(kind=rp) :: reproduction_tolerance
  real(kind=rp) :: control_seconds, reorth_seconds, audit_timer
  real(kind=rp) :: failed_velocity_difference, failed_native_cont_abs, failed_native_cont_rel
  real(kind=rp) :: failed_native_bc, failed_clone_cont_abs, failed_clone_cont_rel, failed_clone_bc
  real(kind=rp) :: failed_grad_abs, failed_grad_rel
  integer :: audit_sample_caps(3)

  argc = command_argument_count()
  if (argc < 1) call neko_error( &
       'Usage: ROM_rugby_ball_projection_preconditioner_diagnostic CASE [CSV] '// &
       '[restart_probe|restart30|restart60|restart120|restart200|restart200_caps|'// &
       'arnoldi_audit|arnoldi_reorth_compare|arnoldi_reorth_mfg|arnoldi_reorth_actual|'// &
       'arnoldi_restart_mfg|arnoldi_restart_actual]')
  call get_command_argument(1, case_file)
  if (argc >= 2) then
     call get_command_argument(2, csv_file)
  else
     csv_file = 'projection_preconditioner_diagnostic.csv'
  end if
  probe_mode = ''
  if (argc >= 3) call get_command_argument(3, probe_mode)
  restart_probe = trim(probe_mode) == 'restart_probe'
  arnoldi_reorth_compare = trim(probe_mode) == 'arnoldi_reorth_compare'
  arnoldi_reorth_mfg = trim(probe_mode) == 'arnoldi_reorth_mfg'
  arnoldi_reorth_actual = trim(probe_mode) == 'arnoldi_reorth_actual'
  arnoldi_restart_mfg = trim(probe_mode) == 'arnoldi_restart_mfg'
  arnoldi_restart_actual = trim(probe_mode) == 'arnoldi_restart_actual'
  arnoldi_audit = trim(probe_mode) == 'arnoldi_audit' .or. arnoldi_reorth_compare .or. &
       arnoldi_reorth_mfg .or. arnoldi_reorth_actual .or. &
       arnoldi_restart_mfg .or. arnoldi_restart_actual
  restart_compare = trim(probe_mode) == 'restart30' .or. &
       trim(probe_mode) == 'restart60' .or. trim(probe_mode) == 'restart120' .or. &
       trim(probe_mode) == 'restart200' .or. &
       trim(probe_mode) == 'restart200_caps' .or. arnoldi_audit
  if (argc >= 3 .and. .not. restart_probe .and. .not. restart_compare) call neko_error( &
       'Unknown projection preconditioner diagnostic mode')
  run_caps = [20, 50, 100, 200, 201]
  n_run_caps = n_caps
  rhs_first = 1; rhs_last = 2; pc_first = 1; pc_last = 2
  if (restart_probe) then
     run_caps(1:3) = [29, 30, 31]
     n_run_caps = 3
     rhs_first = 2; rhs_last = 2
     pc_first = 1; pc_last = 1
  else if (restart_compare) then
     restart_length = 30
     if (trim(probe_mode) == 'restart60') restart_length = 60
     if (trim(probe_mode) == 'restart120') restart_length = 120
     if (trim(probe_mode) == 'restart200' .or. &
          trim(probe_mode) == 'restart200_caps' .or. arnoldi_audit) restart_length = 200
     run_caps(1:3) = [60, 120, 200]
     n_run_caps = 3
     rhs_first = 2; rhs_last = 2
     pc_first = 1; pc_last = 1
     if (trim(probe_mode) == 'restart200') then
        run_caps(1) = 200
        n_run_caps = 1
     else if (trim(probe_mode) == 'restart200_caps') then
        run_caps = [120, 160, 199, 200, 201]
        n_run_caps = 5
     else if (arnoldi_reorth_mfg .or. arnoldi_reorth_actual) then
        run_caps(1:4) = [160, 199, 200, 201]
        n_run_caps = 4
        pc_first = 1; pc_last = 1
        if (arnoldi_reorth_mfg) then
           rhs_first = 1; rhs_last = 1
        else
           rhs_first = 2; rhs_last = 2
        end if
     else if (arnoldi_restart_mfg .or. arnoldi_restart_actual) then
        run_caps(1:3) = arnoldi_restart_caps
        n_run_caps = 3
        pc_first = 1; pc_last = 1
        if (arnoldi_restart_mfg) then
           rhs_first = 1; rhs_last = 1
        else
           rhs_first = 2; rhs_last = 2
        end if
     else if (arnoldi_audit) then
        run_caps(1:3) = [160, 199, 200]
        n_run_caps = 3
     end if
  end if

  call neko_init()
  if (arnoldi_audit .and. pe_size /= 1) call neko_error( &
       'arnoldi_audit requires exactly one MPI rank')
  call neko_top_register_types()
  parameters = json_read_file(trim(case_file))
  call parameters%update('unsteady', .false., found)
  call user_setup(simulation%neko_case%user)
  call simulation%init(parameters)
  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('Projection preconditioner diagnostic requires Pn/Pn fluid')
  end select
  n = fluid%u%dof%size()

  call q%init(fluid%dm_Xh, 'preconditioner_q')
  call rhs_manufactured%init(fluid%dm_Xh, 'preconditioner_rhs_mfg')
  call rhs_actual%init(fluid%dm_Xh, 'preconditioner_rhs_actual')
  call phi%init(fluid%dm_Xh, 'preconditioner_phi')
  call work%init(fluid%dm_Xh, 'preconditioner_work')
  call pzero%init(fluid%dm_Xh, 'preconditioner_pzero')
  call pres%init(fluid%dm_Xh, 'preconditioner_pres')
  call gx%init(fluid%dm_Xh, 'preconditioner_gx')
  call gy%init(fluid%dm_Xh, 'preconditioner_gy')
  call gz%init(fluid%dm_Xh, 'preconditioner_gz')
  call fx%init(fluid%dm_Xh, 'preconditioner_fx')
  call fy%init(fluid%dm_Xh, 'preconditioner_fy')
  call fz%init(fluid%dm_Xh, 'preconditioner_fz')
  call uzero%init(fluid%dm_Xh, 'preconditioner_uzero')
  call vzero%init(fluid%dm_Xh, 'preconditioner_vzero')
  call wzero%init(fluid%dm_Xh, 'preconditioner_wzero')
  call vu%init(fluid%dm_Xh, 'preconditioner_vu')
  call vv%init(fluid%dm_Xh, 'preconditioner_vv')
  call vw%init(fluid%dm_Xh, 'preconditioner_vw')
  call gq_u%init(fluid%dm_Xh, 'preconditioner_gq_u')
  call gq_v%init(fluid%dm_Xh, 'preconditioner_gq_v')
  call gq_w%init(fluid%dm_Xh, 'preconditioner_gq_w')
  call gp_u%init(fluid%dm_Xh, 'preconditioner_gp_u')
  call gp_v%init(fluid%dm_Xh, 'preconditioner_gp_v')
  call gp_w%init(fluid%dm_Xh, 'preconditioner_gp_w')
  call q_constant%init(fluid%dm_Xh, 'preconditioner_q_constant')
  call s_constant%init(fluid%dm_Xh, 'preconditioner_s_constant')
  allocate(rhs(n))
  if (arnoldi_audit) then
     audit_sample_caps = arnoldi_caps
     if (arnoldi_restart_mfg .or. arnoldi_restart_actual) &
          audit_sample_caps = arnoldi_restart_caps
     allocate(native_x(n,4), native_true(4), native_ksp(4), native_true_abs(4))
     allocate(audit_x(n,3), audit_true(3), audit_monitor(3), audit_true_abs(3))
     allocate(audit_grad_diff(3), audit_bc_idempotence(3))
     if (arnoldi_restart_mfg .or. arnoldi_restart_actual) then
        allocate(audit_bc_correction(3), audit_continuity_abs(3), &
             audit_continuity_rel(3), audit_grad_relative(3), audit_cycle(3))
     end if
     allocate(audit_orth_fro(3), audit_orth_max(3), audit_relation(3), &
          audit_repeat(3), audit_vnext(3))
     allocate(reproduction_solution(3), reproduction_true(3), reproduction_monitor(3))
     if (arnoldi_reorth_compare .or. arnoldi_reorth_mfg .or. arnoldi_reorth_actual) then
        allocate(reorth_x(n,3), reorth_true(3), reorth_monitor(3), reorth_true_abs(3))
        allocate(reorth_grad_diff(3), reorth_bc_idempotence(3))
        allocate(reorth_orth_fro(3), reorth_orth_max(3), reorth_relation(3), &
             reorth_repeat(3), reorth_vnext(3))
        allocate(reorth_solution_difference(3), reorth_true_difference(3), &
             reorth_monitor_difference(3))
     end if
     if (arnoldi_restart_mfg .or. arnoldi_restart_actual) then
        allocate(reorth_x(n,3), reorth_true(3), reorth_monitor(3), reorth_true_abs(3))
        allocate(reorth_grad_diff(3), reorth_bc_idempotence(3), reorth_bc_correction(3), &
             reorth_continuity_abs(3), reorth_continuity_rel(3), reorth_grad_relative(3), &
             reorth_cycle(3))
        allocate(reorth_orth_fro(3), reorth_orth_max(3), reorth_relation(3), &
             reorth_repeat(3), reorth_vnext(3))
        allocate(reorth_solution_difference(3), reorth_true_difference(3), &
             reorth_monitor_difference(3))
     end if
     native_x = 0.0_rp; native_true = 0.0_rp; native_ksp = 0.0_rp; native_true_abs = 0.0_rp
  end if
  call field_rzero(pzero); call field_rzero(uzero)
  call field_rzero(vzero); call field_rzero(wzero)

  xmin = glmin(fluid%u%dof%x%x, n); xmax = glmax(fluid%u%dof%x%x, n)
  ymin = glmin(fluid%u%dof%y%x, n); ymax = glmax(fluid%u%dof%y%x, n)
  zmin = glmin(fluid%u%dof%z%x, n); zmax = glmax(fluid%u%dof%z%x, n)
  rho_min = glmin(fluid%rho%x, n); rho_max = glmax(fluid%rho%x, n)
  call schur%init(fluid, q, pzero, pres, gx, gy, gz, fx, fy, fz, &
       uzero, vzero, wzero, mass_weighted_velocity_average = .true.)

  ! Restricted smooth pressure q gives the manufactured compatible rhs b=S q.
  do i = 1, n
     x = normalize(fluid%u%dof%x%x(i,1,1,1), xmin, xmax)
     y = normalize(fluid%u%dof%y%x(i,1,1,1), ymin, ymax)
     z = normalize(fluid%u%dof%z%x(i,1,1,1), zmin, zmax)
     q%x(i,1,1,1) = x + 0.3_rp*y*y - 0.2_rp*z
  end do
  call continuous_scalar(q)
  call fluid%bcs_prs_projector%apply(q%x, n)
  if (.not. fluid%prs_dirichlet) call ortho(q%x, fluid%glb_n_points, n)
  call schur%compute(rhs_manufactured%x, q%x, fluid%c_Xh, fluid%msh, fluid%Xh)
  call assemble_pressure(rhs_manufactured)
  call apply_gradient_map(q, gq_u, gq_v, gq_w)

  ! An independent homogeneous no-slip continuous field gives b=C v.
  call make_velocity_probe()
  call apply_continuity(vu, vv, vw, rhs_actual)

  ! Measure pressure-nullspace and RHS compatibility without a solve.
  q_constant%x = 1.0_rp
  call continuous_scalar(q_constant)
  call fluid%bcs_prs_projector%apply(q_constant%x, n)
  call schur%compute(s_constant%x, q_constant%x, fluid%c_Xh, fluid%msh, fluid%Xh)
  call assemble_pressure(s_constant)
  constant_norm = pressure_norm(q_constant)
  s_constant_norm = pressure_norm(s_constant)

  saved_abs = fluid%ksp_prs%abs_tol
  saved_rel = fluid%ksp_prs%rel_tol
  saved_max = fluid%ksp_prs%max_iter
  saved_monitor = fluid%ksp_prs%monitor
  saved_pc => fluid%ksp_prs%M
  if (restart_compare) then
     select type (solver => fluid%ksp_prs)
     type is (gmres_t)
        solver%lgmres = restart_length
        call solver%init(n, saved_max, M=saved_pc, rel_tol=saved_rel, &
             abs_tol=saved_abs, monitor=saved_monitor)
     class default
        call neko_error('Restart comparison requires the CPU GMRES solver')
     end select
     saved_pc => fluid%ksp_prs%M
  end if
  fluid%ksp_prs%abs_tol = 1.0e-12_rp
  fluid%ksp_prs%rel_tol = 0.0_rp
  call fluid%pc_prs%update()

  if (pe_rank == 0) then
     open(newunit=unit, file=trim(csv_file), status='replace', action='write')
     write(unit, '(A)') 'cap,rhs,preconditioner,rhs_norm,constant_mode_norm,' // &
          'S_constant_norm,rhs_constant_compatibility,ksp_residual,' // &
          'ksp_relative_volume_normalized,true_residual,true_relative,' // &
          'iterations,converged,' // &
          'correction_norm,Gphi_minus_Gq_norm,Gphi_relative_error'
  end if

  do cap_id = 1, n_run_caps
   current_cap = run_caps(cap_id)
   fluid%ksp_prs%max_iter = current_cap
   do rhs_id = rhs_first, rhs_last
     if (rhs_id == 1) then
        work%x = rhs_manufactured%x
     else
        work%x = rhs_actual%x
     end if
     call fluid%bcs_prs_projector%apply(work%x, n)
     if (.not. fluid%prs_dirichlet) call ortho(work%x, fluid%glb_n_points, n)
     rhs_norm = pressure_norm(work)
     rhs = reshape(work%x, [n])
     compatibility = abs(glsc3(q_constant%x, fluid%c_Xh%mult, &
          work%x, n)) / max(constant_norm, tiny(1.0_rp))

     do pc_id = pc_first, pc_last
        if (pc_id == 1) then
           fluid%ksp_prs%M => saved_pc
           call fluid%pc_prs%update()
        else
           fluid%ksp_prs%M => identity_pc
        end if
        call field_rzero(phi)
        result = fluid%ksp_prs%solve(schur, phi, rhs, n, fluid%c_Xh, &
             fluid%bcs_prs_projector, fluid%gs_Xh, niter=current_cap)
        ksp_residual = result%res_final
        iterations = real(result%iter, rp)
        convergence = merge(1.0_rp, 0.0_rp, result%converged)
        call true_schur_residual(phi, rhs, true_residual)
        true_relative = true_residual / max(rhs_norm, tiny(1.0_rp))
        if (arnoldi_audit) then
           native_x(:,cap_id) = reshape(phi%x, [n])
           native_true(cap_id) = true_relative
           native_true_abs(cap_id) = true_residual
           native_ksp(cap_id) = ksp_residual
        end if
        ! GMRES returns its pressure norm divided by sqrt(volume).  Convert
        ! back to a dimensionless RHS-relative residual before comparison.
        ksp_volume_relative = ksp_residual * sqrt(fluid%c_Xh%volume) / &
             max(rhs_norm, tiny(1.0_rp))

        correction_norm = 0.0_rp
        correction_difference = -1.0_rp
        correction_relative = -1.0_rp
        if (rhs_id == 1) then
           call apply_gradient_map(phi, gp_u, gp_v, gp_w)
           correction_norm = velocity_norm(gp_u, gp_v, gp_w)
           gp_u%x = gp_u%x - gq_u%x
           gp_v%x = gp_v%x - gq_v%x
           gp_w%x = gp_w%x - gq_w%x
           correction_difference = velocity_norm(gp_u, gp_v, gp_w)
           correction_relative = correction_difference / &
                max(velocity_norm(gq_u, gq_v, gq_w), tiny(1.0_rp))
        else
           call apply_gradient_map(phi, gp_u, gp_v, gp_w)
           correction_norm = velocity_norm(gp_u, gp_v, gp_w)
        end if
        if (pe_rank == 0) then
           write(unit, '(I0,",",I0,",",I0,13(",",ES24.16))') &
                current_cap, rhs_id, pc_id, &
                rhs_norm, constant_norm, s_constant_norm, compatibility, &
                ksp_residual, ksp_volume_relative, &
                true_residual, true_relative, iterations, convergence, &
                correction_norm, correction_difference, correction_relative
        end if
     end do
   end do
  end do

  if (pe_rank == 0) close(unit)
  if (arnoldi_audit) then
     audit_timer = MPI_Wtime()
     if (arnoldi_restart_mfg .or. arnoldi_restart_actual) then
        call instrumented_restarted_arnoldi(rhs, audit_x, audit_true, audit_monitor, &
             audit_true_abs, audit_grad_diff, audit_grad_relative, audit_bc_correction, &
             audit_continuity_abs, audit_continuity_rel, audit_cycle, audit_orth_fro, &
             audit_orth_max, audit_relation, audit_repeat, audit_vnext, .false., &
             arnoldi_restart_mfg)
     else
        call instrumented_arnoldi(rhs, audit_x, audit_true, audit_monitor, &
             audit_true_abs, audit_grad_diff, audit_bc_idempotence, &
             audit_orth_fro, audit_orth_max, audit_relation, audit_repeat, &
             audit_vnext, .false., arnoldi_reorth_mfg)
     end if
     control_seconds = MPI_Wtime()-audit_timer
     reproduction_tolerance = max(100.0_rp*epsilon(1.0_rp), 5.0e-8_rp)
     do i = 1, 3
        reproduction_solution(i) = weighted_norm(native_x(:,i)-audit_x(:,i)) / &
             max(weighted_norm(native_x(:,i)), tiny(1.0_rp))
        reproduction_true(i) = abs(native_true(i)-audit_true(i))
        reproduction_monitor(i) = abs(native_ksp(i)*sqrt(fluid%c_Xh%volume) / &
             max(rhs_norm,tiny(1.0_rp))-audit_monitor(i))
     end do
     if (maxval(reproduction_solution) > reproduction_tolerance .or. &
          maxval(reproduction_true) > reproduction_tolerance .or. &
          maxval(reproduction_monitor) > reproduction_tolerance) then
        if (pe_rank == 0) then
           write(*,'(A,ES12.4)') 'Arnoldi clone/native scale-aware gate=', &
                reproduction_tolerance
           do i = 1, 3
              write(*,'(A,I0,7(1X,ES16.8))') 'clone/native comparison cap ', &
                audit_sample_caps(i), reproduction_solution(i), reproduction_true(i), &
                   reproduction_monitor(i), native_true(i), audit_true(i), &
                   native_ksp(i)*sqrt(fluid%c_Xh%volume)/max(rhs_norm,tiny(1.0_rp)), &
                   audit_monitor(i)
           end do
        end if
        if (arnoldi_restart_actual .and. reproduction_solution(3) > reproduction_tolerance) then
           call projected_velocity_difference(native_x(:,3),audit_x(:,3), &
                failed_velocity_difference)
           call restarted_velocity_metrics(native_x(:,3),rhs,.false.,failed_grad_abs, &
                failed_grad_rel,failed_native_bc,failed_native_cont_abs,failed_native_cont_rel)
           call restarted_velocity_metrics(audit_x(:,3),rhs,.false.,failed_grad_abs, &
                failed_grad_rel,failed_clone_bc,failed_clone_cont_abs,failed_clone_cont_rel)
           if (pe_rank == 0) then
              write(*,'(A,*(1X,ES16.8))') 'cap400 failed-control diagnostics: '// &
                   'pressure_solution_relative,projected_velocity_gradient_relative, '// &
                   'native_true_abs/native_rel,clone_true_abs/clone_rel, '// &
                   'native_monitor_rel,clone_monitor_rel,native_cont_abs/rel, '// &
                   'clone_cont_abs/rel,native_BC_correction,clone_BC_correction', &
                   reproduction_solution(3), failed_velocity_difference, &
                   native_true_abs(3),native_true(3),audit_true_abs(3),audit_true(3), &
                   native_ksp(3)*sqrt(fluid%c_Xh%volume)/max(rhs_norm,tiny(1.0_rp)), &
                   audit_monitor(3),failed_native_cont_abs,failed_native_cont_rel, &
                   failed_clone_cont_abs,failed_clone_cont_rel,failed_native_bc,failed_clone_bc
           end if
        end if
        call neko_error('Arnoldi audit clone differs from native GMRES; internal metrics withheld')
     end if
     if (arnoldi_restart_mfg .or. arnoldi_restart_actual) then
        audit_timer = MPI_Wtime()
        call instrumented_restarted_arnoldi(rhs, reorth_x, reorth_true, reorth_monitor, &
             reorth_true_abs, reorth_grad_diff, reorth_grad_relative, reorth_bc_correction, &
             reorth_continuity_abs, reorth_continuity_rel, reorth_cycle, reorth_orth_fro, &
             reorth_orth_max, reorth_relation, reorth_repeat, reorth_vnext, .true., &
             arnoldi_restart_mfg)
        reorth_seconds = MPI_Wtime()-audit_timer
        do i = 1, 3
           reorth_solution_difference(i) = weighted_norm(native_x(:,i)-reorth_x(:,i)) / &
                max(weighted_norm(native_x(:,i)),tiny(1.0_rp))
           reorth_true_difference(i) = abs(native_true(i)-reorth_true(i))
           reorth_monitor_difference(i) = abs(native_ksp(i)*sqrt(fluid%c_Xh%volume) / &
                max(rhs_norm,tiny(1.0_rp))-reorth_monitor(i))
        end do
     else if (arnoldi_reorth_compare .or. arnoldi_reorth_mfg .or. arnoldi_reorth_actual) then
        audit_timer = MPI_Wtime()
        call instrumented_arnoldi(rhs, reorth_x, reorth_true, reorth_monitor, &
             reorth_true_abs, reorth_grad_diff, reorth_bc_idempotence, &
             reorth_orth_fro, reorth_orth_max, reorth_relation, reorth_repeat, &
             reorth_vnext, .true., arnoldi_reorth_mfg)
        reorth_seconds = MPI_Wtime()-audit_timer
        do i = 1, 3
           reorth_solution_difference(i) = weighted_norm(native_x(:,i)-reorth_x(:,i)) / &
                max(weighted_norm(native_x(:,i)),tiny(1.0_rp))
           reorth_true_difference(i) = abs(native_true(i)-reorth_true(i))
           reorth_monitor_difference(i) = abs(native_ksp(i)*sqrt(fluid%c_Xh%volume) / &
                max(rhs_norm,tiny(1.0_rp))-reorth_monitor(i))
        end do
     end if
     arnoldi_csv = trim(csv_file)//'.arnoldi.csv'
     if (pe_rank == 0) then
        open(newunit=unit, file=trim(arnoldi_csv), status='replace', action='write')
        if (arnoldi_restart_mfg .or. arnoldi_restart_actual) then
           write(unit, '(A)') 'variant,cap,cycle,solution_relative_error_vs_native,'// &
                'native_true_absolute,native_true_relative,clone_true_absolute,clone_true_relative,'// &
                'true_relative_difference,native_monitor_relative,clone_monitor_relative,'// &
                'monitor_minus_true_relative,cycle_weighted_gram_defect_fro,'// &
                'cycle_weighted_gram_max_offdiag,arnoldi_relation_relative_defect,'// &
                'S_repeatability_relative,weak_continuity_absolute,weak_continuity_relative,'// &
                'manufactured_Gphi_minus_Gq_over_Gq,projected_velocity_BC_correction_norm,clone_wall_seconds'
        else
        write(unit, '(A)') 'variant,iteration,solution_relative_error_vs_native,'// &
             'native_true_absolute,native_true_relative,instrumented_true_absolute,'// &
             'instrumented_true_relative,true_relative_difference,'// &
             'native_monitor_volume_relative,instrumented_monitor_volume_relative,'// &
             'monitor_relative_difference,monitor_minus_true_relative,'// &
             'orthogonality_frobenius,orthogonality_max_off_diagonal,'// &
             'sampled_arnoldi_relation_relative_defect,S_repeatability_relative,'// &
             'v_next_weighted_norm,manufactured_gradient_difference,velocity_BC_projector_idempotence,'// &
             'clone_wall_seconds'
        end if
        do i = 1, 3
           clone_solution_error = reproduction_solution(i)
           clone_true_error = reproduction_true(i)
           clone_monitor_error = reproduction_monitor(i)
           if (arnoldi_restart_mfg .or. arnoldi_restart_actual) then
              write(unit, '(A,",",I0,",",I0,18(",",ES24.16))') 'classical_gs_control', &
                   audit_sample_caps(i), audit_cycle(i), clone_solution_error, &
                   native_true_abs(i), native_true(i), audit_true_abs(i), audit_true(i), &
                   clone_true_error, native_ksp(i)*sqrt(fluid%c_Xh%volume)/ &
                   max(rhs_norm,tiny(1.0_rp)), audit_monitor(i), &
                   audit_monitor(i)-audit_true(i), audit_orth_fro(i), audit_orth_max(i), &
                   audit_relation(i), audit_repeat(i), audit_continuity_abs(i), &
                   audit_continuity_rel(i), audit_grad_relative(i), audit_bc_correction(i), &
                   control_seconds
           else
           write(unit, '(A,",",I0,18(",",ES24.16))') 'classical_gs_control', arnoldi_caps(i), &
                clone_solution_error, native_true_abs(i), native_true(i), audit_true_abs(i), &
                audit_true(i), clone_true_error, &
                native_ksp(i)*sqrt(fluid%c_Xh%volume)/max(rhs_norm,tiny(1.0_rp)), &
                audit_monitor(i), clone_monitor_error, &
                audit_monitor(i)-audit_true(i), audit_orth_fro(i), audit_orth_max(i), &
                audit_relation(i), audit_repeat(i), audit_vnext(i), audit_grad_diff(i), &
                audit_bc_idempotence(i), control_seconds
           end if
           if (arnoldi_restart_mfg .or. arnoldi_restart_actual) then
              write(unit, '(A,",",I0,",",I0,18(",",ES24.16))') &
                   'second_pass_reorthogonalized', audit_sample_caps(i), reorth_cycle(i), &
                   reorth_solution_difference(i), native_true_abs(i), native_true(i), &
                   reorth_true_abs(i), reorth_true(i), reorth_true_difference(i), &
                   native_ksp(i)*sqrt(fluid%c_Xh%volume)/max(rhs_norm,tiny(1.0_rp)), &
                   reorth_monitor(i), reorth_monitor(i)-reorth_true(i), reorth_orth_fro(i), &
                   reorth_orth_max(i), reorth_relation(i), reorth_repeat(i), &
                   reorth_continuity_abs(i), reorth_continuity_rel(i), reorth_grad_relative(i), &
                   reorth_bc_correction(i), reorth_seconds
           else if (arnoldi_reorth_compare .or. arnoldi_reorth_mfg .or. arnoldi_reorth_actual) then
              write(unit, '(A,",",I0,18(",",ES24.16))') &
                   'second_pass_reorthogonalized', arnoldi_caps(i), &
                   reorth_solution_difference(i), native_true_abs(i), native_true(i), &
                   reorth_true_abs(i), reorth_true(i), &
                   reorth_true_difference(i), &
                   native_ksp(i)*sqrt(fluid%c_Xh%volume)/max(rhs_norm,tiny(1.0_rp)), &
                   reorth_monitor(i), reorth_monitor_difference(i), &
                   reorth_monitor(i)-reorth_true(i), reorth_orth_fro(i), &
                   reorth_orth_max(i), reorth_relation(i), reorth_repeat(i), &
                   reorth_vnext(i), reorth_grad_diff(i), reorth_bc_idempotence(i), reorth_seconds
           end if
        end do
        close(unit)
     end if
     call MPI_Barrier(NEKO_COMM, mpi_ierr)
  end if
  fluid%ksp_prs%M => saved_pc
  fluid%ksp_prs%abs_tol = saved_abs
  fluid%ksp_prs%rel_tol = saved_rel
  fluid%ksp_prs%max_iter = saved_max
  call schur%free()
  call q%free(); call rhs_manufactured%free(); call rhs_actual%free()
  call phi%free(); call work%free(); call pzero%free(); call pres%free()
  call gx%free(); call gy%free(); call gz%free(); call fx%free(); call fy%free()
  call fz%free(); call uzero%free(); call vzero%free(); call wzero%free()
  call vu%free(); call vv%free(); call vw%free(); call gq_u%free()
  call gq_v%free(); call gq_w%free(); call gp_u%free(); call gp_v%free()
  call gp_w%free(); call q_constant%free(); call s_constant%free()
  deallocate(rhs)
  if (pe_rank == 0) then
     write(*, '(A,L1)') 'pressure Dirichlet restriction=', fluid%prs_dirichlet
     write(*, '(A,*(1X,I0))') 'matched iteration caps=', run_caps(1:n_run_caps)
     if (restart_compare) write(*, '(A,I0)') &
          'GMRES restart dimension=', restart_length
     write(*, '(A,ES12.4)') 'density min=', rho_min
     write(*, '(A,ES12.4)') 'density max=', rho_max
     write(*, '(A,A)') 'preconditioner CSV: ', trim(csv_file)
  end if
  call neko_finalize()

contains

  !> Test-only right-preconditioned Arnoldi clone.  It keeps the raw
  !! Hessenberg coefficients and one extra basis vector that the production
  !! full-cycle GMRES path deliberately does not retain.
  subroutine instrumented_arnoldi(b, x_samples, true_samples, monitor_samples, &
       true_abs_samples, grad_diff_samples, bc_idempotence_samples, &
       orth_fro_samples, orth_max_samples, relation_samples, repeat_samples, &
       vnext_samples, reorthogonalize, manufactured)
    real(kind=rp), intent(in) :: b(:)
    real(kind=rp), intent(out) :: x_samples(:, :), true_samples(:), monitor_samples(:)
    real(kind=rp), intent(out) :: true_abs_samples(:), grad_diff_samples(:), bc_idempotence_samples(:)
    real(kind=rp), intent(out) :: orth_fro_samples(:), orth_max_samples(:)
    real(kind=rp), intent(out) :: relation_samples(:), repeat_samples(:), vnext_samples(:)
    logical, intent(in) :: reorthogonalize
    logical, intent(in) :: manufactured
    real(kind=rp), allocatable :: v(:, :), z(:, :), w(:), sol(:), sz(:), sz_again(:)
    real(kind=xp), allocatable :: wplus(:)
    real(kind=xp), allocatable :: hbar(:, :), hr(:, :), cosine(:), sine(:), gamma(:), y(:)
    real(kind=xp), allocatable :: gram(:, :)
    real(kind=xp) :: alpha, norm0, lr, tmp, dot_value, relation_norm, max_offdiag
    real(kind=xp) :: local_sum
    real(kind=xp) :: hlocal(200)
    real(kind=xp) :: hsecond(200)
    real(kind=rp) :: rhs_norm
    integer :: j, i, k, sample_index, dof_index
    integer, parameter :: sample_steps(3) = [160, 199, 200]

    allocate(v(n,201), z(n,200), w(n), wplus(n), sol(n), sz(n), sz_again(n))
    allocate(hbar(201,200), hr(201,200), cosine(200), sine(200), gamma(201), y(200))
    v = 0.0_rp; z = 0.0_rp; hbar = 0.0_xp; hr = 0.0_xp
    cosine = 1.0_xp; sine = 1.0_xp; gamma = 0.0_xp
    rhs_norm = weighted_norm(b)
    norm0 = sqrt(weighted_dot(b,b))
    if (norm0 <= tiny(1.0_xp)) call neko_error('Arnoldi audit received zero RHS')
    gamma(1) = norm0
    v(:,1) = b * real(1.0_xp/norm0, rp)
    x_samples = 0.0_rp; true_samples = 0.0_rp; monitor_samples = 0.0_rp
    true_abs_samples = 0.0_rp; grad_diff_samples = -1.0_rp; bc_idempotence_samples = 0.0_rp
    orth_fro_samples = 0.0_rp; orth_max_samples = 0.0_rp
    relation_samples = 0.0_rp; repeat_samples = 0.0_rp; vnext_samples = 0.0_rp

    do j = 1, 200
       call fluid%ksp_prs%M%solve(z(:,j), v(:,j), n)
       call apply_s_array(z(:,j), w)
       wplus = 0.0_xp
       hlocal = 0.0_xp
       do i = 1, j
          do k = 1, n
             hlocal(i) = hlocal(i) + real(w(k),xp)*real(v(k,i),xp)* &
                  real(fluid%c_Xh%mult(k,1,1,1),xp)
          end do
       end do
       call MPI_Allreduce(MPI_IN_PLACE,hlocal,j,MPI_EXTRA_PRECISION,MPI_SUM,NEKO_COMM,mpi_ierr)
       do i = 1, j
          hbar(i,j) = hlocal(i)
          wplus = wplus - hlocal(i)*real(v(:,i),xp)
       end do
       w = real(real(w,xp)+wplus,rp)
       if (reorthogonalize) then
          hsecond = 0.0_xp
          do i = 1, j
             do k = 1, n
                hsecond(i) = hsecond(i) + real(w(k),xp)*real(v(k,i),xp)* &
                     real(fluid%c_Xh%mult(k,1,1,1),xp)
             end do
          end do
          call MPI_Allreduce(MPI_IN_PLACE,hsecond,j,MPI_EXTRA_PRECISION,MPI_SUM,NEKO_COMM,mpi_ierr)
          wplus = 0.0_xp
          do i = 1, j
             hbar(i,j) = hbar(i,j) + hsecond(i)
             wplus = wplus - hsecond(i)*real(v(:,i),xp)
          end do
          w = real(real(w,xp)+wplus,rp)
       end if
       alpha = real(weighted_norm(w), xp)
       hbar(j+1,j) = alpha
       if (alpha > tiny(1.0_xp)) v(:,j+1) = w * real(1.0_xp/alpha, rp)

       hr(1:j+1,j) = hbar(1:j+1,j)
       do i = 1, j-1
          tmp = hr(i,j)
          hr(i,j) = cosine(i)*tmp + sine(i)*hr(i+1,j)
          hr(i+1,j) = -sine(i)*tmp + cosine(i)*hr(i+1,j)
       end do
       lr = sqrt(hr(j,j)*hr(j,j) + alpha*alpha)
       if (lr <= tiny(1.0_xp)) exit
       cosine(j) = hr(j,j)/lr
       sine(j) = alpha/lr
       hr(j,j) = lr
       gamma(j+1) = -sine(j)*gamma(j)
       gamma(j) = cosine(j)*gamma(j)

       sample_index = 0
       do i = 1, size(sample_steps)
          if (j == sample_steps(i)) sample_index = i
       end do
       if (sample_index > 0) then
          y(1:j) = gamma(1:j)
          do k = j, 1, -1
             tmp = y(k)
             do i = j, k+1, -1
                tmp = tmp - hr(k,i)*y(i)
             end do
             y(k) = tmp/hr(k,k)
          end do
          sol = 0.0_rp
          do i = 1, n
             dot_value = 0.0_xp
             do k = 1, j
                dot_value = dot_value + y(k)*real(z(i,k),xp)
             end do
             sol(i) = real(dot_value,rp)
          end do
          x_samples(:,sample_index) = sol
          call true_residual_metrics(sol,b,true_abs_samples(sample_index), &
               true_samples(sample_index))
          monitor_samples(sample_index) = abs(gamma(j+1)) / max(rhs_norm,tiny(1.0_rp))
          vnext_samples(sample_index) = weighted_norm(v(:,j+1))
          phi%x = reshape(sol,shape(phi%x))
          call apply_gradient_map(phi,gp_u,gp_v,gp_w)
          fx%x = gp_u%x; fy%x = gp_v%x; fz%x = gp_w%x
          call fluid%bcs_vel_projector%apply(gp_u%x,gp_v%x,gp_w%x,n)
          gp_u%x = gp_u%x - fx%x
          gp_v%x = gp_v%x - fy%x
          gp_w%x = gp_w%x - fz%x
          bc_idempotence_samples(sample_index) = velocity_norm(gp_u,gp_v,gp_w)
          if (manufactured) then
             gp_u%x = fx%x-gq_u%x
             gp_v%x = fy%x-gq_v%x
             gp_w%x = fz%x-gq_w%x
             grad_diff_samples(sample_index) = velocity_norm(gp_u,gp_v,gp_w)
          end if

          allocate(gram(j+1,j+1)); gram = 0.0_xp
          do k = 1, j+1
             do i = 1, j+1
                local_sum = 0.0_xp
                do dof_index = 1, n
                   local_sum = local_sum + real(v(dof_index,i),xp)* &
                        real(v(dof_index,k),xp)* &
                        real(fluid%c_Xh%mult(dof_index,1,1,1),xp)
                end do
                gram(i,k) = local_sum
             end do
          end do
          call MPI_Allreduce(MPI_IN_PLACE,gram(1,1),(j+1)*(j+1), &
               MPI_EXTRA_PRECISION,MPI_SUM,NEKO_COMM,mpi_ierr)
          orth_fro_samples(sample_index) = 0.0_rp
          max_offdiag = 0.0_xp
          do k = 1, j+1
             do i = 1, j+1
                dot_value = gram(i,k)
                if (i == k) dot_value = dot_value - 1.0_xp
                orth_fro_samples(sample_index) = orth_fro_samples(sample_index) + &
                     real(dot_value*dot_value,rp)
                if (i /= k) max_offdiag = max(max_offdiag,abs(dot_value))
             end do
          end do
          orth_fro_samples(sample_index) = sqrt(orth_fro_samples(sample_index))
          orth_max_samples(sample_index) = real(max_offdiag,rp)
          deallocate(gram)

          call apply_s_array(z(:,j),sz)
          call apply_s_array(z(:,j),sz_again)
          repeat_samples(sample_index) = weighted_norm(sz-sz_again) / &
               max(weighted_norm(sz),tiny(1.0_rp))
          w = sz
          do k = 1, j+1
             w = w - real(hbar(k,j),rp)*v(:,k)
          end do
          relation_norm = real(weighted_norm(w),xp)
          relation_samples(sample_index) = real(relation_norm / &
               max(real(weighted_norm(sz),xp),tiny(1.0_xp)),rp)
       end if
       if (alpha <= tiny(1.0_xp)) exit
    end do
    deallocate(v,z,w,wplus,sol,sz,sz_again,hbar,hr,cosine,sine,gamma,y)
  end subroutine instrumented_arnoldi

  !> Test-only restarted GMRES clone for cap-200/201/400 comparisons.
  !! The residual is recomputed as b-Sx at every 200-step cycle boundary.
  subroutine instrumented_restarted_arnoldi(b, x_samples, true_samples, monitor_samples, &
       true_abs_samples, grad_diff_samples, grad_relative_samples, bc_correction_samples, &
       continuity_abs_samples, continuity_relative_samples, cycle_samples, orth_fro_samples, &
       orth_max_samples, relation_samples, repeat_samples, vnext_samples, reorthogonalize, &
       manufactured)
    real(kind=rp), intent(in) :: b(:)
    real(kind=rp), intent(out) :: x_samples(:, :), true_samples(:), monitor_samples(:)
    real(kind=rp), intent(out) :: true_abs_samples(:), grad_diff_samples(:)
    real(kind=rp), intent(out) :: grad_relative_samples(:), bc_correction_samples(:)
    real(kind=rp), intent(out) :: continuity_abs_samples(:), continuity_relative_samples(:)
    real(kind=rp), intent(out) :: orth_fro_samples(:), orth_max_samples(:)
    real(kind=rp), intent(out) :: relation_samples(:), repeat_samples(:), vnext_samples(:)
    integer, intent(out) :: cycle_samples(:)
    logical, intent(in) :: reorthogonalize, manufactured
    real(kind=rp), allocatable :: v(:, :), z(:, :), w(:), xbase(:), residual(:)
    real(kind=rp), allocatable :: candidate(:), applied(:), sz(:), sz_again(:)
    real(kind=xp), allocatable :: wplus(:), hbar(:, :), hr(:, :), cosine(:), sine(:)
    real(kind=xp), allocatable :: gamma(:), y(:), gram(:, :)
    real(kind=xp) :: alpha, beta, lr, tmp, dot_value, max_offdiag, local_sum
    real(kind=xp) :: hlocal(200), hsecond(200)
    real(kind=rp) :: rhs_norm, gq_norm
    integer :: total, cycle, j, i, k, sample_index, dof_index, steps_in_cycle
    integer, parameter :: sample_caps(3) = [200, 201, 400]

    allocate(v(n,201), z(n,200), w(n), xbase(n), residual(n), candidate(n), applied(n))
    allocate(sz(n), sz_again(n), wplus(n))
    allocate(hbar(201,200), hr(201,200), cosine(200), sine(200), gamma(201), y(200))
    x_samples = 0.0_rp; true_samples = 0.0_rp; monitor_samples = 0.0_rp
    true_abs_samples = 0.0_rp; grad_diff_samples = -1.0_rp
    grad_relative_samples = -1.0_rp; bc_correction_samples = 0.0_rp
    continuity_abs_samples = 0.0_rp; continuity_relative_samples = 0.0_rp
    cycle_samples = 0; orth_fro_samples = 0.0_rp; orth_max_samples = 0.0_rp
    relation_samples = 0.0_rp; repeat_samples = 0.0_rp; vnext_samples = 0.0_rp
    rhs_norm = weighted_norm(b)
    gq_norm = velocity_norm(gq_u,gq_v,gq_w)
    xbase = 0.0_rp
    total = 0
    cycle = 0

    do while (total < maxval(sample_caps))
       cycle = cycle + 1
       if (total == 0) then
          residual = b
       else
          call apply_s_array(xbase, applied)
          residual = b - applied
       end if
       beta = sqrt(weighted_dot(residual,residual))
       if (beta <= tiny(1.0_xp)) exit
       v = 0.0_rp; z = 0.0_rp; hbar = 0.0_xp; hr = 0.0_xp
       cosine = 1.0_xp; sine = 0.0_xp; gamma = 0.0_xp
       gamma(1) = beta
       v(:,1) = residual * real(1.0_xp/beta,rp)
       steps_in_cycle = 0

       do j = 1, 200
          steps_in_cycle = j
          total = total + 1
          call fluid%ksp_prs%M%solve(z(:,j),v(:,j),n)
          call apply_s_array(z(:,j),w)
          hlocal = 0.0_xp
          do i = 1, j
             do k = 1, n
                hlocal(i) = hlocal(i) + real(w(k),xp)*real(v(k,i),xp)* &
                     real(fluid%c_Xh%mult(k,1,1,1),xp)
             end do
          end do
          call MPI_Allreduce(MPI_IN_PLACE,hlocal,j,MPI_EXTRA_PRECISION,MPI_SUM,NEKO_COMM,mpi_ierr)
          wplus = 0.0_xp
          do i = 1, j
             hbar(i,j) = hlocal(i)
             wplus = wplus - hlocal(i)*real(v(:,i),xp)
          end do
          w = real(real(w,xp)+wplus,rp)
          if (reorthogonalize) then
             hsecond = 0.0_xp
             do i = 1, j
                do k = 1, n
                   hsecond(i) = hsecond(i) + real(w(k),xp)*real(v(k,i),xp)* &
                        real(fluid%c_Xh%mult(k,1,1,1),xp)
                end do
             end do
             call MPI_Allreduce(MPI_IN_PLACE,hsecond,j,MPI_EXTRA_PRECISION,MPI_SUM,NEKO_COMM,mpi_ierr)
             wplus = 0.0_xp
             do i = 1, j
                hbar(i,j) = hbar(i,j) + hsecond(i)
                wplus = wplus - hsecond(i)*real(v(:,i),xp)
             end do
             w = real(real(w,xp)+wplus,rp)
          end if
          alpha = real(weighted_norm(w),xp)
          hbar(j+1,j) = alpha
          if (alpha > tiny(1.0_xp)) v(:,j+1) = w*real(1.0_xp/alpha,rp)
          hr(1:j+1,j) = hbar(1:j+1,j)
          do i = 1, j-1
             tmp = hr(i,j)
             hr(i,j) = cosine(i)*tmp+sine(i)*hr(i+1,j)
             hr(i+1,j) = -sine(i)*tmp+cosine(i)*hr(i+1,j)
          end do
          lr = sqrt(hr(j,j)*hr(j,j)+alpha*alpha)
          if (lr <= tiny(1.0_xp)) exit
          cosine(j) = hr(j,j)/lr; sine(j) = alpha/lr; hr(j,j) = lr
          gamma(j+1) = -sine(j)*gamma(j); gamma(j) = cosine(j)*gamma(j)

          sample_index = 0
          do i = 1, size(sample_caps)
             if (total == sample_caps(i)) sample_index = i
          end do
          if (sample_index > 0) then
             call form_restarted_candidate(j,n,gamma,hr,z,xbase,candidate)
             x_samples(:,sample_index) = candidate
             call true_residual_metrics(candidate,b,true_abs_samples(sample_index), &
                  true_samples(sample_index))
             monitor_samples(sample_index) = abs(gamma(j+1))/max(rhs_norm,tiny(1.0_rp))
             cycle_samples(sample_index) = cycle
             vnext_samples(sample_index) = weighted_norm(v(:,j+1))
             call restarted_velocity_metrics(candidate,b,manufactured, &
                  grad_diff_samples(sample_index),grad_relative_samples(sample_index), &
                  bc_correction_samples(sample_index),continuity_abs_samples(sample_index), &
                  continuity_relative_samples(sample_index))

             allocate(gram(j+1,j+1)); gram = 0.0_xp
             do k = 1,j+1
                do i = 1,j+1
                   local_sum = 0.0_xp
                   do dof_index = 1,n
                      local_sum = local_sum + real(v(dof_index,i),xp)*real(v(dof_index,k),xp)* &
                           real(fluid%c_Xh%mult(dof_index,1,1,1),xp)
                   end do
                   gram(i,k) = local_sum
                end do
             end do
             call MPI_Allreduce(MPI_IN_PLACE,gram(1,1),(j+1)*(j+1), &
                  MPI_EXTRA_PRECISION,MPI_SUM,NEKO_COMM,mpi_ierr)
             max_offdiag = 0.0_xp; orth_fro_samples(sample_index) = 0.0_rp
             do k = 1,j+1
                do i = 1,j+1
                   dot_value = gram(i,k)
                   if (i == k) dot_value = dot_value-1.0_xp
                   orth_fro_samples(sample_index) = orth_fro_samples(sample_index) + &
                        real(dot_value*dot_value,rp)
                   if (i /= k) max_offdiag = max(max_offdiag,abs(dot_value))
                end do
             end do
             orth_fro_samples(sample_index) = sqrt(orth_fro_samples(sample_index))
             orth_max_samples(sample_index) = real(max_offdiag,rp)
             deallocate(gram)
             call apply_s_array(z(:,j),sz); call apply_s_array(z(:,j),sz_again)
             repeat_samples(sample_index) = weighted_norm(sz-sz_again) / &
                  max(weighted_norm(sz),tiny(1.0_rp))
             w = sz
             do k = 1,j+1
                w = w-real(hbar(k,j),rp)*v(:,k)
             end do
             relation_samples(sample_index) = weighted_norm(w)/ &
                  max(weighted_norm(sz),tiny(1.0_rp))
          end if
          if (alpha <= tiny(1.0_xp) .or. total >= maxval(sample_caps)) exit
       end do

       if (steps_in_cycle > 0) then
          call form_restarted_candidate(steps_in_cycle,n,gamma,hr,z,xbase,candidate)
          xbase = candidate
       end if
    end do
    deallocate(v,z,w,xbase,residual,candidate,applied,sz,sz_again,wplus, &
         hbar,hr,cosine,sine,gamma,y)
  end subroutine instrumented_restarted_arnoldi

  subroutine form_restarted_candidate(m,nlocal,gamma,hr,z,x_old,x_new)
    integer, intent(in) :: m,nlocal
    real(kind=xp), intent(in) :: gamma(:),hr(:,:)
    real(kind=rp), intent(in) :: z(:,:)
    real(kind=rp), intent(in) :: x_old(:)
    real(kind=rp), intent(out) :: x_new(:)
    real(kind=xp) :: coeff(200),accum,tmp_local
    integer :: row,col
    coeff(1:m) = gamma(1:m)
    do row = m,1,-1
       tmp_local = coeff(row)
       do col = m,row+1,-1
          tmp_local = tmp_local-hr(row,col)*coeff(col)
       end do
       coeff(row) = tmp_local/hr(row,row)
    end do
    x_new = x_old
    do row = 1,nlocal
       accum = 0.0_xp
       do col = 1,m
          accum = accum+coeff(col)*real(z(row,col),xp)
       end do
       x_new(row) = x_new(row)+real(accum,rp)
    end do
  end subroutine form_restarted_candidate

  subroutine restarted_velocity_metrics(pressure,b,manufactured,grad_difference, &
       grad_relative,bc_correction,continuity_absolute,continuity_relative)
    real(kind=rp), intent(in) :: pressure(:),b(:)
    logical, intent(in) :: manufactured
    real(kind=rp), intent(out) :: grad_difference,grad_relative,bc_correction
    real(kind=rp), intent(out) :: continuity_absolute,continuity_relative
    phi%x = reshape(pressure,shape(phi%x))
    call apply_gradient_map_unprojected(phi,gp_u,gp_v,gp_w)
    fx%x = gp_u%x; fy%x = gp_v%x; fz%x = gp_w%x
    call fluid%bcs_vel_projector%apply(gp_u%x,gp_v%x,gp_w%x,n)
    fx%x = fx%x-gp_u%x; fy%x = fy%x-gp_v%x; fz%x = fz%x-gp_w%x
    bc_correction = velocity_norm(fx,fy,fz)
    grad_difference = -1.0_rp; grad_relative = -1.0_rp
    if (manufactured) then
       fx%x = gp_u%x-gq_u%x; fy%x = gp_v%x-gq_v%x; fz%x = gp_w%x-gq_w%x
       grad_difference = velocity_norm(fx,fy,fz)
       grad_relative = grad_difference/max(velocity_norm(gq_u,gq_v,gq_w),tiny(1.0_rp))
    end if
    call apply_continuity(gp_u,gp_v,gp_w,pres)
    pres%x = pres%x-reshape(b,shape(pres%x))
    if (.not. fluid%prs_dirichlet) call ortho(pres%x,fluid%glb_n_points,n)
    call fluid%bcs_prs_projector%apply(pres%x,n)
    continuity_absolute = pressure_norm(pres)
    continuity_relative = continuity_absolute/max(weighted_norm(b),tiny(1.0_rp))
  end subroutine restarted_velocity_metrics

  subroutine projected_velocity_difference(pressure_native,pressure_clone,relative_difference)
    real(kind=rp), intent(in) :: pressure_native(:),pressure_clone(:)
    real(kind=rp), intent(out) :: relative_difference
    real(kind=rp) :: native_gradient_norm
    phi%x = reshape(pressure_native,shape(phi%x))
    call apply_gradient_map(phi,gp_u,gp_v,gp_w)
    fx%x = gp_u%x; fy%x = gp_v%x; fz%x = gp_w%x
    native_gradient_norm = velocity_norm(fx,fy,fz)
    phi%x = reshape(pressure_clone,shape(phi%x))
    call apply_gradient_map(phi,gp_u,gp_v,gp_w)
    gp_u%x = gp_u%x-fx%x; gp_v%x = gp_v%x-fy%x; gp_w%x = gp_w%x-fz%x
    relative_difference = velocity_norm(gp_u,gp_v,gp_w) / &
         max(native_gradient_norm,tiny(1.0_rp))
  end subroutine projected_velocity_difference

  subroutine apply_s_array(input, output)
    real(kind=rp), intent(in) :: input(:)
    real(kind=rp), intent(out) :: output(:)
    event = C_NULL_PTR
    q%x = reshape(input,shape(q%x))
    call schur%compute(work%x,q%x,fluid%c_Xh,fluid%msh,fluid%Xh)
    call fluid%gs_Xh%op(work%x,n,GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(work%x,n)
    output = reshape(work%x,[n])
  end subroutine apply_s_array

  subroutine true_residual_metrics(solution, b, absolute, relative)
    real(kind=rp), intent(in) :: solution(:), b(:)
    real(kind=rp), intent(out) :: absolute, relative
    real(kind=rp), allocatable :: applied(:)
    allocate(applied(n))
    call apply_s_array(solution,applied)
    absolute = weighted_norm(applied-b)
    relative = absolute / max(weighted_norm(b),tiny(1.0_rp))
    deallocate(applied)
  end subroutine true_residual_metrics

  function weighted_norm(a) result(value)
    real(kind=rp), intent(in) :: a(:)
    real(kind=rp) :: value
    value = real(sqrt(max(0.0_xp,weighted_dot(a,a))),rp)
  end function weighted_norm

  function weighted_dot(a,b) result(value)
    real(kind=rp), intent(in) :: a(:),b(:)
    real(kind=xp) :: value
    integer :: k, ierr
    value = 0.0_xp
    do k = 1, n
       value = value + real(a(k),xp)*real(b(k),xp)*real(fluid%c_Xh%mult(k,1,1,1),xp)
    end do
    call MPI_Allreduce(MPI_IN_PLACE,value,1,MPI_EXTRA_PRECISION,MPI_SUM,NEKO_COMM,ierr)
  end function weighted_dot

  subroutine continuous_scalar(a)
    type(field_t), intent(inout) :: a
    call fluid%gs_Xh%op(a%x, n, GS_OP_ADD)
    call col2(a%x, fluid%c_Xh%mult, n)
  end subroutine continuous_scalar

  subroutine continuous_vector(a, b, c)
    type(field_t), intent(inout) :: a, b, c
    call fluid%gs_Xh%op(a%x, b%x, c%x, n, GS_OP_ADD)
    call col2(a%x, fluid%c_Xh%mult, n)
    call col2(b%x, fluid%c_Xh%mult, n)
    call col2(c%x, fluid%c_Xh%mult, n)
  end subroutine continuous_vector

  subroutine assemble_pressure(a)
    type(field_t), intent(inout) :: a
    if (.not. fluid%prs_dirichlet) call ortho(a%x, fluid%glb_n_points, n)
    call fluid%gs_Xh%op(a%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(a%x, n)
  end subroutine assemble_pressure

  subroutine apply_continuity(u, v, w, out)
    type(field_t), intent(inout) :: u, v, w, out
    event = C_NULL_PTR
    call field_col3(fx, u, fluid%rho)
    call field_col3(fy, v, fluid%rho)
    call field_col3(fz, w, fluid%rho)
    call col2(fx%x, fluid%c_Xh%B, n)
    call col2(fy%x, fluid%c_Xh%B, n)
    call col2(fz%x, fluid%c_Xh%B, n)
    call field_rzero(pzero)
    call fluid%prs_res%compute(pzero, pres, u, v, w, uzero, vzero, wzero, &
         fx, fy, fz, fluid%c_Xh, fluid%gs_Xh, fluid%bc_prs_surface, &
         fluid%bc_sym_surface, fluid%Ax_prs, 0.0_rp, 1.0_rp, uzero, &
         fluid%rho, event)
    if (.not. fluid%prs_dirichlet) call ortho(pres%x, fluid%glb_n_points, n)
    call fluid%gs_Xh%op(pres%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(pres%x, n)
    out%x = pres%x
  end subroutine apply_continuity

  subroutine apply_gradient_map(pressure, u, v, w)
    type(field_t), intent(inout) :: pressure
    type(field_t), intent(inout) :: u, v, w
    call apply_gradient_map_unprojected(pressure,u,v,w)
    call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
  end subroutine apply_gradient_map

  subroutine apply_gradient_map_unprojected(pressure, u, v, w)
    type(field_t), intent(inout) :: pressure
    type(field_t), intent(inout) :: u, v, w
    if (.not. fluid%prs_dirichlet) call ortho(pressure%x, fluid%glb_n_points, n)
    call fluid%bcs_prs_projector%apply(pressure%x, n)
    call grad(u%x, v%x, w%x, pressure%x, fluid%c_Xh)
    call col2(u%x, fluid%c_Xh%B, n)
    call col2(v%x, fluid%c_Xh%B, n)
    call col2(w%x, fluid%c_Xh%B, n)
    call fluid%gs_Xh%op(u%x, v%x, w%x, n, GS_OP_ADD)
    call col2(u%x, fluid%c_Xh%Binv, n)
    call col2(v%x, fluid%c_Xh%Binv, n)
    call col2(w%x, fluid%c_Xh%Binv, n)
  end subroutine apply_gradient_map_unprojected

  subroutine make_velocity_probe()
    real(kind=rp) :: a, b, c, bubble
    integer :: j
    do j = 1, n
       a = normalize(fluid%u%dof%x%x(j,1,1,1), xmin, xmax)
       b = normalize(fluid%u%dof%y%x(j,1,1,1), ymin, ymax)
       c = normalize(fluid%u%dof%z%x(j,1,1,1), zmin, zmax)
       bubble = (a*(1.0_rp-a)*b*(1.0_rp-b)*c*(1.0_rp-c))**2
       vu%x(j,1,1,1) = bubble*(1.0_rp+a)
       vv%x(j,1,1,1) = bubble*(1.0_rp+b*b)
       vw%x(j,1,1,1) = bubble*(1.0_rp+c*c)
    end do
    call continuous_vector(vu, vv, vw)
    call fluid%bcs_vel_projector%apply(vu%x, vv%x, vw%x, n)
  end subroutine make_velocity_probe

  subroutine true_schur_residual(solution, rhs, value)
    type(field_t), intent(in) :: solution
    real(kind=rp), intent(in) :: rhs(:)
    real(kind=rp), intent(out) :: value
    call schur%compute(work%x, solution%x, fluid%c_Xh, fluid%msh, fluid%Xh)
    call fluid%gs_Xh%op(work%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(work%x, n)
    work%x = work%x - reshape(rhs, shape(work%x))
    call fluid%bcs_prs_projector%apply(work%x, n)
    value = pressure_norm(work)
  end subroutine true_schur_residual

  function pressure_norm(a) result(value)
    type(field_t), intent(in) :: a
    real(kind=rp) :: value
    value = sqrt(max(0.0_rp, glsc3(a%x, fluid%c_Xh%mult, a%x, n)))
  end function pressure_norm

  function velocity_norm(a, b, c) result(value)
    type(field_t), intent(in) :: a, b, c
    real(kind=rp) :: value
    value = sqrt(max(0.0_rp, glsc3(a%x, fluid%c_Xh%B, a%x, n) + &
         glsc3(b%x, fluid%c_Xh%B, b%x, n) + &
         glsc3(c%x, fluid%c_Xh%B, c%x, n)))
  end function velocity_norm

  function normalize(value, low, high) result(res)
    real(kind=rp), intent(in) :: value, low, high
    real(kind=rp) :: res
    res = (value-low) / max(high-low, tiny(1.0_rp))
  end function normalize

end program rom_projection_preconditioner_diagnostic
