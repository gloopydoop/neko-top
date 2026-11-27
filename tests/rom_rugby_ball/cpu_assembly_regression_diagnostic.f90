program rom_cpu_assembly_regression_diagnostic
  ! Compare the packed CPU production contractions with their historical
  ! per-entry global vector_glsc3 reductions on identical live mesh fields.
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use fluid_pnpn, only: fluid_pnpn_t
  use coefs, only: coef_t
  use field, only: field_t
  use operators, only: grad
  use num_types, only: rp
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: HOST_TO_DEVICE, DEVICE_TO_HOST
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use comm, only: NEKO_COMM, pe_rank
  use mpi_f08, only: MPI_Comm_size, MPI_Allreduce, MPI_MAX, MPI_DOUBLE_PRECISION, MPI_Wtime
  use utils, only: neko_error
  use rom_array_test_adapters, only: test_assemble_mass, &
       test_assemble_mass_legacy, test_assemble_diffusion, &
       test_assemble_diffusion_legacy, test_assemble_brinkman, &
       test_assemble_brinkman_legacy, test_assemble_convection, &
       test_assemble_convection_legacy, test_assemble_scalar_transport, &
       test_assemble_scalar_transport_legacy, test_assemble_cross_mass_petrov, &
       test_assemble_cross_mass_petrov_legacy, test_assemble_brinkman_petrov, &
       test_assemble_brinkman_petrov_legacy, test_assemble_convection_dgemm
  implicit none

  integer, parameter :: rank_basis = 3
  real(kind=rp), parameter :: comparison_tolerance = 5.0e-13_rp
  type(simulation_t) :: simulation
  type(json_file) :: parameters
  type(fluid_pnpn_t), pointer :: fluid
  type(field_t) :: u, v, w, gx, gy, gz
  real(kind=rp), allocatable :: weights(:), chi(:), lift(:, :)
  real(kind=rp), allocatable :: modes(:, :, :), gradients(:, :, :, :)
  real(kind=rp), allocatable :: test_modes(:, :, :)
  real(kind=rp), allocatable :: scalar_modes(:, :), scalar_gradients(:, :, :)
  real(kind=rp), allocatable :: scalar_lift(:), lift_scalar_gradient(:, :), diffusivity(:)
  real(kind=rp), allocatable :: mass_packed(:, :), mass_legacy(:, :)
  real(kind=rp), allocatable :: diffusion_packed(:, :), diffusion_legacy(:, :)
  real(kind=rp), allocatable :: brinkman_packed(:, :), brinkman_legacy(:, :)
  real(kind=rp), allocatable :: brinkman_lift_packed(:), brinkman_lift_legacy(:)
  real(kind=rp), allocatable :: convection_packed(:, :, :), convection_legacy(:, :, :)
  real(kind=rp), allocatable :: convection_dgemm(:, :, :), packed_times(:), dgemm_times(:)
  real(kind=rp), allocatable :: scalar_linear_packed(:, :), scalar_linear_legacy(:, :)
  real(kind=rp), allocatable :: scalar_constant_packed(:), scalar_constant_legacy(:)
  real(kind=rp), allocatable :: scalar_quadratic_packed(:, :, :), scalar_quadratic_legacy(:, :, :)
  real(kind=rp), allocatable :: scalar_implicit_packed(:, :), scalar_implicit_legacy(:, :)
  real(kind=rp), allocatable :: petrov_mass_packed(:, :), petrov_mass_legacy(:, :)
  real(kind=rp), allocatable :: petrov_brinkman_packed(:, :), petrov_brinkman_legacy(:, :)
  real(kind=rp), allocatable :: petrov_lift_packed(:), petrov_lift_legacy(:)
  real(kind=rp) :: x, y, z, rho, mu, nu
  real(kind=rp) :: mass_error, diffusion_error, brinkman_error, brinkman_lift_error
  real(kind=rp) :: convection_error, mass_relative, diffusion_relative
  real(kind=rp) :: brinkman_relative, brinkman_lift_relative, convection_relative
  real(kind=rp) :: scalar_linear_error, scalar_constant_error, scalar_quadratic_error
  real(kind=rp) :: scalar_implicit_error, scalar_linear_relative, scalar_constant_relative
  real(kind=rp) :: scalar_quadratic_relative, scalar_implicit_relative
  real(kind=rp) :: petrov_mass_error, petrov_mass_relative
  real(kind=rp) :: petrov_brinkman_error, petrov_brinkman_relative
  real(kind=rp) :: petrov_lift_error, petrov_lift_relative
  real(kind=rp) :: dgemm_packed_error, dgemm_packed_relative
  real(kind=rp) :: dgemm_legacy_error, dgemm_legacy_relative
  real(kind=rp) :: timer_start, elapsed, global_elapsed, packed_mean, packed_spread
  real(kind=rp) :: dgemm_mean, dgemm_spread
  real(kind=rp) :: scalar_weight
  real(kind=rp) :: scale
  integer :: argc, n, i, mode, nranks, ierr, csv_unit, io_status
  integer :: benchmark_repeats, bench_iter, gemm_block_size
  logical :: found
  logical :: write_csv
  logical :: run_gemm_benchmark
  character(len=1024) :: case_file, csv_file, repeat_arg

  argc = command_argument_count()
  if (argc < 1) call neko_error( &
       'Usage: ROM_rugby_ball_cpu_assembly_regression_diagnostic CASE [OUTPUT.csv [GEMM_REPEATS]]')
  call get_command_argument(1, case_file)
  write_csv = argc >= 2
  if (write_csv) call get_command_argument(2, csv_file)
  run_gemm_benchmark = argc >= 3
  benchmark_repeats = 0
  if (run_gemm_benchmark) then
     call get_command_argument(3, repeat_arg, status=io_status)
     if (io_status /= 0) call neko_error('Could not read GEMM repeat count')
     read(repeat_arg,*,iostat=io_status) benchmark_repeats
     if (io_status /= 0 .or. benchmark_repeats < 2) &
          call neko_error('GEMM benchmark repeat count must be an integer >= 2')
  end if
  gemm_block_size = 8

  call neko_init()
  call neko_top_register_types()
  parameters = json_read_file(trim(case_file))
  ! The simulation defaults to steady mode; avoid mutating the loaded case
  ! object here so its nested case.mesh_file entry remains intact.
  call user_setup(simulation%neko_case%user)
  call simulation%init(parameters)
  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('CPU assembly regression requires a Pn/Pn fluid')
  end select

  n = fluid%u%dof%size()
  allocate(weights(n), chi(n), lift(n,3), modes(n,3,rank_basis), &
       gradients(n,3,rank_basis,3))
  allocate(scalar_modes(n,rank_basis), scalar_gradients(n,rank_basis,3))
  allocate(test_modes(n,3,rank_basis))
  allocate(scalar_lift(n), lift_scalar_gradient(n,3), diffusivity(n))
  allocate(mass_packed(rank_basis,rank_basis), mass_legacy(rank_basis,rank_basis))
  allocate(diffusion_packed(rank_basis,rank_basis), diffusion_legacy(rank_basis,rank_basis))
  allocate(brinkman_packed(rank_basis,rank_basis), brinkman_legacy(rank_basis,rank_basis))
  allocate(brinkman_lift_packed(rank_basis), brinkman_lift_legacy(rank_basis))
  allocate(convection_packed(rank_basis,rank_basis,rank_basis), &
       convection_legacy(rank_basis,rank_basis,rank_basis))
  allocate(scalar_linear_packed(rank_basis,rank_basis), &
       scalar_linear_legacy(rank_basis,rank_basis))
  allocate(scalar_constant_packed(rank_basis), scalar_constant_legacy(rank_basis))
  allocate(scalar_quadratic_packed(rank_basis,rank_basis,rank_basis), &
       scalar_quadratic_legacy(rank_basis,rank_basis,rank_basis))
  allocate(scalar_implicit_packed(rank_basis,rank_basis), &
       scalar_implicit_legacy(rank_basis,rank_basis))
  allocate(petrov_mass_packed(rank_basis,rank_basis), petrov_mass_legacy(rank_basis,rank_basis))
  allocate(petrov_brinkman_packed(rank_basis,rank_basis), &
       petrov_brinkman_legacy(rank_basis,rank_basis))
  allocate(petrov_lift_packed(rank_basis), petrov_lift_legacy(rank_basis))
  weights = reshape(fluid%c_Xh%B, [n])
  call u%init(fluid%dm_Xh, 'cpu_assembly_regression_u')
  call v%init(fluid%dm_Xh, 'cpu_assembly_regression_v')
  call w%init(fluid%dm_Xh, 'cpu_assembly_regression_w')
  call gx%init(fluid%dm_Xh, 'cpu_assembly_regression_gx')
  call gy%init(fluid%dm_Xh, 'cpu_assembly_regression_gy')
  call gz%init(fluid%dm_Xh, 'cpu_assembly_regression_gz')

  ! Three deterministic smooth vector modes are sampled on the actual case
  ! GLL coordinates. Both implementations receive these exact values and
  ! the same Neko-computed gradients.
  do mode = 1, rank_basis
     do i = 1, n
        x = fluid%u%dof%x%x(i,1,1,1)
        y = fluid%u%dof%y%x(i,1,1,1)
        z = fluid%u%dof%z%x(i,1,1,1)
        u%x(i,1,1,1) = sin(real(mode,rp)*x) * cos(y) + 0.07_rp*real(mode,rp)*z
        v%x(i,1,1,1) = cos(real(mode,rp)*y) * sin(z) - 0.05_rp*real(mode,rp)*x
        w%x(i,1,1,1) = sin(z + real(mode,rp)*x) * cos(y)
     end do
     modes(:,1,mode) = reshape(u%x, [n])
     modes(:,2,mode) = reshape(v%x, [n])
     modes(:,3,mode) = reshape(w%x, [n])
     call gradient_from_host(u, gx, gy, gz, fluid%c_Xh)
     gradients(:,1,mode,1) = reshape(gx%x, [n])
     gradients(:,1,mode,2) = reshape(gy%x, [n])
     gradients(:,1,mode,3) = reshape(gz%x, [n])
     call gradient_from_host(v, gx, gy, gz, fluid%c_Xh)
     gradients(:,2,mode,1) = reshape(gx%x, [n])
     gradients(:,2,mode,2) = reshape(gy%x, [n])
     gradients(:,2,mode,3) = reshape(gz%x, [n])
     call gradient_from_host(w, gx, gy, gz, fluid%c_Xh)
     gradients(:,3,mode,1) = reshape(gx%x, [n])
     gradients(:,3,mode,2) = reshape(gy%x, [n])
     gradients(:,3,mode,3) = reshape(gz%x, [n])
  end do

  ! A nonconstant design multiplier and spatial lift exercise both Brinkman
  ! outputs while remaining identical for packed and legacy contractions.
  do i = 1, n
     x = fluid%u%dof%x%x(i,1,1,1)
     y = fluid%u%dof%y%x(i,1,1,1)
     z = fluid%u%dof%z%x(i,1,1,1)
     chi(i) = 0.25_rp + 0.35_rp * (0.5_rp + 0.5_rp * sin(0.19_rp*x) * cos(0.23_rp*y))
     lift(i,1) = 0.10_rp + 0.03_rp * x
     lift(i,2) = -0.08_rp + 0.02_rp * y
     lift(i,3) = 0.04_rp + 0.01_rp * z
     scalar_lift(i) = 0.3_rp + 0.07_rp * sin(x+y+z)
     diffusivity(i) = 0.01_rp + 0.003_rp * (1.0_rp + 0.5_rp*sin(x)*cos(z))
  end do
  u%x = reshape(scalar_lift, shape(u%x))
  call gradient_from_host(u, gx, gy, gz, fluid%c_Xh)
  lift_scalar_gradient(:,1) = reshape(gx%x, [n])
  lift_scalar_gradient(:,2) = reshape(gy%x, [n])
  lift_scalar_gradient(:,3) = reshape(gz%x, [n])
  do mode = 1, rank_basis
     do i = 1, n
        x = fluid%u%dof%x%x(i,1,1,1)
        y = fluid%u%dof%y%x(i,1,1,1)
        z = fluid%u%dof%z%x(i,1,1,1)
        u%x(i,1,1,1) = sin((0.5_rp+real(mode,rp))*x) * cos(1.2_rp*y) + &
             0.1_rp*z*x
     end do
     scalar_modes(:,mode) = reshape(u%x, [n])
     call gradient_from_host(u, gx, gy, gz, fluid%c_Xh)
     scalar_gradients(:,mode,1) = reshape(gx%x, [n])
     scalar_gradients(:,mode,2) = reshape(gy%x, [n])
     scalar_gradients(:,mode,3) = reshape(gz%x, [n])
  end do
  do mode = 1, rank_basis
     do i = 1, n
        x = fluid%u%dof%x%x(i,1,1,1)
        y = fluid%u%dof%y%x(i,1,1,1)
        z = fluid%u%dof%z%x(i,1,1,1)
        u%x(i,1,1,1) = cos((real(mode,rp)+0.7_rp)*x) * sin(y) + 0.13_rp*z
        v%x(i,1,1,1) = sin((real(mode,rp)+0.4_rp)*y) * cos(z)
        w%x(i,1,1,1) = cos(z+(real(mode,rp)+0.6_rp)*x) * sin(y)
     end do
     test_modes(:,1,mode) = reshape(u%x, [n])
     test_modes(:,2,mode) = reshape(v%x, [n])
     test_modes(:,3,mode) = reshape(w%x, [n])
  end do

  call test_assemble_mass(weights, modes, mass_packed)
  call test_assemble_mass_legacy(weights, modes, mass_legacy)
  rho = real(fluid%rho%x(1,1,1,1), rp)
  mu = real(fluid%mu%x(1,1,1,1), rp)
  nu = mu / rho
  call test_assemble_diffusion(weights, nu, gradients, diffusion_packed)
  call test_assemble_diffusion_legacy(weights, nu, gradients, diffusion_legacy)
  call test_assemble_brinkman(weights, chi, modes, lift, &
       brinkman_packed, brinkman_lift_packed)
  call test_assemble_brinkman_legacy(weights, chi, modes, lift, &
       brinkman_legacy, brinkman_lift_legacy)
  call test_assemble_cross_mass_petrov(weights, test_modes, modes, petrov_mass_packed)
  call test_assemble_cross_mass_petrov_legacy(weights, test_modes, modes, petrov_mass_legacy)
  call test_assemble_brinkman_petrov(weights, chi, test_modes, modes, lift, &
       petrov_brinkman_packed, petrov_lift_packed)
  call test_assemble_brinkman_petrov_legacy(weights, chi, test_modes, modes, lift, &
       petrov_brinkman_legacy, petrov_lift_legacy)
  call test_assemble_convection(weights, modes, gradients, convection_packed)
  call test_assemble_convection_legacy(weights, modes, gradients, convection_legacy)
  if (run_gemm_benchmark) then
     allocate(convection_dgemm(rank_basis,rank_basis,rank_basis))
     allocate(packed_times(benchmark_repeats), dgemm_times(benchmark_repeats))
     call test_assemble_convection_dgemm(weights, modes, gradients, gemm_block_size, &
          convection_dgemm)
     dgemm_packed_error = maxval(abs(convection_dgemm-convection_packed))
     scale = max(sqrt(sum(convection_packed**2)), tiny(1.0_rp))
     dgemm_packed_relative = sqrt(sum((convection_dgemm-convection_packed)**2)) / scale
     dgemm_legacy_error = maxval(abs(convection_dgemm-convection_legacy))
     scale = max(sqrt(sum(convection_legacy**2)), tiny(1.0_rp))
     dgemm_legacy_relative = sqrt(sum((convection_dgemm-convection_legacy)**2)) / scale
     if (.not. all(ieee_is_finite(convection_dgemm)) .or. &
          .not. ieee_is_finite(dgemm_packed_error) .or. &
          .not. ieee_is_finite(dgemm_packed_relative) .or. &
          .not. ieee_is_finite(dgemm_legacy_error) .or. &
          .not. ieee_is_finite(dgemm_legacy_relative)) &
          call neko_error('DGEMM benchmark produced a nonfinite result or metric')
     if (dgemm_packed_relative > comparison_tolerance .or. &
          dgemm_legacy_relative > comparison_tolerance) &
          call neko_error('DGEMM convection numerical-equivalence gate failed')
     do bench_iter = 1, benchmark_repeats
        if (mod(bench_iter,2) == 1) then
           timer_start = real(MPI_Wtime(),rp)
           call test_assemble_convection(weights, modes, gradients, convection_packed)
           elapsed = real(MPI_Wtime(),rp) - timer_start
           call MPI_Allreduce(elapsed, global_elapsed, 1, MPI_DOUBLE_PRECISION, &
                MPI_MAX, NEKO_COMM, ierr)
           if (ierr /= 0) call neko_error('MPI_Allreduce failed for packed timing')
           packed_times(bench_iter) = global_elapsed
           timer_start = real(MPI_Wtime(),rp)
           call test_assemble_convection_dgemm(weights, modes, gradients, gemm_block_size, &
                convection_dgemm)
           elapsed = real(MPI_Wtime(),rp) - timer_start
           call MPI_Allreduce(elapsed, global_elapsed, 1, MPI_DOUBLE_PRECISION, &
                MPI_MAX, NEKO_COMM, ierr)
           if (ierr /= 0) call neko_error('MPI_Allreduce failed for DGEMM timing')
           dgemm_times(bench_iter) = global_elapsed
        else
           timer_start = real(MPI_Wtime(),rp)
           call test_assemble_convection_dgemm(weights, modes, gradients, gemm_block_size, &
                convection_dgemm)
           elapsed = real(MPI_Wtime(),rp) - timer_start
           call MPI_Allreduce(elapsed, global_elapsed, 1, MPI_DOUBLE_PRECISION, &
                MPI_MAX, NEKO_COMM, ierr)
           if (ierr /= 0) call neko_error('MPI_Allreduce failed for DGEMM timing')
           dgemm_times(bench_iter) = global_elapsed
           timer_start = real(MPI_Wtime(),rp)
           call test_assemble_convection(weights, modes, gradients, convection_packed)
           elapsed = real(MPI_Wtime(),rp) - timer_start
           call MPI_Allreduce(elapsed, global_elapsed, 1, MPI_DOUBLE_PRECISION, &
                MPI_MAX, NEKO_COMM, ierr)
           if (ierr /= 0) call neko_error('MPI_Allreduce failed for packed timing')
           packed_times(bench_iter) = global_elapsed
        end if
     end do
     packed_mean = sum(packed_times) / real(benchmark_repeats,rp)
     packed_spread = sqrt(sum((packed_times-packed_mean)**2) / &
          real(max(1,benchmark_repeats-1),rp))
     dgemm_mean = sum(dgemm_times) / real(benchmark_repeats,rp)
     dgemm_spread = sqrt(sum((dgemm_times-dgemm_mean)**2) / &
          real(max(1,benchmark_repeats-1),rp))
  end if
  scalar_weight = 1.7_rp
  call test_assemble_scalar_transport(weights, scalar_weight, modes, scalar_modes, &
       gradients, scalar_gradients, lift, scalar_lift, lift_scalar_gradient, diffusivity, &
       scalar_linear_packed, scalar_constant_packed, scalar_quadratic_packed, &
       scalar_implicit_packed)
  call test_assemble_scalar_transport_legacy(weights, scalar_weight, modes, scalar_modes, &
       gradients, scalar_gradients, lift, scalar_lift, lift_scalar_gradient, diffusivity, &
       scalar_linear_legacy, scalar_constant_legacy, scalar_quadratic_legacy, &
       scalar_implicit_legacy)

  if (.not. all(ieee_is_finite(mass_packed)) .or. &
       .not. all(ieee_is_finite(mass_legacy)) .or. &
       .not. all(ieee_is_finite(diffusion_packed)) .or. &
       .not. all(ieee_is_finite(diffusion_legacy)) .or. &
       .not. all(ieee_is_finite(brinkman_packed)) .or. &
       .not. all(ieee_is_finite(brinkman_legacy)) .or. &
       .not. all(ieee_is_finite(brinkman_lift_packed)) .or. &
       .not. all(ieee_is_finite(brinkman_lift_legacy)) .or. &
       .not. all(ieee_is_finite(convection_packed)) .or. &
       .not. all(ieee_is_finite(convection_legacy)) .or. &
       .not. all(ieee_is_finite(scalar_linear_packed)) .or. &
       .not. all(ieee_is_finite(scalar_linear_legacy)) .or. &
       .not. all(ieee_is_finite(scalar_constant_packed)) .or. &
       .not. all(ieee_is_finite(scalar_constant_legacy)) .or. &
       .not. all(ieee_is_finite(scalar_quadratic_packed)) .or. &
       .not. all(ieee_is_finite(scalar_quadratic_legacy)) .or. &
       .not. all(ieee_is_finite(scalar_implicit_packed)) .or. &
       .not. all(ieee_is_finite(scalar_implicit_legacy)) .or. &
       .not. all(ieee_is_finite(petrov_mass_packed)) .or. &
       .not. all(ieee_is_finite(petrov_mass_legacy)) .or. &
       .not. all(ieee_is_finite(petrov_brinkman_packed)) .or. &
       .not. all(ieee_is_finite(petrov_brinkman_legacy)) .or. &
       .not. all(ieee_is_finite(petrov_lift_packed)) .or. &
       .not. all(ieee_is_finite(petrov_lift_legacy))) &
       call neko_error('CPU assembly regression produced a nonfinite operator')

  mass_error = maxval(abs(mass_packed-mass_legacy))
  scale = max(sqrt(sum(mass_legacy**2)), tiny(1.0_rp))
  mass_relative = sqrt(sum((mass_packed-mass_legacy)**2)) / scale
  diffusion_error = maxval(abs(diffusion_packed-diffusion_legacy))
  scale = max(sqrt(sum(diffusion_legacy**2)), tiny(1.0_rp))
  diffusion_relative = sqrt(sum((diffusion_packed-diffusion_legacy)**2)) / scale
  brinkman_error = maxval(abs(brinkman_packed-brinkman_legacy))
  scale = max(sqrt(sum(brinkman_legacy**2)), tiny(1.0_rp))
  brinkman_relative = sqrt(sum((brinkman_packed-brinkman_legacy)**2)) / scale
  brinkman_lift_error = maxval(abs(brinkman_lift_packed-brinkman_lift_legacy))
  scale = max(sqrt(sum(brinkman_lift_legacy**2)), tiny(1.0_rp))
  brinkman_lift_relative = sqrt(sum((brinkman_lift_packed-brinkman_lift_legacy)**2)) / scale
  convection_error = maxval(abs(convection_packed-convection_legacy))
  scale = max(sqrt(sum(convection_legacy**2)), tiny(1.0_rp))
  convection_relative = sqrt(sum((convection_packed-convection_legacy)**2)) / scale
  scalar_linear_error = maxval(abs(scalar_linear_packed-scalar_linear_legacy))
  scale = max(sqrt(sum(scalar_linear_legacy**2)), tiny(1.0_rp))
  scalar_linear_relative = sqrt(sum((scalar_linear_packed-scalar_linear_legacy)**2)) / scale
  scalar_constant_error = maxval(abs(scalar_constant_packed-scalar_constant_legacy))
  scale = max(sqrt(sum(scalar_constant_legacy**2)), tiny(1.0_rp))
  scalar_constant_relative = sqrt(sum((scalar_constant_packed-scalar_constant_legacy)**2)) / scale
  scalar_quadratic_error = maxval(abs(scalar_quadratic_packed-scalar_quadratic_legacy))
  scale = max(sqrt(sum(scalar_quadratic_legacy**2)), tiny(1.0_rp))
  scalar_quadratic_relative = &
       sqrt(sum((scalar_quadratic_packed-scalar_quadratic_legacy)**2)) / scale
  scalar_implicit_error = maxval(abs(scalar_implicit_packed-scalar_implicit_legacy))
  scale = max(sqrt(sum(scalar_implicit_legacy**2)), tiny(1.0_rp))
  scalar_implicit_relative = &
       sqrt(sum((scalar_implicit_packed-scalar_implicit_legacy)**2)) / scale
  petrov_mass_error = maxval(abs(petrov_mass_packed-petrov_mass_legacy))
  scale = max(sqrt(sum(petrov_mass_legacy**2)), tiny(1.0_rp))
  petrov_mass_relative = sqrt(sum((petrov_mass_packed-petrov_mass_legacy)**2)) / scale
  petrov_brinkman_error = maxval(abs(petrov_brinkman_packed-petrov_brinkman_legacy))
  scale = max(sqrt(sum(petrov_brinkman_legacy**2)), tiny(1.0_rp))
  petrov_brinkman_relative = &
       sqrt(sum((petrov_brinkman_packed-petrov_brinkman_legacy)**2)) / scale
  petrov_lift_error = maxval(abs(petrov_lift_packed-petrov_lift_legacy))
  scale = max(sqrt(sum(petrov_lift_legacy**2)), tiny(1.0_rp))
  petrov_lift_relative = &
       sqrt(sum((petrov_lift_packed-petrov_lift_legacy)**2)) / scale
  if (.not. ieee_is_finite(mass_error) .or. &
       .not. ieee_is_finite(diffusion_error) .or. &
       .not. ieee_is_finite(brinkman_error) .or. &
       .not. ieee_is_finite(brinkman_lift_error) .or. &
       .not. ieee_is_finite(convection_error) .or. &
       .not. ieee_is_finite(mass_relative) .or. &
       .not. ieee_is_finite(diffusion_relative) .or. &
       .not. ieee_is_finite(brinkman_relative) .or. &
       .not. ieee_is_finite(brinkman_lift_relative) .or. &
       .not. ieee_is_finite(convection_relative) .or. &
       .not. ieee_is_finite(scalar_linear_error) .or. &
       .not. ieee_is_finite(scalar_constant_error) .or. &
       .not. ieee_is_finite(scalar_quadratic_error) .or. &
       .not. ieee_is_finite(scalar_implicit_error) .or. &
       .not. ieee_is_finite(scalar_linear_relative) .or. &
       .not. ieee_is_finite(scalar_constant_relative) .or. &
       .not. ieee_is_finite(scalar_quadratic_relative) .or. &
       .not. ieee_is_finite(scalar_implicit_relative) .or. &
       .not. ieee_is_finite(petrov_mass_error) .or. &
       .not. ieee_is_finite(petrov_mass_relative) .or. &
       .not. ieee_is_finite(petrov_brinkman_error) .or. &
       .not. ieee_is_finite(petrov_brinkman_relative) .or. &
       .not. ieee_is_finite(petrov_lift_error) .or. &
       .not. ieee_is_finite(petrov_lift_relative)) &
       call neko_error('CPU assembly regression produced a nonfinite error metric')
  call MPI_Comm_size(NEKO_COMM, nranks, ierr)
  if (ierr /= 0) call neko_error('MPI_Comm_size failed in CPU assembly diagnostic')
  if (pe_rank == 0) then
     write(*,'(A,I0,A,I0)') 'CPU packed-vs-legacy assembly on ranks=', nranks, &
          ', local_points=', n
     write(*,'(A,ES24.16,A,ES24.16)') 'mass max_abs=', mass_error, &
          ', relative_frobenius=', mass_relative
     write(*,'(A,ES24.16,A,ES24.16)') 'diffusion max_abs=', diffusion_error, &
          ', relative_frobenius=', diffusion_relative
     write(*,'(A,ES24.16,A,ES24.16)') 'brinkman max_abs=', brinkman_error, &
          ', relative_frobenius=', brinkman_relative
     write(*,'(A,ES24.16,A,ES24.16)') 'brinkman_lift max_abs=', brinkman_lift_error, &
          ', relative_frobenius=', brinkman_lift_relative
     write(*,'(A,ES24.16,A,ES24.16)') 'convection max_abs=', convection_error, &
          ', relative_frobenius=', convection_relative
     if (run_gemm_benchmark) then
        write(*,'(A,I0,A,I0,A,I0)') 'DGEMM benchmark repeats=', benchmark_repeats, &
             ', block=', gemm_block_size, ', local_points=', n
        write(*,'(A,ES24.16,A,ES24.16)') 'DGEMM-vs-packed max_abs=', &
             dgemm_packed_error, ', relative_frobenius=', dgemm_packed_relative
        write(*,'(A,ES24.16,A,ES24.16)') 'DGEMM-vs-legacy max_abs=', &
             dgemm_legacy_error, ', relative_frobenius=', dgemm_legacy_relative
        write(*,'(A,ES24.16,A,ES24.16)') 'packed seconds mean/sample_sd=', &
             packed_mean, ', ', packed_spread
        write(*,'(A,ES24.16,A,ES24.16)') 'DGEMM seconds mean/sample_sd=', &
             dgemm_mean, ', ', dgemm_spread
        write(*,'(A,I0,A,I0,A,I0,A)') 'DGEMM scratch bound: pair_panel=', &
             n*gemm_block_size, ' rp values; weighted_test_panel=', n*rank_basis, &
             ' rp values; accumulator=', rank_basis*gemm_block_size, ' rp values'
     end if
     write(*,'(A,ES24.16,A,ES24.16)') 'scalar_linear max_abs=', scalar_linear_error, &
          ', relative_frobenius=', scalar_linear_relative
     write(*,'(A,ES24.16,A,ES24.16)') 'scalar_constant max_abs=', scalar_constant_error, &
          ', relative_frobenius=', scalar_constant_relative
     write(*,'(A,ES24.16,A,ES24.16)') 'scalar_quadratic max_abs=', scalar_quadratic_error, &
          ', relative_frobenius=', scalar_quadratic_relative
     write(*,'(A,ES24.16,A,ES24.16)') 'scalar_implicit max_abs=', scalar_implicit_error, &
          ', relative_frobenius=', scalar_implicit_relative
     write(*,'(A,ES24.16,A,ES24.16)') 'petrov_cross_mass max_abs=', petrov_mass_error, &
          ', relative_frobenius=', petrov_mass_relative
     write(*,'(A,ES24.16,A,ES24.16)') 'petrov_brinkman max_abs=', petrov_brinkman_error, &
          ', relative_frobenius=', petrov_brinkman_relative
     write(*,'(A,ES24.16,A,ES24.16)') 'petrov_brinkman_lift max_abs=', petrov_lift_error, &
          ', relative_frobenius=', petrov_lift_relative
     write(*,'(A,ES12.4)') 'comparison_tolerance=', comparison_tolerance
     if (write_csv) then
        open(newunit=csv_unit, file=trim(csv_file), status='replace', &
             action='write', iostat=io_status)
        if (io_status /= 0) call neko_error('Could not open CPU assembly CSV output')
        write(csv_unit,'(A)') 'operator,implementation,i,j,k,value'
        call write_matrix_pair('mass',mass_packed,mass_legacy)
        call write_matrix_pair('diffusion',diffusion_packed,diffusion_legacy)
        call write_matrix_pair('brinkman',brinkman_packed,brinkman_legacy)
        call write_vector_pair('brinkman_lift',brinkman_lift_packed,brinkman_lift_legacy)
        call write_tensor_pair('convection',convection_packed,convection_legacy)
        if (run_gemm_benchmark) &
             call write_tensor_single('convection', 'blocked_dgemm', convection_dgemm)
        call write_matrix_pair('scalar_linear',scalar_linear_packed,scalar_linear_legacy)
        call write_vector_pair('scalar_constant',scalar_constant_packed,scalar_constant_legacy)
        call write_tensor_pair('scalar_quadratic',scalar_quadratic_packed,scalar_quadratic_legacy)
        call write_matrix_pair('scalar_implicit',scalar_implicit_packed,scalar_implicit_legacy)
        call write_matrix_pair('petrov_cross_mass',petrov_mass_packed,petrov_mass_legacy)
        call write_matrix_pair('petrov_brinkman',petrov_brinkman_packed,petrov_brinkman_legacy)
        call write_vector_pair('petrov_brinkman_lift',petrov_lift_packed,petrov_lift_legacy)
        close(csv_unit)
     end if
  end if
  if (mass_relative > comparison_tolerance .or. &
       diffusion_relative > comparison_tolerance .or. &
       brinkman_relative > comparison_tolerance .or. &
       brinkman_lift_relative > comparison_tolerance .or. &
       convection_relative > comparison_tolerance .or. &
       scalar_linear_relative > comparison_tolerance .or. &
       scalar_constant_relative > comparison_tolerance .or. &
       scalar_quadratic_relative > comparison_tolerance .or. &
       scalar_implicit_relative > comparison_tolerance .or. &
       petrov_mass_relative > comparison_tolerance .or. &
       petrov_brinkman_relative > comparison_tolerance .or. &
       petrov_lift_relative > comparison_tolerance) &
       call neko_error('CPU packed assembly regression gate failed')

  call u%free(); call v%free(); call w%free()
  call gx%free(); call gy%free(); call gz%free()
  call simulation%free()
  call neko_finalize()

contains

  subroutine gradient_from_host(input, gx, gy, gz, coef)
    type(field_t), intent(inout) :: input, gx, gy, gz
    type(coef_t), intent(in) :: coef
    if (NEKO_BCKND_DEVICE .eq. 1) call input%copy_from(HOST_TO_DEVICE, sync=.true.)
    call grad(gx%x, gy%x, gz%x, input%x, coef)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call gx%copy_from(DEVICE_TO_HOST, sync=.true.)
       call gy%copy_from(DEVICE_TO_HOST, sync=.true.)
       call gz%copy_from(DEVICE_TO_HOST, sync=.true.)
    end if
  end subroutine gradient_from_host

  subroutine write_matrix_pair(label, packed, legacy)
    character(len=*), intent(in) :: label
    real(kind=rp), intent(in) :: packed(:, :), legacy(:, :)
    integer :: ii, jj
    do ii = 1, size(packed,1)
       do jj = 1, size(packed,2)
          call write_value(label,'packed',ii,jj,0,packed(ii,jj))
       end do
    end do
    do ii = 1, size(legacy,1)
       do jj = 1, size(legacy,2)
          call write_value(label,'legacy',ii,jj,0,legacy(ii,jj))
       end do
    end do
  end subroutine write_matrix_pair

  subroutine write_vector_pair(label, packed, legacy)
    character(len=*), intent(in) :: label
    real(kind=rp), intent(in) :: packed(:), legacy(:)
    integer :: ii
    do ii = 1, size(packed)
       call write_value(label,'packed',ii,0,0,packed(ii))
    end do
    do ii = 1, size(legacy)
       call write_value(label,'legacy',ii,0,0,legacy(ii))
    end do
  end subroutine write_vector_pair

  subroutine write_tensor_pair(label, packed, legacy)
    character(len=*), intent(in) :: label
    real(kind=rp), intent(in) :: packed(:, :, :), legacy(:, :, :)
    integer :: ii, jj, kk
    do ii = 1, size(packed,1)
       do jj = 1, size(packed,2)
          do kk = 1, size(packed,3)
             call write_value(label,'packed',ii,jj,kk,packed(ii,jj,kk))
          end do
       end do
    end do
    do ii = 1, size(legacy,1)
       do jj = 1, size(legacy,2)
          do kk = 1, size(legacy,3)
             call write_value(label,'legacy',ii,jj,kk,legacy(ii,jj,kk))
          end do
       end do
    end do
  end subroutine write_tensor_pair

  subroutine write_tensor_single(label, implementation, values)
    character(len=*), intent(in) :: label, implementation
    real(kind=rp), intent(in) :: values(:, :, :)
    integer :: ii, jj, kk
    do ii = 1, size(values,1)
       do jj = 1, size(values,2)
          do kk = 1, size(values,3)
             call write_value(label,implementation,ii,jj,kk,values(ii,jj,kk))
          end do
       end do
    end do
  end subroutine write_tensor_single

  subroutine write_value(label, implementation, ii, jj, kk, value)
    character(len=*), intent(in) :: label, implementation
    integer, intent(in) :: ii, jj, kk
    real(kind=rp), intent(in) :: value
    write(csv_unit,'(A,",",A,",",I0,",",I0,",",I0,",",ES24.16E3)') &
         trim(label), trim(implementation), ii, jj, kk, value
  end subroutine write_value

end program rom_cpu_assembly_regression_diagnostic
