program rom_residual_identity_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use simulation_m, only: simulation_t
  use file, only: file_t
  use fld_file_data, only: fld_file_data_t
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use fluid_pnpn, only: fluid_pnpn_t
  use field, only: field_t
  use field_math, only: field_rzero
  use operators, only: grad
  use gather_scatter, only: GS_OP_ADD
  use math, only: glsc3
  use rom_array_test_adapters, only: test_assemble_mass, test_assemble_brinkman, &
       test_assemble_lift_convection, test_assemble_diffusion, &
       test_assemble_convection
  use comm, only: NEKO_COMM, pe_rank, MPI_REAL_PRECISION
  use num_types, only: rp
  use mpi_f08, only: MPI_Allreduce, MPI_SUM, MPI_MAX, MPI_IN_PLACE
  use utils, only: neko_error
  implicit none

  type(simulation_t) :: simulation
  type(fluid_pnpn_t), pointer :: fluid
  type(file_t) :: modes_file
  type(fld_file_data_t) :: modes_data
  type(json_file) :: parameters
  real(kind=rp), allocatable :: modes(:, :, :), gradients(:, :, :, :)
  real(kind=rp), allocatable :: weights(:), lift(:, :), chi(:)
  real(kind=rp), allocatable :: mass(:, :), brinkman(:, :), brinkman_lift(:)
  real(kind=rp), allocatable :: lift_linear(:, :), diffusion(:, :), convection(:, :, :)
  real(kind=rp), allocatable :: expected(:, :), expected_mass(:, :), &
       expected_brinkman(:, :), expected_lift(:)
  real(kind=rp), allocatable :: u_state(:), velocity(:, :), work(:)
  real(kind=rp), allocatable :: direct(:)
  type(field_t) :: u, v, w, gx, gy, gz, axu, axv, axw, fx, fy, fz
  integer :: argc, n, n_modes, n_fields, r, mode, component, direction
  integer :: i, j, k, unit
  real(kind=rp) :: nu, rho, mu, state_norm, direct_norm, error_max, global_error
  real(kind=rp) :: assembly_error
  real(kind=rp) :: brinkman_error, global_brinkman_error
  logical :: found
  character(len=1024) :: case_file, modes_path

  argc = command_argument_count()
  if (argc < 2) call neko_error('Usage: rom_residual_identity_diagnostic CASE MODES.fld')
  call get_command_argument(1, case_file)
  call get_command_argument(2, modes_path)

  call neko_init()
  call neko_top_register_types()
  parameters = json_read_file(trim(case_file))
  call parameters%update('unsteady', .false., found)
  call simulation%init(parameters)
  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('Residual diagnostic requires a Pn/Pn fluid')
  end select

  call modes_data%init()
  call modes_file%init(trim(modes_path))
  call modes_file%set_counter(0)
  call modes_file%read(modes_data)
  if (modes_data%n_scalars < 3 .or. mod(modes_data%n_scalars, 3) /= 0) &
       call neko_error('Residual diagnostic mode file must contain 3*N fields')
  n = fluid%u%dof%size()
  if (modes_data%s(1)%size() /= n) call neko_error('Mode partition does not match case')
  n_modes = modes_data%n_scalars / 3
  r = min(3, n_modes)
  n_fields = 3

  allocate(weights(n), modes(n, 3, r), gradients(n, 3, r, 3), lift(n,3), chi(n))
  allocate(mass(r,r), brinkman(r,r), brinkman_lift(r), lift_linear(r,r), &
       diffusion(r,r), convection(r,r,r), expected(r,r), expected_mass(r,r), &
       expected_brinkman(r,r), expected_lift(r), u_state(r), velocity(n,3), &
       work(n), direct(r))
  weights = reshape(fluid%c_Xh%B, [n])
  lift = 0.0_rp
  lift(:,1) = 1.0_rp
  do i = 1, n
     chi(i) = 0.25_rp + 0.35_rp * real(mod(i, 17), rp) / 16.0_rp
  end do

  call u%init(fluid%dm_Xh, 'residual_diag_u')
  call v%init(fluid%dm_Xh, 'residual_diag_v')
  call w%init(fluid%dm_Xh, 'residual_diag_w')
  call gx%init(fluid%dm_Xh, 'residual_diag_gx')
  call gy%init(fluid%dm_Xh, 'residual_diag_gy')
  call gz%init(fluid%dm_Xh, 'residual_diag_gz')
  call axu%init(fluid%dm_Xh, 'residual_diag_axu')
  call axv%init(fluid%dm_Xh, 'residual_diag_axv')
  call axw%init(fluid%dm_Xh, 'residual_diag_axw')
  call fx%init(fluid%dm_Xh, 'residual_diag_fx')
  call fy%init(fluid%dm_Xh, 'residual_diag_fy')
  call fz%init(fluid%dm_Xh, 'residual_diag_fz')

  do mode = 1, r
     u%x = reshape(modes_data%s(3*(mode-1)+1)%x, shape(u%x))
     v%x = reshape(modes_data%s(3*(mode-1)+2)%x, shape(v%x))
     w%x = reshape(modes_data%s(3*(mode-1)+3)%x, shape(w%x))
     call fluid%gs_Xh%op(u%x, v%x, w%x, n, GS_OP_ADD)
     call col2_local(u%x, fluid%c_Xh%mult, n)
     call col2_local(v%x, fluid%c_Xh%mult, n)
     call col2_local(w%x, fluid%c_Xh%mult, n)
     call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
     modes(:,1,mode) = reshape(u%x, [n])
     modes(:,2,mode) = reshape(v%x, [n])
     modes(:,3,mode) = reshape(w%x, [n])
     call grad(gx%x, gy%x, gz%x, u%x, fluid%c_Xh)
     gradients(:,1,mode,1) = reshape(gx%x, [n])
     gradients(:,1,mode,2) = reshape(gy%x, [n])
     gradients(:,1,mode,3) = reshape(gz%x, [n])
     call grad(gx%x, gy%x, gz%x, v%x, fluid%c_Xh)
     gradients(:,2,mode,1) = reshape(gx%x, [n])
     gradients(:,2,mode,2) = reshape(gy%x, [n])
     gradients(:,2,mode,3) = reshape(gz%x, [n])
     call grad(gx%x, gy%x, gz%x, w%x, fluid%c_Xh)
     gradients(:,3,mode,1) = reshape(gx%x, [n])
     gradients(:,3,mode,2) = reshape(gy%x, [n])
     gradients(:,3,mode,3) = reshape(gz%x, [n])
  end do
  call modes_data%free()
  call modes_file%free()

  call test_assemble_mass(weights, modes, mass)
  call test_assemble_brinkman(weights, chi, modes, lift, brinkman, brinkman_lift)
  call test_assemble_lift_convection(weights, modes, gradients, lift, lift_linear)
  rho = real(fluid%rho%x(1,1,1,1), rp)
  mu = real(fluid%mu%x(1,1,1,1), rp)
  nu = mu / rho
  call test_assemble_diffusion(weights, nu, gradients, diffusion)
  call test_assemble_convection(weights, modes, gradients, convection)

  ! Direct coefficient contractions are intentionally written independently
  ! here.  The live Neko residual comparison below uses the same modes and
  ! Neko's native advection/Helmholtz kernels, not this direct expression.
  expected = 0.0_rp
  expected_mass = 0.0_rp
  expected_brinkman = 0.0_rp
  expected_lift = 0.0_rp
   do i = 1, r
      do j = 1, r
        do component = 1, 3
           expected_mass(i,j) = expected_mass(i,j) + sum(weights * &
                modes(:,component,i) * modes(:,component,j))
        end do
        do component = 1, 3
           expected(i,j) = expected(i,j) - sum(weights * lift(:,1) * &
                modes(:,component,i) * gradients(:,component,j,1))
           do direction = 2, 3
              expected(i,j) = expected(i,j) - sum(weights * lift(:,direction) * &
                   modes(:,component,i) * gradients(:,component,j,direction))
           end do
           expected(i,j) = expected(i,j) - nu * sum(weights * &
                gradients(:,component,i,1) * gradients(:,component,j,1))
           do direction = 2, 3
              expected(i,j) = expected(i,j) - nu * sum(weights * &
                   gradients(:,component,i,direction) * gradients(:,component,j,direction))
           end do
           expected_brinkman(i,j) = expected_brinkman(i,j) + sum(weights * chi * &
                modes(:,component,i) * modes(:,component,j))
        end do
      end do
      do component = 1, 3
         expected_lift(i) = expected_lift(i) + sum(weights * chi * &
              modes(:,component,i) * lift(:,component))
      end do
   end do
  ! The independent reference contractions above are local array sums,
  ! whereas the ROM assembly uses Neko's globally reduced vector_glsc*.
  call MPI_Allreduce(MPI_IN_PLACE, expected_mass, r*r, MPI_REAL_PRECISION, &
       MPI_SUM, NEKO_COMM)
  call MPI_Allreduce(MPI_IN_PLACE, expected, r*r, MPI_REAL_PRECISION, &
       MPI_SUM, NEKO_COMM)
  call MPI_Allreduce(MPI_IN_PLACE, expected_brinkman, r*r, &
       MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
  call MPI_Allreduce(MPI_IN_PLACE, expected_lift, r, MPI_REAL_PRECISION, &
       MPI_SUM, NEKO_COMM)
  error_max = maxval(abs(mass - expected_mass))
  error_max = max(error_max, maxval(abs(expected - lift_linear - diffusion)))
  error_max = max(error_max, maxval(abs(brinkman - expected_brinkman)))
  brinkman_error = maxval(abs(brinkman_lift - expected_lift))
  call MPI_Allreduce(error_max, global_error, 1, MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM)
  call MPI_Allreduce(brinkman_error, global_brinkman_error, 1, MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM)
  assembly_error = global_error
  if (global_error > 5.0e-13_rp .or. global_brinkman_error > 5.0e-13_rp) &
       call neko_error('ROM assembly direct contraction gate failed')

  ! Use a nonzero state and evaluate the actual Neko spatial residual.  The
  ! field source is Neko's native no-dealias advection kernel; Ax_vel is
  ! temporarily configured for diffusion only.  The physical result is
  ! divided by rho to match the ROM's mass-normalised residual.
  u_state = [0.23_rp, -0.17_rp, 0.11_rp]
  velocity(:,1) = lift(:,1)
  velocity(:,2) = lift(:,2)
  velocity(:,3) = lift(:,3)
  do mode = 1, r
     do component = 1, 3
        velocity(:,component) = velocity(:,component) + modes(:,component,mode) * u_state(mode)
     end do
  end do
  u%x = reshape(velocity(:,1), shape(u%x))
  v%x = reshape(velocity(:,2), shape(v%x))
  w%x = reshape(velocity(:,3), shape(w%x))
  call field_rzero(fx); call field_rzero(fy); call field_rzero(fz)
  call fluid%adv%compute(u, v, w, fx, fy, fz, fluid%Xh, fluid%c_Xh, n)
  ! The production velocity Helmholtz operator also carries the BDF mass
  ! coefficient.  This diagnostic isolates the spatial diffusion part.
  fluid%c_Xh%h1 = mu
  fluid%c_Xh%ifh2 = .false.
  call fluid%Ax_vel%compute_vector(axu%x, axv%x, axw%x, u%x, v%x, w%x, &
       fluid%c_Xh, fluid%msh, fluid%Xh)
  direct = 0.0_rp
  do i = 1, r
     do mode = 1, 3
        select case (mode)
        case (1)
           work = reshape(fx%x - axu%x, [n]) / rho
        case (2)
           work = reshape(fy%x - axv%x, [n]) / rho
        case (3)
           work = reshape(fz%x - axw%x, [n]) / rho
        end select
        direct(i) = direct(i) + sum(modes(:,mode,i) * work)
     end do
     do component = 1, 3
        direct(i) = direct(i) - sum(weights * chi * modes(:,component,i) * &
             velocity(:,component))
     end do
   end do
  call MPI_Allreduce(MPI_IN_PLACE, direct, r, MPI_REAL_PRECISION, MPI_SUM, &
       NEKO_COMM)
  expected_lift = matmul(lift_linear, u_state) + matmul(diffusion, u_state)
  expected_lift = expected_lift - brinkman_lift - matmul(brinkman, u_state)
  do i = 1, r
     do j = 1, r
        do k = 1, r
           expected_lift(i) = expected_lift(i) - convection(i,j,k) * &
                u_state(j) * u_state(k)
        end do
     end do
  end do
  direct_norm = sqrt(max(0.0_rp, sum(direct**2)))
  state_norm = sqrt(max(0.0_rp, sum(expected_lift**2)))
  error_max = maxval(abs(direct - expected_lift))
  call MPI_Allreduce(error_max, global_error, 1, MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM)
  if (pe_rank == 0) then
     write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') &
          'Projected Neko residual max/ROM norm/direct norm: ', global_error, ' / ', &
          state_norm, ' / ', direct_norm
  end if
  if (global_error > 5.0e-10_rp * max(1.0_rp, state_norm)) then
     call neko_error('Projected Neko residual identity gate failed')
  end if

  if (pe_rank == 0) write(*,'(A,ES12.4,A,ES12.4)') &
       'ROM assembly max/brinkman-lift contraction error: ', assembly_error, ' / ', &
       global_brinkman_error

  call u%free(); call v%free(); call w%free(); call gx%free(); call gy%free(); call gz%free()
  call axu%free(); call axv%free(); call axw%free(); call fx%free(); call fy%free(); call fz%free()
  call simulation%free()
  call neko_finalize()

contains

  !> Apply a local coefficient-wise multiplication for the diagnostic.
  subroutine col2_local(a, b, n)
    real(kind=rp), intent(inout) :: a(:,:,:,:)
    real(kind=rp), intent(in) :: b(:,:,:,:)
    integer, intent(in) :: n
    a = a * b
  end subroutine col2_local

end program rom_residual_identity_diagnostic
