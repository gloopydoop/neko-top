program rom_objective_equivalence_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use simulation_m, only: simulation_t
  use fluid_pnpn, only: fluid_pnpn_t
  use file, only: file_t
  use fld_file_data, only: fld_file_data_t
  use field, only: field_t
  use math, only: col2
  use operators, only: grad
  use gather_scatter, only: GS_OP_ADD
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use json_utils, only: json_get
  use design, only: design_t, design_factory
  use brinkman_design, only: brinkman_design_t
  use viscous_dissipation_objective, only: viscous_dissipation_objective_t
  use brinkman_dissipation_objective, only: brinkman_dissipation_objective_t
  use registry, only: neko_registry
  use point_zone, only: point_zone_t
  use point_zone_registry, only: neko_point_zone_registry
  use rom_array_test_adapters, only: test_assemble_viscous_quadratic, &
       test_assemble_brinkman_objective, test_viscous_value, test_brinkman_value, &
       test_brinkman_density
  use vector, only: vector_t
  use comm, only: NEKO_COMM, pe_rank, MPI_REAL_PRECISION
  use mpi_f08, only: MPI_Allreduce, MPI_SUM, MPI_MAX, MPI_IN_PLACE
  use num_types, only: rp
  use utils, only: neko_error
  implicit none

  type(simulation_t) :: simulation
  type(fluid_pnpn_t), pointer :: fluid
  class(design_t), allocatable :: design
  type(viscous_dissipation_objective_t) :: viscous_objective
  type(brinkman_dissipation_objective_t) :: brinkman_objective
  type(file_t) :: modes_file
  type(fld_file_data_t) :: modes_data
  type(json_file) :: parameters, design_parameters
  type(field_t), pointer :: chi_field
  class(point_zone_t), pointer :: viscous_zone, brinkman_zone
  real(kind=rp), allocatable :: modes(:, :, :), gradients(:, :, :, :), weights(:)
  real(kind=rp), allocatable :: mask_v(:), mask_b(:), chi(:), chi_plus(:), chi_minus(:), &
       lift(:, :), state(:), state_direction(:), state_plus(:), state_minus(:), &
       chi_direction(:), chi_base(:), rom_density(:)
  real(kind=rp), allocatable :: viscous_quad(:, :), brinkman_quad(:, :), &
       brinkman_linear(:)
  real(kind=rp) :: brinkman_constant, rom_v, rom_b, neko_v, neko_b
  real(kind=rp) :: neko_v_plus, neko_v_minus, neko_b_plus, neko_b_minus
  real(kind=rp) :: volume_v, volume_b, err_v, err_b, global_err
  real(kind=rp) :: fd_map, analytic_map, eps_map, local_dot, global_dot
  real(kind=rp) :: chi_analytic, state_fd_v, state_fd_b
  real(kind=rp) :: state_analytic_v, state_analytic_b
  real(kind=rp) :: state_error_v, state_error_b, direct_fd, direct_analytic
  real(kind=rp) :: direct_error, direct_field_error, eps_state, eps_chi
  real(kind=rp) :: chi_sens_min, chi_sens_max, mapped_min, mapped_max
  integer :: argc, n, n_modes, r, mode, component, direction, i, j
  logical :: found
  character(len=1024) :: case_file, modes_path
  type(field_t) :: u, v, w, gx, gy, gz
  type(vector_t) :: x_base, x_plus, x_minus, direction_vector, chi_sensitivity, mapped_sensitivity

  argc = command_argument_count()
  if (argc < 2) call neko_error('Usage: rom_objective_equivalence_diagnostic CASE MODES.fld')
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
     call neko_error('Objective diagnostic requires a Pn/Pn fluid')
  end select

  call json_get(parameters, 'optimization.design', design_parameters)
  call design_factory(design, design_parameters, simulation)
  select type (design_in => design)
  type is (brinkman_design_t)
     continue
  class default
     call neko_error('Objective diagnostic requires a Brinkman design')
  end select
  call viscous_objective%init_from_attributes(design, simulation, 1.0_rp, &
       'viscous diagnostic', 'objective_domain', 0.0_rp, huge(0.0_rp))
  call brinkman_objective%init_from_attributes(design, simulation, 1.0_rp, &
       'Brinkman diagnostic', 'optimization_domain', .false., 0.0_rp, huge(0.0_rp))

  call modes_data%init()
  call modes_file%init(trim(modes_path))
  call modes_file%set_counter(0)
  call modes_file%read(modes_data)
  n = fluid%u%dof%size()
  if (modes_data%n_scalars < 3) call neko_error('Objective mode file is empty')
  n_modes = modes_data%n_scalars / 3
  r = min(3, n_modes)
  allocate(weights(n), modes(n,3,r), gradients(n,3,r,3), lift(n,3), state(r))
  allocate(mask_v(n), mask_b(n), chi(n), chi_plus(n), chi_minus(n), &
       viscous_quad(r,r), brinkman_quad(r,r), &
       brinkman_linear(r))
  weights = reshape(fluid%c_Xh%B, [n])
  ! Use a nontrivial physical lift.  It is spatially uniform, as required by
  ! the production forward_lift contract, but all three velocity components
  ! are nonzero so the Brinkman constant and cross terms are exercised.
  lift = 0.0_rp
  lift(:,1) = 0.80_rp
  lift(:,2) = -0.35_rp
  lift(:,3) = 0.20_rp
  state = [0.23_rp, -0.17_rp, 0.11_rp]
  call u%init(fluid%dm_Xh, 'objective_diag_u')
  call v%init(fluid%dm_Xh, 'objective_diag_v')
  call w%init(fluid%dm_Xh, 'objective_diag_w')
  call gx%init(fluid%dm_Xh, 'objective_diag_gx')
  call gy%init(fluid%dm_Xh, 'objective_diag_gy')
  call gz%init(fluid%dm_Xh, 'objective_diag_gz')
  do mode = 1, r
     u%x = reshape(modes_data%s(3*(mode-1)+1)%x, shape(u%x))
     v%x = reshape(modes_data%s(3*(mode-1)+2)%x, shape(v%x))
     w%x = reshape(modes_data%s(3*(mode-1)+3)%x, shape(w%x))
     call fluid%gs_Xh%op(u%x, v%x, w%x, n, GS_OP_ADD)
     call col2(u%x, fluid%c_Xh%mult, n); call col2(v%x, fluid%c_Xh%mult, n)
     call col2(w%x, fluid%c_Xh%mult, n)
     call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
     modes(:,1,mode) = reshape(u%x, [n])
     modes(:,2,mode) = reshape(v%x, [n])
     modes(:,3,mode) = reshape(w%x, [n])
     call grad(gx%x, gy%x, gz%x, u%x, fluid%c_Xh)
     gradients(:,1,mode,:) = reshape([reshape(gx%x,[n]), reshape(gy%x,[n]), &
          reshape(gz%x,[n])], [n,3])
     call grad(gx%x, gy%x, gz%x, v%x, fluid%c_Xh)
     gradients(:,2,mode,:) = reshape([reshape(gx%x,[n]), reshape(gy%x,[n]), &
          reshape(gz%x,[n])], [n,3])
     call grad(gx%x, gy%x, gz%x, w%x, fluid%c_Xh)
     gradients(:,3,mode,:) = reshape([reshape(gx%x,[n]), reshape(gy%x,[n]), &
          reshape(gz%x,[n])], [n,3])
  end do
  call modes_data%free(); call modes_file%free()

  viscous_zone => neko_point_zone_registry%get_point_zone('objective_domain')
  brinkman_zone => neko_point_zone_registry%get_point_zone('optimization_domain')
  mask_v = 0.0_rp; mask_b = 0.0_rp
  mask_v(viscous_zone%mask%get()) = 1.0_rp
  mask_b(brinkman_zone%mask%get()) = 1.0_rp
  chi_field => neko_registry%get_field('brinkman_amplitude')
  chi = reshape(chi_field%x, [n])
  volume_v = sum(weights * mask_v); volume_b = sum(weights * mask_b)
  call MPI_Allreduce(MPI_IN_PLACE, volume_v, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
  call MPI_Allreduce(MPI_IN_PLACE, volume_b, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
  call x_base%init(n)
  do i = 1, n
     ! Nonuniform design indicator used for both the original objective and
     ! the ROM mapped Brinkman coefficient.
     x_base%x(i) = 0.20_rp + 0.60_rp * (0.5_rp + 0.5_rp * sin(0.19_rp * &
          fluid%u%dof%x%x(i,1,1,1)) * cos(0.23_rp * fluid%u%dof%y%x(i,1,1,1)))
  end do
  call design%update_design(x_base)
  chi = reshape(chi_field%x, [n])
  call test_assemble_viscous_quadratic(weights, mask_v, real(fluid%mu%x(1,1,1,1),rp), &
       gradients, viscous_quad)
  call test_assemble_brinkman_objective(weights, mask_b, chi, modes, lift, &
       brinkman_quad, brinkman_linear, brinkman_constant)
  ! The test adapters use vector_glsc*, which are global Neko contractions.
  ! Do not reduce their assembled ROM objective terms a second time.

  u%x = reshape(lift(:,1) + matmul(modes(:,1,:), state), shape(u%x))
  v%x = reshape(lift(:,2) + matmul(modes(:,2,:), state), shape(v%x))
  w%x = reshape(lift(:,3) + matmul(modes(:,3,:), state), shape(w%x))
  fluid%u%x = u%x
  fluid%v%x = v%x
  fluid%w%x = w%x
  call viscous_objective%update_value(design)
  call brinkman_objective%update_value(design)
  neko_v = viscous_objective%get_value()
  neko_b = brinkman_objective%get_value()
  rom_v = test_viscous_value(state, viscous_quad, volume_v)
  rom_b = test_brinkman_value(state, brinkman_constant, brinkman_linear, &
       brinkman_quad, volume_b)
  err_v = abs(rom_v - neko_v); err_b = abs(rom_b - neko_b)
  call MPI_Allreduce(max(err_v,err_b), global_err, 1, MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM)
  if (pe_rank == 0) write(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
       'Objective equivalence viscous ROM/Neko/error = ', rom_v, ' / ', neko_v, ' / ', err_v, &
       '; Brinkman ROM/Neko/error = ', rom_b, ' / ', neko_b, ' / ', err_b
  if (global_err > 5.0e-12_rp) call neko_error('Objective equivalence gate failed')

  ! Compare derivatives with respect to reduced coefficients.  The original
  ! objective classes have no reduced-state API, so their directional
  ! derivatives are evaluated by centred perturbations of the identical
  ! reconstructed physical fields.  The ROM derivative is the analytic
  ! derivative of the quadratic reduced counterpart.
  allocate(state_direction(r), state_plus(r), state_minus(r))
  do i = 1, r
     state_direction(i) = 0.17_rp * sin(0.41_rp * real(i, rp))
  end do
  eps_state = 1.0e-6_rp
  state_plus = state + eps_state * state_direction
  state_minus = state - eps_state * state_direction
  call set_reconstructed_state(u, v, w, fluid, lift, modes, state_plus)
  call viscous_objective%update_value(design)
  call brinkman_objective%update_value(design)
  neko_v_plus = viscous_objective%get_value()
  neko_b_plus = brinkman_objective%get_value()
  call set_reconstructed_state(u, v, w, fluid, lift, modes, state_minus)
  call viscous_objective%update_value(design)
  call brinkman_objective%update_value(design)
  neko_v_minus = viscous_objective%get_value()
  neko_b_minus = brinkman_objective%get_value()
  state_fd_v = (neko_v_plus - neko_v_minus) / (2.0_rp * eps_state)
  state_fd_b = (neko_b_plus - neko_b_minus) / (2.0_rp * eps_state)
  state_analytic_v = dot_product(state_direction, matmul(viscous_quad, state)) / volume_v
  state_analytic_b = dot_product(state_direction, brinkman_linear + &
       matmul(brinkman_quad, state)) / volume_b
  state_error_v = abs(state_fd_v - state_analytic_v)
  state_error_b = abs(state_fd_b - state_analytic_b)
  call set_reconstructed_state(u, v, w, fluid, lift, modes, state)
  call viscous_objective%update_value(design)
  call brinkman_objective%update_value(design)
  if (pe_rank == 0) write(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
       'Reduced-state derivative viscous FD/ROM/error = ', state_fd_v, ' / ', &
       state_analytic_v, ' / ', state_error_v, '; Brinkman FD/ROM/error = ', &
       state_fd_b, ' / ', state_analytic_b, ' / ', state_error_b
  if (max(state_error_v, state_error_b) > 5.0e-7_rp * max(1.0_rp, &
       abs(state_fd_v), abs(state_fd_b))) call neko_error( &
       'Reduced-state objective derivative gate failed')

  call chi_sensitivity%init(n)
  call mapped_sensitivity%init(n)
  allocate(rom_density(n), chi_direction(n), chi_base(n))
  call brinkman_objective%update_sensitivity(design)
  call brinkman_objective%get_sensitivity(chi_sensitivity)
  call test_brinkman_density(weights, mask_b, lift, modes, state, volume_b, rom_density)
  direct_field_error = maxval(abs(rom_density - chi_sensitivity%x))
  call MPI_Allreduce(direct_field_error, global_err, 1, MPI_REAL_PRECISION, MPI_MAX, NEKO_COMM)
  direct_field_error = global_err

  ! Independently perturb the mapped Brinkman field itself.  This does not
  ! invoke the design mapping and therefore tests the direct coefficient
  ! derivative represented by the original objective sensitivity field.
  chi_base = reshape(chi_field%x, [n])
  do i = 1, n
     chi_direction(i) = 0.13_rp * cos(0.17_rp * real(i, rp))
  end do
  eps_chi = 1.0e-6_rp
  chi_field%x = reshape(chi_base + eps_chi * chi_direction, shape(chi_field%x))
  call brinkman_objective%update_value(design)
  neko_b_plus = brinkman_objective%get_value()
  chi_field%x = reshape(chi_base - eps_chi * chi_direction, shape(chi_field%x))
  call brinkman_objective%update_value(design)
  neko_b_minus = brinkman_objective%get_value()
  direct_fd = (neko_b_plus - neko_b_minus) / (2.0_rp * eps_chi)
  local_dot = sum(weights * rom_density * chi_direction)
  call MPI_Allreduce(local_dot, direct_analytic, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
  direct_error = abs(direct_fd - direct_analytic)
  chi_field%x = reshape(chi_base, shape(chi_field%x))
  call brinkman_objective%update_value(design)
  if (pe_rank == 0) write(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
       'Direct mapped Brinkman derivative FD/ROM/error = ', direct_fd, ' / ', &
       direct_analytic, ' / ', direct_error, '; pointwise max error = ', &
       direct_field_error
  if (direct_field_error > 5.0e-13_rp .or. direct_error > 5.0e-7_rp * &
       max(1.0_rp, abs(direct_fd))) call neko_error( &
       'Direct mapped Brinkman derivative gate failed')

  ! Exercise the production Brinkman mapping/sensitivity cascade at a fixed
  ! reconstructed velocity.  The directional derivative is compared with the
  ! original objective's dJ/dchi followed by design%map_backward and the
  ! Brinkman design's exact SEM-diagonal directional conversion.
  call x_plus%init(n); call x_minus%init(n)
  call direction_vector%init(n)
  do i = 1, n
     direction_vector%x(i) = 0.1_rp * sin(real(i, rp) * 0.37_rp)
  end do
  eps_map = 1.0e-6_rp
  x_plus%x = x_base%x + eps_map * direction_vector%x
  x_minus%x = x_base%x - eps_map * direction_vector%x
  call design%update_design(x_plus)
  chi_plus = reshape(chi_field%x, [n])
  call brinkman_objective%update_value(design)
  rom_b = brinkman_objective%get_value()
  call design%update_design(x_minus)
  chi_minus = reshape(chi_field%x, [n])
  call brinkman_objective%update_value(design)
  rom_v = brinkman_objective%get_value()
  fd_map = (rom_b - rom_v) / (2.0_rp * eps_map)
  call design%update_design(x_base)
  call brinkman_objective%update_sensitivity(design)
  call brinkman_objective%get_sensitivity(chi_sensitivity)
  chi_sens_min = minval(chi_sensitivity%x); chi_sens_max = maxval(chi_sensitivity%x)
  local_dot = sum(weights * chi_sensitivity%x * &
       (chi_plus - chi_minus) / (2.0_rp * eps_map))
  call MPI_Allreduce(local_dot, global_dot, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
  chi_analytic = global_dot
  call design%map_backward(chi_sensitivity)
  select type (design_in => design)
  type is (brinkman_design_t)
     call design_in%get_sensitivity(mapped_sensitivity)
     mapped_min = minval(mapped_sensitivity%x); mapped_max = maxval(mapped_sensitivity%x)
     call design_in%convert_to_directional_derivative(mapped_sensitivity)
  class default
     call neko_error('Objective mapping diagnostic lost Brinkman design')
  end select
  local_dot = sum(mapped_sensitivity%x * direction_vector%x)
  call MPI_Allreduce(local_dot, global_dot, 1, MPI_REAL_PRECISION, MPI_SUM, NEKO_COMM)
  analytic_map = global_dot
  if (pe_rank == 0) write(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
       'Mapped Brinkman FD/chi-chain/exact-directional/error = ', fd_map, ' / ', chi_analytic, &
       ' / ', analytic_map, ' / ', abs(fd_map - analytic_map)
  if (pe_rank == 0) write(*,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
       'Sensitivity ranges chi/raw-mass-gradient: ', chi_sens_min, ' / ', chi_sens_max, &
       ' / ', mapped_min, ' / ', mapped_max
  if (abs(fd_map - analytic_map) > 5.0e-7_rp * max(1.0_rp, abs(fd_map))) &
       call neko_error('Mapped Brinkman objective sensitivity gate failed')
  call design%update_design(x_base)
  call x_base%free(); call x_plus%free(); call x_minus%free(); call direction_vector%free()
  call chi_sensitivity%free(); call mapped_sensitivity%free()
  deallocate(state_direction, state_plus, state_minus, rom_density, chi_direction, chi_base)

  call viscous_objective%free(); call brinkman_objective%free(); call design%free()
  call u%free(); call v%free(); call w%free(); call gx%free(); call gy%free(); call gz%free()
  call simulation%free(); call neko_finalize()

contains

  !> Reconstruct a physical velocity from lift and reduced coefficients.
  subroutine set_reconstructed_state(u, v, w, fluid, lift, modes, state)
    type(field_t), intent(inout) :: u, v, w
    type(fluid_pnpn_t), intent(inout) :: fluid
    real(kind=rp), intent(in) :: lift(:, :), modes(:, :, :), state(:)
    integer :: n
    n = fluid%u%dof%size()
    u%x = reshape(lift(:,1) + matmul(modes(:,1,:), state), shape(u%x))
    v%x = reshape(lift(:,2) + matmul(modes(:,2,:), state), shape(v%x))
    w%x = reshape(lift(:,3) + matmul(modes(:,3,:), state), shape(w%x))
    fluid%u%x = u%x
    fluid%v%x = v%x
    fluid%w%x = w%x
  end subroutine set_reconstructed_state

end program rom_objective_equivalence_diagnostic
