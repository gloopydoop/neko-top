program rom_forward_continuity_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use simulation_m, only: simulation_t
  use fluid_pnpn, only: fluid_pnpn_t
  use file, only: file_t
  use fld_file_data, only: fld_file_data_t
  use field, only: field_t
  use field_math, only: field_col3, field_rzero
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2, glsc3
  use operators, only: div, ortho
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use num_types, only: rp
  use comm, only: pe_rank
  use utils, only: neko_error
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
  implicit none

  type(simulation_t) :: simulation
  type(fluid_pnpn_t), pointer :: fluid
  type(file_t) :: snapshot_file, modes_file
  type(fld_file_data_t) :: data
  type(json_file) :: parameters
  type(field_t) :: u, v, w, p, p_res, u_zero, v_zero, w_zero
  type(field_t) :: f_x, f_y, f_z, divergence
  real(kind=rp) :: lift(3) = [1.0_rp, 0.0_rp, 0.0_rp]
  real(kind=rp) :: physical_residual, homogeneous_residual
  real(kind=rp) :: physical_divergence, homogeneous_divergence
  integer :: argc, n, snapshot, n_snapshots, n_modes, mode, unit
  character(len=1024) :: case_file, snapshots_path, modes_path, output_file
  logical :: found
  type(c_ptr) :: event

  argc = command_argument_count()
  if (argc < 3) call neko_error( &
       'Usage: forward_continuity_diagnostic CASE SNAPSHOTS.fld MODES.fld '&
       //'[n_snapshots] [CSV]')
  call get_command_argument(1, case_file)
  call get_command_argument(2, snapshots_path)
  call get_command_argument(3, modes_path)
  n_snapshots = 17
  if (argc >= 4) then
     call get_command_argument(4, output_file)
     read(output_file, *) n_snapshots
  end if
  output_file = 'forward_continuity_diagnostics.csv'
  if (argc >= 5) call get_command_argument(5, output_file)

  call neko_init()
  call neko_top_register_types()
  parameters = json_read_file(trim(case_file))
  call parameters%update('unsteady', .false., found)
  call simulation%init(parameters)
  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('Forward continuity diagnostic requires a Pn/Pn fluid')
  end select

  n = fluid%u%dof%size()
  call u%init(fluid%dm_Xh, 'continuity_u')
  call v%init(fluid%dm_Xh, 'continuity_v')
  call w%init(fluid%dm_Xh, 'continuity_w')
  call p%init(fluid%dm_Xh, 'continuity_p')
  call p_res%init(fluid%dm_Xh, 'continuity_p_res')
  call u_zero%init(fluid%dm_Xh, 'continuity_u_zero')
  call v_zero%init(fluid%dm_Xh, 'continuity_v_zero')
  call w_zero%init(fluid%dm_Xh, 'continuity_w_zero')
  call f_x%init(fluid%dm_Xh, 'continuity_f_x')
  call f_y%init(fluid%dm_Xh, 'continuity_f_y')
  call f_z%init(fluid%dm_Xh, 'continuity_f_z')
  call divergence%init(fluid%dm_Xh, 'continuity_divergence')
  call field_rzero(u_zero); call field_rzero(v_zero); call field_rzero(w_zero)
  event = C_NULL_PTR

  call data%init()
  call modes_file%init(trim(modes_path))
  call modes_file%set_counter(0)
  call modes_file%read(data)
  if (data%n_scalars < 3 .or. mod(data%n_scalars, 3) /= 0) then
     call neko_error('Mode file must contain 3*N velocity fields')
  end if
  n_modes = data%n_scalars / 3

  if (pe_rank == 0) then
     open(newunit=unit, file=trim(output_file), status='replace', action='write')
     write(unit, '(A)') 'kind,index,weak_residual,strong_divergence'
  end if

  call snapshot_file%init(trim(snapshots_path))
  do snapshot = 0, n_snapshots - 1
     call snapshot_file%set_counter(snapshot)
     call data%free()
     call data%init()
     call snapshot_file%read(data)
     if (data%n_scalars < 3) call neko_error('Snapshot file has fewer than 3 fields')
     call load_fields(fluid, data, u, v, w, n)
     call measure_fields('snapshot_physical', snapshot, u, v, w)
     u%x = u%x - lift(1); v%x = v%x - lift(2); w%x = w%x - lift(3)
     call measure_fields('snapshot_homogeneous', snapshot, u, v, w)
     call data%free()
  end do

  call data%free()
  call data%init()
  call modes_file%set_counter(0)
  call modes_file%read(data)
  do mode = 0, n_modes - 1
     call load_fields_at(fluid, data, u, v, w, n, 3 * mode)
     call fluid%bcs_vel_projector%apply(u%x, v%x, w%x, n)
     call measure_fields('pod_mode_restricted', mode, u, v, w)
  end do

  if (pe_rank == 0) close(unit)
  call snapshot_file%free(); call modes_file%free()
  call data%free()
  call u%free(); call v%free(); call w%free(); call p%free(); call p_res%free()
  call u_zero%free(); call v_zero%free(); call w_zero%free()
  call f_x%free(); call f_y%free(); call f_z%free(); call divergence%free()
  call simulation%free()
  call neko_finalize()

contains

  !> Load one vector snapshot from a field-file record.
  subroutine load_fields(fluid, data, u, v, w, n)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(fld_file_data_t), intent(in) :: data
    type(field_t), intent(inout) :: u, v, w
    integer, intent(in) :: n
    call load_fields_at(fluid, data, u, v, w, n, 0)
  end subroutine load_fields

  !> Load one vector snapshot at a specified scalar-field offset.
  subroutine load_fields_at(fluid, data, u, v, w, n, offset)
    type(fluid_pnpn_t), intent(inout) :: fluid
    type(fld_file_data_t), intent(in) :: data
    type(field_t), intent(inout) :: u, v, w
    integer, intent(in) :: n, offset
    u%x = reshape(data%s(offset + 1)%x, shape(u%x))
    v%x = reshape(data%s(offset + 2)%x, shape(v%x))
    w%x = reshape(data%s(offset + 3)%x, shape(w%x))
    call fluid%gs_Xh%op(u%x, v%x, w%x, n, GS_OP_ADD)
    call col2(u%x, fluid%c_Xh%mult, n)
    call col2(v%x, fluid%c_Xh%mult, n)
    call col2(w%x, fluid%c_Xh%mult, n)
  end subroutine load_fields_at

  !> Record strong and weak continuity metrics for one field triple.
  subroutine measure_fields(kind, index, u, v, w)
    character(len=*), intent(in) :: kind
    integer, intent(in) :: index
    type(field_t), intent(inout) :: u, v, w
    call field_col3(f_x, u, fluid%rho)
    call field_col3(f_y, v, fluid%rho)
    call field_col3(f_z, w, fluid%rho)
    call col2(f_x%x, fluid%c_Xh%B, n)
    call col2(f_y%x, fluid%c_Xh%B, n)
    call col2(f_z%x, fluid%c_Xh%B, n)
    call field_rzero(p)
    ! Zero viscosity isolates the assembled continuity part of Neko's native
    ! pressure residual; this is not an independently constructed Laplacian.
    call fluid%prs_res%compute(p, p_res, u, v, w, u_zero, v_zero, w_zero, &
         f_x, f_y, f_z, fluid%c_Xh, fluid%gs_Xh, fluid%bc_prs_surface, &
         fluid%bc_sym_surface, fluid%Ax_prs, 0.0_rp, 1.0_rp, u_zero, &
         fluid%rho, event)
    if (.not. fluid%prs_dirichlet) call ortho(p_res%x, fluid%glb_n_points, n)
    call fluid%gs_Xh%op(p_res%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(p_res%x, n)
    homogeneous_residual = sqrt(max(0.0_rp, glsc3(p_res%x, &
         fluid%c_Xh%mult, p_res%x, n)))
    call div(divergence%x, u%x, v%x, w%x, fluid%c_Xh)
    homogeneous_divergence = sqrt(max(0.0_rp, glsc3(divergence%x, &
         fluid%c_Xh%B, divergence%x, n)))
    if (pe_rank == 0) write(unit, '(A,A,I0,2(",",ES24.16))') trim(kind), &
         ',', index, homogeneous_residual, homogeneous_divergence
  end subroutine measure_fields

end program rom_forward_continuity_diagnostic
