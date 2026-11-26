program rom_affine_continuity_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use fluid_pnpn, only: fluid_pnpn_t
  use file, only: file_t
  use fld_file_data, only: fld_file_data_t
  use field, only: field_t
  use field_math, only: field_rzero, field_col3
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2, glsc3
  use operators, only: ortho
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use num_types, only: rp
  use comm, only: pe_rank
  use utils, only: neko_error
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
  implicit none

  type(simulation_t) :: simulation
  type(fluid_pnpn_t), pointer :: fluid
  type(file_t) :: field_file
  type(fld_file_data_t) :: field_data
  type(json_file) :: parameters
  type(field_t) :: u, v, w, p, zero_u, zero_v, zero_w
  type(field_t) :: f_x, f_y, f_z
  type(field_t) :: lift_residual, column_residual, affine_residual
  type(field_t) :: direct_residual, difference_residual
  type(field_t) :: boundary_direct, boundary_lift, boundary_column
  type(field_t) :: boundary_affine, boundary_difference, boundary_flux
  type(field_t) :: boundary_reconstructed, boundary_reconstruction_difference
  type(field_t) :: zero_residual
  type(field_t) :: u_saved, v_saved, w_saved, u_reconstructed
  type(field_t) :: v_reconstructed, w_reconstructed, du, dv, dw
  type(json_file) :: reduced_metadata
  real(kind=rp), allocatable :: basis(:, :), lift(:), coefficients(:)
  character(len=1024) :: case_path, outer_path, output_path
  character(len=1024) :: reduced_path, basis_path, lift_path, coeff_path
  character(len=1024) :: field_prefix
  integer :: argc, n, rank, metadata_rank, metadata_rows, basis_rows, bytes
  integer :: lift_bytes, unit, frame, mode
  integer :: output_unit, ios
  integer :: coeff_step
  logical :: found, has_scalar, normalized_boundary
  real(kind=rp) :: coefficient_time, frame_time, residual_coeff_time
  real(kind=rp) :: norm_lift, norm_direct, norm_affine, norm_difference
  real(kind=rp) :: component_norm, scaled_norm_sum, boundary_scaled_norm_sum
  real(kind=rp) :: cancellation_ratio, cosine, boundary_cosine
  real(kind=rp) :: denominator, reconstruction_error, boundary_affine_defect
  type(c_ptr) :: event

  argc = command_argument_count()
  if (argc < 2) call neko_error( &
       'Usage: rom_affine_continuity_diagnostic CASE OUTER_DIR [CSV] '// &
       '[normalized_boundary]')
  call get_command_argument(1, case_path)
  call get_command_argument(2, outer_path)
  output_path = 'rom_affine_continuity.csv'
  if (argc >= 3) call get_command_argument(3, output_path)
  normalized_boundary = .false.
  if (argc >= 4) then
     call get_command_argument(4, field_prefix)
     if (trim(field_prefix) /= 'normalized_boundary') call neko_error( &
          'Unknown affine continuity diagnostic option')
     normalized_boundary = .true.
  end if

  reduced_path = trim(outer_path)//'/reduced_system'
  basis_path = trim(reduced_path)//'/trial_basis.bin'
  lift_path = trim(reduced_path)//'/lift_velocity.bin'
  coeff_path = trim(outer_path)//'/forward_coefficients_inner_000.csv'
  field_prefix = trim(outer_path)//'/forward_fields0.fld'

  call neko_init()
  call neko_top_register_types()
  parameters = json_read_file(trim(case_path))
  call parameters%update('unsteady', .false., found)
  call user_setup(simulation%neko_case%user)
  call simulation%init(parameters)
  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('Affine continuity diagnostic requires a Pn/Pn fluid')
  end select

  n = fluid%u%dof%size()
  reduced_metadata = json_read_file(trim(reduced_path)//'/metadata.json')
  call reduced_metadata%get('trial_basis_rank', metadata_rank, found)
  if (.not. found) call neko_error('Reduced-system metadata lacks trial_basis_rank')
  call reduced_metadata%get('trial_basis_rows', metadata_rows, found)
  if (.not. found) call neko_error('Reduced-system metadata lacks trial_basis_rows')
  call reduced_metadata%get('has_scalar', has_scalar, found)
  if (.not. found) call neko_error('Reduced-system metadata lacks has_scalar')
  rank = count_csv_rows(trim(reduced_path)//'/H.csv')
  if (rank < 1 .or. rank /= metadata_rank) call neko_error( &
       'H.csv row count disagrees with trial_basis_rank metadata')
  inquire(file=trim(basis_path), size=bytes, iostat=ios)
  if (ios /= 0 .or. bytes <= 0 .or. mod(bytes, 8 * rank) /= 0) then
     call neko_error('trial_basis.bin size is inconsistent with metadata rank')
  end if
  basis_rows = bytes / (8 * rank)
  if (basis_rows /= metadata_rows) call neko_error( &
       'trial_basis.bin row count disagrees with metadata')
  if (has_scalar) then
     if (basis_rows /= 4 * n) call neko_error( &
          'Scalar trial basis must have 3*n velocity rows followed by n scalar rows')
  else
     if (basis_rows /= 3 * n) call neko_error( &
          'Velocity-only trial basis must have exactly 3*n rows')
  end if
  inquire(file=trim(lift_path), size=lift_bytes, iostat=ios)
  if (ios /= 0 .or. lift_bytes /= 8 * 3 * n) call neko_error( &
       'lift_velocity.bin size is inconsistent with 3*n velocity rows')
  allocate(basis(basis_rows,rank), lift(3*n), coefficients(rank))
  open(newunit=unit, file=trim(basis_path), status='old', action='read', &
       access='stream', form='unformatted')
  read(unit) basis
  close(unit)
  open(newunit=unit, file=trim(lift_path), status='old', action='read', &
       access='stream', form='unformatted')
  read(unit) lift
  close(unit)
  call read_final_coefficients(trim(coeff_path), rank, coefficients, &
       coeff_step, coefficient_time)
  if (coeff_step /= 400) call neko_error( &
       'Expected step 400: output frame 0 is the final state of inner_000')

  call u%init(fluid%dm_Xh, 'affine_u')
  call v%init(fluid%dm_Xh, 'affine_v')
  call w%init(fluid%dm_Xh, 'affine_w')
  call p%init(fluid%dm_Xh, 'affine_p')
  call zero_u%init(fluid%dm_Xh, 'affine_zero_u')
  call zero_v%init(fluid%dm_Xh, 'affine_zero_v')
  call zero_w%init(fluid%dm_Xh, 'affine_zero_w')
  call f_x%init(fluid%dm_Xh, 'affine_f_x')
  call f_y%init(fluid%dm_Xh, 'affine_f_y')
  call f_z%init(fluid%dm_Xh, 'affine_f_z')
  call lift_residual%init(fluid%dm_Xh, 'affine_lift_residual')
  call column_residual%init(fluid%dm_Xh, 'affine_column_residual')
  call affine_residual%init(fluid%dm_Xh, 'affine_sum_residual')
  call direct_residual%init(fluid%dm_Xh, 'affine_direct_residual')
  call difference_residual%init(fluid%dm_Xh, 'affine_difference_residual')
  call boundary_direct%init(fluid%dm_Xh, 'affine_boundary_direct')
  call boundary_lift%init(fluid%dm_Xh, 'affine_boundary_lift')
  call boundary_column%init(fluid%dm_Xh, 'affine_boundary_column')
  call boundary_affine%init(fluid%dm_Xh, 'affine_boundary_sum')
  call boundary_difference%init(fluid%dm_Xh, 'affine_boundary_difference')
  call boundary_flux%init(fluid%dm_Xh, 'affine_boundary_flux')
  call boundary_reconstructed%init(fluid%dm_Xh, 'affine_boundary_reconstructed')
  call boundary_reconstruction_difference%init(fluid%dm_Xh, &
       'affine_boundary_reconstruction_difference')
  call zero_residual%init(fluid%dm_Xh, 'affine_zero_residual')
  call u_saved%init(fluid%dm_Xh, 'affine_saved_u')
  call v_saved%init(fluid%dm_Xh, 'affine_saved_v')
  call w_saved%init(fluid%dm_Xh, 'affine_saved_w')
  call u_reconstructed%init(fluid%dm_Xh, 'affine_reconstructed_u')
  call v_reconstructed%init(fluid%dm_Xh, 'affine_reconstructed_v')
  call w_reconstructed%init(fluid%dm_Xh, 'affine_reconstructed_w')
  call du%init(fluid%dm_Xh, 'affine_du')
  call dv%init(fluid%dm_Xh, 'affine_dv')
  call dw%init(fluid%dm_Xh, 'affine_dw')
  call field_rzero(zero_u); call field_rzero(zero_v); call field_rzero(zero_w)

  call field_file%init(trim(field_prefix))
  call field_data%init()
  event = C_NULL_PTR
  if (pe_rank == 0) then
     open(newunit=output_unit, file=trim(output_path), status='replace', &
          action='write')
     write(output_unit, '(A)') &
          'candidate_index,candidate_frame_label,coefficient_physical_time,'// &
          'metric,mode,coefficient,value'
     if (normalized_boundary) write(*,'(A)') &
          'Affine diagnostic boundary path: normalized synthetic operator '// &
          '(bd=1, dt=1, force=rho*B*u); saved fields are ROM candidate states.'
  end if

  do frame = 0, 1
     call field_file%set_counter(frame)
     call field_data%free()
     call field_data%init()
     call field_file%read(field_data)
     if (field_data%n_scalars < 3) call neko_error( &
          'forward_fields0 frame contains fewer than three velocity fields')
     frame_time = real(field_data%time, rp)
     residual_coeff_time = -1.0_rp
     if (frame == 0) residual_coeff_time = coefficient_time
     call load_saved_velocity()
     call weak_residual(u_saved, v_saved, w_saved, direct_residual, &
          0.0_rp, 1.0_rp, 1.0_rp)
     norm_direct = native_norm(direct_residual)
     if (pe_rank == 0) call write_row(frame, frame_time, -1.0_rp, &
          'direct_saved_state_norm', -1, 0.0_rp, norm_direct)

     if (normalized_boundary) then
        call write_row(frame, frame_time, residual_coeff_time, &
             'normalized_operator_bd', -1, 1.0_rp, 1.0_rp)
        call write_row(frame, frame_time, residual_coeff_time, &
             'normalized_operator_dt', -1, 1.0_rp, 1.0_rp)
        call write_row(frame, frame_time, residual_coeff_time, &
             'normalized_operator_bd_over_dt', -1, 1.0_rp, 1.0_rp)
        call write_row(frame, frame_time, residual_coeff_time, &
             'normalized_operator_force_scale', -1, 1.0_rp, 1.0_rp)
        call weak_residual(u_saved, v_saved, w_saved, boundary_direct, &
             1.0_rp, 1.0_rp, 1.0_rp)
        call write_row(frame, frame_time, residual_coeff_time, &
             'normalized_boundary_direct_saved_rom_state_norm', -1, 1.0_rp, &
             native_norm(boundary_direct))
        boundary_flux%x = boundary_direct%x - direct_residual%x
        call write_row(frame, frame_time, residual_coeff_time, &
             'normalized_pressure_surface_flux_direct_norm', -1, 1.0_rp, &
             native_norm(boundary_flux))
     end if

     if (frame == 0) then
        call reconstruct_velocity(coefficients)
        call report_reconstruction(frame, frame_time, coefficient_time)

        if (normalized_boundary) then
           call weak_residual(u_reconstructed, v_reconstructed, w_reconstructed, &
                boundary_reconstructed, 1.0_rp, 1.0_rp, 1.0_rp)
           boundary_difference%x = boundary_reconstructed%x - boundary_direct%x
           call write_row(frame, frame_time, coefficient_time, &
                'normalized_boundary_reconstructed_residual_norm', -1, 1.0_rp, &
                native_norm(boundary_reconstructed))
           call write_row(frame, frame_time, coefficient_time, &
                'normalized_boundary_reconstructed_direct_relative_error', &
                -1, 1.0_rp, native_norm(boundary_difference) / &
                max(native_norm(boundary_direct), tiny(1.0_rp)))
        end if

        call load_velocity_vector(lift, u, v, w)
        call weak_residual(u, v, w, lift_residual, 0.0_rp, 1.0_rp, 1.0_rp)
        affine_residual%x = lift_residual%x
        if (normalized_boundary) then
           call weak_residual(u, v, w, boundary_lift, 1.0_rp, 1.0_rp, 1.0_rp)
           boundary_affine%x = boundary_lift%x
           boundary_scaled_norm_sum = native_norm(boundary_lift)
        end if
        norm_lift = native_norm(lift_residual)
        scaled_norm_sum = norm_lift
        if (pe_rank == 0) call write_row(frame, frame_time, coefficient_time, &
             'lift_residual_norm', 0, 1.0_rp, norm_lift)

        do mode = 1, rank
           call load_basis_column(mode, u, v, w)
           call weak_residual(u, v, w, column_residual, 0.0_rp, 1.0_rp, 1.0_rp)
           if (normalized_boundary) then
              call weak_residual(u, v, w, boundary_column, &
                   1.0_rp, 1.0_rp, 1.0_rp)
              boundary_affine%x = boundary_affine%x + &
                   coefficients(mode) * boundary_column%x
              boundary_scaled_norm_sum = boundary_scaled_norm_sum + &
                   abs(coefficients(mode)) * native_norm(boundary_column)
           end if
           component_norm = native_norm(column_residual)
           scaled_norm_sum = scaled_norm_sum + &
                abs(coefficients(mode)) * component_norm
           affine_residual%x = affine_residual%x + &
                coefficients(mode) * column_residual%x
           if (pe_rank == 0) then
              call write_row(frame, frame_time, coefficient_time, &
                   'column_residual_norm', mode, coefficients(mode), component_norm)
              call write_row(frame, frame_time, coefficient_time, &
                   'coefficient_scaled_column_norm', mode, coefficients(mode), &
                   abs(coefficients(mode)) * component_norm)
           end if
        end do

        norm_affine = native_norm(affine_residual)
        difference_residual%x = affine_residual%x - direct_residual%x
        norm_difference = native_norm(difference_residual)
        cancellation_ratio = norm_affine / max(scaled_norm_sum, tiny(1.0_rp))
        denominator = norm_affine * norm_direct
        cosine = glsc3(affine_residual%x, fluid%c_Xh%mult, &
             direct_residual%x, n) / max(denominator, tiny(1.0_rp))
        if (normalized_boundary) then
           boundary_difference%x = boundary_affine%x - boundary_direct%x
           boundary_flux%x = boundary_affine%x - affine_residual%x
           call write_row(frame, frame_time, coefficient_time, &
                'normalized_boundary_affine_residual_norm', -1, 1.0_rp, &
                native_norm(boundary_affine))
           call write_row(frame, frame_time, coefficient_time, &
                'normalized_boundary_direct_affine_relative_error', -1, 1.0_rp, &
                native_norm(boundary_difference) / &
                max(native_norm(boundary_direct), &
                native_norm(boundary_affine), tiny(1.0_rp)))
           denominator = native_norm(boundary_affine) * native_norm(boundary_direct)
           boundary_cosine = glsc3(boundary_affine%x, fluid%c_Xh%mult, &
                boundary_direct%x, n) / max(denominator, tiny(1.0_rp))
           call write_row(frame, frame_time, coefficient_time, &
                'normalized_boundary_direct_affine_cosine', -1, 1.0_rp, &
                boundary_cosine)
           call write_row(frame, frame_time, coefficient_time, &
                'normalized_pressure_surface_flux_affine_norm', -1, 1.0_rp, &
                native_norm(boundary_flux))
           boundary_reconstruction_difference%x = boundary_reconstructed%x - &
                boundary_affine%x
           boundary_affine_defect = native_norm(boundary_reconstruction_difference) / &
                max(native_norm(boundary_reconstructed), boundary_scaled_norm_sum, &
                tiny(1.0_rp))
           call write_row(frame, frame_time, coefficient_time, &
                'normalized_boundary_reconstructed_affine_relative_defect', &
                -1, 1.0_rp, boundary_affine_defect)
           call field_rzero(u); call field_rzero(v); call field_rzero(w)
           call weak_residual(u, v, w, zero_residual, 1.0_rp, 1.0_rp, 1.0_rp)
           call write_row(frame, frame_time, coefficient_time, &
                'normalized_boundary_zero_input_residual_norm', -1, 1.0_rp, &
                native_norm(zero_residual))
        end if
        if (pe_rank == 0) then
           call write_row(frame, frame_time, coefficient_time, &
                'affine_sum_residual_norm', -1, 0.0_rp, norm_affine)
           call write_row(frame, frame_time, coefficient_time, &
                'affine_direct_difference_norm', -1, 0.0_rp, norm_difference)
           call write_row(frame, frame_time, coefficient_time, &
                'affine_direct_relative_error', -1, 0.0_rp, &
                norm_difference / max(norm_direct, tiny(1.0_rp)))
           call write_row(frame, frame_time, coefficient_time, &
                'affine_direct_cosine', -1, 0.0_rp, cosine)
           call write_row(frame, frame_time, coefficient_time, &
                'cancellation_ratio', -1, 0.0_rp, cancellation_ratio)
        end if
     else
        if (pe_rank == 0) call write_row(frame, frame_time, -1.0_rp, &
             'affine_decomposition_unavailable_inner_001_basis_not_exported', &
             -1, 0.0_rp, -1.0_rp)
     end if
  end do

  if (pe_rank == 0) close(output_unit)
  call field_file%free()
  call field_data%free()
  call u%free(); call v%free(); call w%free(); call p%free()
  call zero_u%free(); call zero_v%free(); call zero_w%free()
  call f_x%free(); call f_y%free(); call f_z%free()
  call lift_residual%free(); call column_residual%free()
  call affine_residual%free(); call direct_residual%free()
  call difference_residual%free()
  call boundary_direct%free(); call boundary_lift%free(); call boundary_column%free()
  call boundary_affine%free(); call boundary_difference%free(); call boundary_flux%free()
  call boundary_reconstructed%free(); call zero_residual%free()
  call boundary_reconstruction_difference%free()
  call u_saved%free(); call v_saved%free(); call w_saved%free()
  call u_reconstructed%free(); call v_reconstructed%free(); call w_reconstructed%free()
  call du%free(); call dv%free(); call dw%free()
  deallocate(basis, lift, coefficients)
  call simulation%free()
  call neko_finalize()

contains

  subroutine read_final_coefficients(filename, n_modes, values, final_step, final_time)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n_modes
    real(kind=rp), intent(out) :: values(n_modes)
    integer, intent(out) :: final_step
    real(kind=rp), intent(out) :: final_time
    integer :: file_unit, read_status, row_step, row_count
    real(kind=rp) :: row_time
    real(kind=rp) :: row_values(n_modes)
    character(len=2048) :: row
    open(newunit=file_unit, file=trim(filename), status='old', action='read')
    read(file_unit, '(A)', iostat=read_status) row
    if (read_status /= 0) call neko_error('Cannot read coefficient CSV header')
    row_count = 0
    do
       read(file_unit, '(A)', iostat=read_status) row
       if (read_status < 0) exit
       if (read_status /= 0) call neko_error('Error reading coefficient CSV')
       read(row, *, iostat=read_status) row_step, row_time, row_values
       if (read_status /= 0) call neko_error('Cannot parse coefficient CSV row')
       values = row_values
       final_step = row_step
       final_time = row_time
       row_count = row_count + 1
    end do
    close(file_unit)
    if (row_count == 0) call neko_error('Coefficient CSV contains no trajectory rows')
  end subroutine read_final_coefficients

  integer function count_csv_rows(filename)
    character(len=*), intent(in) :: filename
    integer :: file_unit, read_status
    character(len=2048) :: row
    count_csv_rows = 0
    open(newunit=file_unit, file=trim(filename), status='old', action='read')
    do
       read(file_unit, '(A)', iostat=read_status) row
       if (read_status < 0) exit
       if (read_status /= 0) call neko_error('Error counting H.csv rows')
       count_csv_rows = count_csv_rows + 1
    end do
    close(file_unit)
  end function count_csv_rows

  subroutine load_saved_velocity()
    u_saved%x = reshape(field_data%s(1)%x, shape(u_saved%x))
    v_saved%x = reshape(field_data%s(2)%x, shape(v_saved%x))
    w_saved%x = reshape(field_data%s(3)%x, shape(w_saved%x))
    call fluid%gs_Xh%op(u_saved%x, v_saved%x, w_saved%x, n, GS_OP_ADD)
    call col2(u_saved%x, fluid%c_Xh%mult, n)
    call col2(v_saved%x, fluid%c_Xh%mult, n)
    call col2(w_saved%x, fluid%c_Xh%mult, n)
  end subroutine load_saved_velocity

  subroutine reconstruct_velocity(a)
    real(kind=rp), intent(in) :: a(rank)
    u_reconstructed%x = reshape(lift(1:n), shape(u_reconstructed%x))
    v_reconstructed%x = reshape(lift(n+1:2*n), shape(v_reconstructed%x))
    w_reconstructed%x = reshape(lift(2*n+1:3*n), shape(w_reconstructed%x))
    do mode = 1, rank
       u_reconstructed%x = u_reconstructed%x + a(mode) * &
            reshape(basis(1:n,mode), shape(u_reconstructed%x))
       v_reconstructed%x = v_reconstructed%x + a(mode) * &
            reshape(basis(n+1:2*n,mode), shape(v_reconstructed%x))
       w_reconstructed%x = w_reconstructed%x + a(mode) * &
            reshape(basis(2*n+1:3*n,mode), shape(w_reconstructed%x))
    end do
  end subroutine reconstruct_velocity

  subroutine report_reconstruction(frame_id, time_field, time_coeff)
    integer, intent(in) :: frame_id
    real(kind=rp), intent(in) :: time_field, time_coeff
    real(kind=rp) :: numerator, denominator_norm
    du%x = u_reconstructed%x - u_saved%x
    dv%x = v_reconstructed%x - v_saved%x
    dw%x = w_reconstructed%x - w_saved%x
    numerator = sqrt(max(0.0_rp, &
         glsc3(du%x, fluid%c_Xh%B, du%x, n) + &
         glsc3(dv%x, fluid%c_Xh%B, dv%x, n) + &
         glsc3(dw%x, fluid%c_Xh%B, dw%x, n)))
    denominator_norm = sqrt(max(0.0_rp, &
         glsc3(u_saved%x, fluid%c_Xh%B, u_saved%x, n) + &
         glsc3(v_saved%x, fluid%c_Xh%B, v_saved%x, n) + &
         glsc3(w_saved%x, fluid%c_Xh%B, w_saved%x, n)))
    reconstruction_error = numerator / max(denominator_norm, tiny(1.0_rp))
    if (pe_rank == 0) call write_row(frame_id, time_field, time_coeff, &
         'reconstruction_relative_error', -1, 0.0_rp, reconstruction_error)
  end subroutine report_reconstruction

  subroutine load_velocity_vector(vector, a, b, c)
    real(kind=rp), intent(in) :: vector(3*n)
    type(field_t), intent(inout) :: a, b, c
    a%x = reshape(vector(1:n), shape(a%x))
    b%x = reshape(vector(n+1:2*n), shape(b%x))
    c%x = reshape(vector(2*n+1:3*n), shape(c%x))
  end subroutine load_velocity_vector

  subroutine load_basis_column(column, a, b, c)
    integer, intent(in) :: column
    type(field_t), intent(inout) :: a, b, c
    a%x = reshape(basis(1:n,column), shape(a%x))
    b%x = reshape(basis(n+1:2*n,column), shape(b%x))
    c%x = reshape(basis(2*n+1:3*n,column), shape(c%x))
  end subroutine load_basis_column

  subroutine weak_residual(a, b, c, result, bd, dt, force_scale)
    type(field_t), intent(inout) :: a, b, c, result
    real(kind=rp), intent(in) :: bd, dt, force_scale
    call field_col3(f_x, a, fluid%rho)
    call field_col3(f_y, b, fluid%rho)
    call field_col3(f_z, c, fluid%rho)
    call col2(f_x%x, fluid%c_Xh%B, n)
    call col2(f_y%x, fluid%c_Xh%B, n)
    call col2(f_z%x, fluid%c_Xh%B, n)
    f_x%x = force_scale * f_x%x
    f_y%x = force_scale * f_y%x
    f_z%x = force_scale * f_z%x
    call field_rzero(p)
    call field_rzero(result)
    call fluid%prs_res%compute(p, result, a, b, c, zero_u, zero_v, zero_w, &
         f_x, f_y, f_z, fluid%c_Xh, fluid%gs_Xh, fluid%bc_prs_surface, &
         fluid%bc_sym_surface, fluid%Ax_prs, bd, dt, zero_u, &
         fluid%rho, event)
    if (.not. fluid%prs_dirichlet) call ortho(result%x, fluid%glb_n_points, n)
    call fluid%gs_Xh%op(result%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(result%x, n)
  end subroutine weak_residual

  real(kind=rp) function native_norm(a)
    type(field_t), intent(in) :: a
    native_norm = sqrt(max(0.0_rp, glsc3(a%x, fluid%c_Xh%mult, a%x, n)))
  end function native_norm

  subroutine write_row(frame_id, time_field, time_coeff, metric, mode_id, coeff, value)
    integer, intent(in) :: frame_id, mode_id
    real(kind=rp), intent(in) :: time_field, time_coeff, coeff, value
    character(len=*), intent(in) :: metric
    if (pe_rank /= 0) return
    write(output_unit, '(I0,2(",",ES24.16),",",A,",",I0,2(",",ES24.16))') &
         frame_id, time_field, time_coeff, trim(metric), mode_id, coeff, value
  end subroutine write_row

end program rom_affine_continuity_diagnostic
