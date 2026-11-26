program lift_weak_decomposition_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use fluid_pnpn, only: fluid_pnpn_t
  use field, only: field_t
  use field_math, only: field_rzero, field_copy, field_col3
  use facet_normal, only: facet_normal_t
  use pnpn_res_cpu, only: pnpn_prs_res_cpu_t
  use rom_lift, only: rom_lift_t
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use gather_scatter, only: GS_OP_ADD
  use operators, only: ortho
  use math, only: col2, glsc3, glsum
  use num_types, only: rp
  use comm, only: pe_rank
  use utils, only: neko_error
  use, intrinsic :: iso_c_binding, only: c_ptr, C_NULL_PTR
  implicit none

  type(simulation_t) :: simulation
  type(json_file) :: parameters
  type(fluid_pnpn_t), pointer :: fluid
  type(rom_lift_t) :: rom_lift
  type(facet_normal_t) :: empty_symmetry
  type(pnpn_prs_res_cpu_t) :: cpu_residual
  type(field_t) :: u, v, w, p, zero_u, zero_v, zero_w
  type(field_t) :: f_x, f_y, f_z
  type(field_t) :: volume_res, full_res, dispatch_res, symmetry_res
  type(field_t) :: recomposed_res, difference_res
  character(len=1024) :: case_file, output_file
  integer :: argc, n, unit
  logical :: found
  type(c_ptr) :: event

  argc = command_argument_count()
  if (argc < 1) call neko_error( &
       'Usage: lift_weak_decomposition_diagnostic CASE [CSV]')
  call get_command_argument(1, case_file)
  output_file = 'lift_weak_decomposition.csv'
  if (argc >= 2) call get_command_argument(2, output_file)

  call neko_init()
  call neko_top_register_types()
  parameters = json_read_file(trim(case_file))
  call parameters%update('unsteady', .false., found)
  call user_setup(simulation%neko_case%user)
  call simulation%init(parameters)
  select type (fluid_in => simulation%fluid)
  type is (fluid_pnpn_t)
     fluid => fluid_in
  class default
     call neko_error('Lift weak decomposition requires a Pn/Pn fluid')
  end select

  n = fluid%u%dof%size()
  call rom_lift%init(simulation%neko_case)
  if (.not. rom_lift%enabled) call neko_error( &
       'Lift weak decomposition requires the active ROM lift')

  call u%init(fluid%dm_Xh, 'lift_weak_u')
  call v%init(fluid%dm_Xh, 'lift_weak_v')
  call w%init(fluid%dm_Xh, 'lift_weak_w')
  call p%init(fluid%dm_Xh, 'lift_weak_p')
  call zero_u%init(fluid%dm_Xh, 'lift_weak_zero_u')
  call zero_v%init(fluid%dm_Xh, 'lift_weak_zero_v')
  call zero_w%init(fluid%dm_Xh, 'lift_weak_zero_w')
  call f_x%init(fluid%dm_Xh, 'lift_weak_f_x')
  call f_y%init(fluid%dm_Xh, 'lift_weak_f_y')
  call f_z%init(fluid%dm_Xh, 'lift_weak_f_z')
  call volume_res%init(fluid%dm_Xh, 'lift_weak_volume')
  call full_res%init(fluid%dm_Xh, 'lift_weak_full')
  call dispatch_res%init(fluid%dm_Xh, 'lift_weak_dispatch')
  call symmetry_res%init(fluid%dm_Xh, 'lift_weak_symmetry')
  call recomposed_res%init(fluid%dm_Xh, 'lift_weak_recomposed')
  call difference_res%init(fluid%dm_Xh, 'lift_weak_difference')
  call field_rzero(zero_u); call field_rzero(zero_v); call field_rzero(zero_w)

  ! The CPU routine's empty facet-normal contribution reads unique_mask(0).
  ! Use its normal base setup, then give it a valid zero-entry mask. Its surface
  ! loop is therefore the same implementation as the physical boundary object.
  call empty_symmetry%init_from_components(fluid%c_Xh)
  allocate(empty_symmetry%unique_mask(0:0))
  empty_symmetry%unique_mask(0) = 0

  event = C_NULL_PTR
  if (pe_rank == 0) then
     open(newunit=unit, file=trim(output_file), status='replace', action='write')
     write(unit, '(A)') 'probe,component,native_B_norm,native_mult_norm,' // &
          'mult_over_velocity_norm,signed_sum,B_pair_with_volume,' // &
          'relative_vector_error'
  end if

  call field_copy(u, rom_lift%velocity(1))
  call field_copy(v, rom_lift%velocity(2))
  call field_copy(w, rom_lift%velocity(3))
  call field_col3(f_x, u, fluid%rho)
  call field_col3(f_y, v, fluid%rho)
  call field_col3(f_z, w, fluid%rho)
  call col2(f_x%x, fluid%c_Xh%B, n)
  call col2(f_y%x, fluid%c_Xh%B, n)
  call col2(f_z%x, fluid%c_Xh%B, n)
  call decompose_probe('bc_applied_lift')

  call field_rzero(u); call field_rzero(v); call field_rzero(w)
  call field_rzero(f_x); call field_rzero(f_y); call field_rzero(f_z)
  call decompose_probe('zero_velocity')

  if (pe_rank == 0) close(unit)
  call empty_symmetry%free()
  call rom_lift%free()
  call u%free(); call v%free(); call w%free(); call p%free()
  call zero_u%free(); call zero_v%free(); call zero_w%free()
  call f_x%free(); call f_y%free(); call f_z%free()
  call volume_res%free(); call full_res%free(); call dispatch_res%free()
  call symmetry_res%free(); call recomposed_res%free(); call difference_res%free()
  call simulation%free()
  call neko_finalize()

contains

  subroutine decompose_probe(label)
    character(len=*), intent(in) :: label
    call field_rzero(p)
    call evaluate_cpu(.true., volume_res)
    call evaluate_cpu(.false., full_res)
    call postprocess(volume_res)
    call postprocess(full_res)
    call field_copy(symmetry_res, full_res)
    symmetry_res%x = full_res%x - volume_res%x
    call field_copy(recomposed_res, volume_res)
    recomposed_res%x = volume_res%x + symmetry_res%x

    call field_rzero(p)
    call field_rzero(dispatch_res)
    call fluid%prs_res%compute(p, dispatch_res, u, v, w, zero_u, zero_v, &
         zero_w, f_x, f_y, f_z, fluid%c_Xh, fluid%gs_Xh, &
         fluid%bc_prs_surface, fluid%bc_sym_surface, fluid%Ax_prs, &
         0.0_rp, 1.0_rp, zero_u, fluid%rho, event)
    call postprocess(dispatch_res)

    call emit(label, 'volume', volume_res, volume_res, 0.0_rp)
    call emit(label, 'symmetry_surface', symmetry_res, volume_res, 0.0_rp)
    call emit(label, 'recomposed_volume_plus_surface', recomposed_res, &
         volume_res, vector_difference(recomposed_res, full_res))
    call emit(label, 'cpu_full', full_res, volume_res, &
         vector_difference(full_res, dispatch_res))
    call emit(label, 'prs_res_dispatch', dispatch_res, volume_res, 0.0_rp)
  end subroutine decompose_probe

  subroutine evaluate_cpu(use_empty_symmetry, result)
    logical, intent(in) :: use_empty_symmetry
    type(field_t), intent(inout) :: result
    call field_rzero(p)
    call field_rzero(result)
    if (use_empty_symmetry) then
       call cpu_residual%compute(p, result, u, v, w, zero_u, zero_v, zero_w, &
            f_x, f_y, f_z, fluid%c_Xh, fluid%gs_Xh, &
            fluid%bc_prs_surface, empty_symmetry, fluid%Ax_prs, &
            0.0_rp, 1.0_rp, zero_u, fluid%rho, event)
    else
       call cpu_residual%compute(p, result, u, v, w, zero_u, zero_v, zero_w, &
            f_x, f_y, f_z, fluid%c_Xh, fluid%gs_Xh, &
            fluid%bc_prs_surface, fluid%bc_sym_surface, fluid%Ax_prs, &
            0.0_rp, 1.0_rp, zero_u, fluid%rho, event)
    end if
  end subroutine evaluate_cpu

  subroutine postprocess(result)
    type(field_t), intent(inout) :: result
    if (.not. fluid%prs_dirichlet) then
       call ortho(result%x, fluid%glb_n_points, n)
    end if
    call fluid%gs_Xh%op(result%x, n, GS_OP_ADD)
    call fluid%bcs_prs_projector%apply(result%x, n)
  end subroutine postprocess

  subroutine emit(probe, component, result, reference, relative_error)
    character(len=*), intent(in) :: probe, component
    type(field_t), intent(in) :: result, reference
    real(kind=rp), intent(in) :: relative_error
    real(kind=rp) :: norm_value, mult_norm, velocity_norm, signed_sum, pair_value
    norm_value = sqrt(max(0.0_rp, glsc3(result%x, fluid%c_Xh%B, &
         result%x, n)))
    mult_norm = sqrt(max(0.0_rp, glsc3(result%x, fluid%c_Xh%mult, &
         result%x, n)))
    velocity_norm = sqrt(max(0.0_rp, &
         glsc3(u%x, fluid%c_Xh%B, u%x, n) + &
         glsc3(v%x, fluid%c_Xh%B, v%x, n) + &
         glsc3(w%x, fluid%c_Xh%B, w%x, n)))
    signed_sum = glsum(result%x, n)
    pair_value = glsc3(result%x, fluid%c_Xh%B, reference%x, n)
    if (pe_rank == 0) write(unit, '(A,A,A,6(",",ES24.16))') trim(probe), ',', &
         trim(component), norm_value, mult_norm, &
         mult_norm / max(velocity_norm, tiny(1.0_rp)), signed_sum, pair_value, &
         relative_error
  end subroutine emit

  real(kind=rp) function vector_difference(a, b)
    type(field_t), intent(in) :: a, b
    difference_res%x = a%x - b%x
    vector_difference = sqrt(max(0.0_rp, glsc3(difference_res%x, &
         fluid%c_Xh%B, difference_res%x, n))) / &
         max(sqrt(max(0.0_rp, glsc3(b%x, fluid%c_Xh%B, b%x, n))), &
         tiny(1.0_rp))
  end function vector_difference

end program lift_weak_decomposition_diagnostic
