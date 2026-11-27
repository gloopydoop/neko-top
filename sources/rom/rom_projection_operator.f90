!> Matrix-free constrained pressure Schur complement for ROM admissibility.
!!
!! The production Pn/Pn pressure matrix is the collocated operator
!! C(grad_h q).  A ROM velocity column must additionally be continuous and
!! satisfy the homogeneous essential velocity restriction.  This operator
!! applies the collocated strong SEM gradient, assembles it into the
!! continuous velocity representation, restricts it with the native velocity
!! projector, and evaluates the native assembled weak continuity residual:
!!
!!   S q = C P_v E M^{-1}_{\mathrm{avg}} E^T grad_h(q).
!!
!! The implementation is deliberately local to the ROM projection path.  It
!! does not replace Neko's production pressure operator or alter generic KSP
!! execution.
module rom_projection_operator
  use ax_product, only: ax_t
  use fluid_pnpn, only: fluid_pnpn_t
  use vector_bc_projector, only: vector_bc_projector_t
  use field, only: field_t
  use field_math, only: field_col3, field_rzero
  use operators, only: grad, ortho
  use gather_scatter, only: GS_OP_ADD
  use math, only: col2
  use device_math, only: device_col2, device_copy
  use device, only: device_get_ptr
  use neko_config, only: NEKO_BCKND_DEVICE
  use coefs, only: coef_t
  use mesh, only: mesh_t
  use space, only: space_t
  use num_types, only: rp
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
  implicit none
  private

  type, public, extends(ax_t) :: rom_projection_ax_t
     type(fluid_pnpn_t), pointer :: fluid => null()
     type(field_t), pointer :: q => null()
     type(field_t), pointer :: p_zero => null()
     type(field_t), pointer :: p_res => null()
     type(field_t), pointer :: grad_x => null()
     type(field_t), pointer :: grad_y => null()
     type(field_t), pointer :: grad_z => null()
     type(field_t), pointer :: f_x => null()
     type(field_t), pointer :: f_y => null()
     type(field_t), pointer :: f_z => null()
     type(field_t), pointer :: u_zero => null()
     type(field_t), pointer :: v_zero => null()
     type(field_t), pointer :: w_zero => null()
     class(vector_bc_projector_t), pointer :: velocity_projector => null()
     ! The legacy basis handoff uses Neko's multiplicity averaging.  The
     ! adaptive candidate path can request the exact mass-weighted nodal map
     ! paired with its correction, B^{-1} GS(B grad q), without changing the
     ! legacy fixed-basis operator.
     logical :: mass_weighted_velocity_average = .false.
     type(c_ptr) :: event
   contains
     procedure, pass(this) :: init => rom_projection_ax_init
     procedure, pass(this) :: compute => rom_projection_ax_compute
     procedure, pass(this) :: compute_vector => rom_projection_ax_compute_vector
     procedure, pass(this) :: free => rom_projection_ax_free
  end type rom_projection_ax_t

contains

  !> Bind the ROM Schur-complement operator to Neko fluid work fields.
  subroutine rom_projection_ax_init(this, fluid, q, p_zero, p_res, grad_x, &
       grad_y, grad_z, f_x, f_y, f_z, u_zero, v_zero, w_zero, velocity_projector, &
       mass_weighted_velocity_average)
    class(rom_projection_ax_t), intent(inout) :: this
    type(fluid_pnpn_t), target, intent(inout) :: fluid
    type(field_t), target, intent(inout) :: q, p_zero, p_res
    type(field_t), target, intent(inout) :: grad_x, grad_y, grad_z
    type(field_t), target, intent(inout) :: f_x, f_y, f_z
    type(field_t), target, intent(inout) :: u_zero, v_zero, w_zero
    class(vector_bc_projector_t), target, intent(inout), optional :: velocity_projector
    logical, intent(in), optional :: mass_weighted_velocity_average

    this%fluid => fluid
    this%q => q
    this%p_zero => p_zero
    this%p_res => p_res
    this%grad_x => grad_x
    this%grad_y => grad_y
    this%grad_z => grad_z
    this%f_x => f_x
    this%f_y => f_y
    this%f_z => f_z
    this%u_zero => u_zero
    this%v_zero => v_zero
    this%w_zero => w_zero
    nullify(this%velocity_projector)
    if (present(velocity_projector)) this%velocity_projector => velocity_projector
    this%mass_weighted_velocity_average = .false.
    if (present(mass_weighted_velocity_average)) then
       this%mass_weighted_velocity_average = mass_weighted_velocity_average
    end if
    this%event = c_null_ptr
  end subroutine rom_projection_ax_init

  !> Apply the constrained weak-continuity operator to one pressure vector.
  subroutine rom_projection_ax_compute(this, w, u, coef, msh, Xh)
    class(rom_projection_ax_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    real(kind=rp), intent(inout) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    type(coef_t), intent(in) :: coef
    integer :: n
    type(c_ptr) :: event

    ! The Ax interface requires INTENT(IN) for this.  Work-field targets
    ! are mutable, but the pointer associations of this must not change.
    ! Use local aliases to avoid CCE treating target updates as changes to this.
    type(fluid_pnpn_t), pointer :: fluid
    type(field_t), pointer :: q, p_zero, p_res
    type(field_t), pointer :: grad_x, grad_y, grad_z, f_x, f_y, f_z
    type(field_t), pointer :: u_zero, v_zero, w_zero
    class(vector_bc_projector_t), pointer :: velocity_projector

    fluid => this%fluid
    q => this%q
    p_zero => this%p_zero
    p_res => this%p_res
    grad_x => this%grad_x
    grad_y => this%grad_y
    grad_z => this%grad_z
    f_x => this%f_x
    f_y => this%f_y
    f_z => this%f_z
    u_zero => this%u_zero
    v_zero => this%v_zero
    w_zero => this%w_zero
    velocity_projector => this%velocity_projector

    n = coef%dof%size()
    ! The Neko KSP owns the pressure-vector representation.  Its Krylov
    ! search vectors are already in the pressure coefficient space used by
    ! Ax; gathering q here would apply an additional averaging projector that
    ! is absent from the KSP inner products and would change the operator
    ! seen by GMRES.  Gather-scatter is applied below to C^T q, where it is
    ! the required assembly of the velocity mass-adjoint gradient.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(q%x_d, device_get_ptr(u), n)
    else
       q%x = u
    end if
    ! KSP constrains the output residual, but a general right
    ! preconditioner is not required to preserve the pressure Dirichlet
    ! subspace in its input vector.  Apply the same essential pressure
    ! restriction before differentiating q; otherwise S is evaluated on a
    ! different domain from the native pressure solve.
    call fluid%bcs_prs_projector%apply(q%x, n)

    ! Neko's pressure Schur complement pairs C with the collocated strong SEM
    ! gradient grad_h.  The ROM restriction is explicit: gather/average the
    ! element-local correction into the continuous velocity space and then
    ! apply the native essential velocity projector.  The adaptive candidate
    ! path selects the same mass-weighted map as the velocity correction,
    ! B^{-1} GS(B grad_h q).  The default retains the legacy multiplicity
    ! average so fixed-basis projection is not changed by this opt-in path.
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call grad(grad_x%x_d, grad_y%x_d, grad_z%x_d, &
            q%x_d, fluid%c_Xh)
    else
       call grad(grad_x%x, grad_y%x, grad_z%x, q%x, &
            fluid%c_Xh)
    end if
    if (this%mass_weighted_velocity_average) then
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(grad_x%x_d, fluid%c_Xh%B_d, n)
          call device_col2(grad_y%x_d, fluid%c_Xh%B_d, n)
          call device_col2(grad_z%x_d, fluid%c_Xh%B_d, n)
       else
          call col2(grad_x%x, fluid%c_Xh%B, n)
          call col2(grad_y%x, fluid%c_Xh%B, n)
          call col2(grad_z%x, fluid%c_Xh%B, n)
       end if
       call fluid%gs_Xh%op(grad_x%x, grad_y%x, grad_z%x, &
            n, GS_OP_ADD)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(grad_x%x_d, fluid%c_Xh%Binv_d, n)
          call device_col2(grad_y%x_d, fluid%c_Xh%Binv_d, n)
          call device_col2(grad_z%x_d, fluid%c_Xh%Binv_d, n)
       else
          call col2(grad_x%x, fluid%c_Xh%Binv, n)
          call col2(grad_y%x, fluid%c_Xh%Binv, n)
          call col2(grad_z%x, fluid%c_Xh%Binv, n)
       end if
    else
       call fluid%gs_Xh%op(grad_x%x, grad_y%x, grad_z%x, &
            n, GS_OP_ADD)
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_col2(grad_x%x_d, fluid%c_Xh%mult_d, n)
          call device_col2(grad_y%x_d, fluid%c_Xh%mult_d, n)
          call device_col2(grad_z%x_d, fluid%c_Xh%mult_d, n)
       else
          call col2(grad_x%x, fluid%c_Xh%mult, n)
          call col2(grad_y%x, fluid%c_Xh%mult, n)
          call col2(grad_z%x, fluid%c_Xh%mult, n)
       end if
    end if
    if (associated(velocity_projector)) then
       call velocity_projector%apply(grad_x%x, grad_y%x, &
            grad_z%x, n)
    else
       call fluid%bcs_vel_projector%apply(grad_x%x, grad_y%x, &
            grad_z%x, n)
    end if

    ! Apply C to the restricted mass-adjoint gradient.  The result is left
    ! in local elemental form because the surrounding Neko KSP applies its
    ! standard gather-scatter and pressure essential-BC restriction.
    call field_col3(f_x, grad_x, fluid%rho)
    call field_col3(f_y, grad_y, fluid%rho)
    call field_col3(f_z, grad_z, fluid%rho)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(f_x%x_d, fluid%c_Xh%B_d, n)
       call device_col2(f_y%x_d, fluid%c_Xh%B_d, n)
       call device_col2(f_z%x_d, fluid%c_Xh%B_d, n)
    else
       call col2(f_x%x, fluid%c_Xh%B, n)
       call col2(f_y%x, fluid%c_Xh%B, n)
       call col2(f_z%x, fluid%c_Xh%B, n)
    end if
    call field_rzero(p_zero)
    event = c_null_ptr
    call fluid%prs_res%compute(p_zero, p_res, grad_x, &
         grad_y, grad_z, u_zero, v_zero, w_zero, &
         f_x, f_y, f_z, fluid%c_Xh, fluid%gs_Xh, &
         fluid%bc_prs_surface, fluid%bc_sym_surface, &
         fluid%Ax_prs, 0.0_rp, 1.0_rp, u_zero, fluid%rho, &
         event)
    if (.not. fluid%prs_dirichlet) then
       call ortho(p_res%x, fluid%glb_n_points, n)
    end if
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_copy(device_get_ptr(w), p_res%x_d, n)
    else
       w = p_res%x
    end if
  end subroutine rom_projection_ax_compute

  !> Apply the constrained operator to three vector-valued right-hand sides.
  subroutine rom_projection_ax_compute_vector(this, au, av, aw, u, v, w, &
       coef, msh, Xh)
    class(rom_projection_ax_t), intent(in) :: this
    type(mesh_t), intent(in) :: msh
    type(space_t), intent(in) :: Xh
    real(kind=rp), intent(inout) :: au(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: av(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(inout) :: aw(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: u(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: v(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    real(kind=rp), intent(in) :: w(Xh%lx, Xh%ly, Xh%lz, msh%nelv)
    type(coef_t), intent(in) :: coef
    call this%compute(au, u, coef, msh, Xh)
    call this%compute(av, v, coef, msh, Xh)
    call this%compute(aw, w, coef, msh, Xh)
  end subroutine rom_projection_ax_compute_vector

  !> Release the non-owning field associations of the projection operator.
  subroutine rom_projection_ax_free(this)
    class(rom_projection_ax_t), intent(inout) :: this
    nullify(this%fluid, this%q, this%p_zero, this%p_res)
    nullify(this%grad_x, this%grad_y, this%grad_z)
    nullify(this%f_x, this%f_y, this%f_z)
    nullify(this%u_zero, this%v_zero, this%w_zero)
    nullify(this%velocity_projector)
    this%event = c_null_ptr
  end subroutine rom_projection_ax_free

end module rom_projection_operator
