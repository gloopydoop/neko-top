!> Test-only adapters for legacy manufactured diagnostics.
!!
!! The production ROM assembly/objective APIs use Neko vector_t/matrix_t.  A
!! few small historical diagnostics intentionally keep their manufactured
!! arrays for readable reference formulas; this module is the only bridge.
module rom_array_test_adapters
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t
  use rom_assembly, only: rom_assemble_mass, rom_assemble_brinkman, &
       rom_assemble_lift_convection, rom_assemble_diffusion, rom_assemble_convection, &
       rom_assemble_scalar_transport
  use rom_objectives, only: rom_assemble_viscous_quadratic, &
       rom_assemble_brinkman_objective, rom_viscous_value, rom_brinkman_value, &
       rom_brinkman_density, rom_assemble_scalar_mixing_terms
  implicit none
  private
  public :: test_assemble_mass, test_assemble_brinkman
  public :: test_assemble_lift_convection, test_assemble_diffusion
  public :: test_assemble_convection, test_assemble_viscous_quadratic
  public :: test_assemble_brinkman_objective, test_viscous_value
  public :: test_brinkman_value, test_brinkman_density
  public :: test_assemble_scalar_mixing, test_scalar_mixing_value
  public :: test_assemble_scalar_transport

contains

  !> Expose reduced mass assembly to the array-based test harness.
  subroutine test_assemble_mass(weights, modes, mass)
    real(kind=rp), intent(in) :: weights(:), modes(:, :, :)
    real(kind=rp), intent(out) :: mass(:, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: m(:,:)
    type(matrix_t) :: result
    integer :: n, nf, r, i, j
    n=size(weights); nf=size(modes,2); r=size(modes,3)
    call w%init(n); w%x=weights; allocate(m(nf,r))
    do j=1,r; do i=1,nf; call m(i,j)%init(n); m(i,j)%x=modes(:,i,j); end do; end do
    call rom_assemble_mass(w,m,result); mass=result%x
    do j=1,r; do i=1,nf; call m(i,j)%free(); end do; end do
    deallocate(m); call w%free(); call result%free()
  end subroutine test_assemble_mass

  !> Expose reduced Brinkman assembly to the array-based test harness.
  subroutine test_assemble_brinkman(weights, chi, modes, lift, brinkman, lift_vector)
    real(kind=rp), intent(in) :: weights(:), chi(:), modes(:, :, :), lift(:, :)
    real(kind=rp), intent(out) :: brinkman(:, :), lift_vector(:)
    type(vector_t) :: w, c
    type(vector_t), allocatable :: m(:,:), l(:)
    type(matrix_t) :: result
    type(vector_t) :: result_lift
    integer :: n, r, i, j
    n=size(weights); r=size(modes,3)
    call w%init(n); call c%init(n); w%x=weights; c%x=chi
    allocate(m(3,r),l(3))
    do j=1,r; do i=1,3; call m(i,j)%init(n); m(i,j)%x=modes(:,i,j); end do; end do
    do i=1,3; call l(i)%init(n); l(i)%x=lift(:,i); end do
    call rom_assemble_brinkman(w,c,m,l,result,result_lift)
    brinkman=result%x; lift_vector=result_lift%x
    do j=1,r; do i=1,3; call m(i,j)%free(); end do; end do
    do i=1,3; call l(i)%free(); end do
    deallocate(m,l); call w%free(); call c%free(); call result%free(); call result_lift%free()
  end subroutine test_assemble_brinkman

  !> Expose lift-convection assembly to the array-based test harness.
  subroutine test_assemble_lift_convection(weights, modes, gradients, lift, linear)
    real(kind=rp), intent(in) :: weights(:), modes(:, :, :), gradients(:, :, :, :), lift(:, :)
    real(kind=rp), intent(out) :: linear(:, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: m(:,:), g(:,:,:), l(:), lg(:,:)
    type(matrix_t) :: result
    type(vector_t) :: result_constant
    integer :: n,r,i,j,k,d
    n=size(weights); r=size(modes,3); call w%init(n); w%x=weights
    allocate(m(3,r),g(3,r,3),l(3),lg(3,3))
    do j=1,r; do i=1,3
       call m(i,j)%init(n); m(i,j)%x=modes(:,i,j)
       do d=1,3; call g(i,j,d)%init(n); g(i,j,d)%x=gradients(:,i,j,d); end do
    end do; end do
    do i=1,3; call l(i)%init(n); l(i)%x=lift(:,i); end do
    ! The legacy residual diagnostic uses a uniform lift, so grad(lift)=0.
    do i=1,3; do d=1,3; call lg(i,d)%init(n); lg(i,d)=0.0_rp; end do; end do
    call rom_assemble_lift_convection(w,m,g,l,lg,result,result_constant)
    linear=result%x
    do j=1,r; do i=1,3; do d=1,3; call g(i,j,d)%free(); end do; call m(i,j)%free(); end do; end do
    do i=1,3; call l(i)%free(); end do
    do i=1,3; do d=1,3; call lg(i,d)%free(); end do; end do
    deallocate(m,g,l,lg); call w%free(); call result%free(); call result_constant%free()
  end subroutine test_assemble_lift_convection

  !> Expose diffusion assembly to the array-based test harness.
  subroutine test_assemble_diffusion(weights, nu, gradients, diffusion)
    real(kind=rp), intent(in) :: weights(:), nu, gradients(:, :, :, :)
    real(kind=rp), intent(out) :: diffusion(:, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: g(:,:,:)
    type(matrix_t) :: result
    integer :: n,r,i,j,d
    n=size(weights); r=size(gradients,3); call w%init(n); w%x=weights
    allocate(g(3,r,3))
    do j=1,r; do i=1,3; do d=1,3; call g(i,j,d)%init(n); g(i,j,d)%x=gradients(:,i,j,d); end do; end do; end do
    call rom_assemble_diffusion(w,nu,g,result); diffusion=result%x
    do j=1,r; do i=1,3; do d=1,3; call g(i,j,d)%free(); end do; end do; end do
    deallocate(g); call w%free(); call result%free()
  end subroutine test_assemble_diffusion

  !> Expose convection assembly to the array-based test harness.
  subroutine test_assemble_convection(weights, modes, gradients, quadratic)
    real(kind=rp), intent(in) :: weights(:), modes(:, :, :), gradients(:, :, :, :)
    real(kind=rp), intent(out) :: quadratic(:, :, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: m(:,:), g(:,:,:)
    type(matrix_t), allocatable :: result(:)
    integer :: n,r,i,j,d
    n=size(weights); r=size(modes,3); call w%init(n); w%x=weights
    allocate(m(3,r),g(3,r,3),result(r))
    do j=1,r; do i=1,3
       call m(i,j)%init(n); m(i,j)%x=modes(:,i,j)
       do d=1,3; call g(i,j,d)%init(n); g(i,j,d)%x=gradients(:,i,j,d); end do
    end do; end do
    call rom_assemble_convection(w,m,g,result)
    do i=1,r; quadratic(i,:,:)=result(i)%x; end do
    do j=1,r; do i=1,3; do d=1,3; call g(i,j,d)%free(); end do; call m(i,j)%free(); end do; end do
    do i=1,r; call result(i)%free(); end do
    deallocate(m,g,result); call w%free()
  end subroutine test_assemble_convection

  !> Expose viscous objective assembly to the array-based test harness.
  subroutine test_assemble_viscous_quadratic(weights, mask, viscosity, gradients, quadratic)
    real(kind=rp), intent(in) :: weights(:), mask(:), viscosity, gradients(:, :, :, :)
    real(kind=rp), intent(out) :: quadratic(:, :)
    type(vector_t) :: w,m
    type(vector_t), allocatable :: g(:,:,:)
    type(matrix_t) :: result
    integer :: n,r,i,j,d
    n=size(weights); r=size(gradients,3); call w%init(n); call m%init(n); w%x=weights; m%x=mask
    allocate(g(3,r,3))
    do j=1,r; do i=1,3; do d=1,3; call g(i,j,d)%init(n); g(i,j,d)%x=gradients(:,i,j,d); end do; end do; end do
    call rom_assemble_viscous_quadratic(w,m,viscosity,g,result); quadratic=result%x
    do j=1,r; do i=1,3; do d=1,3; call g(i,j,d)%free(); end do; end do; end do
    deallocate(g); call w%free(); call m%free(); call result%free()
  end subroutine test_assemble_viscous_quadratic

  !> Expose Brinkman objective assembly to the array-based test harness.
  subroutine test_assemble_brinkman_objective(weights, mask, chi, modes, lift, quadratic, lift_linear, constant)
    real(kind=rp), intent(in) :: weights(:), mask(:), chi(:), modes(:, :, :), lift(:, :)
    real(kind=rp), intent(out) :: quadratic(:, :), lift_linear(:), constant
    type(vector_t) :: w,m,c
    type(vector_t), allocatable :: v(:,:),l(:)
    type(matrix_t) :: result
    type(vector_t) :: result_lift
    integer :: n,r,i,j
    n=size(weights); r=size(modes,3); call w%init(n); call m%init(n); call c%init(n)
    w%x=weights; m%x=mask; c%x=chi; allocate(v(3,r),l(3))
    do j=1,r; do i=1,3; call v(i,j)%init(n); v(i,j)%x=modes(:,i,j); end do; end do
    do i=1,3; call l(i)%init(n); l(i)%x=lift(:,i); end do
    call rom_assemble_brinkman_objective(w,m,c,v,l,result,result_lift,constant)
    quadratic=result%x; lift_linear=result_lift%x
    do j=1,r; do i=1,3; call v(i,j)%free(); end do; end do
    do i=1,3; call l(i)%free(); end do
    deallocate(v,l); call w%free(); call m%free(); call c%free(); call result%free(); call result_lift%free()
  end subroutine test_assemble_brinkman_objective

  !> Expose viscous objective evaluation to the array-based test harness.
  function test_viscous_value(state, quadratic, domain) result(value)
    real(kind=rp), intent(in) :: state(:), quadratic(:, :), domain
    real(kind=rp) :: value
    real(kind=rp) :: constant
    type(vector_t) :: s, l; type(matrix_t) :: q
    call s%init(size(state)); call l%init(size(state)); call q%init(size(state),size(state))
    s%x=state; l=0.0_rp; q%x=quadratic; constant=0.0_rp
    value=rom_viscous_value(s,constant,l,q,domain)
    call s%free(); call l%free(); call q%free()
  end function test_viscous_value

  !> Expose Brinkman objective evaluation to the array-based test harness.
  function test_brinkman_value(state, constant, lift_linear, quadratic, domain) result(value)
    real(kind=rp), intent(in) :: state(:), constant, lift_linear(:), quadratic(:, :), domain
    real(kind=rp) :: value
    type(vector_t) :: s,l; type(matrix_t) :: q
    call s%init(size(state)); call l%init(size(state)); call q%init(size(state),size(state))
    s%x=state; l%x=lift_linear; q%x=quadratic
    value=rom_brinkman_value(s,constant,l,q,domain); call s%free(); call l%free(); call q%free()
  end function test_brinkman_value

  !> Expose Brinkman density assembly to the array-based test harness.
  subroutine test_brinkman_density(weights, mask, lift, modes, state, domain, density)
    real(kind=rp), intent(in) :: weights(:), mask(:), lift(:, :), modes(:, :, :), state(:), domain
    real(kind=rp), intent(out) :: density(:)
    type(vector_t) :: w,m,s,d
    type(vector_t), allocatable :: v(:,:),l(:)
    integer :: n,r,i,j
    n=size(weights); r=size(state); call w%init(n); call m%init(n); call s%init(r)
    w%x=weights; m%x=mask; s%x=state; allocate(v(3,r),l(3))
    do j=1,r; do i=1,3; call v(i,j)%init(n); v(i,j)%x=modes(:,i,j); end do; end do
    do i=1,3; call l(i)%init(n); l(i)%x=lift(:,i); end do
    call rom_brinkman_density(w,m,l,v,s,domain,d); density=d%x
    do j=1,r; do i=1,3; call v(i,j)%free(); end do; end do
    do i=1,3; call l(i)%free(); end do
    deallocate(v,l); call w%free(); call m%free(); call s%free(); call d%free()
  end subroutine test_brinkman_density

  !> Expose the physical scalar-mixing quadratic form to the test harness.
  subroutine test_assemble_scalar_mixing(weights, mask, lift, modes, reference, &
       constant, linear, quadratic)
    real(kind=rp), intent(in) :: weights(:), mask(:), lift(:), modes(:, :), reference
    real(kind=rp), intent(out) :: constant, linear(:), quadratic(:, :)
    type(vector_t) :: w, m, l
    type(vector_t), allocatable :: s(:)
    type(vector_t) :: result_constant, result_linear
    type(matrix_t) :: result_quadratic
    integer :: n, r, i

    n = size(weights); r = size(modes, 2)
    call w%init(n); call m%init(n); call l%init(n)
    w%x = weights; m%x = mask; l%x = lift
    allocate(s(r))
    do i = 1, r
       call s(i)%init(n); s(i)%x = modes(:, i)
    end do
    call rom_assemble_scalar_mixing_terms(w, m, l, s, reference, &
         result_constant, result_linear, result_quadratic)
    constant = result_constant%x(1)
    linear = result_linear%x
    quadratic = result_quadratic%x
    do i = 1, r; call s(i)%free(); end do
    deallocate(s)
    call w%free(); call m%free(); call l%free()
    call result_constant%free(); call result_linear%free(); call result_quadratic%free()
  end subroutine test_assemble_scalar_mixing

  !> Expose scalar-mixing objective evaluation to the native-objective test.
  function test_scalar_mixing_value(state, constant, linear, quadratic, domain) &
       result(value)
    real(kind=rp), intent(in) :: state(:), constant, linear(:), quadratic(:, :), domain
    real(kind=rp) :: value
    type(vector_t) :: s, l
    type(matrix_t) :: q

    call s%init(size(state)); call l%init(size(state))
    call q%init(size(state), size(state))
    s%x = state; l%x = linear; q%x = quadratic
    value = rom_viscous_value(s, constant, l, q, domain)
    call s%free(); call l%free(); call q%free()
  end function test_scalar_mixing_value

  !> Expose coupled scalar transport assembly to the manufactured test.
  !! Every mesh-sized argument is converted to a Neko vector_t before calling
  !! the production assembler; the adapter contains no alternate ROM formula.
  subroutine test_assemble_scalar_transport(weights, scalar_weight, velocity_modes, &
       scalar_modes, velocity_gradients, scalar_gradients, lift_velocity, scalar_lift, &
       lift_scalar_gradient, diffusivity, linear, constant, quadratic, implicit_linear)
    real(kind=rp), intent(in) :: weights(:), scalar_weight
    real(kind=rp), intent(in) :: velocity_modes(:, :, :), scalar_modes(:, :)
    real(kind=rp), intent(in) :: velocity_gradients(:, :, :, :)
    real(kind=rp), intent(in) :: scalar_gradients(:, :, :)
    real(kind=rp), intent(in) :: lift_velocity(:, :), scalar_lift(:)
    real(kind=rp), intent(in) :: lift_scalar_gradient(:, :), diffusivity(:)
    real(kind=rp), intent(out) :: linear(:, :), constant(:), quadratic(:, :, :), &
         implicit_linear(:, :)
    type(vector_t) :: w, sl, diff
    type(vector_t), allocatable :: vm(:, :), sm(:), vg(:, :, :), sg(:, :), lv(:), lsg(:)
    type(matrix_t) :: result_linear, result_implicit
    type(vector_t) :: result_constant
    type(matrix_t), allocatable :: result_quadratic(:)
    integer :: n, r, i, j, c, d

    n = size(weights); r = size(scalar_modes, 2)
    call w%init(n); w%x = weights
    call sl%init(n); sl%x = scalar_lift
    call diff%init(n); diff%x = diffusivity
    allocate(vm(3,r), sm(r), vg(3,r,3), sg(r,3), lv(3), lsg(3))
    do j = 1, r
       call sm(j)%init(n); sm(j)%x = scalar_modes(:,j)
       do c = 1, 3
          call vm(c,j)%init(n); vm(c,j)%x = velocity_modes(:,c,j)
          do d = 1, 3
             call vg(c,j,d)%init(n); vg(c,j,d)%x = velocity_gradients(:,c,j,d)
          end do
       end do
       do d = 1, 3
          call sg(j,d)%init(n); sg(j,d)%x = scalar_gradients(:,j,d)
       end do
    end do
    do c = 1, 3
       call lv(c)%init(n); lv(c)%x = lift_velocity(:,c)
       call lsg(c)%init(n); lsg(c)%x = lift_scalar_gradient(:,c)
    end do
    allocate(result_quadratic(r))
    call result_linear%init(r,r); result_linear = 0.0_rp
    call result_implicit%init(r,r); result_implicit = 0.0_rp
    call result_constant%init(r); result_constant = 0.0_rp
    do i = 1, r
       call result_quadratic(i)%init(r,r)
       result_quadratic(i) = 0.0_rp
    end do
    call rom_assemble_scalar_transport(w, scalar_weight, vm, sm, vg, sg, lv, sl, lsg, &
         diff, result_linear, result_constant, result_quadratic, result_implicit)
    linear = result_linear%x
    constant = result_constant%x
    implicit_linear = result_implicit%x
    do i = 1, r
       quadratic(i,:,:) = result_quadratic(i)%x
    end do
    do j = 1, r
       call sm(j)%free()
       do c = 1, 3
          call vm(c,j)%free()
          do d = 1, 3; call vg(c,j,d)%free(); end do
       end do
       do d = 1, 3; call sg(j,d)%free(); end do
    end do
    do c = 1, 3
       call lv(c)%free(); call lsg(c)%free()
    end do
    do i = 1, r; call result_quadratic(i)%free(); end do
    deallocate(vm, sm, vg, sg, lv, lsg, result_quadratic)
    call w%free(); call sl%free(); call diff%free()
    call result_linear%free(); call result_constant%free(); call result_implicit%free()
  end subroutine test_assemble_scalar_transport

end module rom_array_test_adapters
