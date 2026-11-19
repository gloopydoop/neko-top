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
       rom_assemble_lift_convection, rom_assemble_diffusion, rom_assemble_convection
  use rom_objectives, only: rom_assemble_viscous_quadratic, &
       rom_assemble_brinkman_objective, rom_viscous_value, rom_brinkman_value, &
       rom_brinkman_density
  implicit none
  private
  public :: test_assemble_mass, test_assemble_brinkman
  public :: test_assemble_lift_convection, test_assemble_diffusion
  public :: test_assemble_convection, test_assemble_viscous_quadratic
  public :: test_assemble_brinkman_objective, test_viscous_value
  public :: test_brinkman_value, test_brinkman_density

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

end module rom_array_test_adapters
