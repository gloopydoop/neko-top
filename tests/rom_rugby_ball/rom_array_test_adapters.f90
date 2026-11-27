!> Test-only adapters for legacy manufactured diagnostics.
!!
!! The production ROM assembly/objective APIs use Neko vector_t/matrix_t.  A
!! few small historical diagnostics intentionally keep their manufactured
!! arrays for readable reference formulas; this module is the only bridge.
module rom_array_test_adapters
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t
  use vector_math, only: vector_glsc3, vector_col3, vector_add2, vector_cmult
  use rom_assembly, only: rom_assemble_mass, rom_assemble_cross_mass, &
       rom_assemble_brinkman, rom_assemble_brinkman_petrov, &
       rom_assemble_lift_convection, rom_assemble_diffusion, rom_assemble_convection, &
       rom_assemble_scalar_transport, rom_reduce_buffer
  use rom_objectives, only: rom_assemble_viscous_quadratic, &
       rom_assemble_brinkman_objective, rom_viscous_value, rom_brinkman_value, &
       rom_scalar_mixing_value, rom_brinkman_density, rom_assemble_scalar_mixing_terms
  implicit none
  private
  public :: test_assemble_mass, test_assemble_brinkman
  public :: test_assemble_mass_legacy
  public :: test_assemble_cross_mass_petrov, test_assemble_cross_mass_petrov_legacy
  public :: test_assemble_brinkman_legacy
  public :: test_assemble_brinkman_petrov, test_assemble_brinkman_petrov_legacy
  public :: test_assemble_lift_convection, test_assemble_diffusion
  public :: test_assemble_convection, test_assemble_viscous_quadratic
  public :: test_assemble_diffusion_legacy
  public :: test_assemble_convection_legacy
  public :: test_assemble_scalar_transport_legacy
  public :: test_assemble_brinkman_objective, test_viscous_value
  public :: test_brinkman_value, test_brinkman_density
  public :: test_assemble_scalar_mixing, test_scalar_mixing_value
  public :: test_assemble_scalar_transport
  public :: test_assemble_convection_dgemm

  interface
     subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
       import :: rp
       character(len=1), intent(in) :: transa, transb
       integer, intent(in) :: m, n, k, lda, ldb, ldc
       real(kind=rp), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
       real(kind=rp), intent(inout) :: c(ldc,*)
     end subroutine dgemm
  end interface

contains

  !> Test-only blocked DGEMM convection path. Each spatial pair field is
  !! formed once, then contracted against a bounded panel of weighted tests.
  !! The local tensor is reduced in the production tensor's packed order.
  subroutine test_assemble_convection_dgemm(weights, modes, gradients, block_size, quadratic)
    real(kind=rp), intent(in) :: weights(:), modes(:, :, :), gradients(:, :, :, :)
    integer, intent(in) :: block_size
    real(kind=rp), intent(out) :: quadratic(:, :, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: m(:, :), g(:, :, :)
    type(vector_t) :: pointwise, product
    real(kind=rp), allocatable :: weighted_tests(:, :), pair_panel(:, :)
    real(kind=rp), allocatable :: gemm_panel(:, :), local_tensor(:, :, :), packed(:)
    integer :: n, r, advect_mode, state_mode, test_mode, field, component
    integer :: first_pair, last_pair, panel_count, pair, pair_advect, pair_state
    integer :: packed_index

    n = size(weights); r = size(modes,3)
    if (block_size < 1 .or. size(modes,1) /= n .or. size(modes,2) /= 3 .or. &
         size(gradients,1) /= n .or. size(gradients,2) /= 3 .or. &
         size(gradients,3) /= r .or. size(gradients,4) /= 3 .or. &
         any(shape(quadratic) /= [r,r,r])) error stop 'DGEMM convection shape mismatch'
    call w%init(n); w = weights
    allocate(m(3,r), g(3,r,3), weighted_tests(n,r))
    do state_mode = 1, r
       do component = 1, 3
          call m(component,state_mode)%init(n)
          m(component,state_mode) = modes(:,component,state_mode)
          do field = 1, 3
             call g(component,state_mode,field)%init(n)
             g(component,state_mode,field) = gradients(:,component,state_mode,field)
          end do
       end do
    end do
    call pointwise%init(n); call product%init(n)
    allocate(pair_panel(n,block_size), gemm_panel(r,block_size))
    allocate(local_tensor(r,r,r)); local_tensor = 0.0_rp

    do field = 1, 3
       do test_mode = 1, r
          weighted_tests(:,test_mode) = weights * modes(:,field,test_mode)
       end do
       do first_pair = 1, r*r, block_size
          last_pair = min(first_pair + block_size - 1, r*r)
          panel_count = last_pair - first_pair + 1
          do pair = first_pair, last_pair
             pair_advect = mod(pair-1,r) + 1
             pair_state = (pair-1)/r + 1
             pointwise = 0.0_rp
             do component = 1, 3
                call vector_col3(product, m(component,pair_advect), &
                     g(field,pair_state,component))
                call vector_add2(pointwise, product)
             end do
             pair_panel(:,pair-first_pair+1) = pointwise%x
          end do
          gemm_panel(:,1:panel_count) = 0.0_rp
          call dgemm('T','N',r,panel_count,n,1.0_rp,weighted_tests,n, &
               pair_panel,n,0.0_rp,gemm_panel,r)
          do pair = first_pair, last_pair
             pair_advect = mod(pair-1,r) + 1
             pair_state = (pair-1)/r + 1
             local_tensor(:,pair_advect,pair_state) = &
                  local_tensor(:,pair_advect,pair_state) + &
                  gemm_panel(:,pair-first_pair+1)
          end do
       end do
    end do
    allocate(packed(size(local_tensor)))
    ! Match rom_reduce_matrix_tensor's order: test matrix outermost, then
    ! each matrix's column-major (advect, state) entries.
    packed_index = 0
    do test_mode = 1, r
       do state_mode = 1, r
          do advect_mode = 1, r
             packed_index = packed_index + 1
             packed(packed_index) = local_tensor(test_mode,advect_mode,state_mode)
          end do
       end do
    end do
    call rom_reduce_buffer(packed)
    packed_index = 0
    do test_mode = 1, r
       do state_mode = 1, r
          do advect_mode = 1, r
             packed_index = packed_index + 1
             local_tensor(test_mode,advect_mode,state_mode) = packed(packed_index)
          end do
       end do
    end do
    quadratic = local_tensor
    do state_mode = 1, r
       do component = 1, 3
          do field = 1, 3; call g(component,state_mode,field)%free(); end do
          call m(component,state_mode)%free()
       end do
    end do
    call w%free(); call pointwise%free(); call product%free()
    deallocate(m,g,weighted_tests,pair_panel,gemm_panel,local_tensor,packed)
  end subroutine test_assemble_convection_dgemm

  !> Legacy per-contraction reference for the current packed mass assembly.
  !! Each vector_glsc3 call performs the historical global reduction.
  subroutine test_assemble_mass_legacy(weights, modes, mass)
    real(kind=rp), intent(in) :: weights(:), modes(:, :, :)
    real(kind=rp), intent(out) :: mass(:, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: m(:,:)
    integer :: n, nf, r, i, j, field

    n = size(weights); nf = size(modes,2); r = size(modes,3)
    call w%init(n); w = weights
    allocate(m(nf,r))
    do j = 1, r
       do field = 1, nf
          call m(field,j)%init(n); m(field,j) = modes(:,field,j)
       end do
    end do
    mass = 0.0_rp
    do i = 1, r
       do j = 1, r
          do field = 1, nf
             mass(i,j) = mass(i,j) + vector_glsc3(m(field,i), w, m(field,j))
          end do
       end do
    end do
    do j = 1, r
       do field = 1, nf
          call m(field,j)%free()
       end do
    end do
    deallocate(m); call w%free()
  end subroutine test_assemble_mass_legacy

  !> Expose reduced mass assembly to the array-based test harness.
  subroutine test_assemble_mass(weights, modes, mass)
    real(kind=rp), intent(in) :: weights(:), modes(:, :, :)
    real(kind=rp), intent(out) :: mass(:, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: m(:,:)
    type(matrix_t) :: result
    integer :: n, nf, r, i, j
    n=size(weights); nf=size(modes,2); r=size(modes,3)
    call w%init(n); w =weights; allocate(m(nf,r))
    do j=1,r; do i=1,nf; call m(i,j)%init(n); m(i,j) =modes(:,i,j); end do; end do
    call rom_assemble_mass(w,m,result); mass=result%x
    do j=1,r; do i=1,nf; call m(i,j)%free(); end do; end do
    deallocate(m); call w%free(); call result%free()
  end subroutine test_assemble_mass

  !> Run the production Petrov cross-mass assembly on distinct test/trial bases.
  subroutine test_assemble_cross_mass_petrov(weights, test_modes, trial_modes, mass)
    real(kind=rp), intent(in) :: weights(:), test_modes(:, :, :), trial_modes(:, :, :)
    real(kind=rp), intent(out) :: mass(:, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: test(:,:), trial(:,:)
    type(matrix_t) :: result
    integer :: n, r, field, j

    n = size(weights); r = size(trial_modes,3)
    call w%init(n); w = weights
    allocate(test(3,r), trial(3,r))
    do j = 1, r
       do field = 1, 3
          call test(field,j)%init(n); test(field,j) = test_modes(:,field,j)
          call trial(field,j)%init(n); trial(field,j) = trial_modes(:,field,j)
       end do
    end do
    call rom_assemble_cross_mass(w, trial, test, result)
    mass = result%x
    do j = 1, r
       do field = 1, 3
          call test(field,j)%free(); call trial(field,j)%free()
       end do
    end do
    deallocate(test,trial); call w%free(); call result%free()
  end subroutine test_assemble_cross_mass_petrov

  !> Legacy global contraction oracle for Petrov velocity cross mass.
  subroutine test_assemble_cross_mass_petrov_legacy(weights, test_modes, trial_modes, mass)
    real(kind=rp), intent(in) :: weights(:), test_modes(:, :, :), trial_modes(:, :, :)
    real(kind=rp), intent(out) :: mass(:, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: test(:,:), trial(:,:)
    integer :: n, r, field, i, j

    n = size(weights); r = size(trial_modes,3)
    call w%init(n); w = weights
    allocate(test(3,r), trial(3,r))
    do j = 1, r
       do field = 1, 3
          call test(field,j)%init(n); test(field,j) = test_modes(:,field,j)
          call trial(field,j)%init(n); trial(field,j) = trial_modes(:,field,j)
       end do
    end do
    mass = 0.0_rp
    do i = 1, r
       do j = 1, r
          do field = 1, 3
             mass(i,j) = mass(i,j) + &
                  vector_glsc3(test(field,i), w, trial(field,j))
          end do
       end do
    end do
    do j = 1, r
       do field = 1, 3
          call test(field,j)%free(); call trial(field,j)%free()
       end do
    end do
    deallocate(test,trial); call w%free()
  end subroutine test_assemble_cross_mass_petrov_legacy

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
    call w%init(n); call c%init(n); w =weights; c =chi
    allocate(m(3,r),l(3))
    do j=1,r; do i=1,3; call m(i,j)%init(n); m(i,j) =modes(:,i,j); end do; end do
    do i=1,3; call l(i)%init(n); l(i) =lift(:,i); end do
    call rom_assemble_brinkman(w,c,m,l,result,result_lift)
    brinkman=result%x; lift_vector=result_lift%x
    do j=1,r; do i=1,3; call m(i,j)%free(); end do; end do
    do i=1,3; call l(i)%free(); end do
    deallocate(m,l); call w%free(); call c%free(); call result%free(); call result_lift%free()
  end subroutine test_assemble_brinkman

  !> Run production Petrov Brinkman assembly with distinct test/trial bases.
  subroutine test_assemble_brinkman_petrov(weights, chi, test_modes, trial_modes, &
       lift, brinkman, lift_vector)
    real(kind=rp), intent(in) :: weights(:), chi(:), test_modes(:, :, :), trial_modes(:, :, :)
    real(kind=rp), intent(in) :: lift(:, :)
    real(kind=rp), intent(out) :: brinkman(:, :), lift_vector(:)
    type(vector_t) :: w, c
    type(vector_t), allocatable :: test(:,:), trial(:,:), l(:)
    type(matrix_t) :: result
    type(vector_t) :: result_lift
    integer :: n, r, field, j

    n = size(weights); r = size(trial_modes,3)
    call w%init(n); call c%init(n); w =weights; c =chi
    allocate(test(3,r), trial(3,r), l(3))
    do j = 1, r
       do field = 1, 3
          call test(field,j)%init(n); test(field,j) = test_modes(:,field,j)
          call trial(field,j)%init(n); trial(field,j) = trial_modes(:,field,j)
       end do
    end do
    do field = 1, 3; call l(field)%init(n); l(field) =lift(:,field); end do
    call rom_assemble_brinkman_petrov(w, c, trial, test, l, result, result_lift)
    brinkman=result%x; lift_vector=result_lift%x
    do j = 1, r
       do field = 1, 3
          call test(field,j)%free(); call trial(field,j)%free()
       end do
    end do
    do field = 1, 3; call l(field)%free(); end do
    deallocate(test,trial,l)
    call w%free(); call c%free(); call result%free(); call result_lift%free()
  end subroutine test_assemble_brinkman_petrov

  !> Legacy per-contraction oracle for Petrov Brinkman matrix and lift vector.
  subroutine test_assemble_brinkman_petrov_legacy(weights, chi, test_modes, &
       trial_modes, lift, brinkman, lift_vector)
    real(kind=rp), intent(in) :: weights(:), chi(:), test_modes(:, :, :), trial_modes(:, :, :)
    real(kind=rp), intent(in) :: lift(:, :)
    real(kind=rp), intent(out) :: brinkman(:, :), lift_vector(:)
    type(vector_t) :: w, c, effective
    type(vector_t), allocatable :: test(:,:), trial(:,:), l(:)
    integer :: n, r, field, i, j

    n = size(weights); r = size(trial_modes,3)
    call w%init(n); call c%init(n); call effective%init(n)
    w =weights; c =chi; call vector_col3(effective,w,c)
    allocate(test(3,r), trial(3,r), l(3))
    do j = 1, r
       do field = 1, 3
          call test(field,j)%init(n); test(field,j) = test_modes(:,field,j)
          call trial(field,j)%init(n); trial(field,j) = trial_modes(:,field,j)
       end do
    end do
    do field = 1, 3; call l(field)%init(n); l(field) =lift(:,field); end do
    brinkman=0.0_rp; lift_vector=0.0_rp
    do i = 1, r
       do j = 1, r
          do field = 1, 3
             brinkman(i,j) = brinkman(i,j) + &
                  vector_glsc3(test(field,i),effective,trial(field,j))
          end do
       end do
       do field = 1, 3
          lift_vector(i) = lift_vector(i) + &
               vector_glsc3(test(field,i),effective,l(field))
       end do
    end do
    do j = 1, r
       do field = 1, 3
          call test(field,j)%free(); call trial(field,j)%free()
       end do
    end do
    do field = 1, 3; call l(field)%free(); end do
    deallocate(test,trial,l)
    call w%free(); call c%free(); call effective%free()
  end subroutine test_assemble_brinkman_petrov_legacy

  !> Legacy per-contraction reference for Brinkman matrix and lift vector.
  !! This matches the committed implementation: form weights*chi once, then
  !! globally reduce each component contraction with vector_glsc3.
  subroutine test_assemble_brinkman_legacy(weights, chi, modes, lift, brinkman, lift_vector)
    real(kind=rp), intent(in) :: weights(:), chi(:), modes(:, :, :), lift(:, :)
    real(kind=rp), intent(out) :: brinkman(:, :), lift_vector(:)
    type(vector_t) :: w, c, effective
    type(vector_t), allocatable :: m(:,:), l(:)
    integer :: n, r, i, j, component

    n = size(weights); r = size(modes,3)
    call w%init(n); call c%init(n); call effective%init(n)
    w = weights; c = chi; call vector_col3(effective, w, c)
    allocate(m(3,r), l(3))
    do j = 1, r
       do component = 1, 3
          call m(component,j)%init(n)
          m(component,j) = modes(:,component,j)
       end do
    end do
    do component = 1, 3
       call l(component)%init(n); l(component) = lift(:,component)
    end do
    brinkman = 0.0_rp; lift_vector = 0.0_rp
    do i = 1, r
       do j = 1, r
          do component = 1, 3
             brinkman(i,j) = brinkman(i,j) + vector_glsc3( &
                  m(component,i), effective, m(component,j))
          end do
       end do
       do component = 1, 3
          lift_vector(i) = lift_vector(i) + vector_glsc3( &
               m(component,i), effective, l(component))
       end do
    end do
    do j = 1, r
       do component = 1, 3; call m(component,j)%free(); end do
    end do
    do component = 1, 3; call l(component)%free(); end do
    deallocate(m,l)
    call w%free(); call c%free(); call effective%free()
  end subroutine test_assemble_brinkman_legacy

  !> Expose lift-convection assembly to the array-based test harness.
  subroutine test_assemble_lift_convection(weights, modes, gradients, lift, linear)
    real(kind=rp), intent(in) :: weights(:), modes(:, :, :), gradients(:, :, :, :), lift(:, :)
    real(kind=rp), intent(out) :: linear(:, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: m(:,:), g(:,:,:), l(:), lg(:,:)
    type(matrix_t) :: result
    type(vector_t) :: result_constant
    integer :: n,r,i,j,k,d
    n=size(weights); r=size(modes,3); call w%init(n); w =weights
    allocate(m(3,r),g(3,r,3),l(3),lg(3,3))
    do j=1,r; do i=1,3
       call m(i,j)%init(n); m(i,j) =modes(:,i,j)
       do d=1,3; call g(i,j,d)%init(n); g(i,j,d) =gradients(:,i,j,d); end do
    end do; end do
    do i=1,3; call l(i)%init(n); l(i) =lift(:,i); end do
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
    n=size(weights); r=size(gradients,3); call w%init(n); w =weights
    allocate(g(3,r,3))
    do j=1,r; do i=1,3; do d=1,3; call g(i,j,d)%init(n); g(i,j,d) =gradients(:,i,j,d); end do; end do; end do
    call rom_assemble_diffusion(w,nu,g,result); diffusion=result%x
    do j=1,r; do i=1,3; do d=1,3; call g(i,j,d)%free(); end do; end do; end do
    deallocate(g); call w%free(); call result%free()
  end subroutine test_assemble_diffusion

  !> Legacy per-contraction reference for current packed velocity diffusion.
  !! The inputs and contraction ordering match rom_assemble_diffusion; only
  !! the reduction strategy differs (global vector_glsc3 for each term).
  subroutine test_assemble_diffusion_legacy(weights, nu, gradients, diffusion)
    real(kind=rp), intent(in) :: weights(:), nu, gradients(:, :, :, :)
    real(kind=rp), intent(out) :: diffusion(:, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: g(:,:,:)
    integer :: n, r, i, j, component, direction

    n = size(weights); r = size(gradients,3)
    call w%init(n); w = weights
    allocate(g(3,r,3))
    do j = 1, r
       do component = 1, 3
          do direction = 1, 3
             call g(component,j,direction)%init(n)
             g(component,j,direction) = gradients(:,component,j,direction)
          end do
       end do
    end do
    diffusion = 0.0_rp
    do i = 1, r
       do j = 1, r
          do component = 1, 3
             do direction = 1, 3
                diffusion(i,j) = diffusion(i,j) - nu * vector_glsc3( &
                     g(component,i,direction), w, g(component,j,direction))
             end do
          end do
       end do
    end do
    do j = 1, r
       do component = 1, 3
          do direction = 1, 3
             call g(component,j,direction)%free()
          end do
       end do
    end do
    deallocate(g); call w%free()
  end subroutine test_assemble_diffusion_legacy

  !> Expose convection assembly to the array-based test harness.
  subroutine test_assemble_convection(weights, modes, gradients, quadratic)
    real(kind=rp), intent(in) :: weights(:), modes(:, :, :), gradients(:, :, :, :)
    real(kind=rp), intent(out) :: quadratic(:, :, :)
    type(vector_t) :: w
    type(vector_t), allocatable :: m(:,:), g(:,:,:)
    type(matrix_t), allocatable :: result(:)
    integer :: n,r,i,j,d
    n=size(weights); r=size(modes,3); call w%init(n); w =weights
    allocate(m(3,r),g(3,r,3),result(r))
    do j=1,r; do i=1,3
       call m(i,j)%init(n); m(i,j) =modes(:,i,j)
       do d=1,3; call g(i,j,d)%init(n); g(i,j,d) =gradients(:,i,j,d); end do
    end do; end do
    call rom_assemble_convection(w,m,g,result)
    do i=1,r; quadratic(i,:,:)=result(i)%x; end do
    do j=1,r; do i=1,3; do d=1,3; call g(i,j,d)%free(); end do; call m(i,j)%free(); end do; end do
    do i=1,r; call result(i)%free(); end do
    deallocate(m,g,result); call w%free()
  end subroutine test_assemble_convection

  !> Legacy per-entry global-reduction reference for velocity convection.
  !! Pointwise products and accumulation follow the committed production
  !! ordering exactly; only vector_glsc3 reduction is retained per entry.
  subroutine test_assemble_convection_legacy(weights, modes, gradients, quadratic)
    real(kind=rp), intent(in) :: weights(:), modes(:, :, :), gradients(:, :, :, :)
    real(kind=rp), intent(out) :: quadratic(:, :, :)
    type(vector_t) :: w, pointwise, product
    type(vector_t), allocatable :: m(:,:), g(:,:,:)
    integer :: n, r, test_mode, advect_mode, state_mode, field, component

    n = size(weights); r = size(modes,3)
    call w%init(n); w =weights
    call pointwise%init(n); call product%init(n)
    allocate(m(3,r),g(3,r,3))
    do state_mode = 1, r
       do component = 1, 3
          call m(component,state_mode)%init(n)
          m(component,state_mode) = modes(:,component,state_mode)
          do field = 1, 3
             call g(component,state_mode,field)%init(n)
             g(component,state_mode,field) = gradients(:,component,state_mode,field)
          end do
       end do
    end do
    quadratic = 0.0_rp
    do test_mode = 1, r
       do advect_mode = 1, r
          do state_mode = 1, r
             do field = 1, 3
                pointwise = 0.0_rp
                do component = 1, 3
                   call vector_col3(product, m(component,advect_mode), &
                        g(field,state_mode,component))
                   call vector_add2(pointwise, product)
                end do
                quadratic(test_mode,advect_mode,state_mode) = &
                     quadratic(test_mode,advect_mode,state_mode) + &
                     vector_glsc3(m(field,test_mode), w, pointwise)
             end do
          end do
       end do
    end do
    do state_mode = 1, r
       do component = 1, 3
          do field = 1, 3; call g(component,state_mode,field)%free(); end do
          call m(component,state_mode)%free()
       end do
    end do
    deallocate(m,g)
    call w%free(); call pointwise%free(); call product%free()
  end subroutine test_assemble_convection_legacy

  !> Expose viscous objective assembly to the array-based test harness.
  subroutine test_assemble_viscous_quadratic(weights, mask, viscosity, gradients, quadratic)
    real(kind=rp), intent(in) :: weights(:), mask(:), viscosity, gradients(:, :, :, :)
    real(kind=rp), intent(out) :: quadratic(:, :)
    type(vector_t) :: w,m
    type(vector_t), allocatable :: g(:,:,:)
    type(matrix_t) :: result
    integer :: n,r,i,j,d
    n=size(weights); r=size(gradients,3); call w%init(n); call m%init(n); w =weights; m =mask
    allocate(g(3,r,3))
    do j=1,r; do i=1,3; do d=1,3; call g(i,j,d)%init(n); g(i,j,d) =gradients(:,i,j,d); end do; end do; end do
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
    w =weights; m =mask; c =chi; allocate(v(3,r),l(3))
    do j=1,r; do i=1,3; call v(i,j)%init(n); v(i,j) =modes(:,i,j); end do; end do
    do i=1,3; call l(i)%init(n); l(i) =lift(:,i); end do
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
    s =state; l=0.0_rp; q%x =quadratic; constant=0.0_rp
    value=rom_viscous_value(s,constant,l,q,domain)
    call s%free(); call l%free(); call q%free()
  end function test_viscous_value

  !> Expose Brinkman objective evaluation to the array-based test harness.
  function test_brinkman_value(state, constant, lift_linear, quadratic, domain) result(value)
    real(kind=rp), intent(in) :: state(:), constant, lift_linear(:), quadratic(:, :), domain
    real(kind=rp) :: value
    type(vector_t) :: s,l; type(matrix_t) :: q
    call s%init(size(state)); call l%init(size(state)); call q%init(size(state),size(state))
    s =state; l =lift_linear; q%x =quadratic
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
    w =weights; m =mask; s =state; allocate(v(3,r),l(3))
    do j=1,r; do i=1,3; call v(i,j)%init(n); v(i,j) =modes(:,i,j); end do; end do
    do i=1,3; call l(i)%init(n); l(i) =lift(:,i); end do
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
    w = weights; m = mask; l = lift
    allocate(s(r))
    do i = 1, r
       call s(i)%init(n); s(i) = modes(:, i)
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
    s = state; l = linear; q%x = quadratic
    value = rom_scalar_mixing_value(s, constant, l, q, domain)
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
    call w%init(n); w = weights
    call sl%init(n); sl = scalar_lift
    call diff%init(n); diff = diffusivity
    allocate(vm(3,r), sm(r), vg(3,r,3), sg(r,3), lv(3), lsg(3))
    do j = 1, r
       call sm(j)%init(n); sm(j) = scalar_modes(:,j)
       do c = 1, 3
          call vm(c,j)%init(n); vm(c,j) = velocity_modes(:,c,j)
          do d = 1, 3
             call vg(c,j,d)%init(n); vg(c,j,d) = velocity_gradients(:,c,j,d)
          end do
       end do
       do d = 1, 3
          call sg(j,d)%init(n); sg(j,d) = scalar_gradients(:,j,d)
       end do
    end do
    do c = 1, 3
       call lv(c)%init(n); lv(c) = lift_velocity(:,c)
       call lsg(c)%init(n); lsg(c) = lift_scalar_gradient(:,c)
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

  !> Legacy per-contraction reference for all scalar-transport outputs.
  !! The pointwise operations and component/direction order match the
  !! committed scalar transport implementation; vector_glsc3 retains one
  !! global reduction per scalar contraction.
  subroutine test_assemble_scalar_transport_legacy(weights, scalar_weight, &
       velocity_modes, scalar_modes, velocity_gradients, scalar_gradients, &
       lift_velocity, scalar_lift, lift_scalar_gradient, diffusivity, &
       linear, constant, quadratic, implicit_linear)
    real(kind=rp), intent(in) :: weights(:), scalar_weight
    real(kind=rp), intent(in) :: velocity_modes(:, :, :), scalar_modes(:, :)
    real(kind=rp), intent(in) :: velocity_gradients(:, :, :, :), scalar_gradients(:, :, :)
    real(kind=rp), intent(in) :: lift_velocity(:, :), scalar_lift(:)
    real(kind=rp), intent(in) :: lift_scalar_gradient(:, :), diffusivity(:)
    real(kind=rp), intent(out) :: linear(:, :), constant(:), quadratic(:, :, :)
    real(kind=rp), intent(out) :: implicit_linear(:, :)
    type(vector_t) :: w, effective, product, pointwise, diff, scalar_lift_vec
    type(vector_t), allocatable :: vm(:,:), sm(:), vg(:,:,:), sg(:,:), lv(:), lsg(:)
    integer :: n, r, i, j, k, direction

    n = size(weights); r = size(scalar_modes,2)
    call w%init(n); w = weights
    call effective%init(n); effective = w
    call vector_cmult(effective, scalar_weight)
    call product%init(n); call pointwise%init(n)
    call diff%init(n); diff = diffusivity
    call scalar_lift_vec%init(n); scalar_lift_vec = scalar_lift
    allocate(vm(3,r), sm(r), vg(3,r,3), sg(r,3), lv(3), lsg(3))
    do j = 1, r
       call sm(j)%init(n); sm(j) = scalar_modes(:,j)
       do i = 1, 3
          call vm(i,j)%init(n); vm(i,j) = velocity_modes(:,i,j)
       end do
       do direction = 1, 3
          call sg(j,direction)%init(n)
          sg(j,direction) = scalar_gradients(:,j,direction)
          do i = 1, 3
             call vg(i,j,direction)%init(n)
             vg(i,j,direction) = velocity_gradients(:,i,j,direction)
          end do
       end do
    end do
    do i = 1, 3
       call lv(i)%init(n); lv(i) = lift_velocity(:,i)
       call lsg(i)%init(n); lsg(i) = lift_scalar_gradient(:,i)
    end do
    linear = 0.0_rp; constant = 0.0_rp
    quadratic = 0.0_rp; implicit_linear = 0.0_rp

    do i = 1, r
       pointwise = 0.0_rp
       do direction = 1, 3
          call vector_col3(product, lv(direction), lsg(direction))
          call vector_add2(pointwise, product)
       end do
       constant(i) = constant(i) - &
            vector_glsc3(sm(i), effective, pointwise)
       do j = 1, r
          pointwise = 0.0_rp
          do direction = 1, 3
             call vector_col3(product, vm(direction,j), lsg(direction))
             call vector_add2(pointwise, product)
             call vector_col3(product, lv(direction), sg(j,direction))
             call vector_add2(pointwise, product)
          end do
          linear(i,j) = linear(i,j) - &
               vector_glsc3(sm(i), effective, pointwise)
       end do
    end do
    do i = 1, r
       do j = 1, r
          do k = 1, r
             pointwise = 0.0_rp
             do direction = 1, 3
                call vector_col3(product, vm(direction,j), sg(k,direction))
                call vector_add2(pointwise, product)
             end do
             quadratic(i,j,k) = quadratic(i,j,k) - &
                  vector_glsc3(sm(i), effective, pointwise)
          end do
       end do
    end do

    call vector_col3(effective, w, diff)
    call vector_cmult(effective, scalar_weight)
    do i = 1, r
       do j = 1, r
          do direction = 1, 3
             implicit_linear(i,j) = implicit_linear(i,j) - &
                  vector_glsc3(sg(i,direction), effective, sg(j,direction))
          end do
       end do
       do direction = 1, 3
          constant(i) = constant(i) - vector_glsc3( &
               sg(i,direction), effective, lsg(direction))
       end do
    end do

    do j = 1, r
       call sm(j)%free()
       do direction = 1, 3
          call sg(j,direction)%free()
          do i = 1, 3; call vg(i,j,direction)%free(); end do
       end do
       do i = 1, 3; call vm(i,j)%free(); end do
    end do
    do i = 1, 3; call lv(i)%free(); call lsg(i)%free(); end do
    deallocate(vm,sm,vg,sg,lv,lsg)
    call w%free(); call effective%free(); call product%free(); call pointwise%free()
    call diff%free(); call scalar_lift_vec%free()
  end subroutine test_assemble_scalar_transport_legacy

end module rom_array_test_adapters
