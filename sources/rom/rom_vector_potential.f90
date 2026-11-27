!> Matrix-free vector-potential basis map.
!! The reference potential uses Neko's Helmholtz/CG infrastructure with h1=1,
!! h2=0. Candidate maps use Neko's native gathered curl. The transpose reverses
!! that exact sequence (Binv, GS_ADD, cdtp, whose kernel applies local B) and
!! the componentwise homogeneous velocity boundary projector. CUDA work uses
!! Neko's device-resident cdtp/curl/GS/Helmholtz paths; only the one-time
!! geometric metric setup is copied from host to device.
module rom_vector_potential
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
  use num_types, only: rp
  use coefs, only: coef_t
  use vector, only: vector_t
  use field, only: field_t
  use ax_product, only: ax_t, ax_helm_allocator
  use krylov, only: ksp_t, ksp_monitor_t, krylov_solver_factory
  use precon, only: pc_t, precon_allocator
  use jacobi, only: jacobi_t
  use device_jacobi, only: device_jacobi_t
  use scalar_bc_projector, only: scalar_bc_projector_t
  use vector_bc_projector, only: vector_bc_projector_t, &
       segregated_vector_bc_projector_t
  use operators, only: curl, cdtp
  use vector_math, only: vector_copy, vector_col2, vector_add2
  use field_math, only: field_rzero, field_copy, field_sub2, field_glsc2
  use device, only: device_map, device_unmap, device_memcpy, device_sync, &
       HOST_TO_DEVICE
  use device_math, only: device_copy, device_rzero, device_add2, &
       device_cmult, device_col2, device_glsc3
  use gs_ops, only: GS_OP_ADD
  use utils, only: neko_error
  use mpi_f08, only: MPI_Allreduce, MPI_MIN, MPI_MAX, MPI_SUM, MPI_INTEGER
  use comm, only: NEKO_COMM, MPI_REAL_PRECISION
  use neko_config, only: NEKO_BCKND_DEVICE, NEKO_BCKND_SX, NEKO_BCKND_XSMM
  implicit none
  private

  type, public :: rom_vector_potential_t
     type(vector_t), allocatable :: potential(:,:) ! (component, mode)
     type(coef_t) :: poisson_coef                  ! private operator state
     class(ax_t), allocatable :: poisson_A
     class(ksp_t), allocatable :: poisson_ksp
     class(pc_t), allocatable :: poisson_pc
     type(scalar_bc_projector_t) :: poisson_bc(3)
     class(vector_bc_projector_t), pointer :: velocity_bc => null()
     integer :: nmode = 0
     integer :: solve_count = 0
     integer, allocatable :: solve_iterations(:,:)
     real(kind=rp), allocatable :: solve_residual(:,:)
     real(kind=rp), allocatable :: solve_true_residual(:,:)
     real(kind=rp), allocatable :: solve_true_residual_normalized(:,:)
     ! Persistent element-layout metric views for native device cdtp calls.
     ! They are copied once at reference-potential construction and reused by
     ! every exact pullback in the outer iteration.
     real(kind=rp), allocatable :: drx(:,:), dsx(:,:), dtx(:,:)
     real(kind=rp), allocatable :: dry(:,:), dsy(:,:), dty(:,:)
     real(kind=rp), allocatable :: drz(:,:), dsz(:,:), dtz(:,:)
     type(c_ptr) :: drx_d = c_null_ptr, dsx_d = c_null_ptr, dtx_d = c_null_ptr
     type(c_ptr) :: dry_d = c_null_ptr, dsy_d = c_null_ptr, dty_d = c_null_ptr
     type(c_ptr) :: drz_d = c_null_ptr, dsz_d = c_null_ptr, dtz_d = c_null_ptr
     logical :: ready = .false.
   contains
     procedure :: construct => rom_vp_construct
     procedure :: apply => rom_vp_apply
     procedure :: pullback => rom_vp_pullback
     procedure :: free => rom_vp_free
  end type rom_vector_potential_t

contains

  subroutine build_box_normal_masks(projector, coef)
    type(scalar_bc_projector_t), intent(inout) :: projector(3)
    type(coef_t), intent(in) :: coef
    real(kind=rp), allocatable :: xyz(:,:)
    real(kind=rp) :: lo(3),hi(3),glo(3),ghi(3),tol
    integer, allocatable :: ids(:)
    integer :: n,c,i,k,j,e,lx,ly,lz,ierr,local_bad,global_bad

    n=coef%dof%size()
    allocate(xyz(3,n),ids(n))
    xyz(1,:)=reshape(coef%dof%x%x,[n])
    xyz(2,:)=reshape(coef%dof%y%x,[n])
    xyz(3,:)=reshape(coef%dof%z%x,[n])
    do c=1,3
       lo(c)=minval(xyz(c,:)); hi(c)=maxval(xyz(c,:))
    end do
    call MPI_Allreduce(lo,glo,3,MPI_REAL_PRECISION,MPI_MIN,NEKO_COMM,ierr)
    call MPI_Allreduce(hi,ghi,3,MPI_REAL_PRECISION,MPI_MAX,NEKO_COMM,ierr)
    ! Coordinate-extrema faces are valid normal faces only for a globally
    ! axis-aligned tensor-product box. Reject rotated/curved elements instead
    ! of silently imposing the wrong component boundary condition.
    lx=coef%Xh%lx; ly=coef%Xh%ly; lz=coef%Xh%lz
    tol=256.0_rp*epsilon(1.0_rp)*max(1.0_rp,maxval(abs(glo)),maxval(abs(ghi)))
    local_bad=0
    do e=1,coef%msh%nelv
       do k=1,lz; do j=1,ly; do i=1,lx
          if (abs(coef%dof%x%x(i,j,k,e)-coef%dof%x%x(i,1,1,e))>tol .or. &
               abs(coef%dof%y%x(i,j,k,e)-coef%dof%y%x(1,j,1,e))>tol .or. &
               abs(coef%dof%z%x(i,j,k,e)-coef%dof%z%x(1,1,k,e))>tol) &
               local_bad=1
       end do; end do; end do
    end do
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,NEKO_COMM,ierr)
    if (global_bad/=0) call neko_error( &
         'Vector-potential Poisson requires an axis-aligned Cartesian box')
    do c=1,3
       tol=128.0_rp*epsilon(1.0_rp)*max(1.0_rp,abs(glo(c)),abs(ghi(c)))
       k=0
       do i=1,n
          if (abs(xyz(c,i)-glo(c))<=tol .or. &
               abs(xyz(c,i)-ghi(c))<=tol) then
             k=k+1; ids(k)=i
          end if
       end do
       if (k>0) call projector(c)%dof_mask%set(ids(1:k),k)
    end do
  end subroutine build_box_normal_masks

  subroutine rom_vp_construct(this, coef, velocity_modes, velocity_bc, &
       tolerance, max_iter)
    class(rom_vector_potential_t), intent(inout) :: this
    type(coef_t), target, intent(in) :: coef
    type(vector_t), intent(in) :: velocity_modes(:,:)
    class(vector_bc_projector_t), target, intent(inout) :: velocity_bc
    real(kind=rp), intent(in) :: tolerance
    integer, intent(in) :: max_iter
    type(field_t) :: rhs, sol, axsol, residual
    type(ksp_monitor_t) :: mon
    real(kind=rp), allocatable :: weak(:,:,:,:)
    real(kind=rp), allocatable :: dtx(:,:), x(:,:), vec(:)
    real(kind=rp), allocatable :: drx(:,:),dsx(:,:),dtxx(:,:), &
         dry(:,:),dsy(:,:),dty(:,:),drz(:,:),dsz(:,:),dtz(:,:)
    integer :: n, nm, m, c

    call this%free()
    if (NEKO_BCKND_SX /= 0 .or. NEKO_BCKND_XSMM /= 0) &
         call neko_error('Vector-potential transpose is not enabled for SX or XSMM')
    if (coef%msh%gdim /= 3 .or. coef%cyclic) &
         call neko_error('Vector potential supports noncyclic 3-D Cartesian meshes only')
    select type (velocity_bc)
    type is (segregated_vector_bc_projector_t)
       continue
    class default
       call neko_error('Vector potential requires a segregated homogeneous velocity projector')
    end select
    if (size(velocity_modes,1) /= 3 .or. size(velocity_modes,2) < 1) &
         call neko_error('Invalid vector-potential reference mode array')
    if (tolerance <= 0.0_rp .or. max_iter < 1) &
         call neko_error('Invalid vector-potential Poisson controls')

    this%velocity_bc => velocity_bc
    this%nmode = size(velocity_modes,2)
    n = coef%dof%size()
    nm = coef%msh%nelv
    allocate(this%potential(3,this%nmode))
    allocate(this%solve_iterations(3,this%nmode), this%solve_residual(3,this%nmode))
    allocate(this%solve_true_residual(3,this%nmode))
    allocate(this%solve_true_residual_normalized(3,this%nmode))
    this%solve_iterations = 0
    this%solve_residual = huge(1.0_rp)
    this%solve_true_residual = huge(1.0_rp)
    this%solve_true_residual_normalized = huge(1.0_rp)
    call build_box_normal_masks(this%poisson_bc, coef)

    ! A distinct coefficient object borrows geometry/GS topology but owns h1/h2.
    ! The live fluid coefficient state is never changed.
    call this%poisson_coef%init(coef%gs_h)
    this%poisson_coef%h1 = 1.0_rp
    this%poisson_coef%h2 = 0.0_rp
    this%poisson_coef%ifh2 = .false.
    call ax_helm_allocator(this%poisson_A, 'standard')
    call krylov_solver_factory(this%poisson_ksp, n, 'cg', max_iter, tolerance)
    call precon_allocator(this%poisson_pc, 'jacobi')
    select type (pcp => this%poisson_pc)
    type is (jacobi_t)
       call pcp%init_from_components(this%poisson_coef, &
            this%poisson_coef%dof, this%poisson_coef%gs_h)
    type is (device_jacobi_t)
       call pcp%init_from_components(this%poisson_coef, &
            this%poisson_coef%dof, this%poisson_coef%gs_h)
    class default
       call neko_error('Unexpected vector-potential Poisson preconditioner')
    end select
    call this%poisson_ksp%set_pc(this%poisson_pc)

    call rhs%init(this%poisson_coef%dof, 'vp_rhs')
    call sol%init(this%poisson_coef%dof, 'vp_sol')
    call axsol%init(this%poisson_coef%dof, 'vp_axsol')
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call rom_vp_prepare_device_metrics(this)
       call rom_vp_construct_device(this, velocity_modes, rhs, sol, axsol, residual, &
            tolerance, n, nm)
       call rhs%free(); call sol%free(); call axsol%free(); call residual%free()
       this%ready = .true.
       return
    end if
    allocate(weak(coef%Xh%lx,coef%Xh%ly,coef%Xh%lz,nm))
    allocate(dtx(coef%Xh%lxyz,nm), x(coef%Xh%lxyz,nm))
    allocate(vec(n))
    allocate(drx(coef%Xh%lxyz,nm),dsx(coef%Xh%lxyz,nm),dtxx(coef%Xh%lxyz,nm))
    allocate(dry(coef%Xh%lxyz,nm),dsy(coef%Xh%lxyz,nm),dty(coef%Xh%lxyz,nm))
    allocate(drz(coef%Xh%lxyz,nm),dsz(coef%Xh%lxyz,nm),dtz(coef%Xh%lxyz,nm))
    drx=reshape(coef%drdx,[coef%Xh%lxyz,nm]); dsx=reshape(coef%dsdx,[coef%Xh%lxyz,nm])
    dtxx=reshape(coef%dtdx,[coef%Xh%lxyz,nm]); dry=reshape(coef%drdy,[coef%Xh%lxyz,nm])
    dsy=reshape(coef%dsdy,[coef%Xh%lxyz,nm]); dty=reshape(coef%dtdy,[coef%Xh%lxyz,nm])
    drz=reshape(coef%drdz,[coef%Xh%lxyz,nm]); dsz=reshape(coef%dsdz,[coef%Xh%lxyz,nm])
    dtz=reshape(coef%dtdz,[coef%Xh%lxyz,nm])

    do m = 1, this%nmode
       do c = 1, 3
          if (velocity_modes(c,m)%size() /= n) &
               call neko_error('Vector-potential mode has inconsistent local size')
       end do
       ! RHS_c = (grad b, e_c x v), using Neko's weighted discrete
       ! derivative transpose. cdtp includes the local SEM mass weighting.
       do c = 1, 3
          weak=0.0_rp
          select case(c)
          case(1)
             x = reshape(-velocity_modes(3,m)%x, shape(x))
             dtx = 0.0_rp
             call cdtp(dtx,x,dry,dsy,dty,coef)
             weak = reshape(dtx,shape(weak))
             x = reshape(velocity_modes(2,m)%x,shape(x))
             dtx = 0.0_rp
             call cdtp(dtx,x,drz,dsz,dtz,coef)
             weak = weak + reshape(dtx,shape(weak))
          case(2)
             x = reshape(velocity_modes(3,m)%x,shape(x))
             dtx = 0.0_rp
             call cdtp(dtx,x,drx,dsx,dtxx,coef)
             weak = reshape(dtx,shape(weak))
             x = reshape(-velocity_modes(1,m)%x,shape(x))
             dtx = 0.0_rp
             call cdtp(dtx,x,drz,dsz,dtz,coef)
             weak = weak + reshape(dtx,shape(weak))
          case(3)
             x = reshape(-velocity_modes(2,m)%x,shape(x))
             dtx = 0.0_rp
             call cdtp(dtx,x,drx,dsx,dtxx,coef)
             weak = reshape(dtx,shape(weak))
             x = reshape(velocity_modes(1,m)%x,shape(x))
             dtx = 0.0_rp
             call cdtp(dtx,x,dry,dsy,dty,coef)
             weak = weak + reshape(dtx,shape(weak))
          end select
          call this%poisson_coef%gs_h%op(weak,n,GS_OP_ADD)
          vec=reshape(weak,[n])
          call this%poisson_bc(c)%apply(vec,n)
          weak=reshape(vec,shape(weak))
          rhs%x = reshape(weak,shape(rhs%x))
          sol%x = 0.0_rp
          mon = this%poisson_ksp%solve(this%poisson_A, sol, rhs%x, n, &
               this%poisson_coef, this%poisson_bc(c), this%poisson_coef%gs_h)
          this%solve_count = this%solve_count + 1
          this%solve_iterations(c,m) = mon%iter
          this%solve_residual(c,m) = mon%res_final
          if (.not. mon%converged) &
               call neko_error('Vector-potential Poisson solve did not converge')
          call this%poisson_A%compute(axsol%x,sol%x,this%poisson_coef, &
               this%poisson_coef%msh,this%poisson_coef%Xh)
          call this%poisson_coef%gs_h%op(axsol,GS_OP_ADD)
          vec=reshape(axsol%x-rhs%x,[n])
          call this%poisson_bc(c)%apply(vec,n)
          this%solve_true_residual(c,m)=relative_algebraic_residual( &
               vec,reshape(rhs%x,[n]))
          this%solve_true_residual_normalized(c,m)= &
               normalized_algebraic_residual(vec,this%poisson_coef)
          if (this%solve_true_residual_normalized(c,m)> &
               max(1.0e-8_rp,100*tolerance)) then
             write(*,'(A,I0,A,I0,A,I0,A,ES14.6,A,ES14.6)') &
                  'Vector-potential Poisson true residual failed: mode=',m, &
                  ' component=',c,' iterations=',mon%iter, &
                  ' true=',this%solve_true_residual(c,m), &
                  ' monitored=',mon%res_final
             call neko_error('Vector-potential Poisson independent residual failed')
          end if
          call this%potential(c,m)%init(n,name='vp_potential')
          this%potential(c,m)%x = reshape(sol%x,[n])
       end do
    end do
    call rhs%free()
    call sol%free()
    call axsol%free()
    this%ready = .true.
  end subroutine rom_vp_construct

  subroutine rom_vp_prepare_device_metrics(this)
    class(rom_vector_potential_t), intent(inout) :: this
    integer :: lxyz, nel, n
    lxyz = this%poisson_coef%Xh%lxyz
    nel = this%poisson_coef%msh%nelv
    n = lxyz * nel
    allocate(this%drx(lxyz,nel), this%dsx(lxyz,nel), this%dtx(lxyz,nel))
    allocate(this%dry(lxyz,nel), this%dsy(lxyz,nel), this%dty(lxyz,nel))
    allocate(this%drz(lxyz,nel), this%dsz(lxyz,nel), this%dtz(lxyz,nel))
    this%drx = reshape(this%poisson_coef%drdx,[lxyz,nel])
    this%dsx = reshape(this%poisson_coef%dsdx,[lxyz,nel])
    this%dtx = reshape(this%poisson_coef%dtdx,[lxyz,nel])
    this%dry = reshape(this%poisson_coef%drdy,[lxyz,nel])
    this%dsy = reshape(this%poisson_coef%dsdy,[lxyz,nel])
    this%dty = reshape(this%poisson_coef%dtdy,[lxyz,nel])
    this%drz = reshape(this%poisson_coef%drdz,[lxyz,nel])
    this%dsz = reshape(this%poisson_coef%dsdz,[lxyz,nel])
    this%dtz = reshape(this%poisson_coef%dtdz,[lxyz,nel])
    call device_map(this%drx,this%drx_d,n); call device_memcpy(this%drx,this%drx_d,n, &
         HOST_TO_DEVICE,.true.)
    call device_map(this%dsx,this%dsx_d,n); call device_memcpy(this%dsx,this%dsx_d,n, &
         HOST_TO_DEVICE,.true.)
    call device_map(this%dtx,this%dtx_d,n); call device_memcpy(this%dtx,this%dtx_d,n, &
         HOST_TO_DEVICE,.true.)
    call device_map(this%dry,this%dry_d,n); call device_memcpy(this%dry,this%dry_d,n, &
         HOST_TO_DEVICE,.true.)
    call device_map(this%dsy,this%dsy_d,n); call device_memcpy(this%dsy,this%dsy_d,n, &
         HOST_TO_DEVICE,.true.)
    call device_map(this%dty,this%dty_d,n); call device_memcpy(this%dty,this%dty_d,n, &
         HOST_TO_DEVICE,.true.)
    call device_map(this%drz,this%drz_d,n); call device_memcpy(this%drz,this%drz_d,n, &
         HOST_TO_DEVICE,.true.)
    call device_map(this%dsz,this%dsz_d,n); call device_memcpy(this%dsz,this%dsz_d,n, &
         HOST_TO_DEVICE,.true.)
    call device_map(this%dtz,this%dtz_d,n); call device_memcpy(this%dtz,this%dtz_d,n, &
         HOST_TO_DEVICE,.true.)
  end subroutine rom_vp_prepare_device_metrics

  subroutine rom_vp_construct_device(this, velocity_modes, rhs, sol, axsol, &
       residual, tolerance, n, nm)
    class(rom_vector_potential_t), intent(inout) :: this
    type(vector_t), intent(in) :: velocity_modes(:,:)
    type(field_t), intent(inout) :: rhs, sol, axsol, residual
    real(kind=rp), intent(in) :: tolerance
    integer, intent(in) :: n, nm
    type(ksp_monitor_t) :: mon
    real(kind=rp), allocatable :: weak(:,:,:,:), dtx(:,:), x(:,:), vec(:)
    type(c_ptr) :: weak_d, dtx_d, x_d, vec_d
    real(kind=rp) :: rhs_norm2, residual_norm2
    real(kind=rp) :: residual_weighted_norm2
    integer :: lxyz, m, c

    lxyz = this%poisson_coef%Xh%lxyz
    allocate(weak(this%poisson_coef%Xh%lx,this%poisson_coef%Xh%ly, &
         this%poisson_coef%Xh%lz,nm))
    allocate(dtx(lxyz,nm),x(lxyz,nm),vec(n))
    weak_d=c_null_ptr; dtx_d=c_null_ptr; x_d=c_null_ptr; vec_d=c_null_ptr
    call device_map(weak,weak_d,n); call device_map(dtx,dtx_d,n)
    call device_map(x,x_d,n); call device_map(vec,vec_d,n)
    call residual%init(this%poisson_coef%dof,'vp_residual')

    do m=1,this%nmode
       do c=1,3
          if (velocity_modes(c,m)%size()/=n) &
               call neko_error('Vector-potential mode has inconsistent local size')
       end do
       do c=1,3
          call device_rzero(weak_d,n)
          select case(c)
          case(1)
             call device_copy(x_d,velocity_modes(3,m)%x_d,n)
             call device_cmult(x_d,-1.0_rp,n)
             call device_rzero(dtx_d,n)
             call cdtp(dtx,x,this%dry,this%dsy,this%dty,this%poisson_coef)
             call device_copy(weak_d,dtx_d,n)
             call device_copy(x_d,velocity_modes(2,m)%x_d,n)
             call device_rzero(dtx_d,n)
             call cdtp(dtx,x,this%drz,this%dsz,this%dtz,this%poisson_coef)
             call device_add2(weak_d,dtx_d,n)
          case(2)
             call device_copy(x_d,velocity_modes(3,m)%x_d,n)
             call device_rzero(dtx_d,n)
             call cdtp(dtx,x,this%drx,this%dsx,this%dtx,this%poisson_coef)
             call device_copy(weak_d,dtx_d,n)
             call device_copy(x_d,velocity_modes(1,m)%x_d,n)
             call device_cmult(x_d,-1.0_rp,n)
             call device_rzero(dtx_d,n)
             call cdtp(dtx,x,this%drz,this%dsz,this%dtz,this%poisson_coef)
             call device_add2(weak_d,dtx_d,n)
          case(3)
             call device_copy(x_d,velocity_modes(2,m)%x_d,n)
             call device_cmult(x_d,-1.0_rp,n)
             call device_rzero(dtx_d,n)
             call cdtp(dtx,x,this%drx,this%dsx,this%dtx,this%poisson_coef)
             call device_copy(weak_d,dtx_d,n)
             call device_copy(x_d,velocity_modes(1,m)%x_d,n)
             call device_rzero(dtx_d,n)
             call cdtp(dtx,x,this%dry,this%dsy,this%dty,this%poisson_coef)
             call device_add2(weak_d,dtx_d,n)
          end select
          call this%poisson_coef%gs_h%op(weak,n,GS_OP_ADD)
          ! Neko's CUDA gather-scatter may use its own stream. The following
          ! vector copy/KSP setup must observe the assembled right-hand side.
          call device_sync()
          call device_copy(vec_d,weak_d,n)
          call this%poisson_bc(c)%apply(vec,n)
          call device_copy(rhs%x_d,vec_d,n)
          call field_rzero(sol,n)
          mon=this%poisson_ksp%solve(this%poisson_A,sol,rhs%x,n, &
               this%poisson_coef,this%poisson_bc(c),this%poisson_coef%gs_h)
          this%solve_count=this%solve_count+1
          this%solve_iterations(c,m)=mon%iter
          this%solve_residual(c,m)=mon%res_final
          if (.not.mon%converged) &
               call neko_error('Vector-potential Poisson solve did not converge')
          call this%poisson_A%compute(axsol%x,sol%x,this%poisson_coef, &
               this%poisson_coef%msh,this%poisson_coef%Xh)
          call this%poisson_coef%gs_h%op(axsol,GS_OP_ADD)
          call device_sync()
          call field_copy(residual,axsol,n)
          call field_sub2(residual,rhs,n)
          call this%poisson_bc(c)%apply(residual%x,n)
          residual_norm2=field_glsc2(residual,residual,n)
          rhs_norm2=field_glsc2(rhs,rhs,n)
          residual_weighted_norm2=device_glsc3(residual%x_d, &
               this%poisson_coef%mult_d,residual%x_d,n)
          this%solve_true_residual(c,m)=sqrt(max(0.0_rp,residual_norm2))/ &
               max(sqrt(max(0.0_rp,rhs_norm2)),tiny(1.0_rp))
          this%solve_true_residual_normalized(c,m)= &
               sqrt(max(0.0_rp,residual_weighted_norm2))/ &
               sqrt(max(this%poisson_coef%volume,tiny(1.0_rp)))
          if (this%solve_true_residual_normalized(c,m)> &
               max(1.0e-8_rp,100*tolerance)) then
             write(*,'(A,I0,A,I0,A,I0,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6)') &
                  'Vector-potential CUDA Poisson true residual failed: mode=',m, &
                  ' component=',c,' iterations=',mon%iter, &
                  ' relative=',this%solve_true_residual(c,m), &
                  ' normalized=',this%solve_true_residual_normalized(c,m), &
                  ' monitored=',mon%res_final, &
                  ' abs residual=',sqrt(max(0.0_rp,residual_norm2))
             call neko_error('Vector-potential Poisson independent residual failed')
          end if
          call this%potential(c,m)%init(n,name='vp_potential')
          call device_copy(this%potential(c,m)%x_d,sol%x_d,n)
       end do
    end do

    call device_unmap(weak,weak_d); call device_unmap(dtx,dtx_d)
    call device_unmap(x,x_d); call device_unmap(vec,vec_d)
  end subroutine rom_vp_construct_device

  function relative_algebraic_residual(residual,rhs) result(rel)
    real(kind=rp),intent(in)::residual(:),rhs(:)
    real(kind=rp)::rel,local(2),global(2)
    integer::ierr
    local=[sum(residual*residual),sum(rhs*rhs)]
    call MPI_Allreduce(local,global,2,MPI_REAL_PRECISION,MPI_SUM,NEKO_COMM,ierr)
    rel=sqrt(global(1))/max(sqrt(global(2)),tiny(1.0_rp))
  end function relative_algebraic_residual

  function normalized_algebraic_residual(residual,coef) result(rel)
    real(kind=rp),intent(in)::residual(:)
    type(coef_t),intent(in)::coef
    real(kind=rp)::rel,local,global
    integer::ierr
    local=sum(residual*reshape(coef%mult,[size(residual)])*residual)
    call MPI_Allreduce(local,global,1,MPI_REAL_PRECISION,MPI_SUM,NEKO_COMM,ierr)
    rel=sqrt(max(0.0_rp,global))/sqrt(max(coef%volume,tiny(1.0_rp)))
  end function normalized_algebraic_residual

  subroutine rom_vp_apply(this, mode_index, factor, correction)
    class(rom_vector_potential_t), intent(in) :: this
    integer, intent(in) :: mode_index
    type(vector_t), intent(in) :: factor
    type(vector_t), intent(inout) :: correction(3)
    type(vector_t) :: product(3)
    type(field_t) :: a1,a2,a3,w1,w2,w3,t1,t2
    integer :: n
    call vp_check(this,mode_index,factor%size())
    n = factor%size()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call product(1)%init(n,'vp_product_x')
       call product(2)%init(n,'vp_product_y')
       call product(3)%init(n,'vp_product_z')
       call vector_copy(product(1),this%potential(1,mode_index),n)
       call vector_copy(product(2),this%potential(2,mode_index),n)
       call vector_copy(product(3),this%potential(3,mode_index),n)
       call vector_col2(product(1),factor,n)
       call vector_col2(product(2),factor,n)
       call vector_col2(product(3),factor,n)
       call a1%init(this%poisson_coef%dof,'vp_a1')
       call a2%init(this%poisson_coef%dof,'vp_a2')
       call a3%init(this%poisson_coef%dof,'vp_a3')
       call w1%init(this%poisson_coef%dof,'vp_curl1')
       call w2%init(this%poisson_coef%dof,'vp_curl2')
       call w3%init(this%poisson_coef%dof,'vp_curl3')
       call t1%init(this%poisson_coef%dof,'vp_work1')
       call t2%init(this%poisson_coef%dof,'vp_work2')
       call device_copy(a1%x_d,product(1)%x_d,n)
       call device_copy(a2%x_d,product(2)%x_d,n)
       call device_copy(a3%x_d,product(3)%x_d,n)
       call curl(w1,w2,w3,a1,a2,a3,t1,t2,this%poisson_coef)
       if (.not. allocated(correction(1)%x)) call correction(1)%init(n,'vp_corr_x')
       if (.not. allocated(correction(2)%x)) call correction(2)%init(n,'vp_corr_y')
       if (.not. allocated(correction(3)%x)) call correction(3)%init(n,'vp_corr_z')
       call device_copy(correction(1)%x_d,w1%x_d,n)
       call device_copy(correction(2)%x_d,w2%x_d,n)
       call device_copy(correction(3)%x_d,w3%x_d,n)
       call this%velocity_bc%apply(correction(1)%x,correction(2)%x, &
            correction(3)%x,n)
       call product(1)%free(); call product(2)%free(); call product(3)%free()
       call a1%free(); call a2%free(); call a3%free()
       call w1%free(); call w2%free(); call w3%free()
       call t1%free(); call t2%free()
       return
    end if
    call a1%init(this%poisson_coef%dof,'vp_a1')
    call a2%init(this%poisson_coef%dof,'vp_a2')
    call a3%init(this%poisson_coef%dof,'vp_a3')
    call w1%init(this%poisson_coef%dof,'vp_curl1')
    call w2%init(this%poisson_coef%dof,'vp_curl2')
    call w3%init(this%poisson_coef%dof,'vp_curl3')
    call t1%init(this%poisson_coef%dof,'vp_work1')
    call t2%init(this%poisson_coef%dof,'vp_work2')
    a1%x = reshape(factor%x*this%potential(1,mode_index)%x,shape(a1%x))
    a2%x = reshape(factor%x*this%potential(2,mode_index)%x,shape(a2%x))
    a3%x = reshape(factor%x*this%potential(3,mode_index)%x,shape(a3%x))
    call curl(w1,w2,w3,a1,a2,a3,t1,t2,this%poisson_coef)
    if (.not. allocated(correction(1)%x)) call correction(1)%init(n,'vp_corr_x')
    if (.not. allocated(correction(2)%x)) call correction(2)%init(n,'vp_corr_y')
    if (.not. allocated(correction(3)%x)) call correction(3)%init(n,'vp_corr_z')
    correction(1)%x=reshape(w1%x,[n]); correction(2)%x=reshape(w2%x,[n])
    correction(3)%x=reshape(w3%x,[n])
    call this%velocity_bc%apply(correction(1)%x,correction(2)%x, &
         correction(3)%x,n)
    call a1%free(); call a2%free(); call a3%free()
    call w1%free(); call w2%free(); call w3%free()
    call t1%free(); call t2%free()
  end subroutine rom_vp_apply

  subroutine rom_vp_pullback(this, mode_index, factor, correction_cotangent, &
       factor_cotangent)
    class(rom_vector_potential_t), intent(in) :: this
    integer, intent(in) :: mode_index
    type(vector_t), intent(in) :: factor
    type(vector_t), intent(in) :: correction_cotangent(3)
    type(vector_t), intent(inout) :: factor_cotangent
    real(kind=rp), allocatable :: q1(:),q2(:),q3(:),g1(:,:),g2(:,:),g3(:,:),x(:,:)
    real(kind=rp), allocatable :: drx(:,:),dsx(:,:),dtxx(:,:), &
         dry(:,:),dsy(:,:),dty(:,:),drz(:,:),dsz(:,:),dtz(:,:)
    real(kind=rp), allocatable :: t(:,:)
    integer :: n,nel,lxyz,i
    call vp_check(this,mode_index,factor%size())
    n=factor%size(); nel=this%poisson_coef%msh%nelv
    lxyz=this%poisson_coef%Xh%lxyz
    do i=1,3
       if (correction_cotangent(i)%size()/=n) &
            call neko_error('Vector-potential cotangent size mismatch')
    end do
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call rom_vp_pullback_device(this,mode_index,correction_cotangent, &
            factor_cotangent,n)
       return
    end if
    allocate(q1(n),q2(n),q3(n))
    q1=correction_cotangent(1)%x; q2=correction_cotangent(2)%x
    q3=correction_cotangent(3)%x
    call this%velocity_bc%apply(q1,q2,q3,n)
    call opcolv_local(q1,q2,q3,reshape(this%poisson_coef%Binv,[n]),n)
    call this%poisson_coef%gs_h%op(q1,q2,q3,n,GS_OP_ADD)
    allocate(g1(lxyz,nel),g2(lxyz,nel),g3(lxyz,nel),x(lxyz,nel))
    allocate(drx(lxyz,nel),dsx(lxyz,nel),dtxx(lxyz,nel))
    allocate(dry(lxyz,nel),dsy(lxyz,nel),dty(lxyz,nel))
    allocate(drz(lxyz,nel),dsz(lxyz,nel),dtz(lxyz,nel))
    allocate(t(lxyz,nel))
    drx=reshape(this%poisson_coef%drdx,[lxyz,nel]); dsx=reshape(this%poisson_coef%dsdx,[lxyz,nel])
    dtxx=reshape(this%poisson_coef%dtdx,[lxyz,nel]); dry=reshape(this%poisson_coef%drdy,[lxyz,nel])
    dsy=reshape(this%poisson_coef%dsdy,[lxyz,nel]); dty=reshape(this%poisson_coef%dtdy,[lxyz,nel])
    drz=reshape(this%poisson_coef%drdz,[lxyz,nel]); dsz=reshape(this%poisson_coef%dsdz,[lxyz,nel])
    dtz=reshape(this%poisson_coef%dtdz,[lxyz,nel])
    g1=0.0_rp; g2=0.0_rp; g3=0.0_rp
    x=reshape(q2,shape(x)); call cdtp(g1,x,drz,dsz,dtz,this%poisson_coef)
    x=reshape(-q3,shape(x)); t=0.0_rp
    call cdtp(t,x,dry,dsy,dty,this%poisson_coef); g1=g1+t
    x=reshape(q3,shape(x)); call cdtp(g2,x,drx,dsx,dtxx,this%poisson_coef)
    x=reshape(-q1,shape(x)); t=0.0_rp
    call cdtp(t,x,drz,dsz,dtz,this%poisson_coef); g2=g2+t
    x=reshape(q1,shape(x)); call cdtp(g3,x,dry,dsy,dty,this%poisson_coef)
    x=reshape(-q2,shape(x)); t=0.0_rp
    call cdtp(t,x,drx,dsx,dtxx,this%poisson_coef); g3=g3+t
    if (.not. allocated(factor_cotangent%x)) &
         call factor_cotangent%init(n,'vp_factor_cotangent')
    factor_cotangent%x = this%potential(1,mode_index)%x*reshape(g1,[n]) + &
         this%potential(2,mode_index)%x*reshape(g2,[n]) + &
         this%potential(3,mode_index)%x*reshape(g3,[n])
  end subroutine rom_vp_pullback

  subroutine rom_vp_pullback_device(this,mode_index,correction_cotangent, &
       factor_cotangent,n)
    class(rom_vector_potential_t), intent(in) :: this
    integer, intent(in) :: mode_index,n
    type(vector_t), intent(in) :: correction_cotangent(3)
    type(vector_t), intent(inout) :: factor_cotangent
    real(kind=rp), allocatable :: q1(:),q2(:),q3(:),x(:,:),dtx(:,:), &
         g1(:,:),g2(:,:),g3(:,:)
    type(c_ptr) :: q1_d,q2_d,q3_d,x_d,dtx_d,g1_d,g2_d,g3_d
    type(vector_t) :: g1v,g2v,g3v,tmp
    integer :: nel,lxyz

    nel=this%poisson_coef%msh%nelv
    lxyz=this%poisson_coef%Xh%lxyz
    allocate(q1(n),q2(n),q3(n))
    allocate(x(lxyz,nel),dtx(lxyz,nel),g1(lxyz,nel),g2(lxyz,nel),g3(lxyz,nel))
    q1_d=c_null_ptr; q2_d=c_null_ptr; q3_d=c_null_ptr
    x_d=c_null_ptr; dtx_d=c_null_ptr; g1_d=c_null_ptr; g2_d=c_null_ptr; g3_d=c_null_ptr
    call device_map(q1,q1_d,n); call device_map(q2,q2_d,n); call device_map(q3,q3_d,n)
    call device_map(x,x_d,n); call device_map(dtx,dtx_d,n)
    call device_map(g1,g1_d,n); call device_map(g2,g2_d,n); call device_map(g3,g3_d,n)
    call device_copy(q1_d,correction_cotangent(1)%x_d,n)
    call device_copy(q2_d,correction_cotangent(2)%x_d,n)
    call device_copy(q3_d,correction_cotangent(3)%x_d,n)
    call this%velocity_bc%apply(q1,q2,q3,n)
    call device_col2(q1_d,this%poisson_coef%Binv_d,n)
    call device_col2(q2_d,this%poisson_coef%Binv_d,n)
    call device_col2(q3_d,this%poisson_coef%Binv_d,n)
    call this%poisson_coef%gs_h%op(q1,q2,q3,n,GS_OP_ADD)

    ! Exact Euclidean transpose of cdtp uses Neko's weighted discrete
    ! derivative transpose and preserves the forward map's component signs.
    call device_copy(x_d,q2_d,n)
    call device_rzero(g1_d,n)
    call cdtp(g1,x,this%drz,this%dsz,this%dtz,this%poisson_coef)
    call device_copy(x_d,q3_d,n); call device_cmult(x_d,-1.0_rp,n)
    call device_rzero(dtx_d,n)
    call cdtp(dtx,x,this%dry,this%dsy,this%dty,this%poisson_coef)
    call device_add2(g1_d,dtx_d,n)

    call device_copy(x_d,q3_d,n)
    call device_rzero(g2_d,n)
    call cdtp(g2,x,this%drx,this%dsx,this%dtx,this%poisson_coef)
    call device_copy(x_d,q1_d,n); call device_cmult(x_d,-1.0_rp,n)
    call device_rzero(dtx_d,n)
    call cdtp(dtx,x,this%drz,this%dsz,this%dtz,this%poisson_coef)
    call device_add2(g2_d,dtx_d,n)

    call device_copy(x_d,q1_d,n)
    call device_rzero(g3_d,n)
    call cdtp(g3,x,this%dry,this%dsy,this%dty,this%poisson_coef)
    call device_copy(x_d,q2_d,n); call device_cmult(x_d,-1.0_rp,n)
    call device_rzero(dtx_d,n)
    call cdtp(dtx,x,this%drx,this%dsx,this%dtx,this%poisson_coef)
    call device_add2(g3_d,dtx_d,n)

    call g1v%init(n,'vp_pullback_x'); call g2v%init(n,'vp_pullback_y')
    call g3v%init(n,'vp_pullback_z'); call tmp%init(n,'vp_pullback_tmp')
    call device_copy(g1v%x_d,g1_d,n); call device_copy(g2v%x_d,g2_d,n)
    call device_copy(g3v%x_d,g3_d,n)
    if (.not. allocated(factor_cotangent%x)) &
         call factor_cotangent%init(n,'vp_factor_cotangent')
    call vector_copy(factor_cotangent,this%potential(1,mode_index),n)
    call vector_col2(factor_cotangent,g1v,n)
    call vector_copy(tmp,this%potential(2,mode_index),n)
    call vector_col2(tmp,g2v,n)
    call vector_add2(factor_cotangent,tmp,n)
    call vector_copy(tmp,this%potential(3,mode_index),n)
    call vector_col2(tmp,g3v,n)
    call vector_add2(factor_cotangent,tmp,n)

    call g1v%free(); call g2v%free(); call g3v%free(); call tmp%free()
    call device_unmap(q1,q1_d); call device_unmap(q2,q2_d); call device_unmap(q3,q3_d)
    call device_unmap(x,x_d); call device_unmap(dtx,dtx_d)
    call device_unmap(g1,g1_d); call device_unmap(g2,g2_d); call device_unmap(g3,g3_d)
  end subroutine rom_vp_pullback_device

  subroutine opcolv_local(x,y,z,c,n)
    integer,intent(in)::n
    real(kind=rp),intent(inout)::x(n),y(n),z(n)
    real(kind=rp),intent(in)::c(n)
    x=x*c; y=y*c; z=z*c
  end subroutine opcolv_local

  subroutine vp_check(this,mode_index,n)
    class(rom_vector_potential_t),intent(in)::this
    integer,intent(in)::mode_index,n
    if (.not.this%ready .or. mode_index<1 .or. mode_index>this%nmode) &
         call neko_error('Vector-potential object or mode index is invalid')
    if (this%potential(1,mode_index)%size()/=n) &
         call neko_error('Vector-potential field size mismatch')
  end subroutine vp_check

  subroutine rom_vp_free(this)
    class(rom_vector_potential_t),intent(inout)::this
    integer::i,j
    if(allocated(this%potential)) then
       do j=1,size(this%potential,2); do i=1,3
          call this%potential(i,j)%free()
       end do; end do
       deallocate(this%potential)
    end if
    if(allocated(this%poisson_A)) then
       call this%poisson_A%free(); deallocate(this%poisson_A)
    end if
    if(allocated(this%poisson_ksp)) then
       call this%poisson_ksp%free(); deallocate(this%poisson_ksp)
    end if
    if(allocated(this%poisson_pc)) then
       call this%poisson_pc%free(); deallocate(this%poisson_pc)
    end if
    if(allocated(this%drx)) then
       call device_unmap(this%drx,this%drx_d); deallocate(this%drx)
    end if
    if(allocated(this%dsx)) then
       call device_unmap(this%dsx,this%dsx_d); deallocate(this%dsx)
    end if
    if(allocated(this%dtx)) then
       call device_unmap(this%dtx,this%dtx_d); deallocate(this%dtx)
    end if
    if(allocated(this%dry)) then
       call device_unmap(this%dry,this%dry_d); deallocate(this%dry)
    end if
    if(allocated(this%dsy)) then
       call device_unmap(this%dsy,this%dsy_d); deallocate(this%dsy)
    end if
    if(allocated(this%dty)) then
       call device_unmap(this%dty,this%dty_d); deallocate(this%dty)
    end if
    if(allocated(this%drz)) then
       call device_unmap(this%drz,this%drz_d); deallocate(this%drz)
    end if
    if(allocated(this%dsz)) then
       call device_unmap(this%dsz,this%dsz_d); deallocate(this%dsz)
    end if
    if(allocated(this%dtz)) then
       call device_unmap(this%dtz,this%dtz_d); deallocate(this%dtz)
    end if
    call this%poisson_coef%free()
    do i=1,3; call this%poisson_bc(i)%free(); end do
    if(allocated(this%solve_iterations)) deallocate(this%solve_iterations)
    if(allocated(this%solve_residual)) deallocate(this%solve_residual)
    if(allocated(this%solve_true_residual)) deallocate(this%solve_true_residual)
    if(allocated(this%solve_true_residual_normalized)) &
         deallocate(this%solve_true_residual_normalized)
    nullify(this%velocity_bc)
    this%nmode=0; this%solve_count=0; this%ready=.false.
  end subroutine rom_vp_free

end module rom_vector_potential
