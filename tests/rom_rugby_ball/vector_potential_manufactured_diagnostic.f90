!> Manufactured Hodge/Poisson construction check on the case mesh.
!! Builds a divergence-free analytic potential with zero normal trace, applies
!! Neko's native curl to make the input mode, constructs its ROM potential, and
!! reports the independently measured Poisson, curl-reconstruction, gauge, and
!! normal-boundary residuals. This is a construction diagnostic, not an ROM or
!! optimization run.
program vector_potential_manufactured_diagnostic
  use neko, only: neko_init, neko_finalize
  use neko_top, only: neko_top_register_types
  use user, only: user_setup
  use simulation_m, only: simulation_t
  use fluid_pnpn, only: fluid_pnpn_t
  use json_module, only: json_file
  use json_utils_ext, only: json_read_file
  use field, only: field_t
  use vector, only: vector_t
  use operators, only: curl, div
  use math, only: glsc3
  use device, only: HOST_TO_DEVICE, DEVICE_TO_HOST
  use rom_vector_potential, only: rom_vector_potential_t
  use num_types, only: rp
  use comm, only: pe_rank
  use utils, only: neko_error
  implicit none

  type(simulation_t) :: simulation
  type(json_file) :: parameters
  type(fluid_pnpn_t), pointer :: fluid
  type(field_t) :: ax,ay,az,vx,vy,vz,cx,cy,cz,w1,w2,div_a
  type(vector_t) :: modes(3,1)
  type(rom_vector_potential_t) :: vp
  real(kind=rp), allocatable :: xx(:),yy(:),zz(:)
  real(kind=rp) :: lo(3),hi(3),lx,ly,lz,pi
  real(kind=rp) :: reconstruction,gauge,boundary,velocity_norm,potential_norm
  integer, pointer :: ids(:)
  integer :: argc,n,i,c,nranks,ierr
  logical :: found
  character(len=1024) :: case_file

  argc=command_argument_count()
  if(argc/=1) call neko_error( &
       'Usage: vector_potential_manufactured_diagnostic CASE')
  call get_command_argument(1,case_file)
  call neko_init()
  call neko_top_register_types()
  parameters=json_read_file(trim(case_file))
  call parameters%update('unsteady',.false.,found)
  call user_setup(simulation%neko_case%user)
  call simulation%init(parameters)
  select type(f=>simulation%fluid)
  type is(fluid_pnpn_t)
     fluid=>f
  class default
     call neko_error('Manufactured vector-potential diagnostic requires Pn/Pn')
  end select
  call MPI_Comm_size_from_neko(nranks)
  if(nranks/=1) call neko_error( &
       'Manufactured vector-potential diagnostic requires one MPI rank')
  n=fluid%c_Xh%dof%size()
  allocate(xx(n),yy(n),zz(n))
  xx=reshape(fluid%c_Xh%dof%x%x,[n])
  yy=reshape(fluid%c_Xh%dof%y%x,[n])
  zz=reshape(fluid%c_Xh%dof%z%x,[n])
  lo=[minval(xx),minval(yy),minval(zz)]
  hi=[maxval(xx),maxval(yy),maxval(zz)]
  lx=hi(1)-lo(1); ly=hi(2)-lo(2); lz=hi(3)-lo(3)
  if(min(lx,ly,lz)<=0.0_rp) call neko_error('Degenerate test box')
  pi=acos(-1.0_rp)

  call ax%init(fluid%u%dof,'vp_mfg_ax'); call ay%init(fluid%u%dof,'vp_mfg_ay')
  call az%init(fluid%u%dof,'vp_mfg_az'); call vx%init(fluid%u%dof,'vp_mfg_vx')
  call vy%init(fluid%u%dof,'vp_mfg_vy'); call vz%init(fluid%u%dof,'vp_mfg_vz')
  call cx%init(fluid%u%dof,'vp_mfg_cx'); call cy%init(fluid%u%dof,'vp_mfg_cy')
  call cz%init(fluid%u%dof,'vp_mfg_cz'); call w1%init(fluid%u%dof,'vp_mfg_w1')
  call w2%init(fluid%u%dof,'vp_mfg_w2'); call div_a%init(fluid%u%dof,'vp_mfg_diva')
  ax%x=reshape(lx*sin(pi*(xx-lo(1))/lx)*cos(pi*(yy-lo(2))/ly)* &
       sin(pi*(zz-lo(3))/lz),shape(ax%x))
  ay%x=reshape(-ly*cos(pi*(xx-lo(1))/lx)*sin(pi*(yy-lo(2))/ly)* &
       sin(pi*(zz-lo(3))/lz),shape(ay%x))
  az%x=0.0_rp
  call ax%copy_from(HOST_TO_DEVICE,.false.)
  call ay%copy_from(HOST_TO_DEVICE,.false.)
  call az%copy_from(HOST_TO_DEVICE,.true.)
  call curl(vx,vy,vz,ax,ay,az,w1,w2,fluid%c_Xh)
  call vx%copy_from(DEVICE_TO_HOST,.false.)
  call vy%copy_from(DEVICE_TO_HOST,.false.)
  call vz%copy_from(DEVICE_TO_HOST,.true.)
  do c=1,3
     call modes(c,1)%init(n,'vp_mfg_velocity')
  end do
  modes(1,1)=reshape(vx%x,[n]); modes(2,1)=reshape(vy%x,[n])
  modes(3,1)=reshape(vz%x,[n])
  call vp%construct(fluid%c_Xh,modes,fluid%bcs_vel_projector,1.0e-10_rp,500)

  call vp%potential(1,1)%copy_from(DEVICE_TO_HOST,.false.)
  call vp%potential(2,1)%copy_from(DEVICE_TO_HOST,.false.)
  call vp%potential(3,1)%copy_from(DEVICE_TO_HOST,.true.)
  cx%x=reshape(vp%potential(1,1)%x,shape(cx%x))
  cy%x=reshape(vp%potential(2,1)%x,shape(cy%x))
  cz%x=reshape(vp%potential(3,1)%x,shape(cz%x))
  call cx%copy_from(HOST_TO_DEVICE,.false.)
  call cy%copy_from(HOST_TO_DEVICE,.false.)
  call cz%copy_from(HOST_TO_DEVICE,.true.)
  call curl(vx,vy,vz,cx,cy,cz,w1,w2,fluid%c_Xh)
  call vx%copy_from(DEVICE_TO_HOST,.false.)
  call vy%copy_from(DEVICE_TO_HOST,.false.)
  call vz%copy_from(DEVICE_TO_HOST,.true.)
  call div(div_a%x,cx%x,cy%x,cz%x,fluid%c_Xh)
  call div_a%copy_from(DEVICE_TO_HOST,.true.)
  potential_norm=sqrt(max(0.0_rp,glsc3(cx%x,fluid%c_Xh%B,cx%x,n)+ &
       glsc3(cy%x,fluid%c_Xh%B,cy%x,n)+glsc3(cz%x,fluid%c_Xh%B,cz%x,n)))
  velocity_norm=sqrt(max(0.0_rp, &
       glsc3(reshape(modes(1,1)%x,shape(vx%x)),fluid%c_Xh%B, &
            reshape(modes(1,1)%x,shape(vx%x)),n)+ &
       glsc3(reshape(modes(2,1)%x,shape(vy%x)),fluid%c_Xh%B, &
            reshape(modes(2,1)%x,shape(vy%x)),n)+ &
       glsc3(reshape(modes(3,1)%x,shape(vz%x)),fluid%c_Xh%B, &
            reshape(modes(3,1)%x,shape(vz%x)),n)))
  reconstruction=sqrt(max(0.0_rp, &
       glsc3(vx%x-reshape(modes(1,1)%x,shape(vx%x)),fluid%c_Xh%B, &
            vx%x-reshape(modes(1,1)%x,shape(vx%x)),n)+ &
       glsc3(vy%x-reshape(modes(2,1)%x,shape(vy%x)),fluid%c_Xh%B, &
            vy%x-reshape(modes(2,1)%x,shape(vy%x)),n)+ &
       glsc3(vz%x-reshape(modes(3,1)%x,shape(vz%x)),fluid%c_Xh%B, &
            vz%x-reshape(modes(3,1)%x,shape(vz%x)),n))) / &
       max(velocity_norm,tiny(1.0_rp))
  gauge=sqrt(max(0.0_rp,glsc3(div_a%x,fluid%c_Xh%B,div_a%x,n))) / &
       max(velocity_norm,tiny(1.0_rp))
  boundary=0.0_rp
  do c=1,3
     if(.not.vp%poisson_bc(c)%dof_mask%is_set()) &
          call neko_error('Manufactured normal boundary mask is empty')
     ids=>vp%poisson_bc(c)%dof_mask%get()
     boundary=max(boundary,maxval(abs(vp%potential(c,1)%x(ids))))
  end do
  if(pe_rank==0) then
     write(*,'(A,3(1X,ES14.6))') 'Manufactured Poisson normalized true residuals:', &
          vp%solve_true_residual_normalized(:,1)
     write(*,'(A,3(1X,ES14.6))') 'Manufactured Poisson relative residuals:', &
          vp%solve_true_residual(:,1)
     write(*,'(A,ES14.6)') 'Native curl reconstruction relative error = ',reconstruction
     write(*,'(A,ES14.6)') 'Potential divergence norm / velocity norm = ',gauge
     write(*,'(A,ES14.6)') 'Potential normal-boundary max / potential norm = ', &
          boundary/max(potential_norm,tiny(1.0_rp))
     write(*,'(A,I0)') 'Poisson solves = ',vp%solve_count
  end if
  ! The solver uses an absolute residual normalized by sqrt(domain volume).
  ! Retain and report relative residuals as well; this manufactured field has
  ! an analytically zero third component, so its relative residual is ill-scaled.
  if(any(vp%solve_true_residual_normalized(:,1)>1.0e-8_rp)) &
       call neko_error('Manufactured Poisson residual check failed')
  call vp%free()
  do c=1,3
     call modes(c,1)%free()
  end do
  call ax%free(); call ay%free(); call az%free()
  call vx%free(); call vy%free(); call vz%free()
  call cx%free(); call cy%free(); call cz%free()
  call w1%free(); call w2%free(); call div_a%free()
  deallocate(xx,yy,zz)
  call simulation%free()
  call neko_finalize()

contains
  subroutine MPI_Comm_size_from_neko(size_out)
    use mpi_f08, only: MPI_Comm_size
    use comm, only: NEKO_COMM
    integer,intent(out)::size_out
    integer::mpierr
    call MPI_Comm_size(NEKO_COMM,size_out,mpierr)
  end subroutine MPI_Comm_size_from_neko
end program vector_potential_manufactured_diagnostic
