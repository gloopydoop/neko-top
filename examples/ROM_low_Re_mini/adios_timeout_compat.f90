! Adapt the affected Neko Fortran-to-C++ SST timeout argument.  The six
! preceding arguments are already pointers on both sides of the interface.
subroutine wrap_adios_timeout(npts, nelv, nelb, nelgv, gdim, comm, timeout) &
     bind(C, name="__wrap_adios2_initialize_")
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_f_pointer
  implicit none
  type(c_ptr), value :: npts, nelv, nelb, nelgv, gdim, comm, timeout
  integer(c_int), pointer :: timeout_value
  interface
     subroutine real_adios_init(npts, nelv, nelb, nelgv, gdim, comm, timeout) &
          bind(C, name="__real_adios2_initialize_")
       import c_ptr, c_int
       type(c_ptr), value :: npts, nelv, nelb, nelgv, gdim, comm
       integer(c_int), value :: timeout
     end subroutine real_adios_init
  end interface
  call c_f_pointer(timeout, timeout_value)
  call real_adios_init(npts, nelv, nelb, nelgv, gdim, comm, timeout_value)
end subroutine wrap_adios_timeout
