program rom_reduced_system_export_diagnostic
  use neko, only: neko_init, neko_finalize
  use num_types, only: rp
  use matrix, only: matrix_t
  use vector, only: vector_t
  use rom_model, only: rom_model_t, rom_bdf1_export_context_t
  use utils, only: neko_error
  implicit none

  integer, parameter :: n = 2, n_steps = 6
  real(kind=rp), parameter :: dt = 0.05_rp
  type(rom_model_t) :: model
  type(rom_bdf1_export_context_t) :: export_context
  type(matrix_t) :: mass, linear, implicit_linear, states
  type(matrix_t), allocatable :: quadratic(:)
  type(vector_t) :: constant, initial_state
  character(len=1024) :: output_directory
  integer :: i

  if (command_argument_count() /= 1) &
       call neko_error('Usage: reduced_system_export_diagnostic OUTPUT_DIRECTORY')
  call get_command_argument(1, output_directory)
  call neko_init()

  allocate(quadratic(n))
  call mass%init(n,n); call linear%init(n,n); call implicit_linear%init(n,n)
  call constant%init(n); call initial_state%init(n)
  do i = 1, n
     call quadratic(i)%init(n,n)
     quadratic(i)%x = 0.0_rp
  end do
  mass%x = reshape([1.2_rp, 0.05_rp, 0.05_rp, 0.9_rp], [n,n])
  linear%x = reshape([-0.12_rp, 0.03_rp, -0.02_rp, -0.08_rp], [n,n])
  implicit_linear%x = reshape([-0.30_rp, 0.02_rp, 0.01_rp, -0.20_rp], [n,n])
  quadratic(1)%x(1,1) = 0.04_rp
  quadratic(2)%x(1,2) = -0.015_rp
  constant%x = [0.20_rp, -0.12_rp]
  initial_state%x = [0.12_rp, -0.21_rp]

  export_context%directory = trim(output_directory)
  export_context%design_id = 'export_diagnostic'
  export_context%basis_forward = 'forward'
  export_context%basis_adjoint = 'enriched'
  export_context%outer_iteration = 2
  export_context%inner_iteration = 4
  export_context%evaluation_index = 7
  export_context%trial_basis_rank = n
  export_context%test_basis_rank = n
  export_context%forward_basis_rank = n
  export_context%adjoint_basis_rank = n + 1
  export_context%enriched_basis_rank = n + 1
  export_context%time_start = 0.1_rp
  call model%init(mass, constant, linear, quadratic, implicit_linear)
  call states%init(n, n_steps + 1)
  call model%forward_bdf1(initial_state, dt, states, export_context)
  write(*,'(A,2(1X,ES16.8))') 'export diagnostic final state:', states%x(:,n_steps+1)

  call model%free(); call states%free()
  call mass%free(); call linear%free(); call implicit_linear%free()
  call constant%free(); call initial_state%free()
  do i = 1, n
     call quadratic(i)%free()
  end do
  deallocate(quadratic)
  call neko_finalize()
end program rom_reduced_system_export_diagnostic
