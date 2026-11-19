program rom_validation_diagnostic
  use num_types, only: rp
  use rom_validation, only: rom_validation_metrics_t, rom_validation_accept
  implicit none

  type(rom_validation_metrics_t) :: metrics
  real(kind=rp), parameter :: value_tolerance = 5.0e-2_rp
  real(kind=rp), parameter :: gradient_tolerance = 5.0e-2_rp

  metrics%objective = 1.0e-3_rp
  metrics%constraint = 2.0e-3_rp
  metrics%objective_gradient = 3.0e-3_rp
  metrics%constraint_gradient = 4.0e-3_rp
  if (.not. rom_validation_accept(metrics, value_tolerance, gradient_tolerance)) &
       error stop 'trust acceptance diagnostic failed'

  metrics%objective_gradient = 6.0e-2_rp
  if (rom_validation_accept(metrics, value_tolerance, gradient_tolerance)) &
       error stop 'trust rejection diagnostic failed'

  write(*,'(A)') 'ROM value/gradient trust acceptance-rejection diagnostic passed'
end program rom_validation_diagnostic
