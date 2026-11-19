!> Value and gradient agreement used by the fixed-basis ROM trust gate.
!!
!! Keeping the comparison decision in a small ROM module makes the policy
!! independently testable without constructing a CFD problem.
module rom_validation
  use num_types, only: rp
  implicit none
  private

  type, public :: rom_validation_metrics_t
     real(kind=rp) :: objective = huge(0.0_rp)
     real(kind=rp) :: constraint = huge(0.0_rp)
     real(kind=rp) :: objective_gradient = huge(0.0_rp)
     real(kind=rp) :: constraint_gradient = huge(0.0_rp)
  end type rom_validation_metrics_t

  public :: rom_validation_accept

contains

  !> Decide whether objective and constraint trust metrics pass their gates.
  pure logical function rom_validation_accept(metrics, value_tolerance, &
       gradient_tolerance) result(accepted)
    type(rom_validation_metrics_t), intent(in) :: metrics
    real(kind=rp), intent(in) :: value_tolerance, gradient_tolerance

    accepted = metrics%objective <= value_tolerance .and. &
         metrics%constraint <= value_tolerance .and. &
         metrics%objective_gradient <= gradient_tolerance .and. &
         metrics%constraint_gradient <= gradient_tolerance
  end function rom_validation_accept

end module rom_validation
