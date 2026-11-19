!> Mass inner products and fail-closed orthonormalisation for ROM columns.
!!
!! This module deliberately has no simulation or fluid dependency.  The
!! caller supplies the already assembled SEM velocity mass diagonal, so the
!! same algebra can be tested on small manufactured column sets and used by
!! the production ROM handoff without embedding basis construction in the
!! generic simulation implementation.
module rom_basis_orthogonalization
  use math, only: glsc3
  use num_types, only: rp
  use vector, only: vector_t
  use scratch_registry, only: neko_scratch_registry
  implicit none
  private

  public :: rom_mass_dot
  public :: rom_mass_orthonormalize

contains

  !> Compute the three-component SEM mass inner product of two columns.
  function rom_mass_dot(a, b, mass, n) result(value)
    real(kind=rp), intent(in) :: a(:), b(:), mass(:)
    integer, intent(in) :: n
    real(kind=rp) :: value
    integer :: component

    if (n < 1 .or. size(a) < 3 * n .or. size(b) < 3 * n .or. &
         size(mass) < n) then
       value = 0.0_rp
       return
    end if

    value = 0.0_rp
    do component = 1, 3
       value = value + glsc3(a((component - 1) * n + 1:component * n), &
            mass(1:n), b((component - 1) * n + 1:component * n), n)
    end do
  end function rom_mass_dot

  !> Modified Gram-Schmidt with a second pass and an explicit SEM mass.
  !!
  !! `projected` is compacted in place when a column is numerically dependent.
  !! `final_mode_index` maps each input column to its accepted column (zero
  !! means rejected).  `ok` is false for malformed input, a non-positive
  !! norm, or an orthogonality failure; callers should then fail closed.
  subroutine rom_mass_orthonormalize(projected, n, mass, tolerance, accepted, &
       final_mode_index, mass_norm_final, mass_cross_final, ok)
    real(kind=rp), intent(inout) :: projected(:, :)
    real(kind=rp), intent(in) :: mass(:)
    real(kind=rp), intent(in) :: tolerance
    integer, intent(in) :: n
    integer, intent(out) :: accepted
    integer, intent(out) :: final_mode_index(:)
    real(kind=rp), intent(out) :: mass_norm_final(:), mass_cross_final(:)
    logical, intent(out) :: ok
    type(vector_t), pointer :: vector
    integer :: vector_index
    real(kind=rp) :: dot_value, norm_squared
    integer :: mode, previous, n_modes

    n_modes = size(projected, 2)
    accepted = 0
    ok = .false.
    final_mode_index = 0
    mass_norm_final = 0.0_rp
    mass_cross_final = 0.0_rp
    if (n < 1 .or. size(projected, 1) /= 3 * n .or. &
         size(mass) < n .or. size(final_mode_index) /= n_modes .or. &
         size(mass_norm_final) /= n_modes .or. &
         size(mass_cross_final) /= n_modes .or. tolerance <= 0.0_rp) return

    call neko_scratch_registry%request_vector(vector, vector_index, 3 * n, .true.)
    do mode = 1, n_modes
       vector%x = projected(:, mode)

       ! Two MGS passes are intentional: the projection can make columns
       ! nearly dependent, and a single pass is not a sufficient numerical
       ! test of the final common basis.
       do previous = 1, accepted
          dot_value = rom_mass_dot(projected(:, previous), vector%x, mass, n)
          vector%x = vector%x - dot_value * projected(:, previous)
       end do
       do previous = 1, accepted
          dot_value = rom_mass_dot(projected(:, previous), vector%x, mass, n)
          vector%x = vector%x - dot_value * projected(:, previous)
       end do

       norm_squared = rom_mass_dot(vector%x, vector%x, mass, n)
       if (norm_squared <= max(1.0e-24_rp, tolerance * tolerance)) cycle
       accepted = accepted + 1
       projected(:, accepted) = vector%x / sqrt(norm_squared)
       final_mode_index(mode) = accepted
       mass_norm_final(mode) = sqrt(max(0.0_rp, rom_mass_dot( &
            projected(:, accepted), projected(:, accepted), mass, n)))
       do previous = 1, accepted - 1
          mass_cross_final(mode) = max(mass_cross_final(mode), abs( &
               rom_mass_dot(projected(:, accepted), projected(:, previous), &
               mass, n)))
       end do
    end do

    if (accepted < 1) then
       call neko_scratch_registry%relinquish_vector(vector_index)
       nullify(vector)
       return
    end if
    do mode = 1, n_modes
       if (final_mode_index(mode) > 0) then
          if (abs(mass_norm_final(mode) - 1.0_rp) > tolerance .or. &
               mass_cross_final(mode) > tolerance) then
             call neko_scratch_registry%relinquish_vector(vector_index)
             nullify(vector)
             return
          end if
       end if
    end do
    ok = .true.
    call neko_scratch_registry%relinquish_vector(vector_index)
    nullify(vector)
  end subroutine rom_mass_orthonormalize

end module rom_basis_orthogonalization
