!> Native-GLL algebra for the supported reduced objectives.
!!
!! Pointwise data use vector_t and all pointwise products/contractions go
!! through vector_math.  The reduced quadratic forms use matrix_t/vector_t.
module rom_objectives
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t
  use vector_math, only: vector_col3, vector_add2s2, vector_vdot3, vector_cmult, &
       vector_cadd, &
       vector_glsc2, vector_glsc3
  use scratch_registry, only: neko_scratch_registry
  implicit none
  private

  public :: rom_assemble_viscous_quadratic
  public :: rom_assemble_viscous_terms
  public :: rom_assemble_scalar_mixing_terms
  public :: rom_assemble_brinkman_objective
  public :: rom_viscous_value
  public :: rom_brinkman_value
  public :: rom_brinkman_density

contains

  !> Assemble all terms of the spatial-lift viscous objective.
  !!
  !! For `u = lift + sum(a_j phi_j)`, the unnormalised objective is
  !! `nu * <grad(u),grad(u)>`.  The returned terms are the constant,
  !! lift-linear and modal quadratic contractions before the full-order
  !! factor one-half and objective-domain normalisation are applied.
  !! @param[in] weights SEM quadrature weights.
  !! @param[in] mask Objective-domain mask.
  !! @param[in] viscosity Dynamic viscosity.
  !! @param[in] gradients Mode gradients.
  !! @param[in] lift_gradients Spatial lift gradients.
  !! @param[out] constant Lift-only contraction.
  !! @param[out] linear Lift/mode contraction.
  !! @param[out] quadratic Mode contraction.
  subroutine rom_assemble_viscous_terms(weights, mask, viscosity, gradients, &
       lift_gradients, constant, linear, quadratic)
    type(vector_t), intent(in) :: weights, mask
    real(kind=rp), intent(in) :: viscosity
    type(vector_t), intent(in) :: gradients(:, :, :), lift_gradients(:, :)
    type(vector_t), intent(inout) :: constant, linear
    type(matrix_t), intent(inout) :: quadratic
    type(vector_t), pointer :: effective
    integer :: effective_index
    integer :: n, r, i, j, component, direction

    n = weights%size(); r = size(gradients,2)
    if (mask%size() /= n .or. size(gradients,1) /= 3 .or. &
         size(gradients,3) /= 3 .or. size(lift_gradients,1) /= 3 .or. &
         size(lift_gradients,2) /= 3) error stop 'ROM viscous lift shape mismatch'
    call neko_scratch_registry%request_vector(effective, effective_index, n, .true.)
    call vector_col3(effective, weights, mask)
    call constant%init(1); call linear%init(r); call quadratic%init(r,r)
    constant%x = 0.0_rp
    linear%x = 0.0_rp
    quadratic%x = 0.0_rp
    do component = 1, 3
       do direction = 1, 3
          constant%x(1) = constant%x(1) + viscosity * vector_glsc3( &
               lift_gradients(component,direction), effective, &
               lift_gradients(component,direction))
          do i = 1, r
             linear%x(i) = linear%x(i) + viscosity * vector_glsc3( &
                  lift_gradients(component,direction), effective, &
                  gradients(component,i,direction))
             do j = 1, r
                quadratic%x(i,j) = quadratic%x(i,j) + viscosity * &
                     vector_glsc3(gradients(component,i,direction), effective, &
                     gradients(component,j,direction))
             end do
          end do
       end do
    end do
    call neko_scratch_registry%relinquish_vector(effective_index)
    nullify(effective)
  end subroutine rom_assemble_viscous_terms

  !> Assemble the scalar-mixing objective for a spatial scalar lift.
  !! The full-order functional is
  !!  0.5 / |Omega| * <mask*(phi-phi_ref)^2, B>.
  !! The returned constant, linear and quadratic terms omit the factor
  !! one-half and domain normalisation, matching the velocity objective API.
  subroutine rom_assemble_scalar_mixing_terms(weights, mask, scalar_lift, &
       scalar_modes, scalar_reference_value, constant, linear, quadratic)
    type(vector_t), intent(in) :: weights, mask, scalar_lift
    type(vector_t), intent(in) :: scalar_modes(:)
    real(kind=rp), intent(in) :: scalar_reference_value
    type(vector_t), intent(inout) :: constant, linear
    type(matrix_t), intent(inout) :: quadratic
    type(vector_t), pointer :: effective, offset
    integer :: effective_index, offset_index
    integer :: n, r, i, j

    n = weights%size(); r = size(scalar_modes)
    if (mask%size() /= n .or. scalar_lift%size() /= n) &
         error stop 'ROM scalar objective shape mismatch'
    call neko_scratch_registry%request_vector(effective, effective_index, n, .true.)
    call neko_scratch_registry%request_vector(offset, offset_index, n, .true.)
    call vector_col3(effective, weights, mask)
    offset = scalar_lift
    call vector_cadd(offset, -scalar_reference_value)
    call constant%init(1); call linear%init(r); call quadratic%init(r,r)
    constant%x(1) = vector_glsc3(offset, effective, offset)
    linear%x = 0.0_rp
    quadratic%x = 0.0_rp
    do i = 1, r
       linear%x(i) = vector_glsc3(offset, effective, scalar_modes(i))
       do j = 1, r
          quadratic%x(i,j) = vector_glsc3(scalar_modes(i), effective, &
               scalar_modes(j))
       end do
    end do
    call neko_scratch_registry%relinquish_vector([effective_index, offset_index])
    nullify(effective); nullify(offset)
  end subroutine rom_assemble_scalar_mixing_terms

  !> Assemble the reduced viscous-dissipation quadratic form.
  subroutine rom_assemble_viscous_quadratic(weights, mask, viscosity, gradients, &
       quadratic)
    type(vector_t), intent(in) :: weights, mask
    real(kind=rp), intent(in) :: viscosity
    type(vector_t), intent(in) :: gradients(:, :, :)
    type(matrix_t), intent(inout) :: quadratic
    type(vector_t), pointer :: effective
    integer :: effective_index
    integer :: n, r, i, j, component, direction

    n = weights%size(); r = size(gradients,2)
    if (mask%size() /= n .or. size(gradients,1) /= 3 .or. &
         size(gradients,3) /= 3) error stop 'ROM viscous objective shape mismatch'
    call neko_scratch_registry%request_vector(effective, effective_index, n, .true.)
    call vector_col3(effective, weights, mask)
    call quadratic%init(r, r)
    do i = 1, r
       do j = 1, r
          do component = 1, 3
             do direction = 1, 3
                quadratic%x(i,j) = quadratic%x(i,j) + viscosity * &
                     vector_glsc3(gradients(component,i,direction), effective, &
                     gradients(component,j,direction))
             end do
          end do
       end do
    end do
    call neko_scratch_registry%relinquish_vector(effective_index); nullify(effective)
  end subroutine rom_assemble_viscous_quadratic

  !> Assemble Brinkman dissipation constant, lift-linear and quadratic terms.
  subroutine rom_assemble_brinkman_objective(weights, mask, chi, velocity_modes, &
       lift_velocity, quadratic, lift_linear, constant)
    type(vector_t), intent(in) :: weights, mask, chi
    type(vector_t), intent(in) :: velocity_modes(:, :), lift_velocity(:)
    type(matrix_t), intent(inout) :: quadratic
    type(vector_t), intent(inout) :: lift_linear
    real(kind=rp), intent(out) :: constant
    type(vector_t), pointer :: effective, masked_weights, lift_square
    integer :: effective_index, masked_weights_index, lift_square_index
    integer :: n, r, i, j, component

    n = weights%size(); r = size(velocity_modes,2)
    if (mask%size() /= n .or. chi%size() /= n .or. size(velocity_modes,1) /= 3 .or. &
         size(lift_velocity) /= 3) error stop 'ROM Brinkman objective shape mismatch'
    call neko_scratch_registry%request_vector(masked_weights, masked_weights_index, n, .true.)
    call neko_scratch_registry%request_vector(effective, effective_index, n, .true.)
    call neko_scratch_registry%request_vector(lift_square, lift_square_index, n, .true.)
    call vector_col3(masked_weights, weights, mask)
    call vector_col3(effective, masked_weights, chi)
    call vector_vdot3(lift_square, lift_velocity(1), lift_velocity(2), lift_velocity(3), &
         lift_velocity(1), lift_velocity(2), lift_velocity(3))
    constant = vector_glsc2(effective, lift_square)
    call quadratic%init(r, r); call lift_linear%init(r)
    quadratic%x = 0.0_rp
    lift_linear%x = 0.0_rp
    do i = 1, r
       do j = 1, r
          do component = 1, 3
             quadratic%x(i,j) = quadratic%x(i,j) + vector_glsc3( &
                  velocity_modes(component,i), effective, velocity_modes(component,j))
          end do
       end do
       do component = 1, 3
          lift_linear%x(i) = lift_linear%x(i) + vector_glsc3( &
               velocity_modes(component,i), effective, lift_velocity(component))
       end do
    end do
    call neko_scratch_registry%relinquish_vector([masked_weights_index, effective_index, lift_square_index])
    nullify(masked_weights); nullify(effective); nullify(lift_square)
  end subroutine rom_assemble_brinkman_objective

  !> Evaluate the reduced viscous-dissipation quadratic objective.
  function rom_viscous_value(state, constant, linear, quadratic, domain) result(value)
    type(vector_t), intent(in) :: state
    real(kind=rp), intent(in) :: constant
    type(vector_t), intent(in) :: linear
    type(matrix_t), intent(in) :: quadratic
    real(kind=rp), intent(in) :: domain
    real(kind=rp) :: value
    if (linear%size() /= state%size() .or. quadratic%get_nrows() /= state%size() .or. &
         quadratic%get_ncols() /= state%size()) error stop 'ROM viscous objective value shape mismatch'
    value = 0.5_rp * (constant + 2.0_rp * dot_product(linear%x, state%x) + &
         dot_product(state%x, matmul(quadratic%x, state%x))) / &
         max(domain, tiny(1.0_rp))
  end function rom_viscous_value

  !> Evaluate the reduced Brinkman-dissipation quadratic objective.
  function rom_brinkman_value(state, constant, lift_linear, quadratic, domain) result(value)
    type(vector_t), intent(in) :: state, lift_linear
    real(kind=rp), intent(in) :: constant, domain
    type(matrix_t), intent(in) :: quadratic
    real(kind=rp) :: value
    if (lift_linear%size() /= state%size() .or. quadratic%get_nrows() /= state%size() .or. &
         quadratic%get_ncols() /= state%size()) error stop 'ROM Brinkman value shape mismatch'
    value = 0.5_rp * (constant + 2.0_rp * dot_product(lift_linear%x, state%x) + &
         dot_product(state%x, matmul(quadratic%x, state%x))) / max(domain, tiny(1.0_rp))
  end function rom_brinkman_value

  !> Assemble the direct Brinkman objective derivative with respect to design.
  subroutine rom_brinkman_density(weights, mask, lift_velocity, velocity_modes, &
       state, domain, density)
    type(vector_t), intent(in) :: weights, mask, lift_velocity(:), state
    type(vector_t), intent(in) :: velocity_modes(:, :)
    real(kind=rp), intent(in) :: domain
    type(vector_t), intent(inout) :: density
    type(vector_t), pointer :: velocity_u, velocity_v, velocity_w
    type(vector_t), pointer :: speed_square, masked_speed
    integer :: velocity_index(3), speed_square_index, masked_speed_index
    integer :: n, r, component, i

    n = weights%size(); r = size(velocity_modes,2)
    if (mask%size() /= n .or. size(lift_velocity) /= 3 .or. state%size() /= r .or. &
         size(velocity_modes,1) /= 3) error stop 'ROM Brinkman density shape mismatch'
    do component = 1, 3
       select case (component)
       case (1)
          call neko_scratch_registry%request_vector(velocity_u, velocity_index(component), n, .true.)
          velocity_u = lift_velocity(component)
       case (2)
          call neko_scratch_registry%request_vector(velocity_v, velocity_index(component), n, .true.)
          velocity_v = lift_velocity(component)
       case (3)
          call neko_scratch_registry%request_vector(velocity_w, velocity_index(component), n, .true.)
          velocity_w = lift_velocity(component)
       end select
       do i = 1, r
          select case (component)
          case (1)
             call vector_add2s2(velocity_u, velocity_modes(component,i), state%x(i))
          case (2)
             call vector_add2s2(velocity_v, velocity_modes(component,i), state%x(i))
          case (3)
             call vector_add2s2(velocity_w, velocity_modes(component,i), state%x(i))
          end select
       end do
    end do
    call neko_scratch_registry%request_vector(speed_square, speed_square_index, n, .true.)
    call neko_scratch_registry%request_vector(masked_speed, masked_speed_index, n, .true.)
    call density%init(n)
    call vector_vdot3(speed_square, velocity_u, velocity_v, velocity_w, &
         velocity_u, velocity_v, velocity_w)
    call vector_col3(masked_speed, mask, speed_square)
    call vector_cmult(masked_speed, 0.5_rp / max(domain, tiny(1.0_rp)))
    density = masked_speed
    do component = 1, 3
       call neko_scratch_registry%relinquish_vector(velocity_index(component))
    end do
    call neko_scratch_registry%relinquish_vector([speed_square_index, masked_speed_index])
    nullify(velocity_u); nullify(velocity_v); nullify(velocity_w)
    nullify(speed_square); nullify(masked_speed)
  end subroutine rom_brinkman_density

end module rom_objectives
