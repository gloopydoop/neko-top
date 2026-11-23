!> SEM-weighted reduced Galerkin assembly.
!!
!! Pointwise data are carried by Neko vector_t containers and reduced
!! operators by matrix_t/vector_t.  The routines use vector_math for all
!! pointwise products and weighted contractions, so the same host/device path
!! is used as in the rest of Neko-TOP.
module rom_assembly
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t
  use vector_math, only: vector_col3, vector_add2, vector_cmult, vector_cmult2, &
       vector_glsc2, vector_glsc3
  use scratch_registry, only: neko_scratch_registry
  implicit none
  private

  public :: rom_assemble_mass
  public :: rom_assemble_mass_coupled
  public :: rom_assemble_cross_mass
  public :: rom_assemble_brinkman
  public :: rom_assemble_brinkman_petrov
  public :: rom_assemble_lift_convection
  public :: rom_assemble_lift_convection_petrov
  public :: rom_assemble_diffusion
  public :: rom_assemble_diffusion_petrov
  public :: rom_assemble_convection
  public :: rom_assemble_convection_petrov
  public :: rom_add_linearized_convection
  public :: rom_add_linearized_scalar_advection
  public :: rom_assemble_scalar_transport
  public :: rom_assemble_scalar_transport_petrov

contains

  !> Assemble the SEM mass matrix in the reduced velocity basis.
  subroutine rom_assemble_mass(weights, modes, mass)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: modes(:, :)
    type(matrix_t), intent(inout) :: mass
    integer :: n, n_fields, n_modes, i, j, field

    n = weights%size(); n_fields = size(modes, 1); n_modes = size(modes, 2)
    if (n < 1 .or. n_fields < 1 .or. any([(modes(field,1)%size(), field=1,n_fields)] /= n)) &
         error stop 'ROM mass assembly shape mismatch'
    call rom_assemble_cross_mass(weights, modes, modes, mass)
  end subroutine rom_assemble_mass

  !> Assemble the Petrov--Galerkin cross mass matrix.
  !! The first mode argument is the trial (right) basis and the second is the
  !! test (left) basis.  Optional scalar modes use the same coupled-state
  !! metric as the existing Galerkin assembly.
  subroutine rom_assemble_cross_mass(weights, trial_modes, test_modes, mass, &
       scalar_trial_modes, scalar_test_modes, scalar_weight)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: trial_modes(:, :), test_modes(:, :)
    type(matrix_t), intent(inout) :: mass
    type(vector_t), intent(in), optional :: scalar_trial_modes(:), scalar_test_modes(:)
    real(kind=rp), intent(in), optional :: scalar_weight
    integer :: n, r_trial, r_test, i, j, field
    real(kind=rp) :: weight

    n = weights%size(); r_trial = size(trial_modes, 2); r_test = size(test_modes, 2)
    if (r_trial < 1 .or. r_trial /= r_test .or. size(trial_modes, 1) /= size(test_modes, 1)) &
         error stop 'ROM cross mass basis shape mismatch'
    do field = 1, size(trial_modes, 1)
       if (trial_modes(field,1)%size() /= n .or. test_modes(field,1)%size() /= n) &
            error stop 'ROM cross mass vector size mismatch'
    end do
    call mass%init(r_test, r_trial)
    mass%x = 0.0_rp
    do i = 1, r_test
       do j = 1, r_trial
          do field = 1, size(trial_modes, 1)
             mass%x(i,j) = mass%x(i,j) + vector_glsc3( &
                  test_modes(field,i), weights, trial_modes(field,j))
          end do
       end do
    end do
    if (present(scalar_trial_modes) .or. present(scalar_test_modes)) then
       if (.not. present(scalar_trial_modes) .or. .not. present(scalar_test_modes) .or. &
            size(scalar_trial_modes) /= r_trial .or. size(scalar_test_modes) /= r_test) &
            error stop 'ROM cross scalar mass shape mismatch'
       weight = 1.0_rp
       if (present(scalar_weight)) weight = scalar_weight
       do i = 1, r_test
          do j = 1, r_trial
             mass%x(i,j) = mass%x(i,j) + weight * vector_glsc3( &
                  scalar_test_modes(i), weights, scalar_trial_modes(j))
          end do
       end do
    end if
  end subroutine rom_assemble_cross_mass

  !> Assemble the coupled velocity/scalar SEM mass matrix.
  !! The scalar multiplier is the declared state metric.  It is applied to
  !! both the scalar mass and scalar residual contractions, so the reduced
  !! equations remain the Galerkin projection of the block-diagonal physical
  !! mass operator.
  subroutine rom_assemble_mass_coupled(weights, velocity_modes, scalar_modes, &
       scalar_weight, mass)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: velocity_modes(:, :)
    type(vector_t), intent(in), optional :: scalar_modes(:)
    real(kind=rp), intent(in), optional :: scalar_weight
    type(matrix_t), intent(inout) :: mass
    integer :: r

    r = size(velocity_modes, 2)
    call rom_assemble_cross_mass(weights, velocity_modes, velocity_modes, mass, &
         scalar_modes, scalar_modes, scalar_weight)
    if (.not. present(scalar_modes)) return
    if (size(scalar_modes) /= r) error stop 'ROM scalar mass shape mismatch'
  end subroutine rom_assemble_mass_coupled

  !> Assemble the reduced Brinkman operator and its lift contribution.
  subroutine rom_assemble_brinkman(weights, chi, velocity_modes, lift_velocity, &
       brinkman, lift_vector)
    type(vector_t), intent(in) :: weights, chi
    type(vector_t), intent(in) :: velocity_modes(:, :), lift_velocity(:)
    type(matrix_t), intent(inout) :: brinkman
    type(vector_t), intent(inout) :: lift_vector
    call rom_assemble_brinkman_petrov(weights, chi, velocity_modes, velocity_modes, &
         lift_velocity, brinkman, lift_vector)
  end subroutine rom_assemble_brinkman

  !> Assemble the Petrov--Galerkin Brinkman operator and lift contribution.
  !! Test modes are used only on the left; trial modes are used on the right.
  subroutine rom_assemble_brinkman_petrov(weights, chi, trial_velocity_modes, &
       test_velocity_modes, lift_velocity, brinkman, lift_vector)
    type(vector_t), intent(in) :: weights, chi
    type(vector_t), intent(in) :: trial_velocity_modes(:, :), test_velocity_modes(:, :)
    type(vector_t), intent(in) :: lift_velocity(:)
    type(matrix_t), intent(inout) :: brinkman
    type(vector_t), intent(inout) :: lift_vector
    type(vector_t), pointer :: effective
    integer :: effective_index
    integer :: n, n_trial, n_test, i, j, component

    n = weights%size(); n_trial = size(trial_velocity_modes, 2); n_test = size(test_velocity_modes, 2)
    if (chi%size() /= n .or. size(trial_velocity_modes, 1) /= 3 .or. &
         size(test_velocity_modes, 1) /= 3 .or. size(lift_velocity) /= 3 .or. &
         any([(trial_velocity_modes(component,1)%size(), component=1,3)] /= n) .or. &
         any([(test_velocity_modes(component,1)%size(), component=1,3)] /= n)) &
         error stop 'ROM Brinkman assembly shape mismatch'
    call neko_scratch_registry%request_vector(effective, effective_index, n, .true.)
    call vector_col3(effective, weights, chi)
    call brinkman%init(n_test, n_trial); call lift_vector%init(n_test)
    ! The contractions return host scalars and these reduced outputs are
    ! deliberately host-owned.
    brinkman%x = 0.0_rp
    lift_vector%x = 0.0_rp
    do i = 1, n_test
       do j = 1, n_trial
          do component = 1, 3
             brinkman%x(i,j) = brinkman%x(i,j) + vector_glsc3( &
                  test_velocity_modes(component,i), effective, trial_velocity_modes(component,j))
          end do
       end do
       do component = 1, 3
          lift_vector%x(i) = lift_vector%x(i) + vector_glsc3( &
               test_velocity_modes(component,i), effective, lift_velocity(component))
       end do
    end do
    call neko_scratch_registry%relinquish_vector(effective_index); nullify(effective)
  end subroutine rom_assemble_brinkman_petrov

  !> Assemble the convection generated by a spatial physical lift.
  !!
  !! For `u = lift + sum(a_j phi_j)`, this routine returns the negative
  !! weak residual contributions from `lift.grad(phi_j)`,
  !! `phi_j.grad(lift)` and `lift.grad(lift)`.  The latter is returned in the
  !! constant vector, while the two cross terms are returned in the matrix.
  !! @param[in] weights SEM quadrature weights.
  !! @param[in] velocity_modes Homogeneous velocity modes.
  !! @param[in] gradients Mode gradients, indexed by component, mode, direction.
  !! @param[in] lift_velocity Spatial lift components.
  !! @param[in] lift_gradients Spatial lift gradients, indexed by component and direction.
  !! @param[out] linear Reduced cross-convection matrix.
  !! @param[out] constant Reduced lift-convection vector.
  subroutine rom_assemble_lift_convection(weights, velocity_modes, gradients, &
       lift_velocity, lift_gradients, linear, constant)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: velocity_modes(:, :), gradients(:, :, :)
    type(vector_t), intent(in) :: lift_velocity(:), lift_gradients(:, :)
    type(matrix_t), intent(inout) :: linear
    type(vector_t), intent(inout) :: constant
    call rom_assemble_lift_convection_petrov(weights, velocity_modes, gradients, &
         velocity_modes, lift_velocity, lift_gradients, linear, constant)
  end subroutine rom_assemble_lift_convection

  !> Assemble lift convection with independent trial and test spaces.
  subroutine rom_assemble_lift_convection_petrov(weights, trial_velocity_modes, &
       trial_gradients, test_velocity_modes, lift_velocity, lift_gradients, &
       linear, constant)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: trial_velocity_modes(:, :), trial_gradients(:, :, :)
    type(vector_t), intent(in) :: test_velocity_modes(:, :)
    type(vector_t), intent(in) :: lift_velocity(:), lift_gradients(:, :)
    type(matrix_t), intent(inout) :: linear
    type(vector_t), intent(inout) :: constant
    type(vector_t), pointer :: pointwise, product
    integer :: pointwise_index, product_index
    integer :: n, r, i, j, component, direction

    n = weights%size(); r = size(trial_velocity_modes, 2)
    if (size(trial_velocity_modes,1) /= 3 .or. size(test_velocity_modes,1) /= 3 .or. &
         size(test_velocity_modes,2) /= r .or. size(trial_gradients,1) /= 3 .or. &
         size(trial_gradients,2) /= r .or. size(trial_gradients,3) /= 3 .or. &
         size(lift_velocity) /= 3 .or. size(lift_gradients,1) /= 3 .or. &
         size(lift_gradients,2) /= 3) error stop 'ROM lift convection shape mismatch'
    call linear%init(r, r)
    call constant%init(r)
    linear%x = 0.0_rp
    constant%x = 0.0_rp
    call neko_scratch_registry%request_vector(pointwise, pointwise_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    do i = 1, r
       do j = 1, r
          do component = 1, 3
             pointwise = 0.0_rp
             do direction = 1, 3
                call vector_col3(product, lift_velocity(direction), &
                     trial_gradients(component,j,direction))
                call vector_add2(pointwise, product)
                call vector_col3(product, trial_velocity_modes(direction,j), &
                     lift_gradients(component,direction))
                call vector_add2(pointwise, product)
             end do
             linear%x(i,j) = linear%x(i,j) - vector_glsc3( &
                  test_velocity_modes(component,i), weights, pointwise)
          end do
       end do
    end do
    do i = 1, r
       do component = 1, 3
          pointwise = 0.0_rp
          do direction = 1, 3
             call vector_col3(product, lift_velocity(direction), &
                  lift_gradients(component,direction))
             call vector_add2(pointwise, product)
          end do
          constant%x(i) = constant%x(i) - vector_glsc3( &
               test_velocity_modes(component,i), weights, pointwise)
       end do
    end do
    call neko_scratch_registry%relinquish_vector([pointwise_index, product_index])
    nullify(pointwise); nullify(product)
  end subroutine rom_assemble_lift_convection_petrov

  !> Assemble the implicit viscous residual matrix and lift diffusion vector.
  !!
  !! The matrix is `-nu <grad(test),grad(state)>`.  When lift gradients and a
  !! constant vector are supplied, the routine additionally accumulates
  !! `-nu <grad(test),grad(lift)>` in that vector.
  subroutine rom_assemble_diffusion(weights, viscosity_over_density, gradients, &
       implicit_linear, lift_gradients, lift_constant)
    type(vector_t), intent(in) :: weights
    real(kind=rp), intent(in) :: viscosity_over_density
    type(vector_t), intent(in) :: gradients(:, :, :)
    type(matrix_t), intent(inout) :: implicit_linear
    type(vector_t), intent(in), optional :: lift_gradients(:, :)
    type(vector_t), intent(inout), optional :: lift_constant
    call rom_assemble_diffusion_petrov(weights, viscosity_over_density, gradients, &
         gradients, implicit_linear, lift_gradients, lift_constant)
  end subroutine rom_assemble_diffusion

  !> Assemble diffusion with independent test and trial gradients.
  subroutine rom_assemble_diffusion_petrov(weights, viscosity_over_density, &
       test_gradients, trial_gradients, implicit_linear, lift_gradients, lift_constant)
    type(vector_t), intent(in) :: weights
    real(kind=rp), intent(in) :: viscosity_over_density
    type(vector_t), intent(in) :: test_gradients(:, :, :), trial_gradients(:, :, :)
    type(matrix_t), intent(inout) :: implicit_linear
    type(vector_t), intent(in), optional :: lift_gradients(:, :)
    type(vector_t), intent(inout), optional :: lift_constant
    integer :: n, r, i, j, component, direction

    n = weights%size(); r = size(trial_gradients,2)
    if (size(test_gradients,1) /= 3 .or. size(test_gradients,2) /= r .or. &
         size(test_gradients,3) /= 3 .or. size(trial_gradients,1) /= 3 .or. &
         size(trial_gradients,3) /= 3) &
         error stop 'ROM diffusion assembly shape mismatch'
    if (present(lift_gradients) .neqv. present(lift_constant)) &
         error stop 'ROM diffusion lift arguments must be supplied together'
    if (present(lift_gradients)) then
       if (size(lift_gradients,1) /= 3 .or. size(lift_gradients,2) /= 3) &
            error stop 'ROM diffusion lift shape mismatch'
    end if
    call implicit_linear%init(r, r)
    implicit_linear%x = 0.0_rp
    do i = 1, r
       if (present(lift_gradients)) then
          do component = 1, 3
             do direction = 1, 3
                lift_constant%x(i) = lift_constant%x(i) - viscosity_over_density * &
                     vector_glsc3(test_gradients(component,i,direction), weights, &
                     lift_gradients(component,direction))
             end do
          end do
       end if
       do j = 1, r
          do component = 1, 3
             do direction = 1, 3
                implicit_linear%x(i,j) = implicit_linear%x(i,j) - &
                     viscosity_over_density * vector_glsc3( &
                     test_gradients(component,i,direction), weights, &
                     trial_gradients(component,j,direction))
             end do
          end do
       end do
    end do
  end subroutine rom_assemble_diffusion_petrov

  !> Assemble the quadratic convective residual tensor.
  subroutine rom_assemble_convection(weights, modes, gradients, quadratic)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: modes(:, :), gradients(:, :, :)
    type(matrix_t), intent(inout) :: quadratic(:)
    call rom_assemble_convection_petrov(weights, modes, modes, gradients, quadratic)
  end subroutine rom_assemble_convection

  !> Assemble the Petrov--Galerkin quadratic velocity convection tensor.
  !! The advecting and transported fields are trial modes; only the test
  !! function is taken from the independent test basis.
  subroutine rom_assemble_convection_petrov(weights, test_modes, trial_modes, &
       trial_gradients, quadratic)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: test_modes(:, :), trial_modes(:, :)
    type(vector_t), intent(in) :: trial_gradients(:, :, :)
    type(matrix_t), intent(inout) :: quadratic(:)
    type(vector_t), pointer :: pointwise, product
    integer :: pointwise_index, product_index
    integer :: n, r, test_mode, advect_mode, state_mode, field, component

    n = weights%size(); r = size(trial_modes,2)
    if (size(test_modes,1) < 3 .or. size(test_modes,2) /= r .or. &
         size(trial_modes,1) < 3 .or. size(trial_gradients,1) < 3 .or. &
         size(trial_gradients,2) /= r .or. size(trial_gradients,3) /= 3 .or. &
         size(quadratic) /= r) error stop 'ROM convection assembly shape mismatch'
    do field = 1, r
       call quadratic(field)%init(r, r)
       quadratic(field)%x = 0.0_rp
    end do
    call neko_scratch_registry%request_vector(pointwise, pointwise_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    ! The velocity-only ROM has three transported fields.  Additional scalar
    ! fields are intentionally unsupported by rom_problem's configuration gate.
    do test_mode = 1, r
       do advect_mode = 1, r
          do state_mode = 1, r
             do field = 1, 3
                pointwise = 0.0_rp
                do component = 1, 3
                   call vector_col3(product, trial_modes(component,advect_mode), &
                        trial_gradients(field,state_mode,component))
                   call vector_add2(pointwise, product)
                end do
                quadratic(test_mode)%x(advect_mode,state_mode) = &
                     quadratic(test_mode)%x(advect_mode,state_mode) + &
                     vector_glsc3(test_modes(field,test_mode), weights, pointwise)
             end do
          end do
       end do
    end do
    call neko_scratch_registry%relinquish_vector([pointwise_index, product_index])
    nullify(pointwise); nullify(product)
  end subroutine rom_assemble_convection_petrov

  !> Add the velocity-convection Jacobian about a forward reduced state.
  !!
  !! The forward state is represented by `coefficients` in the forward trial
  !! basis.  The perturbation and test spaces may be enriched and therefore
  !! have an independent rank.  The two contributions are
  !! `-(delta_u.grad(u_f) + u_f.grad(delta_u))`.
  subroutine rom_add_linearized_convection(weights, forward_modes, forward_gradients, &
       test_modes, perturb_modes, perturb_gradients, coefficients, linearized)
    type(vector_t), intent(in) :: weights
    type(vector_t), intent(in) :: forward_modes(:, :), forward_gradients(:, :, :)
    type(vector_t), intent(in) :: test_modes(:, :), perturb_modes(:, :)
    type(vector_t), intent(in) :: perturb_gradients(:, :, :)
    real(kind=rp), intent(in) :: coefficients(:)
    type(matrix_t), intent(inout) :: linearized
    type(vector_t), pointer :: pointwise, product
    integer :: pointwise_index, product_index
    integer :: n, rf, re, i, j, k, field, direction

    n = weights%size(); rf = size(forward_modes,2); re = size(perturb_modes,2)
    if (size(forward_modes,1) /= 3 .or. size(test_modes,1) /= 3 .or. &
         size(perturb_modes,1) /= 3 .or. size(coefficients) /= rf .or. &
         size(forward_gradients,1) /= 3 .or. size(forward_gradients,2) /= rf .or. &
         size(forward_gradients,3) /= 3 .or. size(perturb_gradients,1) /= 3 .or. &
         size(perturb_gradients,2) /= re .or. size(perturb_gradients,3) /= 3 .or. &
         size(test_modes,2) /= re .or. linearized%get_nrows() /= re .or. &
         linearized%get_ncols() /= re) error stop 'ROM linearized convection shape mismatch'
    do field = 1, 3
       if (forward_modes(field,1)%size() /= n .or. test_modes(field,1)%size() /= n .or. &
            perturb_modes(field,1)%size() /= n) error stop 'ROM linearized convection vector mismatch'
    end do

    call neko_scratch_registry%request_vector(pointwise, pointwise_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    do i = 1, re
       do j = 1, re
          do k = 1, rf
             do field = 1, 3
                pointwise = 0.0_rp
                do direction = 1, 3
                   call vector_col3(product, perturb_modes(direction,j), &
                        forward_gradients(field,k,direction))
                   call vector_add2(pointwise, product)
                   call vector_col3(product, forward_modes(direction,k), &
                        perturb_gradients(field,j,direction))
                   call vector_add2(pointwise, product)
                end do
                linearized%x(i,j) = linearized%x(i,j) - coefficients(k) * &
                     vector_glsc3(test_modes(field,i), weights, pointwise)
             end do
          end do
       end do
    end do
    call neko_scratch_registry%relinquish_vector([pointwise_index, product_index])
    nullify(pointwise); nullify(product)
  end subroutine rom_add_linearized_convection

  !> Add the scalar-advection Jacobian about a forward coupled state.
  !! The derivative is `-(delta_u.grad(phi_f) + u_f.grad(delta_phi))`, with
  !! the configured scalar metric applied to the scalar residual row.
  subroutine rom_add_linearized_scalar_advection(weights, scalar_weight, &
       forward_velocity_modes, forward_scalar_modes, forward_scalar_gradients, &
       test_scalar_modes, perturb_velocity_modes, perturb_scalar_modes, &
       perturb_scalar_gradients, coefficients, linearized)
    type(vector_t), intent(in) :: weights
    real(kind=rp), intent(in) :: scalar_weight
    type(vector_t), intent(in) :: forward_velocity_modes(:, :)
    type(vector_t), intent(in) :: forward_scalar_modes(:), forward_scalar_gradients(:, :)
    type(vector_t), intent(in) :: test_scalar_modes(:)
    type(vector_t), intent(in) :: perturb_velocity_modes(:, :)
    type(vector_t), intent(in) :: perturb_scalar_modes(:), perturb_scalar_gradients(:, :)
    real(kind=rp), intent(in) :: coefficients(:)
    type(matrix_t), intent(inout) :: linearized
    type(vector_t), pointer :: pointwise, product
    integer :: pointwise_index, product_index
    integer :: n, rf, re, i, j, k, direction

    n = weights%size(); rf = size(forward_velocity_modes,2); re = size(perturb_velocity_modes,2)
    if (size(forward_velocity_modes,1) /= 3 .or. size(perturb_velocity_modes,1) /= 3 .or. &
         size(forward_scalar_modes) /= rf .or. size(perturb_scalar_modes) /= re .or. &
         size(forward_scalar_gradients,1) /= rf .or. size(forward_scalar_gradients,2) /= 3 .or. &
         size(perturb_scalar_gradients,1) /= re .or. size(perturb_scalar_gradients,2) /= 3 .or. &
         size(test_scalar_modes) /= re .or. size(coefficients) /= rf .or. &
         linearized%get_nrows() /= re .or. linearized%get_ncols() /= re) &
         error stop 'ROM linearized scalar advection shape mismatch'
    if (scalar_weight <= 0.0_rp) error stop 'ROM scalar metric must be positive'

    call neko_scratch_registry%request_vector(pointwise, pointwise_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    do i = 1, re
       do j = 1, re
          do k = 1, rf
             pointwise = 0.0_rp
             do direction = 1, 3
                call vector_col3(product, perturb_velocity_modes(direction,j), &
                     forward_scalar_gradients(k,direction))
                call vector_add2(pointwise, product)
                call vector_col3(product, forward_velocity_modes(direction,k), &
                     perturb_scalar_gradients(j,direction))
                call vector_add2(pointwise, product)
             end do
             linearized%x(i,j) = linearized%x(i,j) - scalar_weight * coefficients(k) * &
                  vector_glsc3(test_scalar_modes(i), weights, pointwise)
          end do
       end do
    end do
    call neko_scratch_registry%relinquish_vector([pointwise_index, product_index])
    nullify(pointwise); nullify(product)
  end subroutine rom_add_linearized_scalar_advection

  !> Assemble scalar advection and diffusion in the common ROM basis.
  !!
  !! For
  !!   phi = phi_L + sum_j a_j phi_j and
  !!   u   = u_L   + sum_j a_j u_j,
  !! this adds the negative weak residual of u.grad(phi) and the implicit
  !! diffusion matrix to the same reduced rows used by the velocity equation.
  !! The routine deliberately takes spatial diffusivity so nonuniform material
  !! properties follow Neko's scalar scheme rather than a guessed Péclet
  !! constant.
  subroutine rom_assemble_scalar_transport(weights, scalar_weight, velocity_modes, &
       scalar_modes, velocity_gradients, scalar_gradients, lift_velocity, &
       scalar_lift, lift_scalar_gradient, diffusivity, linear, constant, quadratic, &
       implicit_linear)
    type(vector_t), intent(in) :: weights
    real(kind=rp), intent(in) :: scalar_weight
    type(vector_t), intent(in) :: velocity_modes(:, :), scalar_modes(:)
    type(vector_t), intent(in) :: velocity_gradients(:, :, :)
    type(vector_t), intent(in) :: scalar_gradients(:, :)
    type(vector_t), intent(in) :: lift_velocity(:), scalar_lift, lift_scalar_gradient(:)
    type(vector_t), intent(in) :: diffusivity
    type(matrix_t), intent(inout) :: linear, implicit_linear
    type(vector_t), intent(inout) :: constant
    type(matrix_t), intent(inout) :: quadratic(:)
    call rom_assemble_scalar_transport_petrov(weights, scalar_weight, velocity_modes, &
         scalar_modes, velocity_gradients, scalar_gradients, scalar_modes, scalar_gradients, &
         lift_velocity, &
         scalar_lift, lift_scalar_gradient, diffusivity, linear, constant, quadratic, &
         implicit_linear)
  end subroutine rom_assemble_scalar_transport

  !> Assemble scalar transport with trial velocity/scalar states and test
  !! scalar modes.  The scalar equation shares the Petrov test index with the
  !! coupled reduced state.
  subroutine rom_assemble_scalar_transport_petrov(weights, scalar_weight, &
       trial_velocity_modes, trial_scalar_modes, trial_velocity_gradients, &
       trial_scalar_gradients, test_scalar_modes, test_scalar_gradients, lift_velocity, scalar_lift, &
       lift_scalar_gradient, diffusivity, linear, constant, quadratic, implicit_linear)
    type(vector_t), intent(in) :: weights
    real(kind=rp), intent(in) :: scalar_weight
    type(vector_t), intent(in) :: trial_velocity_modes(:, :), trial_scalar_modes(:)
    type(vector_t), intent(in) :: trial_velocity_gradients(:, :, :)
    type(vector_t), intent(in) :: trial_scalar_gradients(:, :)
    type(vector_t), intent(in) :: test_scalar_modes(:)
    type(vector_t), intent(in) :: test_scalar_gradients(:, :)
    type(vector_t), intent(in) :: lift_velocity(:), scalar_lift, lift_scalar_gradient(:)
    type(vector_t), intent(in) :: diffusivity
    type(matrix_t), intent(inout) :: linear, implicit_linear
    type(vector_t), intent(inout) :: constant
    type(matrix_t), intent(inout) :: quadratic(:)
    type(vector_t), pointer :: effective, pointwise, product
    integer :: effective_index, pointwise_index, product_index
    integer :: n, r, i, j, k, direction

    n = weights%size(); r = size(trial_scalar_modes)
    if (size(trial_velocity_modes,1) /= 3 .or. size(trial_velocity_modes,2) /= r .or. &
         size(test_scalar_modes) /= r .or. size(trial_velocity_gradients,1) /= 3 .or. &
         size(trial_velocity_gradients,2) /= r .or. size(trial_velocity_gradients,3) /= 3 .or. &
         size(trial_scalar_gradients,1) /= r .or. size(trial_scalar_gradients,2) /= 3 .or. &
         size(test_scalar_gradients,1) /= r .or. size(test_scalar_gradients,2) /= 3 .or. &
         size(lift_velocity) /= 3 .or. &
         scalar_lift%size() /= n .or. &
         size(lift_scalar_gradient) /= 3 .or. diffusivity%size() /= n) then
       error stop 'ROM scalar transport shape mismatch'
    end if
    if (linear%get_nrows() /= r .or. linear%get_ncols() /= r .or. &
         implicit_linear%get_nrows() /= r .or. implicit_linear%get_ncols() /= r .or. &
         size(quadratic) /= r .or. any([(quadratic(i)%get_nrows(), i=1,r)] /= r) .or. &
         any([(quadratic(i)%get_ncols(), i=1,r)] /= r)) then
       error stop 'ROM scalar transport operator shape mismatch'
    end if
    call neko_scratch_registry%request_vector(effective, effective_index, n, .true.)
    call neko_scratch_registry%request_vector(pointwise, pointwise_index, n, .true.)
    call neko_scratch_registry%request_vector(product, product_index, n, .true.)
    effective = weights
    call vector_cmult(effective, scalar_weight)

    ! The outputs are small host-owned reduced coefficients.  This routine
    ! accumulates the scalar block into the velocity operators assembled by
    ! the preceding routines, so their host values are intentionally retained.

    ! Lift--lift and lift--mode advection terms.  The lift--lift term must use
    ! grad(phi_L); using grad(phi_i) here would silently turn a constant
    ! residual into a state-dependent term and is not the Galerkin expansion
    ! of u.grad(phi).
    do i = 1, r
       pointwise = 0.0_rp
       do direction = 1, 3
          call vector_col3(product, lift_velocity(direction), &
               lift_scalar_gradient(direction))
          call vector_add2(pointwise, product)
       end do
       constant%x(i) = constant%x(i) - vector_glsc3( &
            test_scalar_modes(i), effective, pointwise)
       do j = 1, r
          pointwise = 0.0_rp
          do direction = 1, 3
             call vector_col3(product, trial_velocity_modes(direction,j), &
                  lift_scalar_gradient(direction))
             call vector_add2(pointwise, product)
             call vector_col3(product, lift_velocity(direction), &
                  trial_scalar_gradients(j,direction))
             call vector_add2(pointwise, product)
          end do
          linear%x(i,j) = linear%x(i,j) - vector_glsc3( &
               test_scalar_modes(i), effective, pointwise)
       end do
    end do

    ! Nonlinear scalar advection u_j.grad(phi_k).
    do i = 1, r
       do j = 1, r
          do k = 1, r
             pointwise = 0.0_rp
             do direction = 1, 3
                call vector_col3(product, trial_velocity_modes(direction,j), &
                     trial_scalar_gradients(k,direction))
                call vector_add2(pointwise, product)
             end do
             quadratic(i)%x(j,k) = quadratic(i)%x(j,k) - &
                  vector_glsc3(test_scalar_modes(i), effective, pointwise)
          end do
       end do
    end do

    ! Scalar diffusion is implicit in the BDF1 map, exactly as the velocity
    ! diffusion assembled by rom_assemble_diffusion.
    call vector_col3(effective, weights, diffusivity)
    call vector_cmult(effective, scalar_weight)
    do i = 1, r
       do j = 1, r
          do direction = 1, 3
             implicit_linear%x(i,j) = implicit_linear%x(i,j) - &
                  vector_glsc3(test_scalar_gradients(i,direction), effective, &
                  trial_scalar_gradients(j,direction))
          end do
       end do
       do direction = 1, 3
          ! The lift diffusion is a state-independent residual term.  It
          ! belongs in the explicit constant vector; fixed_implicit is
          ! reserved for the homogeneous scalar fluctuation.
          constant%x(i) = constant%x(i) - &
               vector_glsc3(test_scalar_gradients(i,direction), effective, &
               lift_scalar_gradient(direction))
       end do
    end do
    call neko_scratch_registry%relinquish_vector([effective_index, pointwise_index, &
         product_index])
    nullify(effective); nullify(pointwise); nullify(product)
  end subroutine rom_assemble_scalar_transport_petrov

end module rom_assembly
