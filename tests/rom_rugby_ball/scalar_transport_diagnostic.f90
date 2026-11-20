program rom_rugby_ball_scalar_transport_diagnostic
  use neko, only: neko_init, neko_finalize
  use num_types, only: rp
  use rom_array_test_adapters, only: test_assemble_scalar_transport
  implicit none

  integer, parameter :: n = 9, r = 2
  real(kind=rp) :: weights(n), velocity_modes(n,3,r), scalar_modes(n,r)
  real(kind=rp) :: velocity_gradients(n,3,r,3), scalar_gradients(n,r,3)
  real(kind=rp) :: lift_velocity(n,3), scalar_lift(n), lift_scalar_gradient(n,3)
  real(kind=rp) :: diffusivity(n), linear(r,r), constant(r), quadratic(r,r,r)
  real(kind=rp) :: implicit_linear(r,r)
  real(kind=rp) :: expected_linear(r,r), expected_constant(r), expected_quadratic(r,r,r)
  real(kind=rp) :: expected_implicit(r,r), error_value, scale, pointwise
  real(kind=rp) :: scalar_weight
  integer :: i, j, k, c, d

  call neko_init()
  scalar_weight = 0.7_rp
  do i = 1, n
     weights(i) = 0.4_rp + 0.02_rp * real(i, rp)
     diffusivity(i) = 0.3_rp + 0.01_rp * real(mod(i, 4), rp)
     scalar_lift(i) = 0.2_rp + 0.03_rp * real(i, rp)
     do c = 1, 3
        lift_velocity(i,c) = 0.04_rp * real(c, rp) + 0.005_rp * real(i, rp)
        do d = 1, 3
           lift_scalar_gradient(i,d) = 0.01_rp * real(i + c + d, rp)
        end do
        do j = 1, r
           velocity_modes(i,c,j) = 0.02_rp * real(i + c + j, rp)
           do d = 1, 3
              velocity_gradients(i,c,j,d) = 0.003_rp * real(i + c + j + d, rp)
           end do
        end do
     end do
     do j = 1, r
        scalar_modes(i,j) = 0.015_rp * real(i + 2*j, rp)
        do d = 1, 3
           scalar_gradients(i,j,d) = 0.004_rp * real(i + j + d, rp)
        end do
     end do
  end do

  call test_assemble_scalar_transport(weights, scalar_weight, velocity_modes, &
       scalar_modes, velocity_gradients, scalar_gradients, lift_velocity, scalar_lift, &
       lift_scalar_gradient, diffusivity, linear, constant, quadratic, implicit_linear)

  expected_constant = 0.0_rp
  expected_linear = 0.0_rp
  expected_quadratic = 0.0_rp
  expected_implicit = 0.0_rp
  do i = 1, r
     do d = 1, 3
        expected_constant(i) = expected_constant(i) - scalar_weight * &
             sum(weights * scalar_modes(:,i) * lift_velocity(:,d) * &
             lift_scalar_gradient(:,d))
        expected_constant(i) = expected_constant(i) - scalar_weight * &
             sum(weights * diffusivity * scalar_gradients(:,i,d) * &
             lift_scalar_gradient(:,d))
     end do
     do j = 1, r
        do d = 1, 3
           expected_linear(i,j) = expected_linear(i,j) - scalar_weight * &
                sum(weights * scalar_modes(:,i) * &
                (velocity_modes(:,d,j) * lift_scalar_gradient(:,d) + &
                 lift_velocity(:,d) * scalar_gradients(:,j,d)))
           expected_implicit(i,j) = expected_implicit(i,j) - scalar_weight * &
                sum(weights * diffusivity * scalar_gradients(:,i,d) * &
                scalar_gradients(:,j,d))
        end do
        do k = 1, r
           do d = 1, 3
              pointwise = sum(weights * scalar_modes(:,i) * &
                   velocity_modes(:,d,j) * scalar_gradients(:,k,d))
              expected_quadratic(i,j,k) = expected_quadratic(i,j,k) - &
                   scalar_weight * pointwise
           end do
        end do
     end do
  end do

  error_value = maxval(abs(linear - expected_linear))
  error_value = max(error_value, maxval(abs(constant - expected_constant)))
  error_value = max(error_value, maxval(abs(implicit_linear - expected_implicit)))
  error_value = max(error_value, maxval(abs(quadratic - expected_quadratic)))
  scale = max(1.0_rp, maxval(abs(expected_linear)), maxval(abs(expected_constant)), &
       maxval(abs(expected_implicit)), maxval(abs(expected_quadratic)))
  write(*,'(A,ES12.4)') 'Scalar transport assembly relative error: ', error_value / scale
  if (error_value / scale > 1.0e-13_rp) &
       error stop 'scalar transport assembly mismatch'

  call neko_finalize()
end program rom_rugby_ball_scalar_transport_diagnostic
