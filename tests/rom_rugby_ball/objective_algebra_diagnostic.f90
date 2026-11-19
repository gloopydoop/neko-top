program rom_objective_algebra_diagnostic
  use neko, only: neko_init, neko_finalize
  use num_types, only: rp
  use rom_array_test_adapters, only: test_assemble_viscous_quadratic, &
       test_assemble_brinkman_objective, test_viscous_value, test_brinkman_value, &
       test_brinkman_density
  implicit none

  integer, parameter :: n = 11, r = 3
  real(kind=rp) :: weights(n), mask(n), chi(n), modes(n,3,r), gradients(n,3,r,3)
  real(kind=rp) :: lift(n,3), viscous(r,r), brinkman(r,r), lift_linear(r)
  real(kind=rp) :: state(r), density(n), expected, value, fd, eps, domain
  real(kind=rp) :: chi_plus(n), chi_minus(n), objective_plus, objective_minus
  integer :: i, j, c, d

  call neko_init()

  do i = 1, n
     weights(i) = 0.2_rp + 0.03_rp * real(i, rp)
     mask(i) = merge(1.0_rp, 0.0_rp, mod(i, 3) /= 0)
     chi(i) = 0.4_rp + 0.02_rp * real(i, rp)
     do c = 1, 3
        lift(i,c) = 0.05_rp * real(c, rp) + 0.01_rp * real(i, rp)
        do j = 1, r
           modes(i,c,j) = 0.02_rp * real(i + c + 2*j, rp)
           do d = 1, 3
              gradients(i,c,j,d) = 0.01_rp * real(i + c + j + d, rp)
           end do
        end do
     end do
  end do
  state = [0.3_rp, -0.2_rp, 0.4_rp]
  domain = sum(weights * mask)

  call test_assemble_viscous_quadratic(weights, mask, 0.7_rp, gradients, viscous)
  expected = 0.0_rp
  do i = 1, r
     do j = 1, r
        do c = 1, 3
           do d = 1, 3
              expected = expected + 0.5_rp * state(i) * state(j) * &
                   sum(weights * mask * gradients(:,c,i,d) * gradients(:,c,j,d)) * 0.7_rp
           end do
        end do
     end do
  end do
  value = test_viscous_value(state, viscous, domain)
  if (abs(value - expected / domain) > 1.0e-13_rp) &
       error stop 'viscous objective algebra mismatch'

  call test_assemble_brinkman_objective(weights, mask, chi, modes, lift, &
       brinkman, lift_linear, expected)
  value = test_brinkman_value(state, expected, lift_linear, brinkman, domain)
  expected = 0.0_rp
  do i = 1, n
     expected = expected + 0.5_rp * weights(i) * mask(i) * chi(i) * &
          sum((lift(i,:) + matmul(modes(i,:,:), state))**2)
  end do
  expected = expected / domain
  if (abs(value - expected) > 1.0e-13_rp) &
       error stop 'Brinkman objective algebra mismatch'

  call test_brinkman_density(weights, mask, lift, modes, state, domain, density)
  eps = 1.0e-7_rp
  chi_plus = chi + eps * ([(0.1_rp * real(i, rp), i=1,n)])
  chi_minus = chi - eps * ([(0.1_rp * real(i, rp), i=1,n)])
  objective_plus = 0.0_rp; objective_minus = 0.0_rp
  do i = 1, n
     objective_plus = objective_plus + 0.5_rp * weights(i) * mask(i) * chi_plus(i) * &
          sum((lift(i,:) + matmul(modes(i,:,:), state))**2) / domain
     objective_minus = objective_minus + 0.5_rp * weights(i) * mask(i) * chi_minus(i) * &
          sum((lift(i,:) + matmul(modes(i,:,:), state))**2) / domain
  end do
  fd = (objective_plus - objective_minus) / (2.0_rp * eps)
  value = sum(weights * 0.1_rp * [(real(i, rp), i=1,n)] * density)
  if (abs(fd - value) > 1.0e-8_rp) error stop 'Brinkman density FD mismatch'

  write(*,'(A,ES12.4,A,ES12.4)') 'ROM objective algebra passed: value/FD=', value, ' / ', fd
  call neko_finalize()
end program rom_objective_algebra_diagnostic
