program test_maxwell
   implicit none
   integer, parameter :: dp = selected_real_kind(15, 307)
   real(dp), parameter :: PI = 3.141592653589793d0
   real(dp), parameter :: M_D  = 3.344d-27
   real(dp), parameter :: EV = 1.602d-19

   real(dp) :: T_eV, v_th, vx, vy, vz, E_sum, E
   integer :: i, n_samples
   real(dp) :: r1, r2, r3, r4, r5, r6

   T_eV = 10.0d0
   n_samples = 100000
   E_sum = 0.0d0

   v_th = sqrt(T_eV * EV / M_D)

   do i = 1, n_samples
      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      call random_number(r4)
      call random_number(r5)
      call random_number(r6)

      vx = v_th * sqrt(-2.0d0 * log(max(r1, 1.0d-30))) * cos(2.0d0 * PI * r2)
      vy = v_th * sqrt(-2.0d0 * log(max(r3, 1.0d-30))) * cos(2.0d0 * PI * r4)
      vz = v_th * sqrt(-2.0d0 * log(max(r5, 1.0d-30))) * cos(2.0d0 * PI * r6)

      E = 0.5d0 * M_D * (vx**2 + vy**2 + vz**2) / EV
      E_sum = E_sum + E
   end do

   print *, "T_eV = ", T_eV
   print *, "Expected <E> = 3/2 * T = ", 1.5d0 * T_eV
   print *, "Sampled <E> = ", E_sum / n_samples

end program test_maxwell
