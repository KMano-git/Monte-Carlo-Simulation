program test_modules
   use constants
   use random_utils
   implicit none

   real(dp) :: T_eV, vx, vy, vz, E, E_sum, E_mean
   integer :: i
   integer, parameter :: N = 100000

   print *, "Testing Maxwell sampling using modules..."
   print *, "M_D in module:", M_D
   print *, "EV in module:", EV

   T_eV = 10.0d0
   E_sum = 0.0d0

   do i = 1, N
      call sample_maxwell_velocity(T_eV, vx, vy, vz)
      E = 0.5d0 * M_D * (vx**2 + vy**2 + vz**2) * J_TO_EV
      E_sum = E_sum + E
   end do

   E_mean = E_sum / real(N, dp)

   print *, "Target T_bg:", T_eV
   print *, "Expected E_mean (3/2 T):", 1.5d0 * T_eV
   print *, "Sampled E_mean:", E_mean
   print *, "Ratio:", E_mean / (1.5d0 * T_eV)

end program test_modules
