program test_collision_sampling
   use constants
   use random_utils
   implicit none

   real(dp) :: T_eV, v_th
   real(dp) :: vx_test, vy_test, vz_test
   real(dp) :: vx_bg, vy_bg, vz_bg
   real(dp) :: E_sum, E_mean, v_bg_mag
   real(dp) :: r, ratio
   integer :: i, count
   integer, parameter :: N = 100000

   ! テスト設定
   T_eV = 10.0d0
   v_th = sqrt(T_eV * EV_TO_J / M_D)

   ! テスト粒子（熱速度程度とする）
   vx_test = v_th
   vy_test = 0.0d0
   vz_test = 0.0d0

   print *, "Testing Partner Sampling Logic..."
   print *, "Test Particle v/v_th =", 1.0

   ! 1. 単純なMaxwellサンプリング（現在の実装）
   E_sum = 0.0d0
   do i = 1, N
      call sample_maxwell_velocity(T_eV, vx_bg, vy_bg, vz_bg)
      v_bg_mag = sqrt(vx_bg**2 + vy_bg**2 + vz_bg**2)
      E_sum = E_sum + 0.5d0 * M_D * v_bg_mag**2 * J_TO_EV
   end do
   print *, "Standard Sampling Mean E [eV]:", E_sum/N, " Expected:", 1.5*T_eV

   ! 2. 相対速度を考慮したサンプリング（修正案）
   E_sum = 0.0d0
   count = 0
   do while (count < N)
      ! 候補をサンプリング
      call sample_maxwell_velocity(T_eV, vx_bg, vy_bg, vz_bg)

      ! 相対速度
      ratio = sqrt((vx_test - vx_bg)**2 + (vy_test - vy_bg)**2 + (vz_test - vz_bg)**2)

      ! 棄却法判定 (最大相対速度をおおよそ 5*v_th + v_test と仮定)
      call random_number(r)
      if (r < ratio / (sqrt(vx_test**2 + vy_test**2 + vz_test**2) + 5.0d0*v_th)) then
         v_bg_mag = sqrt(vx_bg**2 + vy_bg**2 + vz_bg**2)
         E_sum = E_sum + 0.5d0 * M_D * v_bg_mag**2 * J_TO_EV
         count = count + 1
      end if
   end do
   print *, "Weighted Sampling Mean E [eV]:", E_sum/N
   print *, "  (Should be higher than standard due to preference for head-on/fast collisions)"

end program test_collision_sampling
