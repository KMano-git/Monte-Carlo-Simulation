!===============================================================================
! Module: random_utils
! 乱数生成とMaxwell分布サンプリング
!===============================================================================
module random_utils
   use constants, only: dp, M_D_kg, EV_TO_J, PI
   implicit none

contains

   !-----------------------------------------------------------------------------
   ! 乱数シードの初期化
   !-----------------------------------------------------------------------------
   subroutine init_random_seed(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: seed_array(:)

      call random_seed(size=n)
      allocate(seed_array(n))

      ! シード値の設定
      do i = 1, n
         seed_array(i) = seed + i - 1
      end do

      call random_seed(put=seed_array)
      deallocate(seed_array)
   end subroutine init_random_seed

   !-----------------------------------------------------------------------------
   ! Maxwell分布から速度をサンプリング（Box-Muller法）
   !-----------------------------------------------------------------------------
   subroutine sample_maxwell_velocity(T_eV, vx, vy, vz)
      real(dp), intent(in)  :: T_eV  ! 温度 [eV]
      real(dp), intent(out) :: vx, vy, vz

      real(dp) :: v_th, r1, r2, r3, r4, r5, r6

      ! 熱速度
      v_th = sqrt(T_eV * EV_TO_J / M_D_kg)

      ! Box-Muller法でガウス分布を生成
      ! 各成分のペアに対して独立した乱数を使用
      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      call random_number(r4)
      call random_number(r5)
      call random_number(r6)

      ! 各成分は独立したガウス分布
      ! (r1, r2) で vx を生成
      vx = v_th * sqrt(-2.0d0 * log(max(r1, 1.0d-30))) * cos(2.0d0 * PI * r2)
      ! (r3, r4) で vy を生成
      vy = v_th * sqrt(-2.0d0 * log(max(r3, 1.0d-30))) * cos(2.0d0 * PI * r4)
      ! (r5, r6) で vz を生成
      vz = v_th * sqrt(-2.0d0 * log(max(r5, 1.0d-30))) * cos(2.0d0 * PI * r6)
   end subroutine sample_maxwell_velocity

   !-----------------------------------------------------------------------------
   ! 単一方向の速度を指定エネルギーで生成
   !-----------------------------------------------------------------------------
   subroutine set_beam_velocity(E_eV, vx, vy, vz)
      real(dp), intent(in)  :: E_eV  ! エネルギー [eV]
      real(dp), intent(out) :: vx, vy, vz

      real(dp) :: v_mag

      ! 速度の大きさ
      v_mag = sqrt(2.0d0 * E_eV * EV_TO_J / M_D_kg)

      ! +x方向のビーム
      vx = v_mag
      vy = 0.0d0
      vz = 0.0d0
   end subroutine set_beam_velocity

end module random_utils
