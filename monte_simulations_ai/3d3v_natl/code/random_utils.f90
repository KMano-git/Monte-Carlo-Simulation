!===============================================================================
! Module: random_utils
! 乱数生成とMaxwell分布サンプリング
!===============================================================================
module random_utils
   use constants, only: dp, M_D_kg, EV_TO_J, PI
   use data_types, only: plasma_params
   implicit none

contains

   !---------------------------------------------------------------------------
   ! 乱数シードの初期化
   !---------------------------------------------------------------------------
   subroutine init_random_seed(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: seed_array(:)

      call random_seed(size=n)
      allocate(seed_array(n))

      do i = 1, n
         seed_array(i) = seed + i - 1
      end do

      call random_seed(put=seed_array)
      deallocate(seed_array)
   end subroutine init_random_seed

   !---------------------------------------------------------------------------
   ! Maxwell分布から速度をサンプリング（Box-Muller法）
   !---------------------------------------------------------------------------
   subroutine sample_maxwell_velocity(T_eV, vx, vy, vz)
      real(dp), intent(in)  :: T_eV
      real(dp), intent(out) :: vx, vy, vz

      real(dp) :: v_th, r1, r2, r3, r4, r5, r6

      v_th = sqrt(T_eV * EV_TO_J / M_D_kg)

      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      call random_number(r4)
      call random_number(r5)
      call random_number(r6)

      vx = v_th * sqrt(-2.0d0 * log(max(r1, 1.0d-30))) * cos(2.0d0 * PI * r2)
      vy = v_th * sqrt(-2.0d0 * log(max(r3, 1.0d-30))) * cos(2.0d0 * PI * r4)
      vz = v_th * sqrt(-2.0d0 * log(max(r5, 1.0d-30))) * cos(2.0d0 * PI * r6)
   end subroutine sample_maxwell_velocity

   !---------------------------------------------------------------------------
   ! Maxwell分布から背景イオン速度をサンプリング（流速オフセット付き）
   !---------------------------------------------------------------------------
   subroutine sample_maxwell_velocity_ion(T_eV, plasma, vx, vy, vz)
      real(dp), intent(in)  :: T_eV
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(out) :: vx, vy, vz

      real(dp) :: v_th, r1, r2, r3, r4, r5, r6

      v_th = sqrt(T_eV * EV_TO_J / M_D_kg)

      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      call random_number(r4)
      call random_number(r5)
      call random_number(r6)

      vx = plasma%u_x + v_th * sqrt(-2.0d0 * log(max(r1, 1.0d-30))) * cos(2.0d0 * PI * r2)
      vy = plasma%u_y + v_th * sqrt(-2.0d0 * log(max(r3, 1.0d-30))) * cos(2.0d0 * PI * r4)
      vz = plasma%u_z + v_th * sqrt(-2.0d0 * log(max(r5, 1.0d-30))) * cos(2.0d0 * PI * r6)
   end subroutine sample_maxwell_velocity_ion

   !---------------------------------------------------------------------------
   ! 単一方向の速度を指定エネルギーで生成（ビーム）
   !---------------------------------------------------------------------------
   subroutine set_beam_velocity(E_eV, vx, vy, vz)
      real(dp), intent(in)  :: E_eV
      real(dp), intent(out) :: vx, vy, vz

      real(dp) :: v_mag

      v_mag = sqrt(2.0d0 * E_eV * EV_TO_J / M_D_kg)
      vx = v_mag
      vy = 0.0d0
      vz = 0.0d0
   end subroutine set_beam_velocity

end module random_utils
