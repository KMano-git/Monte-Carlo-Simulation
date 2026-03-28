!===============================================================================
! Module: random_utils
! xoshiro256+乱数生成器、Maxwell分布サンプリング
! スレッドセーフ: 各粒子が独立したRNG状態を持つ
!===============================================================================
module random_utils
   use constants, only: dp, M_D_kg, EV_TO_J, PI
   use data_types, only: rng_state, plasma_params
   implicit none

   private
   public :: init_rng, random_double
   public :: init_random_seed
   public :: sample_maxwell_velocity, sample_maxwell_velocity_ion
   public :: set_beam_velocity

contains

   !---------------------------------------------------------------------------
   ! xoshiro256+ の内部ヘルパー: 64bitローテーション
   !---------------------------------------------------------------------------
   pure function rotl(x, k) result(res)
      integer(8), intent(in) :: x
      integer, intent(in) :: k
      integer(8) :: res
      res = ior(ishft(x, k), ishft(x, -(64 - k)))
   end function rotl

   !---------------------------------------------------------------------------
   ! xoshiro256+: [0, 1) の倍精度一様乱数を返す
   ! 参照: https://prng.di.unimi.it/xoshiro256plus.c
   !---------------------------------------------------------------------------
   function random_double(rng) result(r)
      type(rng_state), intent(inout) :: rng
      real(dp) :: r
      integer(8) :: result_val, t

      result_val = rng%s(1) + rng%s(4)

      ! [0, 1) への変換: 上位53bitを使用
      ! 2^-53 = 1.1102230246251565d-16
      r = dble(iand(result_val, int(Z'001FFFFFFFFFFFFF', 8))) * (1.0d0 / 9007199254740992.0d0)

      ! 状態遷移
      t = ishft(rng%s(2), 17)
      rng%s(3) = ieor(rng%s(3), rng%s(1))
      rng%s(4) = ieor(rng%s(4), rng%s(2))
      rng%s(2) = ieor(rng%s(2), rng%s(3))
      rng%s(1) = ieor(rng%s(1), rng%s(4))
      rng%s(3) = ieor(rng%s(3), t)
      rng%s(4) = rotl(rng%s(4), 45)
   end function random_double

   !---------------------------------------------------------------------------
   ! RNG状態の初期化（SplitMix64でシードから4つの64bit状態を生成）
   !---------------------------------------------------------------------------
   subroutine init_rng(rng, seed)
      type(rng_state), intent(out) :: rng
      integer(8), intent(in) :: seed
      integer(8) :: z
      integer :: i

      z = seed
      do i = 1, 4
         z = z + int(Z'9E3779B97F4A7C15', 8)
         z = ieor(z, ishft(z, -30))
         z = z * int(Z'BF58476D1CE4E5B9', 8)
         z = ieor(z, ishft(z, -27))
         z = z * int(Z'94D049BB133111EB', 8)
         z = ieor(z, ishft(z, -31))
         rng%s(i) = z
      end do
   end subroutine init_rng

   !---------------------------------------------------------------------------
   ! 従来の乱数シード初期化（非並列用、初期化段階でのみ使用）
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
   ! rng: 粒子固有のRNG状態
   !---------------------------------------------------------------------------
   subroutine sample_maxwell_velocity(rng, T_eV, vx, vy, vz)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in)  :: T_eV
      real(dp), intent(out) :: vx, vy, vz

      real(dp) :: v_th, r1, r2, r3, r4, r5, r6

      v_th = sqrt(T_eV * EV_TO_J / M_D_kg)

      r1 = random_double(rng)
      r2 = random_double(rng)
      r3 = random_double(rng)
      r4 = random_double(rng)
      r5 = random_double(rng)
      r6 = random_double(rng)

      vx = v_th * sqrt(-2.0d0 * log(max(r1, 1.0d-30))) * cos(2.0d0 * PI * r2)
      vy = v_th * sqrt(-2.0d0 * log(max(r3, 1.0d-30))) * cos(2.0d0 * PI * r4)
      vz = v_th * sqrt(-2.0d0 * log(max(r5, 1.0d-30))) * cos(2.0d0 * PI * r6)
   end subroutine sample_maxwell_velocity

   !---------------------------------------------------------------------------
   ! Maxwell分布から背景イオン速度をサンプリング（流速オフセット付き）
   ! rng: 粒子固有のRNG状態
   !---------------------------------------------------------------------------
   subroutine sample_maxwell_velocity_ion(rng, T_eV, plasma, vx, vy, vz)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in)  :: T_eV
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(out) :: vx, vy, vz

      real(dp) :: v_th, r1, r2, r3, r4, r5, r6

      v_th = sqrt(T_eV * EV_TO_J / M_D_kg)

      r1 = random_double(rng)
      r2 = random_double(rng)
      r3 = random_double(rng)
      r4 = random_double(rng)
      r5 = random_double(rng)
      r6 = random_double(rng)

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
