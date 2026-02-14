!===============================================================================
! Module: collisions
! 衝突プロセス（弾性散乱と荷電交換）
!===============================================================================
module collisions
   use constants, only: dp, M_D, EV_TO_J, J_TO_EV, PI
   use particle_data, only: particle_0d, plasma_params
   use random_utils, only: sample_maxwell_velocity
   use cdf_reader, only: sample_scattering_angle
   implicit none

contains

   !-----------------------------------------------------------------------------
   ! 荷電交換処理
   ! 速度が背景イオンの速度に完全に入れ替わる
   !-----------------------------------------------------------------------------
   subroutine collision_cx(p, vx_bg, vy_bg, vz_bg, delta_E)
      type(particle_0d), intent(inout) :: p
      real(dp), intent(in) :: vx_bg, vy_bg, vz_bg
      real(dp), intent(out) :: delta_E

      real(dp) :: E_before, E_after

      ! 衝突前のエネルギー
      E_before = 0.5d0 * M_D * (p%vx**2 + p%vy**2 + p%vz**2) * J_TO_EV

      ! 速度を完全に入れ替え
      p%vx = vx_bg
      p%vy = vy_bg
      p%vz = vz_bg

      ! 衝突後のエネルギー
      E_after = 0.5d0 * M_D * (p%vx**2 + p%vy**2 + p%vz**2) * J_TO_EV

      delta_E = E_after - E_before

   end subroutine collision_cx

   !-----------------------------------------------------------------------------
   ! 弾性散乱処理
   ! 重心系変換を用いた散乱（熱運動するターゲットを考慮）
   !-----------------------------------------------------------------------------
   subroutine collision_elastic(p, vx_bg, vy_bg, vz_bg, use_isotropic, E_collision, delta_E)
      type(particle_0d), intent(inout) :: p
      real(dp), intent(in) :: vx_bg, vy_bg, vz_bg
      logical, intent(in) :: use_isotropic
      real(dp), intent(in) :: E_collision  ! 衝突エネルギー [eV] (CDFサンプリング用)
      real(dp), intent(out) :: delta_E

      real(dp) :: E_before, E_after
      ! vx_bg, vy_bg, vz_bg は引数で受け取る
      real(dp) :: vx_cm, vy_cm, vz_cm      ! 重心速度 [m/s]
      real(dp) :: ux, uy, uz, u_mag        ! 相対速度 [m/s]
      real(dp) :: ux_new, uy_new, uz_new   ! 散乱後の相対速度 [m/s]
      real(dp) :: chi, phi                 ! 散乱角 [rad]
      real(dp) :: r_chi, E_rel             ! CDFサンプリング用

      ! 衝突前のエネルギー
      E_before = 0.5d0 * M_D * (p%vx**2 + p%vy**2 + p%vz**2) * J_TO_EV

      ! 重心速度（等質量なので平均）
      vx_cm = 0.5d0 * (p%vx + vx_bg)
      vy_cm = 0.5d0 * (p%vy + vy_bg)
      vz_cm = 0.5d0 * (p%vz + vz_bg)

      ! 相対速度
      ux = p%vx - vx_bg
      uy = p%vy - vy_bg
      uz = p%vz - vz_bg
      u_mag = sqrt(ux**2 + uy**2 + uz**2)

      ! 相対速度がゼロに近い場合は散乱なし
      if (u_mag < 1.0d-30) then
         delta_E = 0.0d0
         return
      end if

      ! 散乱角の決定
      if (use_isotropic) then
         call sample_scattering_angle_isotropic(chi)
      else
         call random_number(r_chi)
         ! CDFから散乱角をサンプリング
         ! 重心系エネルギー（相対エネルギー）を使用: E_rel = 1/4 * m * v_rel^2
         ! v_relは上で既に計算済み
         E_rel = 0.25d0 * M_D * u_mag**2 * J_TO_EV
         chi = sample_scattering_angle(E_rel, r_chi)
      end if

      ! 方位角
      call random_number(phi)
      phi = 2.0d0 * PI * phi

      ! 相対速度ベクトルを散乱角で回転
      call rotate_vector(ux, uy, uz, chi, phi, ux_new, uy_new, uz_new)

      ! 実験室系に戻す（等質量なので u' = 2*(v_test' - V_cm)）
      ! → v_test' = V_cm + u'/2
      p%vx = vx_cm + 0.5d0 * ux_new
      p%vy = vy_cm + 0.5d0 * uy_new
      p%vz = vz_cm + 0.5d0 * uz_new

      ! 衝突後のエネルギー
      E_after = 0.5d0 * M_D * (p%vx**2 + p%vy**2 + p%vz**2) * J_TO_EV

      delta_E = E_after - E_before

   end subroutine collision_elastic

   !-----------------------------------------------------------------------------
   ! 等方散乱角サンプリング
   !-----------------------------------------------------------------------------
   subroutine sample_scattering_angle_isotropic(chi)
      real(dp), intent(out) :: chi

      real(dp) :: r, cos_chi

      call random_number(r)
      cos_chi = 2.0d0 * r - 1.0d0  ! cos(chi) ∈ [-1, 1]
      chi = acos(cos_chi)

   end subroutine sample_scattering_angle_isotropic

   !-----------------------------------------------------------------------------
   ! 3Dベクトル回転（Rodriguesの回転公式）
   !-----------------------------------------------------------------------------
   subroutine rotate_vector(ux, uy, uz, chi, phi, ux_new, uy_new, uz_new)
      real(dp), intent(in)  :: ux, uy, uz      ! 元のベクトル
      real(dp), intent(in)  :: chi, phi        ! 散乱角と方位角 [rad]
      real(dp), intent(out) :: ux_new, uy_new, uz_new

      real(dp) :: u_mag, ex, ey, ez
      real(dp) :: sin_chi, cos_chi, sin_phi, cos_phi
      real(dp) :: e_perp

      u_mag = sqrt(ux**2 + uy**2 + uz**2)

      if (u_mag < 1.0d-30) then
         ux_new = 0.0d0
         uy_new = 0.0d0
         uz_new = 0.0d0
         return
      end if

      ! 単位ベクトル
      ex = ux / u_mag
      ey = uy / u_mag
      ez = uz / u_mag

      sin_chi = sin(chi)
      cos_chi = cos(chi)
      sin_phi = sin(phi)
      cos_phi = cos(phi)

      ! 垂直方向の基底ベクトルを構築して回転
      ! （Rodriguesの回転公式を使用）
      if (abs(ez) > 0.9d0) then
         ! z方向が支配的な場合、x-y平面で垂直ベクトルを作る
         e_perp = sqrt(ex**2 + ey**2)
         if (e_perp > 1.0d-10) then
            ux_new = u_mag * (sin_chi * (-ey * cos_phi / e_perp) + cos_chi * ex)
            uy_new = u_mag * (sin_chi * (ex * cos_phi / e_perp) + cos_chi * ey)
            uz_new = u_mag * (sin_chi * sin_phi + cos_chi * ez)
         else
            ux_new = u_mag * sin_chi * cos_phi
            uy_new = u_mag * sin_chi * sin_phi
            uz_new = u_mag * cos_chi
         end if
      else
         ! それ以外の場合、x-z平面で垂直ベクトルを作る
         e_perp = sqrt(ex**2 + ez**2)
         if (e_perp > 1.0d-10) then
            ux_new = u_mag * (sin_chi * (-ez * cos_phi / e_perp) + cos_chi * ex)
            uy_new = u_mag * (sin_chi * sin_phi + cos_chi * ey)
            uz_new = u_mag * (sin_chi * (ex * cos_phi / e_perp) + cos_chi * ez)
         else
            ux_new = u_mag * sin_chi * cos_phi
            uy_new = u_mag * cos_chi
            uz_new = u_mag * sin_chi * sin_phi
         end if
      end if

   end subroutine rotate_vector

end module collisions
