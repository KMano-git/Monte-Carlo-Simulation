!===============================================================================
! Module: collisions
! 衝突処理と衝突判定
!===============================================================================
module collisions
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV, PI
   use particle_data, only: particle_t, plasma_params
   use random_utils, only: sample_maxwell_velocity
   use cdf_reader, only: sample_scattering_angle
   use cross_sections, only: sigma_cx, sigma_ei, sigma_el, get_all_cross_sections
   implicit none

contains
   !-----------------------------------------------------------------------------
   ! 荷電交換処理
   ! 速度が背景イオンの速度に完全に入れ替わる
   !-----------------------------------------------------------------------------
   subroutine collision_cx(p, vx_i, vy_i, vz_i, delta_E)
      type(particle_t), intent(inout) :: p
      real(dp), intent(in)  :: vx_i, vy_i, vz_i
      real(dp), intent(out) :: delta_E

      real(dp) :: E_old, E_new

      ! 衝突前のエネルギー
      E_old = 0.5d0 * M_D_kg * (p%vx * p%vx + p%vy * p%vy + p%vz * p%vz)

      ! 速度の入れ替え
      p%vx = vx_i
      p%vy = vy_i
      p%vz = vz_i

      ! 衝突後のエネルギー
      E_new = 0.5d0 * M_D_kg * (p%vx * p%vx + p%vy * p%vy + p%vz * p%vz)

      ! エネルギー変化量
      delta_E = E_new - E_old
   end subroutine collision_cx

   !-----------------------------------------------------------------------------
   ! 弾性散乱処理
   ! 重心系変換を用いた散乱（熱運動するターゲットを考慮）
   !-----------------------------------------------------------------------------
   subroutine collision_el(p, vx_i, vy_i, vz_i, use_isotropic, delta_E)
      type(particle_t), intent(inout) :: p
      real(dp), intent(in)  :: vx_i, vy_i, vz_i
      logical, intent(in)   :: use_isotropic
      real(dp), intent(out) :: delta_E

      real(dp) :: E_old, E_new
      real(dp) :: vx_g, vy_g, vz_g ! 重心速度
      real(dp) :: ux, uy, uz, u_mag ! 相対速度
      real(dp) :: ux_new, uy_new, uz_new ! 散乱後の相対速度
      real(dp) :: chi, phi ! 散乱角と方位角
      real(dp) :: r_chi, E_rel ! CDFサンプリング用

      ! 衝突前のエネルギー
      E_old = 0.5d0 * M_D_kg * (p%vx * p%vx + p%vy * p%vy + p%vz * p%vz)

      ! 重心速度
      vx_g = 0.5d0 * (p%vx + vx_i)
      vy_g = 0.5d0 * (p%vy + vy_i)
      vz_g = 0.5d0 * (p%vz + vz_i)

      ! 相対速度
      ux = p%vx - vx_i
      uy = p%vy - vy_i
      uz = p%vz - vz_i
      u_mag = sqrt(ux*ux + uy*uy + uz*uz)

      if (u_mag < 1.0d-30) return

      if (use_isotropic) then
         call sample_isotropic_scattering(chi)
      else
         ! 衝突の相対エネルギー(eV)
         E_rel = 0.25d0 * M_D_kg * u_mag * u_mag * J_TO_EV
         call random_number(r_chi)
         chi = sample_scattering_angle(E_rel, r_chi)
      end if

      ! 方位角は完全にランダム
      call random_number(phi)
      phi = 2.0d0 * PI * phi

      ! 相対速度ベクトルを散乱角で回転
      call rotate_vector(ux, uy, uz, chi, phi, ux_new, uy_new, uz_new)

      ! 実験室系に速度を戻す
      p%vx = vx_g + 0.5d0 * ux_new
      p%vy = vy_g + 0.5d0 * uy_new
      p%vz = vz_g + 0.5d0 * uz_new

      ! 衝突後のエネルギー
      E_new = 0.5d0 * M_D_kg * (p%vx * p%vx + p%vy * p%vy + p%vz * p%vz)

      ! エネルギー変化量
      delta_E = E_new - E_old
   end subroutine collision_el

   !-----------------------------------------------------------------------------
   ! 等方散乱角サンプリング（use_isotropic = .true. の場合）
   !-----------------------------------------------------------------------------
   subroutine sample_isotropic_scattering(chi)
      real(dp), intent(inout) :: chi

      real(dp) :: r, cos_chi

      call random_number(r)
      cos_chi = 2.0d0 * r - 1.0d0
      chi = acos(cos_chi)
   end subroutine sample_isotropic_scattering

   !-----------------------------------------------------------------------------
   ! 3D�����トル回転（Rodriguesの回転公式）:rotate_vector
   !-----------------------------------------------------------------------------
   subroutine rotate_vector(ux, uy, uz, chi, phi, ux_new, uy_new, uz_new)
      real(dp), intent(in) :: ux, uy, uz
      real(dp), intent(in) :: chi, phi
      real(dp), intent(out):: ux_new, uy_new, uz_new

      real(dp) :: u_mag, ex, ey, ez
      real(dp) :: sin_chi, cos_chi, sin_phi, cos_phi
      real(dp) :: e_perp

      u_mag = sqrt(ux*ux + uy*uy + uz*uz)

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
         e_perp = sqrt(ex*ex + ey*ey)
         if (e_perp > 1.0d-10) then
            ux_new = u_mag * (sin_chi * (-ey * cos_phi / e_perp) + cos_chi * ex)
            uy_new = u_mag * (sin_chi * (ex * cos_phi / e_perp) + cos_chi * ey)
            uz_new = u_mag * (sin_chi * sin_phi + cos_chi * ez)
         else
            ux_new = u_mag * sin_chi
            uy_new = 0.0d0
            uz_new = u_mag * cos_chi
         end if
      else
         ! それ以外の場合、x-z平面で垂直ベクトルを作る
         e_perp = sqrt(1.0d0 - ez*ez)
         ux_new = u_mag * (sin_chi * (-ey * cos_phi / e_perp) + cos_chi * ex)
         uy_new = u_mag * (sin_chi * (ex * cos_phi / e_perp) + cos_chi * ey)
         uz_new = u_mag * (sin_chi * sin_phi + cos_chi * ez)
      end if

   end subroutine

end module collisions
