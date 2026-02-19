!===============================================================================
! Module: dynamics
! 粒子力学: 衝突判定（Non-Analog Null-Collision法）、衝突処理、位置更新
!===============================================================================
module dynamics
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV, PI
   use data_types, only: particle_t, plasma_params, sim_params
   use random_utils, only: sample_maxwell_velocity_ion
   use cross_sections, only: sigma_cx, sigma_el, sigma_v_max, &
      ionization_rate_coeff
   use cdf_reader, only: sample_scattering_angle
   implicit none

   private
   public :: advance_particle, check_collision_nonanalog
   public :: collision_cx, collision_el
   public :: update_weight

   !衝突タイプ定数
   integer, parameter, public :: COLL_NONE = 0
   integer, parameter, public :: COLL_CX   = 1
   integer, parameter, public :: COLL_EL   = 2

contains

   !---------------------------------------------------------------------------
   ! 位置の更新
   !---------------------------------------------------------------------------
   subroutine advance_particle(p, dt)
      type(particle_t), intent(inout) :: p
      real(dp), intent(in) :: dt

      p%x = p%x + p%vx * dt
      p%y = p%y + p%vy * dt
      p%z = p%z + p%vz * dt
   end subroutine advance_particle

   !---------------------------------------------------------------------------
   ! Non-Analog重み更新
   ! w_new = w_old * exp(-R_a * dt)
   !---------------------------------------------------------------------------
   subroutine update_weight(p, plasma, dt, weight_min)
      type(particle_t), intent(inout) :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: weight_min

      real(dp) :: R_a

      R_a = plasma%n_e * ionization_rate_coeff(plasma%T_e)

      if (R_a > 1.0d-30) then
         p%weight = p%weight * exp(-R_a * dt)
      end if

      !最小重み以下なら粒子消滅
      if (p%weight < weight_min) then
         p%alive = .false.
      end if
   end subroutine update_weight

   !---------------------------------------------------------------------------
   ! Non-Analog Null-Collision判定
   ! CX+ELのみでNull-Collision判定（電離は重み減衰で扱う）
   ! ν_max = n_i * sigma_v_max（事前計算済み）
   !---------------------------------------------------------------------------
   subroutine check_collision_nonanalog(p, plasma, sim, dt, &
      vx_i, vy_i, vz_i, &
      v_rel, E_rel, coll_type)
      type(particle_t), intent(in)  :: p
      type(plasma_params), intent(in) :: plasma
      type(sim_params), intent(in) :: sim
      real(dp), intent(in)  :: dt
      real(dp), intent(out) :: vx_i, vy_i, vz_i  !背景イオン速度
      real(dp), intent(out) :: v_rel              !相対速度
      real(dp), intent(out) :: E_rel              !相対エネルギー [eV]
      integer, intent(out)  :: coll_type

      real(dp) :: nu_max, P_coll
      real(dp) :: sig_cx_val, sig_el_val, sig_s
      real(dp) :: P_accept
      real(dp) :: r, r_accept, r_type

      coll_type = COLL_NONE
      vx_i = 0.0d0; vy_i = 0.0d0; vz_i = 0.0d0
      v_rel = 0.0d0; E_rel = 0.0d0

      !最大衝突頻度（事前計算済み）
      nu_max = plasma%n_i * sigma_v_max
      if (nu_max < 1.0d-30) return

      !衝突確率
      P_coll = 1.0d0 - exp(-nu_max * dt)

      !衝突判定
      call random_number(r)
      if (r > P_coll) return

      !背景イオン速度のサンプリング
      call sample_maxwell_velocity_ion(plasma%T_i, plasma, vx_i, vy_i, vz_i)

      !相対速度
      v_rel = sqrt((p%vx - vx_i)**2 + (p%vy - vy_i)**2 + (p%vz - vz_i)**2)

      !相対エネルギー（重心系）
      E_rel = 0.25d0 * M_D_kg * v_rel * v_rel * J_TO_EV

      !実断面積（CX+ELのみ）
      sig_cx_val = sigma_cx(E_rel)
      sig_el_val = sigma_el(E_rel)
      sig_s = sig_cx_val + sig_el_val

      !Rejection判定: P_accept = σ_s * v_rel / (sigma_v_max)
      P_accept = sig_s * v_rel / sigma_v_max

      !安全策: P_accept > 1 の場合（理論上は起きないはず）
      if (P_accept > 1.0d0) P_accept = 1.0d0

      call random_number(r_accept)
      if (r_accept > P_accept) return

      !CX/EL衝突タイプの選択
      call random_number(r_type)
      if (sig_s > 1.0d-30) then
         if (r_type < sig_cx_val / sig_s .and. sim%enable_cx) then
            coll_type = COLL_CX
         else if (sim%enable_el) then
            coll_type = COLL_EL
         end if
      end if

   end subroutine check_collision_nonanalog

   !---------------------------------------------------------------------------
   ! 荷電交換衝突: 中性粒子速度を背景イオン速度で置き換え
   !---------------------------------------------------------------------------
   subroutine collision_cx(p, vx_i, vy_i, vz_i, delta_E)
      type(particle_t), intent(inout) :: p
      real(dp), intent(in)  :: vx_i, vy_i, vz_i
      real(dp), intent(out) :: delta_E

      real(dp) :: E_old, E_new

      E_old = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)

      p%vx = vx_i
      p%vy = vy_i
      p%vz = vz_i

      E_new = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)
      delta_E = E_new - E_old
   end subroutine collision_cx

   !---------------------------------------------------------------------------
   ! 弾性散乱衝突: 重心系で散乱角を適用
   !---------------------------------------------------------------------------
   subroutine collision_el(p, vx_i, vy_i, vz_i, use_isotropic, delta_E)
      type(particle_t), intent(inout) :: p
      real(dp), intent(in)  :: vx_i, vy_i, vz_i
      logical, intent(in)   :: use_isotropic
      real(dp), intent(out) :: delta_E

      real(dp) :: E_old, E_new
      real(dp) :: vx_g, vy_g, vz_g  !重心速度
      real(dp) :: ux, uy, uz, u_mag  !相対速度
      real(dp) :: ux_n, uy_n, uz_n   !散乱後の相対速度
      real(dp) :: chi, phi           !散乱角
      real(dp) :: r_chi, E_rel

      E_old = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)

      !重心速度
      vx_g = 0.5d0 * (p%vx + vx_i)
      vy_g = 0.5d0 * (p%vy + vy_i)
      vz_g = 0.5d0 * (p%vz + vz_i)

      !相対速度
      ux = p%vx - vx_i
      uy = p%vy - vy_i
      uz = p%vz - vz_i
      u_mag = sqrt(ux*ux + uy*uy + uz*uz)

      if (u_mag < 1.0d-30) then
         delta_E = 0.0d0
         return
      end if

      !散乱角
      if (use_isotropic) then
         call random_number(r_chi)
         chi = acos(1.0d0 - 2.0d0 * r_chi)
      else
         E_rel = 0.25d0 * M_D_kg * u_mag * u_mag * J_TO_EV
         call random_number(r_chi)
         chi = sample_scattering_angle(E_rel, r_chi)
      end if

      call random_number(phi)
      phi = 2.0d0 * PI * phi

      !ベクトル回転
      call rotate_vector(ux, uy, uz, chi, phi, ux_n, uy_n, uz_n)

      !Lab系に戻す
      p%vx = vx_g + 0.5d0 * ux_n
      p%vy = vy_g + 0.5d0 * uy_n
      p%vz = vz_g + 0.5d0 * uz_n

      E_new = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)
      delta_E = E_new - E_old
   end subroutine collision_el

   !---------------------------------------------------------------------------
   ! ベクトル回転
   !---------------------------------------------------------------------------
   subroutine rotate_vector(ux, uy, uz, chi, phi, ux_n, uy_n, uz_n)
      real(dp), intent(in)  :: ux, uy, uz
      real(dp), intent(in)  :: chi, phi
      real(dp), intent(out) :: ux_n, uy_n, uz_n

      real(dp) :: u_mag, u_perp
      real(dp) :: cos_chi, sin_chi, cos_phi, sin_phi
      real(dp) :: nx, ny, nz

      u_mag = sqrt(ux*ux + uy*uy + uz*uz)
      if (u_mag < 1.0d-30) then
         ux_n = ux; uy_n = uy; uz_n = uz
         return
      end if

      cos_chi = cos(chi)
      sin_chi = sin(chi)
      cos_phi = cos(phi)
      sin_phi = sin(phi)

      nx = ux / u_mag
      ny = uy / u_mag
      nz = uz / u_mag

      u_perp = sqrt(nx*nx + ny*ny)

      if (u_perp > 1.0d-10) then
         ux_n = u_mag * (nx * cos_chi + &
            (nx * nz * cos_phi - ny * sin_phi) * sin_chi / u_perp)
         uy_n = u_mag * (ny * cos_chi + &
            (ny * nz * cos_phi + nx * sin_phi) * sin_chi / u_perp)
         uz_n = u_mag * (nz * cos_chi - u_perp * cos_phi * sin_chi)
      else
         ux_n = u_mag * sin_chi * cos_phi
         uy_n = u_mag * sin_chi * sin_phi
         uz_n = u_mag * cos_chi * sign(1.0d0, nz)
      end if
   end subroutine rotate_vector

end module dynamics
