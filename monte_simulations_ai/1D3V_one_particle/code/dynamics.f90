!===============================================================================
! Module: dynamics
! 粒子の推進と衝突判定
!===============================================================================
module dynamics
   use constants, only: dp, M_D, EV_TO_J, J_TO_EV, PI
   use types, only: particle, sim_params, plasma_params
   use cross_sections, only: get_all_cross_sections, sigma_cx, sigma_el, sigma_ei
   use cdf_reader, only: sample_scattering_angle
   implicit none

   ! 衝突タイプの定義
   integer, parameter :: COLL_NONE = 0
   integer, parameter :: COLL_CX = 1
   integer, parameter :: COLL_EL = 2
   integer, parameter :: COLL_EI = 3

contains

   !-----------------------------------------------------------------------------
   ! 粒子の推進
   !-----------------------------------------------------------------------------
   subroutine advance_particle(p, dt, params)
      type(particle), intent(inout) :: p
      real(dp), intent(in) :: dt
      type(sim_params), intent(in) :: params

      ! 位置の更新
      p%x = p%x + p%vx * dt

      ! 境界チェック（領域外に出たら消滅）
      if (p%x < params%x_min .or. p%x > params%x_max) then
         p%alive = .false.
      end if

   end subroutine advance_particle

   !-----------------------------------------------------------------------------
   ! 衝突判定
   !-----------------------------------------------------------------------------
   subroutine check_collision(p, plasma, dt, coll_type, delta_E)
      type(particle), intent(inout) :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in) :: dt
      integer, intent(out) :: coll_type
      real(dp), intent(out) :: delta_E

      real(dp) :: v, E_particle
      real(dp) :: sig_cx, sig_el, sig_ei, sig_total
      real(dp) :: nu_total, P_coll
      real(dp) :: rand_coll, rand_type, rand_angle, rand_phi
      real(dp) :: cum_prob

      coll_type = COLL_NONE
      delta_E = 0.0d0

      if (.not. p%alive) return

      ! 粒子速度と運動エネルギー
      v = sqrt(p%vx**2 + p%vy**2 + p%vz**2)
      if (v < 1.0d-10) return

      E_particle = 0.5d0 * M_D * v**2 * J_TO_EV

      ! 断面積の取得
      call get_all_cross_sections(E_particle, sig_cx, sig_el, sig_ei, sig_total)

      ! 衝突頻度
      nu_total = v * plasma%n_bg * sig_total

      ! 衝突確率
      P_coll = 1.0d0 - exp(-nu_total * dt)

      ! 衝突判定
      call random_number(rand_coll)
      if (rand_coll > P_coll) return  ! 衝突なし

      ! 衝突タイプの選択
      call random_number(rand_type)
      cum_prob = sig_cx / sig_total

      if (rand_type < cum_prob) then
         ! 荷電交換
         coll_type = COLL_CX
         call process_cx(p, plasma, delta_E)
      else
         cum_prob = cum_prob + sig_el / sig_total
         if (rand_type < cum_prob) then
            ! 弾性散乱
            coll_type = COLL_EL
            call random_number(rand_angle)
            call random_number(rand_phi)
            call process_elastic(p, plasma, E_particle, rand_angle, rand_phi, delta_E)
         else
            ! 電離
            coll_type = COLL_EI
            call process_ionization(p, delta_E)
         end if
      end if

   end subroutine check_collision

   !-----------------------------------------------------------------------------
   ! 荷電交換処理
   ! 速度が背景イオンの熱分布に入れ替わる
   !-----------------------------------------------------------------------------
   subroutine process_cx(p, plasma, delta_E)
      type(particle), intent(inout) :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(out) :: delta_E

      real(dp) :: E_before, E_after
      real(dp) :: v_th, vx_new, vy_new, vz_new
      real(dp) :: r1, r2, r3, r4

      ! 衝突前のエネルギー
      E_before = 0.5d0 * M_D * (p%vx**2 + p%vy**2 + p%vz**2) * J_TO_EV

      ! 熱速度
      v_th = sqrt(plasma%T_bg * EV_TO_J / M_D)

      ! Box-Muller法でガウス分布を生成
      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      call random_number(r4)

      vx_new = v_th * sqrt(-2.0d0 * log(max(r1, 1.0d-30))) * cos(2.0d0 * PI * r2)
      vy_new = v_th * sqrt(-2.0d0 * log(max(r1, 1.0d-30))) * sin(2.0d0 * PI * r2)
      vz_new = v_th * sqrt(-2.0d0 * log(max(r3, 1.0d-30))) * cos(2.0d0 * PI * r4)

      p%vx = vx_new
      p%vy = vy_new
      p%vz = vz_new

      ! 衝突後のエネルギー
      E_after = 0.5d0 * M_D * (p%vx**2 + p%vy**2 + p%vz**2) * J_TO_EV

      delta_E = E_after - E_before

   end subroutine process_cx

   !-----------------------------------------------------------------------------
   ! 弾性散乱処理
   ! 散乱角に基づいて速度方向を変更
   !-----------------------------------------------------------------------------
   subroutine process_elastic(p, plasma, E_collision, rand_p, rand_phi, delta_E)
      type(particle), intent(inout) :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in) :: E_collision
      real(dp), intent(in) :: rand_p, rand_phi
      real(dp), intent(out) :: delta_E

      real(dp) :: E_before, E_after
      real(dp) :: chi, phi
      real(dp) :: vx_cm, vy_cm, vz_cm  ! 重心系速度
      real(dp) :: vx_rel, vy_rel, vz_rel, v_rel  ! 相対速度
      real(dp) :: sin_chi, cos_chi, sin_phi, cos_phi
      real(dp) :: ex, ey, ez, e_perp
      real(dp) :: vx_new, vy_new, vz_new
      real(dp) :: v_bg_th, vx_bg, vy_bg, vz_bg
      real(dp) :: r1, r2, r3, r4

      E_before = 0.5d0 * M_D * (p%vx**2 + p%vy**2 + p%vz**2) * J_TO_EV

      ! 背景イオンの熱速度
      v_bg_th = sqrt(plasma%T_bg * EV_TO_J / M_D)

      ! 背景イオンの速度をサンプリング
      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      call random_number(r4)
      vx_bg = v_bg_th * sqrt(-2.0d0 * log(max(r1, 1.0d-30))) * cos(2.0d0 * PI * r2)
      vy_bg = v_bg_th * sqrt(-2.0d0 * log(max(r1, 1.0d-30))) * sin(2.0d0 * PI * r2)
      vz_bg = v_bg_th * sqrt(-2.0d0 * log(max(r3, 1.0d-30))) * cos(2.0d0 * PI * r4)

      ! 重心速度（等質量なので平均）
      vx_cm = 0.5d0 * (p%vx + vx_bg)
      vy_cm = 0.5d0 * (p%vy + vy_bg)
      vz_cm = 0.5d0 * (p%vz + vz_bg)

      ! 相対速度
      vx_rel = p%vx - vx_bg
      vy_rel = p%vy - vy_bg
      vz_rel = p%vz - vz_bg
      v_rel = sqrt(vx_rel**2 + vy_rel**2 + vz_rel**2)

      if (v_rel < 1.0d-30) then
         delta_E = 0.0d0
         return
      end if

      ! 散乱角のサンプリング
      chi = sample_scattering_angle(E_collision, rand_p)
      phi = 2.0d0 * PI * rand_phi

      sin_chi = sin(chi)
      cos_chi = cos(chi)
      sin_phi = sin(phi)
      cos_phi = cos(phi)

      ! 相対速度の単位ベクトル
      ex = vx_rel / v_rel
      ey = vy_rel / v_rel
      ez = vz_rel / v_rel

      ! 垂直方向の基底ベクトル
      if (abs(ez) > 0.9d0) then
         e_perp = sqrt(ex**2 + ey**2)
         if (e_perp > 1.0d-10) then
            ! ロドリゲスの回転公式を適用
            vx_new = v_rel * (sin_chi * (-ey * cos_phi / e_perp) + cos_chi * ex)
            vy_new = v_rel * (sin_chi * (ex * cos_phi / e_perp) + cos_chi * ey)
            vz_new = v_rel * (sin_chi * sin_phi + cos_chi * ez)
         else
            vx_new = v_rel * sin_chi * cos_phi
            vy_new = v_rel * sin_chi * sin_phi
            vz_new = v_rel * cos_chi
         end if
      else
         e_perp = sqrt(ex**2 + ez**2)
         if (e_perp > 1.0d-10) then
            vx_new = v_rel * (sin_chi * (-ez * cos_phi / e_perp) + cos_chi * ex)
            vy_new = v_rel * (sin_chi * sin_phi + cos_chi * ey)
            vz_new = v_rel * (sin_chi * (ex * cos_phi / e_perp) + cos_chi * ez)
         else
            vx_new = v_rel * sin_chi * cos_phi
            vy_new = v_rel * cos_chi
            vz_new = v_rel * sin_chi * sin_phi
         end if
      end if

      ! 実験室系に戻す
      p%vx = vx_cm + 0.5d0 * vx_new
      p%vy = vy_cm + 0.5d0 * vy_new
      p%vz = vz_cm + 0.5d0 * vz_new

      E_after = 0.5d0 * M_D * (p%vx**2 + p%vy**2 + p%vz**2) * J_TO_EV
      delta_E = E_after - E_before

   end subroutine process_elastic

   !-----------------------------------------------------------------------------
   ! 電離処理
   ! 粒子は電離後イオン化され追跡終了
   !-----------------------------------------------------------------------------
   subroutine process_ionization(p, delta_E)
      type(particle), intent(inout) :: p
      real(dp), intent(out) :: delta_E

      real(dp), parameter :: E_ionize = 13.6d0  ! 電離エネルギー [eV]

      ! 電離エネルギー分を失う
      delta_E = -E_ionize

      ! 粒子は電離されて消滅（イオン化）
      p%alive = .false.

   end subroutine process_ionization

   !-----------------------------------------------------------------------------
   ! Maxwell分布から速度を生成
   !-----------------------------------------------------------------------------
   subroutine sample_maxwell_velocity(T_eV, vx, vy, vz)
      real(dp), intent(in) :: T_eV  ! 温度 [eV]
      real(dp), intent(out) :: vx, vy, vz

      real(dp) :: v_th, r1, r2, r3, r4

      v_th = sqrt(T_eV * EV_TO_J / M_D)

      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      call random_number(r4)

      vx = v_th * sqrt(-2.0d0 * log(max(r1, 1.0d-30))) * cos(2.0d0 * PI * r2)
      vy = v_th * sqrt(-2.0d0 * log(max(r1, 1.0d-30))) * sin(2.0d0 * PI * r2)
      vz = v_th * sqrt(-2.0d0 * log(max(r3, 1.0d-30))) * cos(2.0d0 * PI * r4)

   end subroutine sample_maxwell_velocity

   !-----------------------------------------------------------------------------
   ! 単一方向の速度を指定エネルギーで生成
   !-----------------------------------------------------------------------------
   subroutine set_mono_velocity(E_eV, direction, vx, vy, vz)
      real(dp), intent(in) :: E_eV       ! エネルギー [eV]
      integer, intent(in) :: direction   ! 1: +x, -1: -x
      real(dp), intent(out) :: vx, vy, vz

      real(dp) :: v_mag

      v_mag = sqrt(2.0d0 * E_eV * EV_TO_J / M_D)

      vx = v_mag * direction
      vy = 0.0d0
      vz = 0.0d0

   end subroutine set_mono_velocity

end module dynamics
