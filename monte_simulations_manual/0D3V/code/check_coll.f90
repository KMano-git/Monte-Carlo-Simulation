!===============================================================================
! Module: check_collision
! 衝突判定
!===============================================================================
module check_coll
   use constants, only: dp, M_D_kg, PI, E_IONIZE_THRESHOLD, EV_TO_J, J_TO_EV
   use particle_data, only: particle_t, sim_params, plasma_params
   use random_utils, only: sample_maxwell_velocity_ion
   use cross_sections, only: get_all_cross_sections
   implicit none

   ! 衝突タイプの定義
   integer, parameter :: COLL_NONE = 0
   integer, parameter :: COLL_CX = 1
   integer, parameter :: COLL_EL = 2
   integer, parameter :: COLL_EI = 3

contains

   !-----------------------------------------------------------------------------
   ! 衝突判定(Null Collision Method)
   !-----------------------------------------------------------------------------
   subroutine check_collision(p, plasma, sim, dt, vx_i, vy_i, vz_i, coll_type)
      type(particle_t), intent(in) :: p
      type(plasma_params), intent(in) :: plasma
      type(sim_params), intent(in) :: sim
      real(dp), intent(in) :: dt
      real(dp), intent(out) :: vx_i, vy_i, vz_i
      integer, intent(out) :: coll_type

      real(dp) :: v_mag, E_particle
      real(dp) :: v_th_i, v_rel_max, E_rel_max, sigma_max
      real(dp) :: sigma_cx, sigma_el, sigma_ei, sigma_total
      real(dp) :: nu_cx, nu_el, nu_ei, nu_total
      real(dp) :: nu_max, P_coll, P_accept
      real(dp) :: v_rel, E_rel ! サンプリング後の相対速度と相対エネルギー
      real(dp) :: r, r_accept, r_type

      ! 衝突タイプ初期化
      coll_type = COLL_NONE

      ! 粒子の速度とエネルギー
      v_mag = sqrt(p%vx*p%vx + p%vy*p%vy + p%vz*p%vz)
      if (v_mag < 1.0d-20) return
      E_particle = 0.5d0 * M_D_kg * (p%vx*p%vx + p%vy*p%vy + p%vz*p%vz) * J_TO_EV

      ! 相対速度の最大値を推定する
      v_th_i = sqrt(plasma%T_i * EV_TO_J / M_D_kg)
      v_rel_max = v_mag + sqrt(plasma%u_x**2 + plasma%u_y**2 + plasma%u_z**2) + 20.0d0 * v_th_i
      E_rel_max = 0.25d0 * M_D_kg * v_rel_max * v_rel_max * J_TO_EV
      ! v_rel_max における断面積の最大値
      call get_all_cross_sections(E_rel_max, sigma_cx, sigma_el, sigma_ei, sigma_max)

      ! 最大衝突頻度
      nu_max = plasma%n_i * sigma_max * v_rel_max

      ! 衝突確率
      P_coll = 1.0d0 - exp(-nu_max * dt)

      ! 衝突判定
      call random_number(r)
      if (r > P_coll) return

      ! 背景粒子のサンプリング
      call sample_maxwell_velocity_ion(plasma%T_i, plasma, vx_i, vy_i, vz_i)

      ! 相対速度
      v_rel = sqrt((p%vx-vx_i)**2 + (p%vy-vy_i)**2 + (p%vz-vz_i)**2)

      E_rel = 0.25d0 * M_D_kg * v_rel*v_rel * J_TO_EV

      call get_all_cross_sections(E_rel, sigma_cx, sigma_el, sigma_ei, sigma_total)

      ! Rejection判定
      ! sigma(v_rel)*v_rel / (sigma_max * v_rel_max)

      ! 安全策: v_rel > v_rel_max の場合は確率1で採用（ただし頻度が過小評価されていることになる）
      if (v_rel > v_rel_max) then
         P_accept = 1.0d0
      else
         P_accept = sigma_total*v_rel / (sigma_max * v_rel_max)
      end if

      call random_number(r_accept)
      if (r_accept > P_accept) return

      call random_number(r_type)

      if (r_type < sigma_cx / sigma_total .and. sim%enable_cx) then
         coll_type = COLL_CX
      else if (r_type < (sigma_cx + sigma_el) / sigma_total .and. sim%enable_el) then
         coll_type = COLL_EL
      else if (sim%enable_ei) then
         coll_type = COLL_EI
      end if

   end subroutine check_collision

end module check_coll
