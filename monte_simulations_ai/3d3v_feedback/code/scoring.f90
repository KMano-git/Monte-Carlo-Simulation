!===============================================================================
! Module: scoring
! Collision Estimator (CL) と Track-Length Estimator (TL) の実装
!===============================================================================
module scoring
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV, E_IONIZE_THRESHOLD
   use data_types, only: particle_t, plasma_params, score_data
   use cross_sections, only: sigma_cx, sigma_el, ionization_rate_coeff
   use random_utils, only: sample_maxwell_velocity_ion, random_double
   use cdf_reader, only: sample_scattering_angle
   implicit none

   private
   public :: score_collision_estimator, score_track_length_estimator

contains

   !---------------------------------------------------------------------------
   ! Collision Estimator (CL)
   ! 実衝突（Rejection通過後）時のみ呼び出す
   ! s_CL = w_k * (R_a/R_s * s_a + s_s)
   !---------------------------------------------------------------------------
   subroutine score_collision_estimator(p, plasma, vx_i, vy_i, vz_i, &
      v_rel, E_rel, coll_type, delta_E_el, enable_ei, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)  :: vx_i, vy_i, vz_i !背景イオン速度 [m/s]
      real(dp), intent(in)  :: v_rel       !相対速度 [m/s]
      real(dp), intent(in)  :: E_rel       !相対エネルギー [eV]
      integer, intent(in)   :: coll_type   !衝突タイプ (1=CX, 2=EL)
      real(dp), intent(in)  :: delta_E_el  !弾性散乱のエネルギー変化 [J]
      logical, intent(in)   :: enable_ei   !電離有効フラグ
      type(score_data), intent(inout) :: score

      real(dp) :: E_kin          !中性粒子の運動エネルギー [J]
      real(dp) :: R_cx, R_el, R_s, R_a
      real(dp) :: s_a, s_cx, s_el
      real(dp) :: u_flow2
      real(dp) :: sig_cx_val, sig_el_val

      !粒子の運動エネルギー
      E_kin = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)

      !断面積
      sig_cx_val = sigma_cx(E_rel)
      sig_el_val = sigma_el(E_rel)

      !レート計算
      R_cx = plasma%n_i * sig_cx_val * v_rel
      R_el = plasma%n_i * sig_el_val * v_rel
      R_s  = R_cx + R_el

      !電離レート（enable_ei=.false.のとき0）
      if (enable_ei) then
         R_a = plasma%n_e * ionization_rate_coeff(plasma%T_e)
      else
         R_a = 0.0d0
      end if

      !背景イオン流速の2乗
      u_flow2 = plasma%u_x**2 + plasma%u_y**2 + plasma%u_z**2

      !各プロセスのスコア [J]
      !電離: 中性粒子の運動エネルギーがイオンに渡る
      s_a = E_kin

      !荷電交換: E_kin - E_ion を[J]で計算 (E_ionはサンプリングされたイオンのエネルギー)
      s_cx = E_kin - 0.5d0 * M_D_kg * (vx_i**2 + vy_i**2 + vz_i**2)

      !弾性散乱: -ΔE_EL
      s_el = -delta_E_el

      !CL公式: s_CL = w_k * (R_a/R_s * s_a + s_s)
      !s_sは実際に起きた衝突タイプに対応するスコア
      if (R_s > 1.0d-30) then
         !電離の寄与（R_a/R_sで補正）
         score%cl_ei = score%cl_ei + p%weight * (R_a / R_s) * s_a

         !散乱の寄与（実衝突タイプのスコア）
         if (coll_type == 1) then  !CX
            score%cl_cx = score%cl_cx + p%weight * s_cx
         else if (coll_type == 2) then  !EL
            score%cl_el = score%cl_el + p%weight * s_el
         end if
      end if

   end subroutine score_collision_estimator

   !---------------------------------------------------------------------------
   ! Track-Length Estimator (TL)
   ! 毎タイムステップ呼び出す（衝突の有無によらず）
   ! 仮想サンプリング（Dummy Collision）を内部で実行
   ! s_TL,j = s_c,j * (w_k / R_a) * (1 - exp(-R_a * dt))
   !---------------------------------------------------------------------------
   subroutine score_track_length_estimator(p, plasma, dt, &
      use_isotropic, enable_ei, score)
      type(particle_t), intent(inout)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)  :: dt           !タイムステップ [s]
      logical, intent(in)   :: use_isotropic !等方散乱フラグ
      logical, intent(in)   :: enable_ei    !電離有効フラグ
      type(score_data), intent(inout) :: score

      real(dp) :: E_kin      !中性粒子の運動エネルギー [J]
      real(dp) :: R_a, R_cx_d, R_el_d
      real(dp) :: s_c_ei, s_c_cx, s_c_el, s_c_ei_e, s_c_ei_n
      real(dp) :: eff_time   !(w_k / R_a) * (1 - exp(-R_a * dt))
      real(dp) :: u_flow2
      real(dp) :: delta_E_el_dummy  !仮想ELエネルギー変化 [J]

      !仮想サンプリング用変数
      real(dp) :: vx_d, vy_d, vz_d  !仮想イオン速度
      real(dp) :: v_rel_d, E_rel_d  !仮想相対速度・エネルギー
      real(dp) :: sig_cx_d, sig_el_d
      real(dp) :: chi_d, phi_d, r_chi, r_phi  !仮想散乱角
      real(dp) :: vx_g, vy_g, vz_g  !重心速度
      real(dp) :: ux, uy, uz, u_mag  !相対速度ベクトル
      real(dp) :: ux_n, uy_n, uz_n   !散乱後の相対速度
      real(dp) :: vx_new, vy_new, vz_new  !仮想散乱後の速度
      real(dp) :: E_old, E_new

      !粒子の運動エネルギー
      E_kin = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)

      !--- 仮想サンプリング（Dummy Collision） ---
      !1. 仮想イオン速度をサンプリング
      call sample_maxwell_velocity_ion(p%rng, plasma%T_i, plasma, vx_d, vy_d, vz_d)

      !2. 仮想相対速度
      v_rel_d = sqrt((p%vx - vx_d)**2 + (p%vy - vy_d)**2 + (p%vz - vz_d)**2)

      !相対エネルギー（重心系）
      E_rel_d = 0.25d0 * M_D_kg * v_rel_d * v_rel_d * J_TO_EV

      !3. 仮想断面積
      sig_cx_d = sigma_cx(E_rel_d)
      sig_el_d = sigma_el(E_rel_d)

      !4. 仮想散乱角のサンプリングとΔE_ELの計算
      !   重心系で散乱角を適用し、エネルギー変化を計算
      !   粒子速度は変更しない
      E_old = E_kin  !散乱前のエネルギー

      !重心速度
      vx_g = 0.5d0 * (p%vx + vx_d)
      vy_g = 0.5d0 * (p%vy + vy_d)
      vz_g = 0.5d0 * (p%vz + vz_d)

      !相対速度ベクトル
      ux = p%vx - vx_d
      uy = p%vy - vy_d
      uz = p%vz - vz_d
      u_mag = sqrt(ux*ux + uy*uy + uz*uz)

      delta_E_el_dummy = 0.0d0
      if (u_mag > 1.0d-30) then
         !散乱角サンプリング
         if (use_isotropic) then
            r_chi = random_double(p%rng)
            chi_d = acos(1.0d0 - 2.0d0 * r_chi)
         else
            r_chi = random_double(p%rng)
            chi_d = sample_scattering_angle(E_rel_d, r_chi)
         end if

         r_phi = random_double(p%rng)
         phi_d = 2.0d0 * 3.141592653589793d0 * r_phi

         !相対速度を散乱角で回転
         call rotate_vector_tl(ux, uy, uz, chi_d, phi_d, ux_n, uy_n, uz_n)

         !仮想散乱後の速度（lab系）
         vx_new = vx_g + 0.5d0 * ux_n
         vy_new = vy_g + 0.5d0 * uy_n
         vz_new = vz_g + 0.5d0 * uz_n

         !仮想散乱後のエネルギー
         E_new = 0.5d0 * M_D_kg * (vx_new**2 + vy_new**2 + vz_new**2)

         delta_E_el_dummy = E_new - E_old
      end if

      !--- レート計算 ---
      if (enable_ei) then
         R_a = plasma%n_e * ionization_rate_coeff(plasma%T_e)
      else
         R_a = 0.0d0
      end if
      R_cx_d = plasma%n_i * sig_cx_d * v_rel_d
      R_el_d = plasma%n_i * sig_el_d * v_rel_d

      !背景イオン流速の2乗
      u_flow2 = plasma%u_x**2 + plasma%u_y**2 + plasma%u_z**2

      !--- 各プロセスの s_c（単位時間あたりの期待スコア）[J/s] ---
      s_c_ei = R_a * E_kin
      s_c_cx = R_cx_d * (E_kin - 0.5d0 * M_D_kg * (vx_d**2 + vy_d**2 + vz_d**2))
      s_c_el = R_el_d * (-delta_E_el_dummy)
      s_c_ei_e = R_a * (-E_IONIZE_THRESHOLD * EV_TO_J)
      s_c_ei_n = R_a

      !--- 有効飛行時間の計算 ---
      !eff_time = (w_k / R_a) * (1 - exp(-R_a * dt))
      !R_a ≈ 0 の場合: eff_time → w_k * dt
      if (R_a > 1.0d-30) then
         eff_time = (p%weight / R_a) * (1.0d0 - exp(-R_a * dt))
      else
         eff_time = p%weight * dt
      end if

      !--- スコア蓄積 [J] ---
      score%tl_ei = score%tl_ei + s_c_ei * eff_time
      score%tl_cx = score%tl_cx + s_c_cx * eff_time
      score%tl_el = score%tl_el + s_c_el * eff_time
      score%tl_ei_e = score%tl_ei_e + s_c_ei_e * eff_time
      score%tl_ei_n = score%tl_ei_n + s_c_ei_n * eff_time

   end subroutine score_track_length_estimator

   !---------------------------------------------------------------------------
   ! ベクトル回転（TLスコアリング内部用）
   !---------------------------------------------------------------------------
   subroutine rotate_vector_tl(ux, uy, uz, chi, phi, ux_n, uy_n, uz_n)
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

      !正規化された入射方向
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
         !z軸に平行な場合の特殊処理
         ux_n = u_mag * sin_chi * cos_phi
         uy_n = u_mag * sin_chi * sin_phi
         uz_n = u_mag * cos_chi * sign(1.0d0, nz)
      end if

   end subroutine rotate_vector_tl

end module scoring
