!===============================================================================
! Module: scoring
! Elastic analog / collision / track-length scoring, plus comparison variants
!===============================================================================
module scoring
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV, PI
   use data_types, only: particle_t, plasma_params, score_data, rng_state
   use cross_sections, only: sigma_cx, sigma_el, ionization_rate_coeff
   use random_utils, only: sample_maxwell_velocity_ion, random_double
   use cdf_reader, only: sample_scattering_angle, get_reaction_rate, get_I_1_0, &
      get_I_1_1_up, get_I_1_2_up2
   use tl_el_table_reader, only: get_tl_el_rate_from_table
   implicit none

   real(dp), parameter :: TLMT_EL = 0.9d0
   real(dp), parameter :: SCORE_TINY = 1.0d-30

   private
   public :: score_analog_elastic, score_collision_estimator
   public :: score_inner_multi_collision_el, score_track_length_estimator
   public :: score_inner_multi_track_length, score_pretabulated_track_length
   public :: score_table_lookup_track_length
   public :: tl_el_sample_once, tl_el_sample_avg, tl_el_sample_stats

contains

   !---------------------------------------------------------------------------
   ! Analog estimator for elastic collisions
   ! 実際に起きた EL 衝突の実現値をそのまま score する
   !---------------------------------------------------------------------------
   subroutine score_analog_elastic(p, delta_E_el, score)
      type(particle_t), intent(in)    :: p
      real(dp), intent(in)            :: delta_E_el
      type(score_data), intent(inout) :: score

      score%a_el = score%a_el + p%weight * (-delta_E_el)
   end subroutine score_analog_elastic

   !---------------------------------------------------------------------------
   ! Collision Estimator (CL)
   ! アクセプトされた scatter collision ごとに、
   ! CX/EL の期待寄与を両方 score する
   !---------------------------------------------------------------------------
   subroutine score_collision_estimator(p, plasma, vx_i, vy_i, vz_i, &
      v_rel, E_rel, enable_cx, enable_el, enable_ionization, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: vx_i, vy_i, vz_i
      real(dp), intent(in)            :: v_rel
      real(dp), intent(in)            :: E_rel
      logical, intent(in)             :: enable_cx
      logical, intent(in)             :: enable_el
      logical, intent(in)             :: enable_ionization
      type(score_data), intent(inout) :: score

      real(dp) :: E_kin
      real(dp) :: R_cx, R_el, R_s, R_a
      real(dp) :: s_a, s_cx
      real(dp) :: sig_cx_val, sig_el_val
      real(dp) :: s_el_lookup, rate_el_dummy, reaction_rate_dummy

      E_kin = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)

      sig_cx_val = sigma_cx(E_rel)
      sig_el_val = sigma_el(E_rel)

      if (enable_cx) then
         R_cx = plasma%ion_density * sig_cx_val * v_rel
      else
         R_cx = 0.0d0
      end if

      if (enable_el) then
         R_el = plasma%ion_density * sig_el_val * v_rel
      else
         R_el = 0.0d0
      end if
      R_s = R_cx + R_el

      if (enable_ionization) then
         R_a = plasma%electron_density * &
            ionization_rate_coeff(plasma%electron_temperature_eV)
      else
         R_a = 0.0d0
      end if

      s_a = E_kin
      s_cx = E_kin - 0.5d0 * M_D_kg * (vx_i**2 + vy_i**2 + vz_i**2)

      if (R_s > SCORE_TINY) then
         score%cl_ei = score%cl_ei + p%weight * (R_a / R_s) * s_a
         score%cl_cx = score%cl_cx + p%weight * (R_cx / R_s) * s_cx
         if (enable_el) then
            call compute_el_lookup_rate_and_score(p, plasma, rate_el_dummy, &
               reaction_rate_dummy, s_el_lookup)
            score%cl_el = score%cl_el + p%weight * (R_el / R_s) * s_el_lookup
         end if
      end if
   end subroutine score_collision_estimator

   !---------------------------------------------------------------------------
   ! EL inner multi-sampling CL
   ! 背景イオン速度と散乱角を内側サンプリングして、1 collision あたり期待寄与を作る
   !---------------------------------------------------------------------------
   subroutine score_inner_multi_collision_el(p, plasma, E_rel, enable_cx, enable_el, &
      use_isotropic, n_inner, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: E_rel
      logical, intent(in)             :: enable_cx, enable_el
      logical, intent(in)             :: use_isotropic
      integer, intent(in)             :: n_inner
      type(score_data), intent(inout) :: score

      real(dp) :: s_el_avg, sig_cx_val, sig_el_val, weight_el

      if (.not. enable_el) return

      if (enable_cx) then
         sig_cx_val = sigma_cx(E_rel)
      else
         sig_cx_val = 0.0d0
      end if
      sig_el_val = sigma_el(E_rel)

      s_el_avg = cl_el_sample_avg(p, plasma, use_isotropic, n_inner)
      if (sig_cx_val + sig_el_val > SCORE_TINY) then
         weight_el = sig_el_val / (sig_cx_val + sig_el_val)
      else
         weight_el = 0.0d0
      end if

      score%cl_inner_el = score%cl_inner_el + p%weight * weight_el * s_el_avg
   end subroutine score_inner_multi_collision_el

   !---------------------------------------------------------------------------
   ! Track-Length Estimator (TR, lookup mainline)
   ! 実衝突までに蓄積した pending_eff_time に、lookup ベースの期待 rate を掛ける
   !---------------------------------------------------------------------------
   subroutine score_track_length_estimator(p, plasma, pending_eff_time, &
      enable_cx, enable_el, enable_ionization, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: pending_eff_time
      logical, intent(in)             :: enable_cx, enable_el, enable_ionization
      type(score_data), intent(inout) :: score

      real(dp) :: E_kin
      real(dp) :: s_c_ei, s_c_cx, s_c_el
      real(dp) :: vx_d, vy_d, vz_d
      real(dp) :: v_rel_d, E_rel_cx, sig_cx_d
      real(dp) :: rate_el_dummy, collision_score_dummy
      type(rng_state) :: rng_work

      if (pending_eff_time <= 0.0d0) return

      E_kin = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)

      if (enable_ionization) then
         s_c_ei = plasma%electron_density * &
            ionization_rate_coeff(plasma%electron_temperature_eV) * E_kin
      else
         s_c_ei = 0.0d0
      end if

      s_c_cx = 0.0d0
      if (enable_cx) then
         rng_work = p%rng
         call sample_maxwell_velocity_ion(rng_work, plasma%ion_temperature_eV, &
            plasma, vx_d, vy_d, vz_d)
         v_rel_d = sqrt((p%vx - vx_d)**2 + (p%vy - vy_d)**2 + (p%vz - vz_d)**2)
         E_rel_cx = 0.25d0 * M_D_kg * v_rel_d * v_rel_d * J_TO_EV
         sig_cx_d = sigma_cx(E_rel_cx)
         s_c_cx = plasma%ion_density * sig_cx_d * v_rel_d * &
            (E_kin - 0.5d0 * M_D_kg * (vx_d**2 + vy_d**2 + vz_d**2))
      end if

      s_c_el = 0.0d0
      if (enable_el) then
         ! Lookup TR is angle-averaged already, so it does not depend on
         ! the runtime elastic angle sampling model.
         call compute_el_lookup_rate_and_score(p, plasma, s_c_el, rate_el_dummy, &
            collision_score_dummy)
      end if

      score%tr_ei = score%tr_ei + s_c_ei * pending_eff_time
      score%tr_cx = score%tr_cx + s_c_cx * pending_eff_time
      score%tr_el = score%tr_el + s_c_el * pending_eff_time
   end subroutine score_track_length_estimator

   !---------------------------------------------------------------------------
   ! Backward-compatible wrapper
   !---------------------------------------------------------------------------
   subroutine score_table_lookup_track_length(p, plasma, pending_eff_time, &
      enable_cx, enable_el, enable_ionization, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: pending_eff_time
      logical, intent(in)             :: enable_cx, enable_el, enable_ionization
      type(score_data), intent(inout) :: score

      call score_track_length_estimator(p, plasma, pending_eff_time, &
         enable_cx, enable_el, enable_ionization, score)
   end subroutine score_table_lookup_track_length

   !---------------------------------------------------------------------------
   ! EL 1サンプル分の rate を、既にサンプリング済みの背景イオン速度から評価する
   ! optional: reaction_rate_el = R_el, collision_score_el = -ΔE_actual
   !---------------------------------------------------------------------------
   function tl_el_rate_from_ion_sample(p, plasma, use_isotropic, rng_work, &
      vx_d, vy_d, vz_d, reaction_rate_el, collision_score_el) result(s_c_el)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      logical, intent(in)             :: use_isotropic
      type(rng_state), intent(inout)  :: rng_work
      real(dp), intent(in)            :: vx_d, vy_d, vz_d
      real(dp), intent(out), optional :: reaction_rate_el, collision_score_el
      real(dp) :: s_c_el

      real(dp) :: v_rel_d, E_rel_d
      real(dp) :: sig_el_d, R_el_d
      real(dp) :: chi_d, phi_d, r_chi, r_phi
      real(dp) :: vx_g, vy_g, vz_g
      real(dp) :: ux, uy, uz, u_mag
      real(dp) :: ux_n, uy_n, uz_n
      real(dp) :: vx_new, vy_new, vz_new
      real(dp) :: E_old, E_new, delta_E_el_dummy

      v_rel_d = sqrt((p%vx - vx_d)**2 + (p%vy - vy_d)**2 + (p%vz - vz_d)**2)
      E_rel_d = 0.25d0 * M_D_kg * v_rel_d * v_rel_d * J_TO_EV
      sig_el_d = sigma_el(E_rel_d)
      R_el_d = plasma%ion_density * sig_el_d * v_rel_d

      E_old = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)
      vx_g = 0.5d0 * (p%vx + vx_d)
      vy_g = 0.5d0 * (p%vy + vy_d)
      vz_g = 0.5d0 * (p%vz + vz_d)

      ux = p%vx - vx_d
      uy = p%vy - vy_d
      uz = p%vz - vz_d
      u_mag = sqrt(ux*ux + uy*uy + uz*uz)

      delta_E_el_dummy = 0.0d0
      if (u_mag > SCORE_TINY) then
         r_chi = random_double(rng_work)
         if (use_isotropic) then
            chi_d = acos(1.0d0 - 2.0d0 * r_chi)
         else
               chi_d = sample_scattering_angle(E_rel_d, r_chi)
         end if

         r_phi = random_double(rng_work)
         phi_d = 2.0d0 * PI * r_phi

         call rotate_vector_tl(ux, uy, uz, chi_d, phi_d, ux_n, uy_n, uz_n)

         vx_new = vx_g + 0.5d0 * ux_n
         vy_new = vy_g + 0.5d0 * uy_n
         vz_new = vz_g + 0.5d0 * uz_n
         E_new = 0.5d0 * M_D_kg * (vx_new**2 + vy_new**2 + vz_new**2)
         delta_E_el_dummy = E_new - E_old
      end if

      if (present(reaction_rate_el)) reaction_rate_el = R_el_d
      if (present(collision_score_el)) collision_score_el = -delta_E_el_dummy

      s_c_el = R_el_d * (-delta_E_el_dummy)
   end function tl_el_rate_from_ion_sample

   !---------------------------------------------------------------------------
   ! EL inner multi-sampling 1サンプル版
   !---------------------------------------------------------------------------
   function tl_el_sample_once(p, plasma, use_isotropic) result(s_c_el)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      logical, intent(in)             :: use_isotropic
      real(dp) :: s_c_el

      type(rng_state) :: rng_work

      rng_work = p%rng
      s_c_el = tl_el_sample_once_with_rng(p, plasma, use_isotropic, rng_work)
   end function tl_el_sample_once

   !---------------------------------------------------------------------------
   ! EL inner multi-sampling 平均版（TR 用）
   !---------------------------------------------------------------------------
   function tl_el_sample_avg(p, plasma, use_isotropic, n_inner) result(s_c_el_avg)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      logical, intent(in)             :: use_isotropic
      integer, intent(in)             :: n_inner
      real(dp) :: s_c_el_avg

      real(dp) :: variance_dummy, stddev_dummy, stderr_dummy

      call tl_el_sample_stats(p, plasma, use_isotropic, n_inner, s_c_el_avg, &
         variance_dummy, stddev_dummy, stderr_dummy)
   end function tl_el_sample_avg

   subroutine tl_el_sample_stats(p, plasma, use_isotropic, n_samples_in, mean_rate, &
      variance_rate, stddev_rate, stderr_rate)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      logical, intent(in)             :: use_isotropic
      integer, intent(in)             :: n_samples_in
      real(dp), intent(out)           :: mean_rate
      real(dp), intent(out)           :: variance_rate
      real(dp), intent(out)           :: stddev_rate
      real(dp), intent(out)           :: stderr_rate

      integer :: m, n_samples
      real(dp) :: sample_rate, delta, delta2, m2
      type(rng_state) :: rng_work

      n_samples = max(1, n_samples_in)
      rng_work = p%rng
      mean_rate = 0.0d0
      m2 = 0.0d0

      do m = 1, n_samples
         sample_rate = tl_el_sample_once_with_rng(p, plasma, use_isotropic, rng_work)
         delta = sample_rate - mean_rate
         mean_rate = mean_rate + delta / real(m, dp)
         delta2 = sample_rate - mean_rate
         m2 = m2 + delta * delta2
      end do

      if (n_samples >= 2) then
         variance_rate = m2 / real(n_samples - 1, dp)
      else
         variance_rate = 0.0d0
      end if
      variance_rate = max(0.0d0, variance_rate)
      stddev_rate = sqrt(variance_rate)
      stderr_rate = stddev_rate / sqrt(real(n_samples, dp))
   end subroutine tl_el_sample_stats

   !---------------------------------------------------------------------------
   ! EL inner multi-sampling 平均版（CL 用）
   ! mean(R_el * score) / mean(R_el)
   !---------------------------------------------------------------------------
   function cl_el_sample_avg(p, plasma, use_isotropic, n_inner) result(s_el_avg)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      logical, intent(in)             :: use_isotropic
      integer, intent(in)             :: n_inner
      real(dp) :: s_el_avg

      integer :: m, n_samples
      real(dp) :: vx_d, vy_d, vz_d
      real(dp) :: sample_rate, sample_reaction_rate
      real(dp) :: sum_rate, sum_reaction_rate, score_dummy
      type(rng_state) :: rng_work

      n_samples = max(1, n_inner)
      rng_work = p%rng
      sum_rate = 0.0d0
      sum_reaction_rate = 0.0d0

      do m = 1, n_samples
         call sample_maxwell_velocity_ion(rng_work, plasma%ion_temperature_eV, &
            plasma, vx_d, vy_d, vz_d)
         sample_rate = tl_el_rate_from_ion_sample(p, plasma, use_isotropic, rng_work, &
            vx_d, vy_d, vz_d, sample_reaction_rate, score_dummy)
         sum_rate = sum_rate + sample_rate
         sum_reaction_rate = sum_reaction_rate + sample_reaction_rate
      end do

      if (sum_reaction_rate > SCORE_TINY) then
         s_el_avg = sum_rate / sum_reaction_rate
      else
         s_el_avg = 0.0d0
      end if
   end function cl_el_sample_avg

   function tl_el_sample_once_with_rng(p, plasma, use_isotropic, rng_work) result(s_c_el)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      logical, intent(in)             :: use_isotropic
      type(rng_state), intent(inout)  :: rng_work
      real(dp) :: s_c_el

      real(dp) :: vx_d, vy_d, vz_d

      call sample_maxwell_velocity_ion(rng_work, plasma%ion_temperature_eV, &
         plasma, vx_d, vy_d, vz_d)
      s_c_el = tl_el_rate_from_ion_sample(p, plasma, use_isotropic, rng_work, &
         vx_d, vy_d, vz_d)
   end function tl_el_sample_once_with_rng

   !---------------------------------------------------------------------------
   ! ベクトル回転（EL scoring 内部用）
   !---------------------------------------------------------------------------
   subroutine rotate_vector_tl(ux, uy, uz, chi, phi, ux_n, uy_n, uz_n)
      real(dp), intent(in)  :: ux, uy, uz
      real(dp), intent(in)  :: chi, phi
      real(dp), intent(out) :: ux_n, uy_n, uz_n

      real(dp) :: u_mag, u_perp
      real(dp) :: cos_chi, sin_chi, cos_phi, sin_phi
      real(dp) :: nx, ny, nz

      u_mag = sqrt(ux*ux + uy*uy + uz*uz)
      if (u_mag < SCORE_TINY) then
         ux_n = ux
         uy_n = uy
         uz_n = uz
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
   end subroutine rotate_vector_tl

   !---------------------------------------------------------------------------
   ! EL inner multi-sampling TR
   !---------------------------------------------------------------------------
   subroutine score_inner_multi_track_length(p, plasma, pending_eff_time, &
      use_isotropic, enable_el, n_inner, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: pending_eff_time
      logical, intent(in)             :: use_isotropic, enable_el
      integer, intent(in)             :: n_inner
      type(score_data), intent(inout) :: score

      real(dp) :: s_c_el

      if (pending_eff_time <= 0.0d0 .or. .not. enable_el) return

      s_c_el = tl_el_sample_avg(p, plasma, use_isotropic, n_inner)
      score%tr_inner_el = score%tr_inner_el + s_c_el * pending_eff_time
   end subroutine score_inner_multi_track_length

   !---------------------------------------------------------------------------
   ! EL pre-tabulated TR
   !---------------------------------------------------------------------------
   subroutine score_pretabulated_track_length(p, plasma, pending_eff_time, &
      enable_el, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: pending_eff_time
      logical, intent(in)             :: enable_el
      type(score_data), intent(inout) :: score

      real(dp) :: etm, tim, s_c_el

      if (pending_eff_time <= 0.0d0 .or. .not. enable_el) return

      call compute_el_lookup_coords(p, plasma, etm, tim)
      s_c_el = get_tl_el_rate_from_table(etm, tim)
      score%tr_pretab_el = score%tr_pretab_el + s_c_el * pending_eff_time
   end subroutine score_pretabulated_track_length

   pure subroutine compute_el_lookup_coords(p, plasma, etm, tim)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(out)           :: etm, tim

      real(dp) :: utx, uty, utz, ut2
      real(dp) :: m_ref, m_i

      utx = p%vx - plasma%ion_flow_vx
      uty = p%vy - plasma%ion_flow_vy
      utz = p%vz - plasma%ion_flow_vz
      ut2 = utx**2 + uty**2 + utz**2

      m_ref = 0.5d0 * M_D_kg
      m_i = M_D_kg

      ! DEGAS/EIRENE の Etm 軸は specific_energy [eV/amu] で、
      ! D + D+ では reduced mass = 1 amu のため数値的に E_rel [eV] と一致する。
      etm = 0.5d0 * m_ref * ut2 * J_TO_EV
      tim = (m_ref / m_i) * plasma%ion_temperature_eV
   end subroutine compute_el_lookup_coords

   !---------------------------------------------------------------------------
   ! EL lookup kernel
   ! el_rate          : unit-time expected source [J/s]
   ! el_reaction_rate : EL collision frequency [1/s]
   ! el_collision_score : expected source per EL collision [J]
   !---------------------------------------------------------------------------
   subroutine compute_el_lookup_rate_and_score(p, plasma, el_rate, el_reaction_rate, &
      el_collision_score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(out)           :: el_rate
      real(dp), intent(out)           :: el_reaction_rate
      real(dp), intent(out)           :: el_collision_score

      real(dp) :: utx, uty, utz, ut2, ut, zuv, va2
      real(dp) :: etm, tim, tim2
      real(dp) :: val_I_1_0, val_I_1_1, val_I_1_2
      real(dp) :: sp0, dtr_eng
      real(dp) :: m_ref, m_t, m_i, ame

      utx = p%vx - plasma%ion_flow_vx
      uty = p%vy - plasma%ion_flow_vy
      utz = p%vz - plasma%ion_flow_vz
      ut2 = utx**2 + uty**2 + utz**2
      ut = sqrt(ut2)
      zuv = utx * p%vx + uty * p%vy + utz * p%vz

      m_t = M_D_kg
      m_i = M_D_kg
      m_ref = 0.5d0 * M_D_kg
      ame = 2.0d0 * m_t * m_i / (m_t + m_i)**2
      va2 = 2.0d0 * plasma%ion_temperature_eV * EV_TO_J / m_i

      call compute_el_lookup_coords(p, plasma, etm, tim)
      tim2 = max(tim, TLMT_EL)

      val_I_1_0 = get_I_1_0(etm, tim2)
      val_I_1_1 = get_I_1_1_up(etm, tim2)
      val_I_1_2 = get_I_1_2_up2(etm, tim)
      el_reaction_rate = plasma%ion_density * get_reaction_rate(etm, tim)

      if (ut > 1.0d-10 .and. ut2 > 1.0d-20) then
         sp0 = val_I_1_1 / ut - 0.5d0 * va2 / ut2 * val_I_1_0
         dtr_eng = -ame * (0.5d0 * (m_t + m_i) * sp0 * zuv - 0.5d0 * m_i * val_I_1_2)
      else
         dtr_eng = 0.0d0
      end if

      el_rate = plasma%ion_density * (-dtr_eng)
      if (el_reaction_rate > SCORE_TINY) then
         el_collision_score = el_rate / el_reaction_rate
      else
         el_collision_score = 0.0d0
      end if
   end subroutine compute_el_lookup_rate_and_score

end module scoring
