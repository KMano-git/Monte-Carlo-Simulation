!===============================================================================
! Module: scoring
! Collision Estimator (CL) と Track-Length Estimator (TL) の実装
!===============================================================================
module scoring
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV, PI
   use data_types, only: particle_t, plasma_params, score_data, rng_state
   use cross_sections, only: sigma_cx, sigma_el, ionization_rate_coeff
   use random_utils, only: sample_maxwell_velocity_ion, random_double
   use cdf_reader, only: sample_scattering_angle, get_I_1_0, get_I_1_1_up, &
      get_I_1_2_up2
   use tl_el_table_reader, only: get_tl_el_rate_from_table
   implicit none

   real(dp), parameter :: TLMT_EL = 0.9d0

   private
   public :: score_collision_estimator, score_track_length_estimator
   public :: score_inner_multi_track_length, score_pretabulated_track_length
   public :: score_table_lookup_track_length
   public :: tl_el_sample_once, tl_el_sample_avg, tl_el_sample_stats

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
      real(dp), intent(in)  :: vx_i, vy_i, vz_i
      real(dp), intent(in)  :: v_rel
      real(dp), intent(in)  :: E_rel
      integer, intent(in)   :: coll_type
      real(dp), intent(in)  :: delta_E_el
      logical, intent(in)   :: enable_ei
      type(score_data), intent(inout) :: score

      real(dp) :: E_kin
      real(dp) :: R_cx, R_el, R_s, R_a
      real(dp) :: s_a, s_cx, s_el
      real(dp) :: sig_cx_val, sig_el_val

      E_kin = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)

      sig_cx_val = sigma_cx(E_rel)
      sig_el_val = sigma_el(E_rel)

      R_cx = plasma%n_i * sig_cx_val * v_rel
      R_el = plasma%n_i * sig_el_val * v_rel
      R_s  = R_cx + R_el

      if (enable_ei) then
         R_a = plasma%n_e * ionization_rate_coeff(plasma%T_e)
      else
         R_a = 0.0d0
      end if

      s_a = E_kin
      s_cx = E_kin - 0.5d0 * M_D_kg * (vx_i**2 + vy_i**2 + vz_i**2)
      s_el = -delta_E_el

      if (R_s > 1.0d-30) then
         score%cl_ei = score%cl_ei + p%weight * (R_a / R_s) * s_a

         if (coll_type == 1) then
            score%cl_cx = score%cl_cx + p%weight * s_cx
         else if (coll_type == 2) then
            score%cl_el = score%cl_el + p%weight * s_el
         end if
      end if

   end subroutine score_collision_estimator

   !---------------------------------------------------------------------------
   ! Track-Length Estimator (TL)
   ! 実衝突までに蓄積した pending_eff_time をまとめてスコアする
   !---------------------------------------------------------------------------
   subroutine score_track_length_estimator(p, plasma, pending_eff_time, &
      use_isotropic, enable_cx, enable_el, enable_ei, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: pending_eff_time
      logical, intent(in)             :: use_isotropic
      logical, intent(in)             :: enable_cx, enable_el, enable_ei
      type(score_data), intent(inout) :: score

      real(dp) :: E_kin
      real(dp) :: R_cx_d
      real(dp) :: s_c_ei, s_c_cx, s_c_el
      real(dp) :: vx_d, vy_d, vz_d
      real(dp) :: v_rel_d, E_rel_d
      real(dp) :: sig_cx_d
      type(rng_state) :: rng_work

      if (pending_eff_time <= 0.0d0) return

      E_kin = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)

      if (enable_ei) then
         s_c_ei = plasma%n_e * ionization_rate_coeff(plasma%T_e) * E_kin
      else
         s_c_ei = 0.0d0
      end if

      s_c_cx = 0.0d0
      s_c_el = 0.0d0
      if (.not. enable_cx .and. .not. enable_el) then
         score%tl_ei = score%tl_ei + s_c_ei * pending_eff_time
         return
      end if

      rng_work = p%rng
      call sample_maxwell_velocity_ion(rng_work, plasma%T_i, plasma, vx_d, vy_d, vz_d)

      v_rel_d = sqrt((p%vx - vx_d)**2 + (p%vy - vy_d)**2 + (p%vz - vz_d)**2)
      E_rel_d = 0.25d0 * M_D_kg * v_rel_d * v_rel_d * J_TO_EV

      if (enable_cx) then
         sig_cx_d = sigma_cx(E_rel_d)
         R_cx_d = plasma%n_i * sig_cx_d * v_rel_d
         s_c_cx = R_cx_d * (E_kin - 0.5d0 * M_D_kg * (vx_d**2 + vy_d**2 + vz_d**2))
      end if

      if (enable_el) then
         s_c_el = tl_el_rate_from_ion_sample(p, plasma, use_isotropic, rng_work, &
            vx_d, vy_d, vz_d)
      end if

      score%tl_ei = score%tl_ei + s_c_ei * pending_eff_time
      score%tl_cx = score%tl_cx + s_c_cx * pending_eff_time
      score%tl_el = score%tl_el + s_c_el * pending_eff_time

   end subroutine score_track_length_estimator

   !---------------------------------------------------------------------------
   ! EL 1サンプル分の rate を、既にサンプリング済みの背景イオン速度から評価する
   !---------------------------------------------------------------------------
   function tl_el_rate_from_ion_sample(p, plasma, use_isotropic, rng_work, &
      vx_d, vy_d, vz_d) result(s_c_el)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      logical, intent(in)             :: use_isotropic
      type(rng_state), intent(inout)  :: rng_work
      real(dp), intent(in)            :: vx_d, vy_d, vz_d
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
      R_el_d = plasma%n_i * sig_el_d * v_rel_d

      E_old = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)
      vx_g = 0.5d0 * (p%vx + vx_d)
      vy_g = 0.5d0 * (p%vy + vy_d)
      vz_g = 0.5d0 * (p%vz + vz_d)

      ux = p%vx - vx_d
      uy = p%vy - vy_d
      uz = p%vz - vz_d
      u_mag = sqrt(ux*ux + uy*uy + uz*uz)

      delta_E_el_dummy = 0.0d0
      if (u_mag > 1.0d-30) then
         r_chi = random_double(rng_work)
         if (use_isotropic) then
            chi_d = acos(1.0d0 - 2.0d0 * r_chi)
         else
            chi_d = sample_scattering_angle(0.5d0 * E_rel_d, r_chi)
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

      s_c_el = R_el_d * (-delta_E_el_dummy)
   end function tl_el_rate_from_ion_sample

   !---------------------------------------------------------------------------
   ! EL naive TL 1サンプル版
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
   ! EL inner multi-sampling 平均版
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

   function tl_el_sample_once_with_rng(p, plasma, use_isotropic, rng_work) result(s_c_el)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      logical, intent(in)             :: use_isotropic
      type(rng_state), intent(inout)  :: rng_work
      real(dp) :: s_c_el

      real(dp) :: vx_d, vy_d, vz_d

      call sample_maxwell_velocity_ion(rng_work, plasma%T_i, plasma, vx_d, vy_d, vz_d)
      s_c_el = tl_el_rate_from_ion_sample(p, plasma, use_isotropic, rng_work, vx_d, vy_d, vz_d)
   end function tl_el_sample_once_with_rng

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
   ! EL inner multi-sampling TL
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
      score%tl_inner_el = score%tl_inner_el + s_c_el * pending_eff_time
   end subroutine score_inner_multi_track_length

   !---------------------------------------------------------------------------
   ! EL pre-tabulated TL
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

      call compute_tl_el_lookup_coords(p, plasma, etm, tim)
      s_c_el = get_tl_el_rate_from_table(etm, tim)
      score%tl_pretab_el = score%tl_pretab_el + s_c_el * pending_eff_time
   end subroutine score_pretabulated_track_length

   pure subroutine compute_tl_el_lookup_coords(p, plasma, etm, tim)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(out)           :: etm, tim

      real(dp) :: utx, uty, utz, ut2
      real(dp) :: m_ref, m_i

      utx = p%vx - plasma%u_x
      uty = p%vy - plasma%u_y
      utz = p%vz - plasma%u_z
      ut2 = utx**2 + uty**2 + utz**2

      m_ref = 0.5d0 * M_D_kg
      m_i = M_D_kg

      etm = 0.5d0 * m_ref * ut2 * J_TO_EV
      tim = (m_ref / m_i) * plasma%T_i
   end subroutine compute_tl_el_lookup_coords

   !---------------------------------------------------------------------------
   ! Table-Lookup Track-Length Estimator (TL Lookup)
   ! 実衝突までに蓄積した pending_eff_time をまとめてスコアする
   !---------------------------------------------------------------------------
   subroutine score_table_lookup_track_length(p, plasma, pending_eff_time, &
      enable_cx, enable_el, enable_ei, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: pending_eff_time
      logical, intent(in)             :: enable_cx, enable_el, enable_ei
      type(score_data), intent(inout) :: score

      real(dp) :: E_kin
      real(dp) :: s_c_ei, s_c_cx, s_c_el
      real(dp) :: utx, uty, utz, ut2, ut, zuv, va2
      real(dp) :: etm, tim, tim2
      real(dp) :: val_I_1_0, val_I_1_1, val_I_1_2
      real(dp) :: sp0, dtr_eng
      real(dp) :: vx_d, vy_d, vz_d, v_rel_d, E_rel_cx, sig_cx_d
      real(dp) :: m_ref, m_t, m_i, ame
      type(rng_state) :: rng_work

      if (pending_eff_time <= 0.0d0) return

      E_kin = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)

      if (enable_ei) then
         s_c_ei = plasma%n_e * ionization_rate_coeff(plasma%T_e) * E_kin
      else
         s_c_ei = 0.0d0
      end if

      s_c_cx = 0.0d0
      if (enable_cx) then
         rng_work = p%rng
         call sample_maxwell_velocity_ion(rng_work, plasma%T_i, plasma, vx_d, vy_d, vz_d)
         v_rel_d = sqrt((p%vx - vx_d)**2 + (p%vy - vy_d)**2 + (p%vz - vz_d)**2)
         E_rel_cx = 0.25d0 * M_D_kg * v_rel_d * v_rel_d * J_TO_EV
         sig_cx_d = sigma_cx(E_rel_cx)
         s_c_cx = plasma%n_i * sig_cx_d * v_rel_d * &
            (E_kin - 0.5d0 * M_D_kg * (vx_d**2 + vy_d**2 + vz_d**2))
      end if

      s_c_el = 0.0d0
      if (enable_el) then
         utx = p%vx - plasma%u_x
         uty = p%vy - plasma%u_y
         utz = p%vz - plasma%u_z
         ut2 = utx**2 + uty**2 + utz**2
         ut = sqrt(ut2)
         zuv = utx*p%vx + uty*p%vy + utz*p%vz

         m_t = M_D_kg
         m_i = M_D_kg
         m_ref = 0.5d0 * M_D_kg
         ame = 2.0d0 * m_t * m_i / (m_t + m_i)**2
         va2 = 2.0d0 * plasma%T_i * EV_TO_J / m_i

         etm = 0.5d0 * m_ref * ut2 * J_TO_EV
         tim = (m_ref / m_i) * plasma%T_i
         tim2 = max(tim, TLMT_EL)

         val_I_1_0 = get_I_1_0(etm, tim2)
         val_I_1_1 = get_I_1_1_up(etm, tim2)
         val_I_1_2 = get_I_1_2_up2(etm, tim)

         if (ut > 1.0d-10 .and. ut2 > 1.0d-20) then
            sp0 = val_I_1_1 / ut - 0.5d0 * va2 / ut2 * val_I_1_0
            dtr_eng = -ame * (0.5d0 * (m_t + m_i) * sp0 * zuv - 0.5d0 * m_i * val_I_1_2)
         else
            dtr_eng = 0.0d0
         end if

         s_c_el = plasma%n_i * (-dtr_eng)
      end if

      score%tl_lookup_ei = score%tl_lookup_ei + s_c_ei * pending_eff_time
      score%tl_lookup_cx = score%tl_lookup_cx + s_c_cx * pending_eff_time
      score%tl_lookup_el = score%tl_lookup_el + s_c_el * pending_eff_time

   end subroutine score_table_lookup_track_length

end module scoring
