!===============================================================================
! Module: scoring
! 1D 空間分布を持つ CL/TL スコアリング
!===============================================================================
module scoring
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV, PI
   use data_types, only: particle_t, plasma_params, sim_params, score_data, &
      profile_data, rng_state
   use cross_sections, only: sigma_cx, sigma_el, ionization_rate_coeff
   use random_utils, only: sample_maxwell_velocity_ion, random_double
   use cdf_reader, only: sample_scattering_angle, get_I_1_0, get_I_1_1_up, &
      get_I_1_2_up2
   implicit none

   real(dp), parameter :: TLMT_EL = 0.9d0
   real(dp), parameter :: EPS_X = 1.0d-12

   private
   public :: zero_score, add_score
   public :: init_profile_data, free_profile_data
   public :: score_collision_estimator
   public :: score_track_length_estimator, score_table_lookup_track_length
   public :: flush_pending_track_scores, accumulate_collision_profile

contains

   subroutine zero_score(score)
      type(score_data), intent(out) :: score

      score%cl_ei = 0.0d0
      score%cl_cx = 0.0d0
      score%cl_el = 0.0d0
      score%tl_ei = 0.0d0
      score%tl_cx = 0.0d0
      score%tl_el = 0.0d0
      score%tl_lookup_ei = 0.0d0
      score%tl_lookup_cx = 0.0d0
      score%tl_lookup_el = 0.0d0
   end subroutine zero_score

   subroutine add_score(target, increment)
      type(score_data), intent(inout) :: target
      type(score_data), intent(in) :: increment

      target%cl_ei = target%cl_ei + increment%cl_ei
      target%cl_cx = target%cl_cx + increment%cl_cx
      target%cl_el = target%cl_el + increment%cl_el
      target%tl_ei = target%tl_ei + increment%tl_ei
      target%tl_cx = target%tl_cx + increment%tl_cx
      target%tl_el = target%tl_el + increment%tl_el
      target%tl_lookup_ei = target%tl_lookup_ei + increment%tl_lookup_ei
      target%tl_lookup_cx = target%tl_lookup_cx + increment%tl_lookup_cx
      target%tl_lookup_el = target%tl_lookup_el + increment%tl_lookup_el
   end subroutine add_score

   subroutine init_profile_data(profile, n_bins)
      type(profile_data), intent(out) :: profile
      integer, intent(in) :: n_bins

      allocate(profile%residence_eff(n_bins))
      allocate(profile%cl_ei(n_bins))
      allocate(profile%cl_cx(n_bins))
      allocate(profile%cl_el(n_bins))
      allocate(profile%tl_ei(n_bins))
      allocate(profile%tl_cx(n_bins))
      allocate(profile%tl_el(n_bins))
      allocate(profile%tl_lookup_ei(n_bins))
      allocate(profile%tl_lookup_cx(n_bins))
      allocate(profile%tl_lookup_el(n_bins))
      allocate(profile%collision_count_cx(n_bins))
      allocate(profile%collision_count_el(n_bins))
      allocate(profile%collision_count_total(n_bins))

      profile%residence_eff = 0.0d0
      profile%cl_ei = 0.0d0
      profile%cl_cx = 0.0d0
      profile%cl_el = 0.0d0
      profile%tl_ei = 0.0d0
      profile%tl_cx = 0.0d0
      profile%tl_el = 0.0d0
      profile%tl_lookup_ei = 0.0d0
      profile%tl_lookup_cx = 0.0d0
      profile%tl_lookup_el = 0.0d0
      profile%collision_count_cx = 0.0d0
      profile%collision_count_el = 0.0d0
      profile%collision_count_total = 0.0d0
   end subroutine init_profile_data

   subroutine free_profile_data(profile)
      type(profile_data), intent(inout) :: profile

      if (allocated(profile%residence_eff)) deallocate(profile%residence_eff)
      if (allocated(profile%cl_ei)) deallocate(profile%cl_ei)
      if (allocated(profile%cl_cx)) deallocate(profile%cl_cx)
      if (allocated(profile%cl_el)) deallocate(profile%cl_el)
      if (allocated(profile%tl_ei)) deallocate(profile%tl_ei)
      if (allocated(profile%tl_cx)) deallocate(profile%tl_cx)
      if (allocated(profile%tl_el)) deallocate(profile%tl_el)
      if (allocated(profile%tl_lookup_ei)) deallocate(profile%tl_lookup_ei)
      if (allocated(profile%tl_lookup_cx)) deallocate(profile%tl_lookup_cx)
      if (allocated(profile%tl_lookup_el)) deallocate(profile%tl_lookup_el)
      if (allocated(profile%collision_count_cx)) deallocate(profile%collision_count_cx)
      if (allocated(profile%collision_count_el)) deallocate(profile%collision_count_el)
      if (allocated(profile%collision_count_total)) deallocate(profile%collision_count_total)
   end subroutine free_profile_data

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

   subroutine score_track_length_estimator(p, plasma, pending_eff_time, &
      use_isotropic, enable_cx, enable_el, enable_ei, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: pending_eff_time
      logical, intent(in)             :: use_isotropic
      logical, intent(in)             :: enable_cx, enable_el, enable_ei
      type(score_data), intent(inout) :: score

      real(dp) :: E_kin
      real(dp) :: R_cx_d, R_el_d
      real(dp) :: s_c_ei, s_c_cx, s_c_el
      real(dp) :: delta_E_el_dummy
      real(dp) :: vx_d, vy_d, vz_d
      real(dp) :: v_rel_d, E_rel_d
      real(dp) :: sig_cx_d, sig_el_d
      real(dp) :: chi_d, phi_d, r_chi, r_phi
      real(dp) :: vx_g, vy_g, vz_g
      real(dp) :: ux, uy, uz, u_mag
      real(dp) :: ux_n, uy_n, uz_n
      real(dp) :: vx_new, vy_new, vz_new
      real(dp) :: E_old, E_new
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
         sig_el_d = sigma_el(E_rel_d)
         R_el_d = plasma%n_i * sig_el_d * v_rel_d

         E_old = E_kin
         vx_g = 0.5d0 * (p%vx + vx_d)
         vy_g = 0.5d0 * (p%vy + vy_d)
         vz_g = 0.5d0 * (p%vz + vz_d)

         ux = p%vx - vx_d
         uy = p%vy - vy_d
         uz = p%vz - vz_d
         u_mag = sqrt(ux*ux + uy*uy + uz*uz)

         delta_E_el_dummy = 0.0d0
         if (u_mag > 1.0d-30) then
            if (use_isotropic) then
               r_chi = random_double(rng_work)
               chi_d = acos(1.0d0 - 2.0d0 * r_chi)
            else
               r_chi = random_double(rng_work)
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
      end if

      score%tl_ei = score%tl_ei + s_c_ei * pending_eff_time
      score%tl_cx = score%tl_cx + s_c_cx * pending_eff_time
      score%tl_el = score%tl_el + s_c_el * pending_eff_time
   end subroutine score_track_length_estimator

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

   subroutine flush_pending_track_scores(p, plasma, sim, profile, step_score)
      type(particle_t), intent(inout) :: p
      type(plasma_params), intent(in) :: plasma
      type(sim_params), intent(in)    :: sim
      type(profile_data), intent(inout) :: profile
      type(score_data), intent(inout) :: step_score

      type(score_data) :: segment_score
      real(dp) :: loss_rate, total_eff_time

      if (p%tl_pending_time <= 0.0d0) return

      if (sim%enable_ei) then
         loss_rate = plasma%n_e * ionization_rate_coeff(plasma%T_e)
      else
         loss_rate = 0.0d0
      end if

      total_eff_time = compute_effective_track_time(p%tl_pending_weight_start, &
         loss_rate, p%tl_pending_time)

      call zero_score(segment_score)
      if (total_eff_time > 0.0d0) then
         call score_track_length_estimator(p, plasma, total_eff_time, &
            sim%use_isotropic, sim%enable_cx, sim%enable_el, sim%enable_ei, segment_score)
         if (sim%enable_tl_lookup) then
            call score_table_lookup_track_length(p, plasma, total_eff_time, &
               sim%enable_cx, sim%enable_el, sim%enable_ei, segment_score)
         end if

         call distribute_segment_profile(profile, sim, p%tl_pending_x_start, p%vx, &
            p%tl_pending_time, p%tl_pending_weight_start, loss_rate, total_eff_time, &
            segment_score)
      end if

      call add_score(step_score, segment_score)
      p%tl_pending_time = 0.0d0
      p%tl_pending_x_start = p%x
      p%tl_pending_weight_start = p%weight
   end subroutine flush_pending_track_scores

   subroutine accumulate_collision_profile(profile, sim, x, event_score, coll_type, event_weight)
      type(profile_data), intent(inout) :: profile
      type(sim_params), intent(in) :: sim
      real(dp), intent(in) :: x
      type(score_data), intent(in) :: event_score
      integer, intent(in) :: coll_type
      real(dp), intent(in) :: event_weight

      integer :: ix

      ix = get_x_bin(x, sim)

      profile%cl_ei(ix) = profile%cl_ei(ix) + event_score%cl_ei
      profile%cl_cx(ix) = profile%cl_cx(ix) + event_score%cl_cx
      profile%cl_el(ix) = profile%cl_el(ix) + event_score%cl_el
      profile%collision_count_total(ix) = profile%collision_count_total(ix) + event_weight

      if (coll_type == 1) then
         profile%collision_count_cx(ix) = profile%collision_count_cx(ix) + event_weight
      else if (coll_type == 2) then
         profile%collision_count_el(ix) = profile%collision_count_el(ix) + event_weight
      end if
   end subroutine accumulate_collision_profile

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

   function compute_effective_track_time(weight_start, loss_rate, dt) result(eff_time)
      real(dp), intent(in) :: weight_start, loss_rate, dt
      real(dp) :: eff_time

      if (dt <= 0.0d0 .or. weight_start <= 0.0d0) then
         eff_time = 0.0d0
      else if (loss_rate > 1.0d-30) then
         eff_time = (weight_start / loss_rate) * (1.0d0 - exp(-loss_rate * dt))
      else
         eff_time = weight_start * dt
      end if
   end function compute_effective_track_time

   function effective_time_interval(weight_start, loss_rate, t0, t1) result(eff_time)
      real(dp), intent(in) :: weight_start, loss_rate, t0, t1
      real(dp) :: eff_time

      if (t1 <= t0 .or. weight_start <= 0.0d0) then
         eff_time = 0.0d0
      else if (loss_rate > 1.0d-30) then
         eff_time = (weight_start / loss_rate) * &
            (exp(-loss_rate * t0) - exp(-loss_rate * t1))
      else
         eff_time = weight_start * (t1 - t0)
      end if
   end function effective_time_interval

   integer function get_x_bin(x, sim) result(ix)
      real(dp), intent(in) :: x
      type(sim_params), intent(in) :: sim

      real(dp) :: x_use

      if (x <= sim%x_min) then
         ix = 1
         return
      else if (x >= sim%x_max) then
         ix = sim%n_x_bins
         return
      end if

      x_use = x
      ix = int((x_use - sim%x_min) / sim%dx) + 1
      ix = max(1, min(ix, sim%n_x_bins))
   end function get_x_bin

   subroutine distribute_segment_profile(profile, sim, x0, vx, dt_total, weight_start, &
      loss_rate, total_eff_time, segment_score)
      type(profile_data), intent(inout) :: profile
      type(sim_params), intent(in) :: sim
      real(dp), intent(in) :: x0, vx, dt_total, weight_start, loss_rate, total_eff_time
      type(score_data), intent(in) :: segment_score

      integer :: ix
      real(dp) :: t_cursor, x_cursor, sub_dt, eff_sub, x_edge, dir_sign

      if (dt_total <= 0.0d0 .or. total_eff_time <= 0.0d0) return

      if (abs(vx) < 1.0d-20) then
         ix = get_x_bin(x0, sim)
         call add_segment_to_bin(profile, ix, segment_score, total_eff_time, total_eff_time)
         return
      end if

      t_cursor = 0.0d0
      x_cursor = x0
      dir_sign = sign(1.0d0, vx)

      do while (t_cursor < dt_total - 1.0d-14)
         if (vx > 0.0d0) then
            ix = get_x_bin(min(x_cursor, sim%x_max - EPS_X), sim)
            x_edge = sim%x_min + dble(ix) * sim%dx
         else
            ix = get_x_bin(max(x_cursor, sim%x_min + EPS_X), sim)
            x_edge = sim%x_min + dble(ix - 1) * sim%dx
         end if

         sub_dt = (x_edge - x_cursor) / vx
         if (sub_dt <= 1.0d-14) then
            x_cursor = x_cursor + dir_sign * EPS_X
            cycle
         end if
         sub_dt = min(sub_dt, dt_total - t_cursor)

         eff_sub = effective_time_interval(weight_start, loss_rate, t_cursor, t_cursor + sub_dt)
         call add_segment_to_bin(profile, ix, segment_score, eff_sub, total_eff_time)

         t_cursor = t_cursor + sub_dt
         x_cursor = x_cursor + vx * sub_dt
         if (t_cursor < dt_total - 1.0d-14) x_cursor = x_cursor + dir_sign * EPS_X
      end do
   end subroutine distribute_segment_profile

   subroutine add_segment_to_bin(profile, ix, segment_score, eff_sub, total_eff_time)
      type(profile_data), intent(inout) :: profile
      integer, intent(in) :: ix
      type(score_data), intent(in) :: segment_score
      real(dp), intent(in) :: eff_sub, total_eff_time

      real(dp) :: frac

      if (eff_sub <= 0.0d0 .or. total_eff_time <= 0.0d0) return
      frac = eff_sub / total_eff_time

      profile%residence_eff(ix) = profile%residence_eff(ix) + eff_sub
      profile%tl_ei(ix) = profile%tl_ei(ix) + segment_score%tl_ei * frac
      profile%tl_cx(ix) = profile%tl_cx(ix) + segment_score%tl_cx * frac
      profile%tl_el(ix) = profile%tl_el(ix) + segment_score%tl_el * frac
      profile%tl_lookup_ei(ix) = profile%tl_lookup_ei(ix) + segment_score%tl_lookup_ei * frac
      profile%tl_lookup_cx(ix) = profile%tl_lookup_cx(ix) + segment_score%tl_lookup_cx * frac
      profile%tl_lookup_el(ix) = profile%tl_lookup_el(ix) + segment_score%tl_lookup_el * frac
   end subroutine add_segment_to_bin

end module scoring
