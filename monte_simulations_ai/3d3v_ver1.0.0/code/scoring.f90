!===============================================================================
! Module: scoring
! CL, CL(avg), and TR source estimators.
!===============================================================================
module scoring
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV
   use data_types, only: particle_t, plasma_params, score_data, rng_state, &
      REACT_EI, REACT_CX, REACT_EL
   use cross_sections, only: sigma_cx, sigma_el, ionization_rate_coeff
   use random_utils, only: sample_maxwell_velocity_ion
   use cdf_reader, only: get_reaction_rate, get_I_1_0, get_I_1_1_up, get_I_1_2_up2
   use tl_el_table_reader, only: get_tl_el_rate_from_table
   use source_terms, only: add_ei_source, add_cx_source, add_el_actual_source, &
      add_el_average_source, add_el_rate_source
   implicit none

   real(dp), parameter :: TLMT_EL = 0.9d0
   real(dp), parameter :: SCORE_TINY = 1.0d-30

   private
   public :: score_collision_estimators
   public :: score_el_actual_collision
   public :: score_track_length_estimator
   public :: score_pretabulated_track_length
   public :: compute_el_lookup_sources

contains

   !---------------------------------------------------------------------------
   ! Accepted scatter-collision estimator.
   !
   ! EI/CX are added to both CL and CL(avg), unchanged. EL is added only to
   ! CL(avg) here. The actual EL contribution to CL is scored after the EL
   ! collision updates the particle velocity.
   !---------------------------------------------------------------------------
   subroutine score_collision_estimators(p, plasma, vx_i, vy_i, vz_i, v_rel, E_rel, &
      enable_cx, enable_el, enable_ionization, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: vx_i, vy_i, vz_i
      real(dp), intent(in)            :: v_rel, E_rel
      logical, intent(in)             :: enable_cx, enable_el, enable_ionization
      type(score_data), intent(inout) :: score

      real(dp) :: sig_cx_val, sig_el_val
      real(dp) :: R_cx, R_el, R_s, R_a
      real(dp) :: dN_ei, dN_cx, dN_el
      real(dp) :: el_energy_rate, el_reaction_rate, el_energy_score
      real(dp) :: el_momentum_rate(3), el_momentum_score(3)

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
      if (R_s <= SCORE_TINY) return

      if (enable_ionization) then
         R_a = plasma%electron_density * &
            ionization_rate_coeff(plasma%electron_temperature_eV)
      else
         R_a = 0.0d0
      end if

      dN_ei = p%weight * (R_a / R_s)
      dN_cx = p%weight * (R_cx / R_s)

      call add_ei_source(score%cl%reaction(REACT_EI), dN_ei, p)
      call add_ei_source(score%cl_avg%reaction(REACT_EI), dN_ei, p)

      call add_cx_source(score%cl%reaction(REACT_CX), dN_cx, p, plasma, &
         vx_i, vy_i, vz_i)
      call add_cx_source(score%cl_avg%reaction(REACT_CX), dN_cx, p, plasma, &
         vx_i, vy_i, vz_i)

      if (enable_el .and. R_el > SCORE_TINY) then
         dN_el = p%weight * (R_el / R_s)
         call compute_el_lookup_sources(p, plasma, el_energy_rate, &
            el_momentum_rate, el_reaction_rate, el_energy_score, &
            el_momentum_score)
         call add_el_average_source(score%cl_avg%reaction(REACT_EL), dN_el, &
            el_energy_score, el_momentum_score)
      end if
   end subroutine score_collision_estimators

   !---------------------------------------------------------------------------
   ! Actual EL contribution to CL.
   !---------------------------------------------------------------------------
   subroutine score_el_actual_collision(p_before, p_after, score)
      type(particle_t), intent(in) :: p_before, p_after
      type(score_data), intent(inout) :: score

      call add_el_actual_source(score%cl%reaction(REACT_EL), p_before, p_after)
   end subroutine score_el_actual_collision

   !---------------------------------------------------------------------------
   ! Track-length estimator for EI/CX/EL.
   !---------------------------------------------------------------------------
   subroutine score_track_length_estimator(p, plasma, pending_eff_time, &
      enable_cx, enable_el, enable_ionization, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: pending_eff_time
      logical, intent(in)             :: enable_cx, enable_el, enable_ionization
      type(score_data), intent(inout) :: score

      real(dp) :: R_a, dN_ei
      real(dp) :: vx_d, vy_d, vz_d, v_rel_d, E_rel_cx
      real(dp) :: R_cx_sample, dN_cx
      real(dp) :: el_energy_rate, el_reaction_rate, el_energy_score
      real(dp) :: el_momentum_rate(3), el_momentum_score(3)
      type(rng_state) :: rng_work

      if (pending_eff_time <= 0.0d0) return

      if (enable_ionization) then
         R_a = plasma%electron_density * &
            ionization_rate_coeff(plasma%electron_temperature_eV)
         dN_ei = R_a * pending_eff_time
         call add_ei_source(score%tr%reaction(REACT_EI), dN_ei, p)
      end if

      if (enable_cx) then
         rng_work = p%rng
         call sample_maxwell_velocity_ion(rng_work, plasma%ion_temperature_eV, &
            plasma, vx_d, vy_d, vz_d)
         v_rel_d = sqrt((p%vx - vx_d)**2 + (p%vy - vy_d)**2 + (p%vz - vz_d)**2)
         E_rel_cx = 0.25d0 * M_D_kg * v_rel_d * v_rel_d * J_TO_EV
         R_cx_sample = plasma%ion_density * sigma_cx(E_rel_cx) * v_rel_d
         dN_cx = R_cx_sample * pending_eff_time
         call add_cx_source(score%tr%reaction(REACT_CX), dN_cx, p, plasma, &
            vx_d, vy_d, vz_d)
      end if

      if (enable_el) then
         call compute_el_lookup_sources(p, plasma, el_energy_rate, &
            el_momentum_rate, el_reaction_rate, el_energy_score, &
            el_momentum_score)
         call add_el_rate_source(score%tr%reaction(REACT_EL), pending_eff_time, &
            el_energy_rate, el_momentum_rate)
      end if
   end subroutine score_track_length_estimator

   !---------------------------------------------------------------------------
   ! Optional pre-tabulated TR comparison. The current table contains only the
   ! EL energy rate, so momentum is intentionally left zero.
   !---------------------------------------------------------------------------
   subroutine score_pretabulated_track_length(p, plasma, pending_eff_time, &
      enable_el, score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in)            :: pending_eff_time
      logical, intent(in)             :: enable_el
      type(score_data), intent(inout) :: score

      real(dp) :: etm, tim, energy_rate
      real(dp) :: momentum_rate(3)

      if (pending_eff_time <= 0.0d0 .or. .not. enable_el) return

      call compute_el_lookup_coords(p, plasma, etm, tim)
      energy_rate = get_tl_el_rate_from_table(etm, tim)
      momentum_rate = 0.0d0
      call add_el_rate_source(score%tr_pretab%reaction(REACT_EL), pending_eff_time, &
         energy_rate, momentum_rate)
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

      etm = 0.5d0 * m_ref * ut2 * J_TO_EV
      tim = (m_ref / m_i) * plasma%ion_temperature_eV
   end subroutine compute_el_lookup_coords

   !---------------------------------------------------------------------------
   ! EL lookup kernel.
   !
   ! energy_rate and momentum_rate are plasma-side source rates.
   ! energy_score and momentum_score are plasma-side source per EL collision.
   !---------------------------------------------------------------------------
   subroutine compute_el_lookup_sources(p, plasma, energy_rate, momentum_rate, &
      reaction_rate, energy_score, momentum_score)
      type(particle_t), intent(in)    :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(out)           :: energy_rate
      real(dp), intent(out)           :: momentum_rate(3)
      real(dp), intent(out)           :: reaction_rate
      real(dp), intent(out)           :: energy_score
      real(dp), intent(out)           :: momentum_score(3)

      real(dp) :: utx, uty, utz, ut2, ut, zuv, va2
      real(dp) :: etm, tim, tim2
      real(dp) :: val_I_1_0, val_I_1_1, val_I_1_2
      real(dp) :: sp0, dtr_eng
      real(dp) :: m_ref, m_t, m_i, ame
      real(dp) :: dtr_mv(3)

      utx = p%vx - plasma%ion_flow_vx
      uty = p%vy - plasma%ion_flow_vy
      utz = p%vz - plasma%ion_flow_vz
      ut2 = utx**2 + uty**2 + utz**2
      ut = sqrt(ut2)
      zuv = utx * p%vx + uty * p%vy + utz * p%vz

      m_t = M_D_kg
      m_i = M_D_kg
      m_ref = m_t * m_i / (m_t + m_i)
      ame = 2.0d0 * m_t * m_i / (m_t + m_i)**2
      va2 = 2.0d0 * plasma%ion_temperature_eV * EV_TO_J / m_i

      call compute_el_lookup_coords(p, plasma, etm, tim)
      tim2 = max(tim, TLMT_EL)

      val_I_1_0 = get_I_1_0(etm, tim2)
      val_I_1_1 = get_I_1_1_up(etm, tim2)
      val_I_1_2 = get_I_1_2_up2(etm, tim)
      reaction_rate = plasma%ion_density * get_reaction_rate(etm, tim)

      if (ut > 1.0d-10 .and. ut2 > 1.0d-20) then
         sp0 = val_I_1_1 / ut - 0.5d0 * va2 / ut2 * val_I_1_0
         dtr_mv(1) = -utx * m_ref * sp0
         dtr_mv(2) = -uty * m_ref * sp0
         dtr_mv(3) = -utz * m_ref * sp0
         dtr_eng = -ame * (0.5d0 * (m_t + m_i) * sp0 * zuv - &
            0.5d0 * m_i * val_I_1_2)
         momentum_rate = -plasma%ion_density * dtr_mv
         energy_rate = -plasma%ion_density * dtr_eng
      else
         momentum_rate = 0.0d0
         energy_rate = 0.0d0
      end if

      if (reaction_rate > SCORE_TINY) then
         energy_score = energy_rate / reaction_rate
         momentum_score = momentum_rate / reaction_rate
      else
         energy_score = 0.0d0
         momentum_score = 0.0d0
      end if
   end subroutine compute_el_lookup_sources

end module scoring
