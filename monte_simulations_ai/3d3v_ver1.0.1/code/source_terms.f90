!===============================================================================
! Module: source_terms
! Plasma-side Sn/Sp/We/Wi source algebra.
!===============================================================================
module source_terms
   use constants, only: dp, M_D_kg, EV_TO_J, E_IONIZE_THRESHOLD
   use data_types, only: particle_t, plasma_params, source_terms_t, &
      estimator_score_t, score_data, system_totals_t, N_REACTIONS, &
      REACT_EI, REACT_CX, REACT_EL, REACT_TOTAL
   implicit none

   private
   public :: zero_source_terms, zero_estimator_score, zero_score_data
   public :: add_source_terms, add_estimator_score, add_score_data
   public :: add_score_data_active
   public :: finalize_estimator_totals, finalize_score_totals
   public :: add_ei_source, add_cx_source
   public :: add_el_actual_source, add_el_average_source, add_el_rate_source
   public :: compute_system_totals

contains

   subroutine zero_source_terms(src)
      type(source_terms_t), intent(out) :: src

      src%sn_plus = 0.0d0
      src%sn_minus = 0.0d0
      src%sn_net = 0.0d0
      src%sp = 0.0d0
      src%we = 0.0d0
      src%wi = 0.0d0
   end subroutine zero_source_terms

   subroutine zero_estimator_score(est)
      type(estimator_score_t), intent(out) :: est
      integer :: ir

      do ir = 1, N_REACTIONS
         call zero_source_terms(est%reaction(ir))
      end do
   end subroutine zero_estimator_score

   subroutine zero_score_data(score)
      type(score_data), intent(out) :: score

      call zero_estimator_score(score%cl)
      call zero_estimator_score(score%cl_avg)
      call zero_estimator_score(score%tr)
      call zero_estimator_score(score%tr_pretab)
   end subroutine zero_score_data

   subroutine add_source_terms(dst, src)
      type(source_terms_t), intent(inout) :: dst
      type(source_terms_t), intent(in)    :: src

      dst%sn_plus = dst%sn_plus + src%sn_plus
      dst%sn_minus = dst%sn_minus + src%sn_minus
      dst%sn_net = dst%sn_net + src%sn_net
      dst%sp = dst%sp + src%sp
      dst%we = dst%we + src%we
      dst%wi = dst%wi + src%wi
   end subroutine add_source_terms

   subroutine add_estimator_score(dst, src)
      type(estimator_score_t), intent(inout) :: dst
      type(estimator_score_t), intent(in)    :: src
      integer :: ir

      do ir = 1, N_REACTIONS
         call add_source_terms(dst%reaction(ir), src%reaction(ir))
      end do
   end subroutine add_estimator_score

   subroutine add_score_data(dst, src)
      type(score_data), intent(inout) :: dst
      type(score_data), intent(in)    :: src

      call add_estimator_score(dst%cl, src%cl)
      call add_estimator_score(dst%cl_avg, src%cl_avg)
      call add_estimator_score(dst%tr, src%tr)
      call add_estimator_score(dst%tr_pretab, src%tr_pretab)
   end subroutine add_score_data

   subroutine add_score_data_active(dst, src, include_pretab)
      type(score_data), intent(inout) :: dst
      type(score_data), intent(in)    :: src
      logical, intent(in) :: include_pretab

      call add_estimator_score(dst%cl, src%cl)
      call add_estimator_score(dst%cl_avg, src%cl_avg)
      call add_estimator_score(dst%tr, src%tr)
      if (include_pretab) call add_estimator_score(dst%tr_pretab, src%tr_pretab)
   end subroutine add_score_data_active

   subroutine finalize_estimator_totals(est)
      type(estimator_score_t), intent(inout) :: est

      call zero_source_terms(est%reaction(REACT_TOTAL))
      call add_source_terms(est%reaction(REACT_TOTAL), est%reaction(REACT_EI))
      call add_source_terms(est%reaction(REACT_TOTAL), est%reaction(REACT_CX))
      call add_source_terms(est%reaction(REACT_TOTAL), est%reaction(REACT_EL))
   end subroutine finalize_estimator_totals

   subroutine finalize_score_totals(score, include_pretab)
      type(score_data), intent(inout) :: score
      logical, intent(in) :: include_pretab

      call finalize_estimator_totals(score%cl)
      call finalize_estimator_totals(score%cl_avg)
      call finalize_estimator_totals(score%tr)
      if (include_pretab) call finalize_estimator_totals(score%tr_pretab)
   end subroutine finalize_score_totals

   subroutine add_ei_source(src, dN, p)
      type(source_terms_t), intent(inout) :: src
      real(dp), intent(in) :: dN
      type(particle_t), intent(in) :: p

      real(dp) :: v2, kinetic, electron_loss

      if (dN <= 0.0d0) return

      v2 = p%vx*p%vx + p%vy*p%vy + p%vz*p%vz
      kinetic = 0.5d0 * M_D_kg * v2
      electron_loss = E_IONIZE_THRESHOLD * EV_TO_J

      src%sn_plus = src%sn_plus + dN
      src%sn_net = src%sn_net + dN
      src%sp(1) = src%sp(1) + dN * M_D_kg * p%vx
      src%sp(2) = src%sp(2) + dN * M_D_kg * p%vy
      src%sp(3) = src%sp(3) + dN * M_D_kg * p%vz
      src%we = src%we - dN * electron_loss
      src%wi = src%wi + dN * kinetic
   end subroutine add_ei_source

   subroutine add_cx_source(src, dN, p, plasma, vx_i, vy_i, vz_i)
      type(source_terms_t), intent(inout) :: src
      real(dp), intent(in) :: dN
      type(particle_t), intent(in) :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in), optional :: vx_i, vy_i, vz_i

      real(dp) :: v2, flow2, kinetic_neutral, kinetic_ion
      real(dp) :: vix, viy, viz
      logical :: use_sampled_ion

      if (dN <= 0.0d0) return

      use_sampled_ion = present(vx_i) .and. present(vy_i) .and. present(vz_i)
      if (use_sampled_ion) then
         vix = vx_i
         viy = vy_i
         viz = vz_i
      else
         vix = plasma%ion_flow_vx
         viy = plasma%ion_flow_vy
         viz = plasma%ion_flow_vz
      end if

      v2 = p%vx*p%vx + p%vy*p%vy + p%vz*p%vz
      kinetic_neutral = 0.5d0 * M_D_kg * v2
      if (use_sampled_ion) then
         kinetic_ion = 0.5d0 * M_D_kg * (vix*vix + viy*viy + viz*viz)
      else
         flow2 = plasma%ion_flow_vx**2 + plasma%ion_flow_vy**2 + plasma%ion_flow_vz**2
         kinetic_ion = 1.5d0 * plasma%ion_temperature_eV * EV_TO_J + &
            0.5d0 * M_D_kg * flow2
      end if

      src%sn_plus = src%sn_plus + dN
      src%sn_minus = src%sn_minus + dN
      src%sp(1) = src%sp(1) + dN * M_D_kg * (p%vx - vix)
      src%sp(2) = src%sp(2) + dN * M_D_kg * (p%vy - viy)
      src%sp(3) = src%sp(3) + dN * M_D_kg * (p%vz - viz)
      src%wi = src%wi + dN * (kinetic_neutral - kinetic_ion)
   end subroutine add_cx_source

   subroutine add_el_actual_source(src, p_before, p_after)
      type(source_terms_t), intent(inout) :: src
      type(particle_t), intent(in) :: p_before, p_after

      real(dp) :: dpx, dpy, dpz, e_before, e_after, dE

      dpx = M_D_kg * (p_after%vx - p_before%vx)
      dpy = M_D_kg * (p_after%vy - p_before%vy)
      dpz = M_D_kg * (p_after%vz - p_before%vz)
      e_before = 0.5d0 * M_D_kg * &
         (p_before%vx**2 + p_before%vy**2 + p_before%vz**2)
      e_after = 0.5d0 * M_D_kg * &
         (p_after%vx**2 + p_after%vy**2 + p_after%vz**2)
      dE = e_after - e_before

      src%sp(1) = src%sp(1) - p_before%weight * dpx
      src%sp(2) = src%sp(2) - p_before%weight * dpy
      src%sp(3) = src%sp(3) - p_before%weight * dpz
      src%wi = src%wi - p_before%weight * dE
   end subroutine add_el_actual_source

   subroutine add_el_average_source(src, dN, energy_score, momentum_score)
      type(source_terms_t), intent(inout) :: src
      real(dp), intent(in) :: dN, energy_score
      real(dp), intent(in) :: momentum_score(3)

      if (dN <= 0.0d0) return

      src%sp = src%sp + dN * momentum_score
      src%wi = src%wi + dN * energy_score
   end subroutine add_el_average_source

   subroutine add_el_rate_source(src, tau_eff, energy_rate, momentum_rate)
      type(source_terms_t), intent(inout) :: src
      real(dp), intent(in) :: tau_eff, energy_rate
      real(dp), intent(in) :: momentum_rate(3)

      if (tau_eff <= 0.0d0) return

      src%sp = src%sp + tau_eff * momentum_rate
      src%wi = src%wi + tau_eff * energy_rate
   end subroutine add_el_rate_source

   subroutine compute_system_totals(particles, totals)
      type(particle_t), intent(in) :: particles(:)
      type(system_totals_t), intent(out) :: totals

      integer :: i
      real(dp) :: v2

      totals%particles = 0.0d0
      totals%momentum = 0.0d0
      totals%energy = 0.0d0

      do i = 1, size(particles)
         if (.not. particles(i)%alive) cycle
         totals%particles = totals%particles + particles(i)%weight
         totals%momentum(1) = totals%momentum(1) + particles(i)%weight * M_D_kg * particles(i)%vx
         totals%momentum(2) = totals%momentum(2) + particles(i)%weight * M_D_kg * particles(i)%vy
         totals%momentum(3) = totals%momentum(3) + particles(i)%weight * M_D_kg * particles(i)%vz
         v2 = particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2
         totals%energy = totals%energy + particles(i)%weight * 0.5d0 * M_D_kg * v2
      end do
   end subroutine compute_system_totals

end module source_terms
