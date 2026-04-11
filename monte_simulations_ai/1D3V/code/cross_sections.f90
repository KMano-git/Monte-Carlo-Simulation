!===============================================================================
! Module: cross_sections
! 断面積計算と ν_max 評価
!===============================================================================
module cross_sections
   use constants, only: dp, M_D_kg, J_TO_EV, E_IONIZE_THRESHOLD
   use cdf_reader, only: get_sigma_elastic
   implicit none

   private
   public :: sigma_cx, sigma_el, sigma_ei, get_scatter_cross_sections
   public :: ionization_rate_coeff, compute_nu_max

   real(dp), public :: sigma_v_max = 0.0d0

contains

   function sigma_cx(E_eV) result(sigma)
      real(dp), intent(in) :: E_eV
      real(dp) :: sigma
      real(dp) :: E_safe

      E_safe = max(E_eV, 0.01d0)
      sigma = 0.6937d-18 * (1.0d0 - 0.155d0 * log10(E_safe))**2 &
         / (1.0d0 + 0.1112d-14 * E_safe**3.3d0)
      sigma = max(sigma, 1.0d-22)
   end function sigma_cx

   function sigma_el(E_eV) result(sigma)
      real(dp), intent(in) :: E_eV
      real(dp) :: sigma

      ! DEGAS/EIRENE 側の elastic grid は specific_energy [eV/amu] だが、
      ! D + D+ では reduced mass = 1 amu なので数値的に E_rel [eV] と一致する。
      sigma = get_sigma_elastic(E_eV)
   end function sigma_el

   function sigma_ei(E_eV) result(sigma)
      real(dp), intent(in) :: E_eV
      real(dp) :: sigma

      real(dp) :: u, A_coeff

      sigma = 0.0d0
      if (E_eV <= E_IONIZE_THRESHOLD) return

      u = E_eV / E_IONIZE_THRESHOLD
      A_coeff = 0.68d-20
      sigma = A_coeff / (E_eV * E_IONIZE_THRESHOLD) * (1.0d0 - 1.0d0 / u) * log(u)
      sigma = max(sigma, 0.0d0)
   end function sigma_ei

   function ionization_rate_coeff(T_e) result(rate)
      real(dp), intent(in) :: T_e
      real(dp) :: rate

      real(dp) :: x

      rate = 0.0d0
      if (T_e <= 0.0d0) return

      x = E_IONIZE_THRESHOLD / T_e
      if (x > 80.0d0) return

      rate = 2.34d-14 * exp(-x) * sqrt(x) / (1.0d0 + 0.25d0 * x)
      rate = max(rate, 0.0d0)
   end function ionization_rate_coeff

   subroutine get_scatter_cross_sections(E_rel, sig_cx, sig_el, sig_s)
      real(dp), intent(in)  :: E_rel
      real(dp), intent(out) :: sig_cx
      real(dp), intent(out) :: sig_el
      real(dp), intent(out) :: sig_s

      sig_cx = sigma_cx(E_rel)
      sig_el = sigma_el(E_rel)
      sig_s = sig_cx + sig_el
   end subroutine get_scatter_cross_sections

   subroutine compute_nu_max(n_i)
      real(dp), intent(in) :: n_i

      integer :: i
      integer, parameter :: n_scan = 1000
      real(dp) :: v_rel, E_rel, sig_cx_val, sig_el_val, sig_s_val
      real(dp) :: sigma_v, sv_max
      real(dp) :: v_min, v_max, log_v_min, log_v_max, log_v

      v_min = 1.0d1
      v_max = 1.0d7
      log_v_min = log(v_min)
      log_v_max = log(v_max)

      sv_max = 0.0d0
      do i = 0, n_scan
         log_v = log_v_min + (log_v_max - log_v_min) * dble(i) / dble(n_scan)
         v_rel = exp(log_v)
         E_rel = 0.25d0 * M_D_kg * v_rel * v_rel * J_TO_EV

         call get_scatter_cross_sections(E_rel, sig_cx_val, sig_el_val, sig_s_val)
         sigma_v = sig_s_val * v_rel
         if (sigma_v > sv_max) sv_max = sigma_v
      end do

      sigma_v_max = sv_max
      write(*,'(A,ES12.4,A)') ' nu_max/n_i = sigma_v_max = ', sv_max, ' [m^3/s]'
      write(*,'(A,ES12.4,A)') ' nu_max = ', n_i * sv_max, ' [1/s]'
   end subroutine compute_nu_max

end module cross_sections
