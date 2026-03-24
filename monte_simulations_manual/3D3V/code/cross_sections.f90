!===============================================================================
! Module: cross_sections
! cross section calculation (CX, EL, EI), ionization rate coefficient, ν_max pre-calculation
!===============================================================================
module cross_sections
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV, CM2_TO_M2, PI, E_IONIZE_THRESHOLD
   use cdf_reader, only: get_sigma_elastic
   implicit none

   private
   public :: sigma_cx, sigma_el, sigma_ei, get_scatter_cross_sections
   public :: ionization_rate_coeff, compute_nu_max

   ! pre-calculated v_max/n_i
   real(dp), public :: sigma_v_max = 0.0d0     !max(σ_s * v_rel) [m^3/s]

contains

   !---------------------------------------------------------------------------
   ! Charge exchange cross section (Janev approximation)
   !---------------------------------------------------------------------------
   function sigma_cx(E_eV) result(sigma)
      real(dp), intent(in) :: E_eV
      real(dp) :: sigma
      real(dp) :: E_safe

      ! Prevent divergence at low energy
      E_safe = max(E_eV, 1.0d-3)

      ! Janev approximation (same as ntcros.f)
      ! σ_CX(E) [m^2], E is relative energy in center-of-mass frame [eV]
      sigma = 0.6937d-18 * (1.0d0 - 0.155d0 * log10(E_safe))**2 &
         / (1.0d0 + 0.1112d-14 * E_safe**3.3d0)

      ! Minimum value
      sigma = max(sigma, 1.0d-22)

   end function sigma_cx

   !---------------------------------------------------------------------------
   ! Elastic scattering cross section (from CDF data)
   !---------------------------------------------------------------------------
   function sigma_el(E_eV) result(sigma)
      real(dp), intent(in) :: E_eV
      real(dp) :: sigma

      sigma = get_sigma_elastic(E_eV)
   end function sigma_el

   !---------------------------------------------------------------------------
   ! Ionization cross section (simple Bethe model)
   !---------------------------------------------------------------------------
   function sigma_ei(E_eV) result(sigma)
      real(dp), intent(in) :: E_eV
      real(dp) :: sigma

      real(dp) :: u, A_coeff

      sigma = 0.0d0
      if (E_eV <= E_IONIZE_THRESHOLD) return

      u = E_eV / E_IONIZE_THRESHOLD
      A_coeff = 0.68d-20  ![m^2 * eV^2]

      sigma = A_coeff / (E_eV * E_IONIZE_THRESHOLD) * &
         (1.0d0 - 1.0d0/u) * log(u)  ![m^2]
      sigma = max(sigma, 0.0d0)
   end function sigma_ei

   !---------------------------------------------------------------------------
   ! Ionization rate coefficient <σv>_EI [m^3/s] (depends on background electron temp T_e)
   ! Approximation based on Lotz's semi-empirical formula
   !---------------------------------------------------------------------------
   function ionization_rate_coeff(T_e) result(rate)
      real(dp), intent(in) :: T_e  ![eV]
      real(dp) :: rate

      real(dp) :: x, sqrt_Te

      rate = 0.0d0
      if (T_e <= 0.0d0) return

      x = E_IONIZE_THRESHOLD / T_e

      if (x > 80.0d0) return  !exp(-x) ≈ 0

      sqrt_Te = sqrt(T_e * EV_TO_J)

      ! Lotz formula approximation: <σv> ≈ C * sqrt(T_e) * x * E1(x)
      ! Simplified approximation: E1(x) ≈ exp(-x) * ln(1 + 1/x) / (1 + x)
      rate = 2.34d-14 * exp(-x) * sqrt(x) / (1.0d0 + 0.25d0 * x)  ![m^3/s]
      rate = max(rate, 0.0d0)
   end function ionization_rate_coeff

   !---------------------------------------------------------------------------
   ! Get both CX and EL cross sections simultaneously
   !---------------------------------------------------------------------------
   subroutine get_scatter_cross_sections(E_rel, sig_cx, sig_el, sig_s)
      real(dp), intent(in)  :: E_rel    !relative energy [eV]
      real(dp), intent(out) :: sig_cx   !CX cross section [m^2]
      real(dp), intent(out) :: sig_el   !EL cross section [m^2]
      real(dp), intent(out) :: sig_s    !total scattering cross section [m^2]

      sig_cx = sigma_cx(E_rel)
      sig_el = sigma_el(E_rel)
      sig_s = sig_cx + sig_el
   end subroutine get_scatter_cross_sections

   !---------------------------------------------------------------------------
   ! Pre-calculate ν_max = n_i * max(σ_s * v_rel)
   ! Scan velocity range to find max(σ_s(E_rel)*v_rel)
   !---------------------------------------------------------------------------
   subroutine compute_nu_max(n_i)
      real(dp), intent(in) :: n_i  !background ion density [m^-3]

      integer :: i
      integer, parameter :: n_scan = 1000
      real(dp) :: v_rel, E_rel, sig_cx_val, sig_el_val, sig_s_val
      real(dp) :: sigma_v, sv_max
      real(dp) :: v_min, v_max, log_v_min, log_v_max, log_v

      ! Velocity range: 10 m/s to 10^7 m/s
      v_min = 1.0d1
      v_max = 1.0d7
      log_v_min = log(v_min)
      log_v_max = log(v_max)

      sv_max = 0.0d0

      do i = 0, n_scan
         log_v = log_v_min + (log_v_max - log_v_min) * dble(i) / dble(n_scan)
         v_rel = exp(log_v)

         ! Relative energy (center-of-mass): E_rel = (1/2)*μ*v_rel^2, μ=m_D/2
         E_rel = 0.25d0 * M_D_kg * v_rel * v_rel * J_TO_EV

         call get_scatter_cross_sections(E_rel, sig_cx_val, sig_el_val, sig_s_val)
         sigma_v = sig_s_val * v_rel

         if (sigma_v > sv_max) then
            sv_max = sigma_v
         end if
      end do

      sigma_v_max = sv_max
      write(*,'(A,ES12.4,A)') ' nu_max/n_i = sigma_v_max = ', sv_max, ' [m^3/s]'
      write(*,'(A,ES12.4,A)') ' nu_max = ', n_i * sv_max, ' [1/s]'

   end subroutine compute_nu_max

end module cross_sections
