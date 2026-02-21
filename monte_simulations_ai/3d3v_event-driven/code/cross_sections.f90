!===============================================================================
! Module: cross_sections
! 断面積の計算（CX, EL, EI）、電離レート係数、ν_max事前計算
!===============================================================================
module cross_sections
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV, CM2_TO_M2, PI, &
      E_IONIZE_THRESHOLD
   use cdf_reader, only: get_sigma_elastic
   implicit none

   private
   public :: sigma_cx, sigma_el, sigma_ei, get_scatter_cross_sections
   public :: ionization_rate_coeff, compute_nu_max

   !ν_max/n_i の事前計算値
   real(dp), public :: sigma_v_max = 0.0d0     !max(σ_s * v_rel) [m^3/s]

contains

   !---------------------------------------------------------------------------
   ! 荷電交換断面積（Janev近似式）
   ! σ_CX(E) [m^2]、E は重心系相対エネルギー [eV]
   !---------------------------------------------------------------------------
   function sigma_cx(E_eV) result(sigma)
      real(dp), intent(in) :: E_eV
      real(dp) :: sigma
      real(dp) :: E_safe, log10E, a1, a2

      E_safe = max(E_eV, 1.0d-3)
      log10E = log10(E_safe)

      !Janev式:  σ = (a1 - a2*log10(E))^2 [cm^2]
      a1 = 3.245d0
      a2 = 0.4060d0

      sigma = (a1 - a2 * log10E)**2 * 1.0d-16 * CM2_TO_M2  ![m^2]
      sigma = max(sigma, 0.0d0)
   end function sigma_cx

   !---------------------------------------------------------------------------
   ! 弾性散乱断面積（CDFデータから）
   !---------------------------------------------------------------------------
   function sigma_el(E_eV) result(sigma)
      real(dp), intent(in) :: E_eV
      real(dp) :: sigma

      sigma = get_sigma_elastic(E_eV)
   end function sigma_el

   !---------------------------------------------------------------------------
   ! 電離断面積（簡易ベーテ式モデル）
   ! σ_EI(E) [m^2]、E は電子エネルギー [eV]
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
   ! 電離レート係数 <σv>_EI [m^3/s]（背景電子温度T_eに依存）
   ! Lotzの半経験式に基づく近似
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

      !Lotz式近似: <σv> ≈ C * sqrt(T_e) * x * E1(x)
      !簡易近似: E1(x) ≈ exp(-x) * ln(1 + 1/x) / (1 + x)
      rate = 2.34d-14 * exp(-x) * sqrt(x) / (1.0d0 + 0.25d0 * x)  ![m^3/s]
      rate = max(rate, 0.0d0)
   end function ionization_rate_coeff

   !---------------------------------------------------------------------------
   ! CX+EL断面積を同時に取得
   !---------------------------------------------------------------------------
   subroutine get_scatter_cross_sections(E_rel, sig_cx, sig_el, sig_s)
      real(dp), intent(in)  :: E_rel    !相対エネルギー [eV]
      real(dp), intent(out) :: sig_cx   !CX断面積 [m^2]
      real(dp), intent(out) :: sig_el   !EL断面積 [m^2]
      real(dp), intent(out) :: sig_s    !散乱断面積合計 [m^2]

      sig_cx = sigma_cx(E_rel)
      sig_el = sigma_el(E_rel)
      sig_s = sig_cx + sig_el
   end subroutine get_scatter_cross_sections

   !---------------------------------------------------------------------------
   ! ν_max = n_i * max(σ_s * v_rel) を事前計算
   ! 速度範囲をスキャンして σ_s(E_rel)*v_rel の最大値を探す
   !---------------------------------------------------------------------------
   subroutine compute_nu_max(n_i)
      real(dp), intent(in) :: n_i  !背景イオン密度 [m^-3]

      integer :: i
      integer, parameter :: n_scan = 1000
      real(dp) :: v_rel, E_rel, sig_cx_val, sig_el_val, sig_s_val
      real(dp) :: sigma_v, sv_max
      real(dp) :: v_min, v_max, log_v_min, log_v_max, log_v

      !速度範囲: 10 m/s ～ 10^7 m/s
      v_min = 1.0d1
      v_max = 1.0d7
      log_v_min = log(v_min)
      log_v_max = log(v_max)

      sv_max = 0.0d0

      do i = 0, n_scan
         log_v = log_v_min + (log_v_max - log_v_min) * dble(i) / dble(n_scan)
         v_rel = exp(log_v)

         !相対エネルギー（重心系）: E_rel = (1/2)*μ*v_rel^2, μ=m_D/2
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
