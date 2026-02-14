!===============================================================================
! Module: cross_sections
! 衝突断面積の計算（荷電交換、弾性衝突、電離衝突）
!===============================================================================
module cross_sections
   use constants, only: dp, EV, E_IONIZE_THRESHOLD
   use cdf_reader, only: get_sigma_elastic
   implicit none

contains

   !-----------------------------------------------------------------------------
   ! 荷電交換断面積（Janev/Rabinovitch近似）
   ! D⁰ + D⁺ → D⁺ + D⁰
   !-----------------------------------------------------------------------------
   function sigma_cx(energy) result(sigma)
      real(dp), intent(in) :: energy  ! 衝突エネルギー [eV]
      real(dp) :: sigma               ! 断面積 [m²]

      real(dp) :: E_clip

      ! 低エネルギーでの発散を防ぐ
      E_clip = max(energy, 0.01d0)

      ! Janev近似式
      ! σ_cx ≈ 3.0×10^-19 × (1 - 0.05×log10(E))^2 [m²]
      ! sigma = 3.0d-19 * (1.0d0 - 0.05d0 * log10(E_clip))**2
      sigma = 0.6937d-18 * (1.0d0 - 0.155d0 * log10(E_clip))**2 / (1.0d0 + 0.1112d-14 * E_clip**3.3d0)

      ! 最小値を設定
      sigma = max(sigma, 1.0d-22)

   end function sigma_cx

   !-----------------------------------------------------------------------------
   ! 弾性散乱断面積（外部CDFデータから）
   ! D⁰ + D⁺ → D⁰ + D⁺（散乱）
   !-----------------------------------------------------------------------------
   function sigma_el(energy) result(sigma)
      real(dp), intent(in) :: energy ! 衝突エネルギー（相対エネルギー）[eV]
      real(dp) :: sigma              ! 断面積 [m²]

      sigma = get_sigma_elastic(energy)

   end function sigma_el

   !-----------------------------------------------------------------------------
   ! 電離断面積（電子衝突電離）
   ! D⁰ + e⁻ → D⁺ + 2e⁻
   ! 電子エネルギーで判定（電子-中性粒子の相対エネルギー ≈ 電子エネルギー）
   !-----------------------------------------------------------------------------
   function sigma_ei(E_electron) result(sigma)
      real(dp), intent(in) :: E_electron  ! 電子エネルギー [eV]
      real(dp) :: sigma                   ! 断面積 [m²]

      real(dp), parameter :: sigma_0 = 1.0d-20  ! スケーリング係数 [m²・eV]

      if (E_electron <= E_IONIZE_THRESHOLD) then
         sigma = 0.0d0
      else
         ! ベーテ式に基づく簡易モデル
         sigma = sigma_0 * log(E_electron / E_IONIZE_THRESHOLD) / E_electron
         sigma = max(sigma, 0.0d0)
      end if

   end function sigma_ei

   !-----------------------------------------------------------------------------
   ! 電離レート係数 <σv>_ei の計算
   ! Maxwell分布の電子に対する平均レート係数
   ! 入力: T_e [eV]
   ! 出力: <σv>_ei [m³/s]
   !-----------------------------------------------------------------------------
   function ionization_rate_coeff(T_e) result(sigmav)
      real(dp), intent(in) :: T_e    ! 電子温度 [eV]
      real(dp) :: sigmav             ! <σv> [m³/s]

      real(dp), parameter :: m_e = 9.109d-31    ! 電子質量 [kg]
      real(dp), parameter :: eV_to_J = 1.60218d-19
      real(dp) :: v_th, E_mean

      ! 電子の熱速度
      v_th = sqrt(2.0d0 * T_e * eV_to_J / m_e)

      ! 平均電子エネルギー ≈ (3/2) * T_e
      E_mean = 1.5d0 * T_e

      ! <σv> ≈ σ(E_mean) * v_th
      sigmav = sigma_ei(E_mean) * v_th

   end function ionization_rate_coeff

   !-----------------------------------------------------------------------------
   ! 全断面積の計算（イオン衝突用）
   !-----------------------------------------------------------------------------
   subroutine get_all_cross_sections(energy, sig_cx, sig_el, sig_ei, sig_total)
      real(dp), intent(in)  :: energy ! 衝突エネルギー（相対エネルギー） [eV]
      real(dp), intent(out) :: sig_cx, sig_el, sig_ei, sig_total

      sig_cx = sigma_cx(energy)
      sig_el = sigma_el(energy)
      ! 電離周りは未実装として先延ばし
      sig_ei = 0.0d0
      sig_total = sig_cx + sig_el + sig_ei

   end subroutine get_all_cross_sections

end module cross_sections
