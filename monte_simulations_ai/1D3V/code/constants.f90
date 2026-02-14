!===============================================================================
! Module: constants
! 物理定数の定義
!===============================================================================
module constants
  implicit none
  
  ! 精度定義
  integer, parameter :: dp = selected_real_kind(15, 307)
  
  ! 物理定数
  real(dp), parameter :: AMU = 1.66054d-27        ! 原子質量単位 [kg]
  real(dp), parameter :: M_D = 2.014d0 * AMU      ! 重水素質量 [kg]
  real(dp), parameter :: EV_TO_J = 1.60218d-19    ! eV → J 変換
  real(dp), parameter :: J_TO_EV = 1.0d0 / EV_TO_J
  real(dp), parameter :: CM2_TO_M2 = 1.0d-4       ! cm² → m²
  real(dp), parameter :: KB = 1.38065d-23         ! ボルツマン定数 [J/K]
  real(dp), parameter :: PI = 3.14159265358979d0
  
  ! 電離閾値エネルギー
  real(dp), parameter :: E_IONIZE_THRESHOLD = 13.6d0  ! [eV] 水素電離エネルギー

end module constants
