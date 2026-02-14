!===============================================================================
! Module: constants
! 物理定数と単位変換
!===============================================================================
module constants
   implicit none

   ! 倍精度定義
   integer, parameter :: dp = selected_real_kind(15, 307)

   ! 数学定数
   real(dp), parameter :: PI = 3.141592653589793d0

   ! 物理定数
   real(dp), parameter :: M_D  = 3.344d-27      ! 重水素質量 [kg]
   real(dp), parameter :: M_E  = 9.109d-31      ! 電子質量 [kg]
   real(dp), parameter :: EV   = 1.602d-19      ! 電子ボルト [J]
   real(dp), parameter :: K_B  = 1.381d-23      ! ボルツマン定数 [J/K]

   ! 単位変換係数
   real(dp), parameter :: EV_TO_J = EV          ! [eV] → [J]
   real(dp), parameter :: J_TO_EV = 1.0d0 / EV  ! [J] → [eV]
   real(dp), parameter :: CM2_TO_M2 = 1.0d-4    ! [cm^2] → [m^2]

end module constants
