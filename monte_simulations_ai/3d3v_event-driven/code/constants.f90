!===============================================================================
! Module: constants
! 物理定数と単位変換
!===============================================================================
module constants
   implicit none

   !倍精度定義
   integer, parameter :: dp = selected_real_kind(15, 307)

   !数学定数
   real(dp), parameter :: PI = 3.141592653589793d0

   !物理定数
   real(dp), parameter :: M_u = 1.660538921d-27    !原子質量単位 [kg]
   real(dp), parameter :: M_P = 1.00782503223d0    !水素原子質量比 [u]
   real(dp), parameter :: M_D = 2.01410177812d0    !重水素原子質量比 [u]
   real(dp), parameter :: M_D_kg = M_D * M_u       !重水素原子質量 [kg]
   real(dp), parameter :: M_E = 9.1093837139d-31   !電子質量 [kg]
   real(dp), parameter :: EV  = 1.602176634d-19    !電子ボルト [J],電気素量 [C]
   real(dp), parameter :: K_B = 1.380649d-23       !ボルツマン定数 [J/K]
   real(dp), parameter :: E_IONIZE_THRESHOLD = 13.598434605d0 !電離閾値エネルギー [eV]

   !単位変換係数
   real(dp), parameter :: EV_TO_J = EV             ![eV] → [J]
   real(dp), parameter :: J_TO_EV = 1.0d0 / EV     ![J] → [eV]
   real(dp), parameter :: CM2_TO_M2 = 1.0d-4       ![cm^2] → [m^2]

end module constants
