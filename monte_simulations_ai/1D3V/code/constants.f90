!===============================================================================
! Module: constants
! 物理定数と単位変換
!===============================================================================
module constants
   implicit none

   integer, parameter :: dp = selected_real_kind(15, 307)

   real(dp), parameter :: PI = 3.141592653589793d0

   real(dp), parameter :: M_u = 1.660538921d-27
   real(dp), parameter :: M_P = 1.00782503223d0
   real(dp), parameter :: M_D = 2.01410177812d0
   real(dp), parameter :: M_D_kg = M_D * M_u
   real(dp), parameter :: M_E = 9.1093837139d-31
   real(dp), parameter :: EV  = 1.602176634d-19
   real(dp), parameter :: K_B = 1.380649d-23
   real(dp), parameter :: E_IONIZE_THRESHOLD = 13.598434605d0

   real(dp), parameter :: EV_TO_J = EV
   real(dp), parameter :: J_TO_EV = 1.0d0 / EV
   real(dp), parameter :: CM2_TO_M2 = 1.0d-4

end module constants
