!===============================================================================
! Module: constants
! physical constants and conversion factors
!===============================================================================
module constants
   implicit none

   ! double precision
   integer, parameter :: dp = selected_real_kind(15, 307)

   ! math constants
   real(dp), parameter :: PI = 3.141592653589793d0

   ! physical constants
   real(dp), parameter :: M_u = 1.660538921d-27    !atomic mass unit [kg]
   real(dp), parameter :: M_P = 1.00782503223d0    !hydrogen atomic mass ratio [u]
   real(dp), parameter :: M_D = 2.01410177812d0    !deuterium atomic mass ratio [u]
   real(dp), parameter :: M_D_kg = M_D * M_u       !deuterium atomic mass [kg]
   real(dp), parameter :: M_E = 9.1093837139d-31   !electron mass [kg]
   real(dp), parameter :: EV  = 1.602176634d-19    !electron volt [J], elementary charge [C]
   real(dp), parameter :: K_B = 1.380649d-23       !boltzman constant [J/K]
   real(dp), parameter :: E_IONIZE_THRESHOLD = 13.598434605d0 !ionization threshold energy [eV]

   ! conversion factors
   real(dp), parameter :: EV_TO_J = EV             ![eV] → [J]
   real(dp), parameter :: J_TO_EV = 1.0d0 / EV     ![J] → [eV]
   real(dp), parameter :: CM2_TO_M2 = 1.0d-4       ![cm^2] → [m^2]

end module constants
