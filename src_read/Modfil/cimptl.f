! added replace all include files with module files by kamata 2021/08/18
      module cimptl
      implicit none

      integer :: kptl = 0, kc = 0, kz = 0
      real(8) :: wrr = 0.0_8, wzz = 0.0_8, wwg = 0.0_8, wvz = 0.0_8
     >    , wv2 = 0.0_8, wvr2 = 0.0_8, wdt = 0.0_8, wtim = 0.0_8
      real(8) :: wrr0 = 0.0_8, wzz0 = 0.0_8, wwt0 = 0.0_8, wroh0 = 0.0_8
      integer :: wic0 = 0, wis0 = 0, wic1 = 0, wis1 = 0

!::point source
      integer :: ic_00 = 0, is_00 = 0
      real(8) :: ri_00 = 0.0_8, zi_00 = 0.0_8, vv_00 = 0.0_8
     >    , vz_00 = 0.0_8

      end module cimptl
