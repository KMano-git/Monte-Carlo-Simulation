! added replace all include files with module files by kamata 2021/08/18
      module cimpuf
      use csize, only : ndgt
      use cimcom, only : impmc_max
      implicit none

!::puff
      integer, parameter :: ndpfz = 5
      integer :: pf_nty = 0, pf_nwl(ndpfz) = 0
      real(8) :: pf_rate(ndpfz) = 0.0_8, pf_mag(ndpfz) = 0.0_8
     >    , pf_px1(ndpfz) = 0.0_8, pf_px2(ndpfz) = 0.0_8
     >    , pf_py1(ndpfz) = 0.0_8, pf_py2(ndpfz) = 0.0_8

      integer :: pf_imx = 0
      integer, allocatable :: pf_iiw(:), pf_ity(:)
      real(8), allocatable :: pf_are(:), pf_deg(:), pf_flx(:)

!::time depending puff rate setting
      logical :: use_gaspufftime_IMP = .false.
      integer, parameter :: pufftime_max_IMP = 20 ! max size of pufftime_IMP
      real*8 :: pufftime_IMP(pufftime_max_IMP) = 0.0_8 ! time when change puff rate
      ! time depending puff rate 
      ! index: pf_mag_time(m,i_time,i_IMPMC)
      real*8 :: pf_mag_time(ndpfz,pufftime_max_IMP,impmc_max) = 0.0_8
      real*8 :: pf_rate_time(ndpfz,pufftime_max_IMP,impmc_max) = 0.0_8
      integer :: pufftime_size_IMP = 0

!::backflow from gate   2012/06/04
      integer :: bk_nty = 0, bk_npt = 0, bk_nwl(ndgt) = 0
      real(8) :: bk_mag(ndgt) = 0.0_8, bk_px1(ndgt) = 0.0_8
     >    , bk_px2(ndgt) = 0.0_8, bk_py1(ndgt) = 0.0_8
     >    , bk_py2(ndgt) = 0.0_8

      integer :: bk_imx = 0
      integer, allocatable :: bk_iiw(:), bk_ity(:)
      real(8), allocatable :: bk_are(:), bk_deg(:), bk_flx(:)

      integer :: lbkstw = 0

      end module cimpuf
