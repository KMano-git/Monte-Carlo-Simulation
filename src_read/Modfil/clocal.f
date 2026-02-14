! added replace all include files with module files by kamata 2021/08/18
! local common in IMPMC and IMPMC_TD
      module clocal
      implicit none

! /cobi/ in gmequ/bicub.f, spln2d.f
      real(8) :: xim1 = 0.0_8, yjm1 = 0.0_8

! /com_bicub/ in gmequ/spln2d.f
      integer, parameter :: ndspx = 400, ndspy = 513
     >    , ndspc = ndspx * ndspy
      real(8) :: dbx = 0.0_8, dby = 0.0_8, tbcf(ndspc,16) = 0.0_8
     >    , tbx(ndspx) = 0.0_8, tby(ndspy) = 0.0_8
      integer :: ntx = 0, nty = 0, tbic(ndspx,ndspy) = 0

! /com_set_roh/ in gmequ/set_roh.f
      real(8), allocatable :: tpsi(:), troh(:), tvol(:)

! /com_imhist/ in IMPMC/imhist.f
      real(8) :: xcpum = 0.0_8, xmype = 0.0_8, xpmax = 0.0_8
     >    , xstep = 0.0_8, xtout = 0.0_8
      real(8) :: xcnpt(0:10) = 0.0_8, xcnwt(0:10) = 0.0_8
     >    , xrgpa(0:10) = 0.0_8, xrgpc(0:10) = 0.0_8

! /com_sput/ in IMPMC_com/set_tomp.f, sputN.f
      integer, parameter ::  msptin = 11
      character(2) :: csptin(msptin) = ' '
      real(8) :: es = 0.0_8
      real(8), dimension(msptin) :: etf = 0.0_8, eth = 0.0_8, fy = 0.0_8
     >    , qy = 0.0_8
! es : surface binding energy

! added 8 linex dynamic allocation of arrays by kamata 2022/05/29
! local variables
! in ntpump1.f
      real(8), allocatable :: vden0(:,:), vpre0(:,:), vdeng(:,:)
     >    , vpreg(:,:)

! in plmanprf.f
      real(8), allocatable :: man_na(:,:), man_ne(:), man_ni(:)
     >    , man_te(:), man_ti(:)

      end module clocal
