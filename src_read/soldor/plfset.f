!**********************************************************************
      subroutine plfset
!**********************************************************************
      use cplcom, only : cfdbg, cfopt
      use csonic, only : lfdbg, lfopt, lpst_bprf, nfdbg, nfopt
      use cunit,  only : n5, n6
      implicit none
!
!::local variables
! modified 1/3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  n, lenx
      integer  n
! function
      integer  lenx
!
!KH111108
      namelist /uplflg/ lfopt,lfdbg,lpst_bprf
!
      write(n6,'(/2x,"*** plfset ***")')
!
!----------------------------------------------------------------------
!::option flags
!----------------------------------------------------------------------
      cfopt( 1) = "pxstep on(1)/off(0)"
      lfopt( 1) = 1
!
      cfopt( 2) = "pystep on(1)/off(0)"
      lfopt( 2) = 1
!
      cfopt( 3) = "pycpnt grdient(1)/plate(0)"
      lfopt( 3) = 1
!
      nfopt = 3
!
!----------------------------------------------------------------------
!::debug flags
!----------------------------------------------------------------------
      cfdbg( 1) = "metlst on(1):list/on(2):cont/off(0)"
      lfdbg( 1) = 0
!
      cfdbg( 2) = "ploutp at start time  on(1)/off(0)"
      lfdbg( 2) = 0
!
      cfdbg( 3) = "ntgrid (neut=>soldor) on(1)/off(0)"
      lfdbg( 3) = 0
!
      cfdbg( 4) = "ntwall (wall-point)   on(1)/off(0)"
      lfdbg( 4) = 0
!
      cfdbg( 5) = "ntplas  init(1)/every time(2)/off(0)"
      lfdbg( 5) = 1
!
      cfdbg( 6) = "ntfldp  init(1)/every time(2)/off(0)"
      lfdbg( 6) = 1
!
      cfdbg( 7) = "ntoutp-tots init(1)/every time(2)/off(0)"
      lfdbg( 7) = 1
!
      cfdbg( 8) = "ntoutp-weit init(1)/every time(2)/off(0)"
      lfdbg( 8) = 2
!
      cfdbg( 9) = "ntoutp-den0 init(1)/every time(2)/off(0)"
      lfdbg( 9) = 1
!
      cfdbg(10) = "plmflx  init(1)/every time(2)/off(0)"
      lfdbg(10) = 2
!
      nfdbg = 10
!
!
      lpst_bprf(1:20) = 0
!
!----------------------------------------------------------------------
!::input /uplflg/
!----------------------------------------------------------------------
      rewind (n5)
      read(n5,uplflg)
!xx   write(n6,uplflg)
!
!----------------------------------------------------------------------
!::output
!----------------------------------------------------------------------
      write(n6,'(2x,"flag               meaning ")')
      do 310 n = 1, nfopt
      write(n6,'(2x,"lfopt(",i2,") = ",i3,2x,a)')
     >  n,lfopt(n),cfopt(n)(1:lenx(cfopt(n)))
 310  continue
      do 320 n = 1, nfdbg
      write(n6,'(2x,"lfdbg(",i2,") = ",i3,2x,a)')
     >   n,lfdbg(n),cfdbg(n)(1:lenx(cfdbg(n)))
 320  continue
!
!KH110912
      write(n6,'(2x,"lpst_bprf =",20i3)') lpst_bprf
!
      return
      end
