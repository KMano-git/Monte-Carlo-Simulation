!%%mkspt +
!**********************************************************************
      subroutine ranini(ktyp)
!**********************************************************************
      use cunit, only : lmype, mype, n6
      implicit none
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer ktyp
      integer, intent(in) :: ktyp
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer iseed; real*8 xran,random
      integer iseed
      character cmsg*20
!
      if( ktyp.le.0 ) then
      iseed = 0
      cmsg  = "no seed"
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   xran  = 0.0d0
!
      elseif( ktyp.eq.1 ) then
      iseed = 761
      cmsg  = "same seed"
!
      elseif( ktyp.eq.2 ) then
      iseed =  mype*761 + 971
      cmsg  = "mype*761 + 971"
!
      elseif( ktyp.eq.3 ) then
      iseed = max0( mype-1, 0 )*761 + 971
      cmsg  = "(mype-1)*761 + 971"
!
      elseif( ktyp.eq.4 ) then
      iseed =  lmype*761 + 971
      cmsg  = "lmype*761 + 971"
!
      else
      iseed =  mype*mype + 1223
      cmsg  = "mype*mype + 1223"
      endif
!
      write(n6,'(4x,"*** ranini ***   mype =",i5,"  ktyp =",i2,
     > "  seed =",i10,2x,a)') mype, ktyp, iseed, cmsg
!
      if( ktyp.gt.0 ) then
        call srandom( iseed )
      endif
!
      return
      end
