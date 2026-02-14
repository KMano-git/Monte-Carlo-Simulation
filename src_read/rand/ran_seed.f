!***********************************************************************
      subroutine ran_init
!***********************************************************************
      use csonic, only : lrand
      use cunit,  only : lmype, n6
      implicit none

      integer :: iseed
      real(8) :: xran
! function
      real(8)    random

      if( lrand /= 0 ) return

      iseed = lmype*17 + 27
      call srandom(iseed)

      xran = random(0)

      write(n6,'(/2x,"*** ran_init  ***  lmype =",i5,2x,a,i8,
     >  "  xran =",f12.8)') lmype, "iseed=ipe*17+27 =", iseed, xran

      return
      end
