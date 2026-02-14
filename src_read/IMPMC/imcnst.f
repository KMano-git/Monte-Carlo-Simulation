!**********************************************************************
      subroutine imcnst
!**********************************************************************
      use cntcom, only : iransd
      use cunit,  only : n6
      implicit none
!
      integer  i
      real*8   x
! function
      real(8)    random, rnlast
!
      call trmark("imcnst","start")
!
!----------------------------------------------------------------------
!::pi,cev,cmp
!----------------------------------------------------------------------
      call phcnst
      call phycnst
!
!----------------------------------------------------------------------
!::random
!----------------------------------------------------------------------
      write(n6,'(2x,"imcnst  lransd =",i10)') iransd
      if( iransd > 0 ) then
        do i = 1, iransd
          x = random(0)
        enddo
      else
        x = rnlast(0)
      endif
!
      write(n6,'(4x,"initial random number iransd =",2f20.12)')
     >   x, rnlast(0)
!
!----------------------------------------------------------------------
!::table
!----------------------------------------------------------------------
      call mtbnrm
      call mtbtri
      write(n6,'(4x,"call mtbnrm : normal random number")')
      write(n6,'(4x,"call mtbtri : tri function")')
!
!----------------------------------------------------------------------
!::atomic hudrogen cross section
!----------------------------------------------------------------------
      write(n6,'(4x)')
      write(n6,'(4x,"call atomhy : ionization & recombination",
     >  " data of Hydrogen  degas2  00/08/29")')
      call atomhy

      call trmark("imcnst","return")

      return
      end
