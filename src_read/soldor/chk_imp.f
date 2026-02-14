!**********************************************************************
      subroutine chk_imp
!**********************************************************************
      use catcom, only : catmz
      use cplimp, only : azmasl, ismaxl, wmc_nty
      use cunit,  only : n6
      implicit none

      integer nty

      write(n6,'(/2x,"impurity species: wmc_nty = ",i3)') wmc_nty

      do nty = 1, wmc_nty
! modified 2/2 lines cpadas to catcom by kamata 2022/04/15
!ik     write(n6,'(4x,i2," impurity : ",a2,"  ismax =",i3,"  azmas =",
!ik  >     f10.3)') nty, catmzL(nty), ismaxL(nty), azmasL(nty)
        write(n6,'(4x,i2," impurity : ",a2,"  ismax =",i3,"  azmas =",
     >     f10.3)') nty, catmz(nty), ismaxL(nty), azmasL(nty)
      enddo

      return
      end
