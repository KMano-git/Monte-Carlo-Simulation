!**********************************************************************
      subroutine ntl_recd
!**********************************************************************
!
!        recd : record of ntl-calculation
!
!----------------------------------------------------------------------
      use cntctl, only : itmnt, itmntb, nclnt
      use cputim, only : cntcal, cntpls, cntrev, cntsnd, cntsum
      use csize,  only : ndmsr
      use csonic, only : itim, lntl, time
      use cunit,  only : n6
      implicit none

!::local varaiables
      integer  i, iclnt

!::clear cpu
      cntsnd = 0.0d0
      cntrev = 0.0d0
      cntpls = 0.0d0
      do i = 0, ndmsr
      cntcal(i) = 0.0d0
      cntsum(i) = 0.0d0
      enddo

      if( lntl.eq.0 ) goto 100
      if( lntl.ne.3 ) goto 100

!::caluclation step
! modified 2/1 lines organize local variables and include files by kamata 2021/06/16
!ik   iclnt = nclnt
!ik   iclnt = iclnt + 1
      iclnt = nclnt + 1
!
!::debug write
      write(n6,'(/30("="),"  START  ntl_cal  itim =",i6,"  time =",
     > 1pe14.6,"  iclnt,itmnt =",3i4)') itim,time,iclnt,itmnt,itmntb
!
!::number of calculation
      itmntb = itmnt
      itmnt  = 1
!
 100  continue
      return
      end
