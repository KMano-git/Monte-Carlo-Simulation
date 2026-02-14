!::@@A  dbg_snflx
!***********************************************************************
      subroutine dbg_snflx(cmsg)
!***********************************************************************
      use cplcom, only : snflx, TN0, tT0
      use csonic, only : itim
      use cunit,  only : n6
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   character :: cmsg*(*)
      character, intent(in) :: cmsg*(*)

      integer :: mj, ic
      integer, parameter :: icN0 = 67

      mj = index(cmsg,"-")
      ic = icN0

! modified 1/1 lines with TOPICS by kamata 2021/12/22
!ik   write(n6,'(1x,a3,1x,a8,1x,a12,i4,i5,1p8e14.6, 1p2e14.6)')
      write(n6,'(1x,a3,1x,a8,1x,a12,i11,i5,1p8e14.6, 1p2e14.6)')
     >   "@@A",cmsg(1:mj-1), cmsg(mj+1:), itim, ic, snflx(1:8),
     >    tN0(ic,1), tT0(ic,1)

      return
      end
