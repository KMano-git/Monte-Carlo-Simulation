!**********************************************************************
      subroutine plmrgn
!**********************************************************************
!
!       it = itmax
!       kreg(j,i) = 6,  mrgn(ic) = 6,  mrgnp(ic) = 7
!                   6 : core edge                  7 : core
!
!       integrate in plasma volume
!        do it = 2, itmax-1
!
!           Note.  tube (itmax: irg2=7) is out of system
!
!----------------------------------------------------------------------
      use cntcom, only : iplx, iply, mcel, mrgn
      use cplmet, only : icel, itmpe, jcel, jtmax, jtmin
      use cplqcn, only : mrgnp
      use csize,  only : ndmc
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer it, jt, j, i, ir, ir2, ic
      integer it, jt, j, i, ir2, ic
!
      write(n6,'(/2x,"*** plmrgn ***")')
!
      do ic = 1, ndmc
      mrgnp(ic) = mrgn(ic)
      enddo
!
      write(n6,'(2x,"  jt  it    j   i    ic  ix   iy ir => ir2")')
      it = itmpe
      do jt = jtmin(it), jtmax(it)
      j = jcel(jt,it)
      i = icel(jt,it)
      ic = mcel(j,i)
      if( ic.le.0 ) cycle
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   ir  = mrgn(ic)
      ir2 = 7
!
      mrgnp(ic) = ir2
      if( mod(jt,20).ne.0 ) cycle
      write(n6,'(2x,2i4,1x,2i4,1x,i5,1x,2i4,i3,2x,i3)')
     >  jt, it, j, i, ic, iplx(ic), iply(ic), mrgn(ic), mrgnp(ic)
      enddo
!
      return
      end
