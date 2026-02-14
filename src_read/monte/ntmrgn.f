!**********************************************************************
      subroutine ntmrgn
!**********************************************************************
!
!       it = itmax     
!       kreg(j,i) = 6,  mrgn(ic) = 6,  mrgnp(ic) = 7
!                   6 : core edge                  7 : core 
!
!       integrate in plasma volume    
!        do it = 2, itmax-1  
!
!           Note.  tube (itmax) is out of system
!
!----------------------------------------------------------------------
      use cntcom, only : iplx, iply, mcel, mrgn
      use cplmet, only : icel, itmpe, jcel, jtmin, jtmax
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer it, jt, j, i, ir, ir0, ic
!
      write(n6,'(2x,"*** ntmrgn ***")')
      write(n6,'(2x,"  jt  it    j   i    ic  ix   iy ir => ir2")')
!
      it = itmpe
      do jt = jtmin(it), jtmax(it)
      j = jcel(jt,it)
      i = icel(jt,it)
      ic = mcel(j,i)
      if( ic.le.0 ) cycle
      ir0 = mrgn(ic)
      ir  = 7
!
      mrgn(ic) = ir
      if( mod(jt,20).ne.0 ) cycle
      write(n6,'(2x,2i4,1x,2i4,1x,i5,1x,2i4,i3,2x,i3)')
     >  jt, it, j, i, ic, iplx(ic), iply(ic), ir0, mrgn(ic)
      enddo
!
      return
      end
