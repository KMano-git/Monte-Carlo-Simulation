!**********************************************************************
      subroutine dspmtn(ip,ic,ln,ri,zi)
!**********************************************************************
!
!     random displacement to include neutral hydorcarbon
!
!----------------------------------------------------------------------
      use cimcom, only : dlmtn
      use cntcom, only : fnfi, tcfi, tsfi
      use cunit,  only : n6
      implicit none
!
!::argument
      integer, intent(in)    :: ip ! dummy
      integer, intent(inout) :: ic
      integer, intent(out)   :: ln
      real(8), intent(inout) :: ri, zi
!
!::local variables for fndway
      integer ndmv; parameter (ndmv=201)
      integer nmv, icmv(ndmv)
      real*8  dtmv(ndmv)
!
!::local variables
      integer    ifi, ic2, ln2, lp
      real(8)    ri2, zi2
! function
      integer    imox, imoy
      real(8)    random
!
      lp = 0
 100  continue
      lp = lp + 1
      if( lp.gt.50 ) goto 910
      ifi = int(fnfi*random(0) + 1.0d0)
      ri2 = ri + dlmtn*tcfi(ifi)
      zi2 = zi + dlmtn*tsfi(ifi)
      ln2 = 0
!
      call fndway(ri,zi,ri2,zi2,ic,ic2, ndmv,nmv,icmv,dtmv)
      if( ic2.le.0 ) goto 100
!
      ri = ri2
      zi = zi2
      ic = ic2
      ln = ln2
      return
!
 910  continue
      write(n6,'(2x,"Warning  sub. dspmtn  lp.gt.50 ",
     >   2f9.4,i8,2i5)') ri,zi,ic,imox(ic),imoy(ic)
      return
      end
