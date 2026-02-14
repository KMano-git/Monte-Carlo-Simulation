!***********************************************************************
      subroutine cdf_attr(ckey,ctyp,nrnk)
!***********************************************************************
!
!       ckey  : (i)  "xs_units"
!       ctyp  : (o)  "c"   ( char-typ )
!       nrnk  : (o)   2    0:y  1:x  2:w    -1:error
!
!-----------------------------------------------------------------------
      use cdfcom, only : hvar, hvrk, hvty, mhvr, nhvr
      implicit none
!
!::arguments
! modified 2/3 lines organize local variables and include files by kamata 2021/06/16
!ik   character ckey*(*), ctyp*1
!ik   integer   nrnk, icon
      character, intent(in)  :: ckey*(*)
      character, intent(out) :: ctyp*1
      integer,   intent(out) :: nrnk

!::local varaibales
! modified 1/3 lines organize local variables and include files by kamata 2021/06/16
!ik   integer i, mj, lenx, ino
      integer  i, mj
! function
      integer  lenx
!
      mj = lenx(ckey)
!
      ctyp = "x"
      nrnk = -1
!
      do i = 1, nhvr
      if( hvar(i)(1:mhvr(i)).eq.ckey(1:mj) ) then
        ctyp = hvty(i)(1:1)
        nrnk = hvrk(i)
        goto 100
      endif
      enddo
!
 100  continue
      return
      end
