!***********************************************************************
      subroutine pxbset(nt)
!***********************************************************************
!
!       boundary condition at divertor plate
!             j = 1,  jmax
!
!-----------------------------------------------------------------------
      use cplcom, only : c00, c11, cc, dd, ee, ff, nion
      use cplmet, only : jtmax
      use csize,  only :
      implicit none
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nt
      integer, intent(in) :: nt
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, i, jmax, nequ, jw, n, m
      integer  jmax, nequ, jw, n, m
!
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   it   = nt
!ik   i    = icel(1,it)
!ik   jmax = jtmax(it)
      jmax = jtmax(nt)
      nequ = 2*nion + 2
!
      jw = 1
      do 210 n = 1, nequ
      do 215 m = 1, nequ
      cc(n,m,jw) = c00
      dd(n,m,jw) = c00
      ee(n,m,jw) = c00
  215 continue
      dd(n,n,jw) = c11
      ff(n,jw)   = c00
  210 continue
!
!-----------------------------------------------------------------------
!::boundary condition j = jmax
!-----------------------------------------------------------------------
      jw = jmax
      do 230 n = 1, nequ
      do 235 m = 1, nequ
      cc(n,m,jw) = c00
      dd(n,m,jw) = c00
      ee(n,m,jw) = c00
  235 continue
      dd(n,n,jw) = c11
      ff(n,jw)   = c00
  230 continue
!
      return
      end
