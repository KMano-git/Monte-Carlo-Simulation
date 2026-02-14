!***********************************************************************
      subroutine pybcon2(j)
!***********************************************************************
!
!   Yacobian at boundary
!
!   boundary condition at divertor plate
!      i=1  (sol wall) / i=imax (prv wall)
!
!  (1) na   na = na_wall
!  (2) va   va = 0
!  (3) Ti   Ti = Ti_wall
!  (4) Te   Te = Te_wall
!
!-----------------------------------------------------------------------
      use cplcom, only : dd, ff, nion, q1a, q2a, q3, q4, qb1a, qb2a
     >    , qb3, qb4, qw1a, qw2a, qw3, qw4
      use cplmet, only : icmax, jcxp1, jcxp2, kcn
      use cunit,  only : n6
      implicit none
!
      integer, intent(in) :: j
!
!::local variables
      integer  imax, iw, kbc
      integer  m1a, m2a, m3, m4, ia
!
      integer  lp; data lp/0/; save lp
!
      if( lp.eq.0 ) then
        write(n6,'(/2x,"***  pybcon2  ***  constant ")')
        lp = 1
      endif
!
!::index
      imax  = icmax(j)
!
      m3   = 2*nion + 1
      m4   = 2*nion + 2
!
!-----------------------------------------------------------------------
!::boundary condition i = 1 (sol-wall)
!-----------------------------------------------------------------------
      iw  = 1
      kbc = 1
!
      do ia = 1, nion
        m1a = 2*ia - 1
        m2a = 2*ia
        dd(m1a,m1a,iw) = 1.0d0
        dd(m2a,m2a,iw) = 1.0d0
        ff(m1a,iw) = qw1a(ia,kbc) - q1a(j,iw,ia)
        ff(m2a,iw) = qw2a(ia,kbc) - q2a(j,iw,ia)
      enddo
      dd(m3,m3,  iw) = 1.0d0
      dd(m4,m4,  iw) = 1.0d0
      ff(m3,iw) = qw3(kbc) - q3(j,iw)
      ff(m4,iw) = qw4(kbc) - q4(j,iw)
!
!-----------------------------------------------------------------------
!::boundary condition i = imax (prv-wall)
!-----------------------------------------------------------------------
      if( j.gt.jcxp1 .and. j.lt.jcxp2 ) goto 300
!
      iw  = imax
      kbc = 2
!
      do ia = 1, nion
        m1a = 2*ia - 1
        m2a = 2*ia
        dd(m1a,m1a,iw) = 1.0d0
        dd(m2a,m2a,iw) = 1.0d0
        ff(m1a,iw) = qw1a(ia,kbc) - q1a(j,iw,ia)
        ff(m2a,iw) = qw2a(ia,kbc) - q2a(j,iw,ia)
      enddo
      dd(m3,m3,  iw) = 1.0d0
      dd(m4,m4,  iw) = 1.0d0
      ff(m3,iw) = qw3(kbc) - q3(j,iw)
      ff(m4,iw) = qw4(kbc) - q4(j,iw)
!
      return
!
!-----------------------------------------------------------------------
!::boundary condition i = imax (edge plasma)
!-----------------------------------------------------------------------
 300  continue
      kbc = kcn
      iw  = imax
      do ia = 1, nion
        m1a = 2*ia - 1
        m2a = 2*ia
        dd(m1a,m1a,iw) = 1.0d0
        dd(m2a,m2a,iw) = 1.0d0
        ff(m1a,iw) = qb1a(j,ia,kbc)-q1a(j,iw,ia)
        ff(m2a,iw) = qb2a(j,ia,kbc)-q2a(j,iw,ia)
      enddo
      dd(m3,m3,  iw) = 1.0d0
      dd(m4,m4,  iw) = 1.0d0
      ff(m3,iw) = qb3(j,kbc)-q3(j,iw)
      ff(m4,iw) = qb4(j,kbc)-q4(j,iw)
!
      return
      end
