!***********************************************************************
      subroutine codrad(jc,rst,ren,rhf,n,nd,ic1,ic2)
!***********************************************************************
!
!   cell name in metric calculation
!
!                       +icaxs--------------------+
!                       !                         !
!                       !icmpe~~~~~~~icmax(jc)~~~~!
!      icwl2 /+^^^^^^^^^+                         +^^^^^^^^^+/
!            /!  (4)    !icmps  Main plasma (6)   !    (5)  !/
!            /!---------!-------------------------!---------!/
!      icspx /!  OUT    !                         !  IN     !/
!            /!   DIV   !      Scrape-Off         !   DIV   !/
!            /!  (1)    !                  (2)    !    (3)  !/
!   A  icwl1 /j+^^^^^^^j^^^^^^^^^^^^^icmin(jc)^^^^^j^^^^^^^^j j
! i !         c        c                           c        c c
! c !         d        x                           x        d m
!   !-->    1=p        p                           p        p=a
!     jc      1        1                           2        2 x
!
!                do jc = 1, jcmax                ir = kreg(jc,ic)
!                do ic = icmin(jc), icmax(jc)
!
!----------------------------------------------------------------------
      use cgdcom, only : grdx, grdy
      use cplmet, only : icspx, jcxp1, jcxp2
      implicit none
!
!::argument
      integer, intent(in)  :: jc, nd, ic1, ic2
      integer, intent(out) :: n
      real(8), intent(out) :: rst(nd), ren(nd), rhf(nd)
!
!::local variables
      integer  ic, kin
      integer  mx1, mx2, mx3, mx4, my1, my2, my3, my4
      real*8   sum, ra, rb, rc, rd, za, zb, zc, zd, rp, zp, rq, zq, dl
!
!
      do ic = 1, nd
      rst(ic) = 0.0d0
      ren(ic) = 0.0d0
      rhf(ic) = 0.0d0
      enddo
!
      sum = 0.0d0
      do ic = ic1, ic2
      call mcpnt(jc,ic,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
!
      ra = grdx(mx1,my1)
      rb = grdx(mx2,my2)
      rc = grdx(mx3,my3)
      rd = grdx(mx4,my4)
      za = grdy(mx1,my1)
      zb = grdy(mx2,my2)
      zc = grdy(mx3,my3)
      zd = grdy(mx4,my4)
!
      rp = 0.5d0*(ra+rb)
      zp = 0.5d0*(za+zb)
      rq = 0.5d0*(rc+rd)
      zq = 0.5d0*(zc+zd)
!
      dl = dsqrt( (rq-rp)**2+(zq-zp)**2)
      rst(ic) = sum
      ren(ic) = sum + dl
      rhf(ic) = sum + 0.5d0*dl
      sum = sum + dl
      enddo
!
      kin = 0
      if( jc.ge.(jcxp1+jcxp2)/2 ) kin = 1
      do ic = ic1, ic2
      rhf(ic) = ren(icspx) - rhf(ic)
      if( kin.eq.1 ) rhf(ic) = -rhf(ic)
      enddo
      n = ic2
!
      return
      end