!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   subroutine pxsorc(nt)
      subroutine pxsorc(it)
!***********************************************************************
!     dq/dt + div(F) + Hp = S
!
!     dlt( dq/dt + divF + Hp - S ) = -( dq/dt + div(F) + Hp - S )^l
!
!       S = Sc + Sv*q
!
!      Cj*q(j-1) + Dj*q(j) + Ej*q(j+1) = fj
!
!-----------------------------------------------------------------------
      use cplcom, only : ama, clpa, ff, nion, q2a, q3, q4, ss, ssnc
     >    , ssnv, sspc, sspv, swec, swev, swic, swiv, vna
      use cplmet, only : hvol, icel, jcel, jtmax
      use cplqcn, only : qvl_sc
      use csize, only :
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nt
      integer, intent(in) :: it
!
!::local variables
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, i, jmax, jmax1, ia, jw, j, m1a, m2a, m3, m4
!ik   real*8   zfac, zss1, zdd1, zws1, zss4, zdd4, zws4
      integer  i, jmax, jmax1, ia, jw, j, m1a, m2a, m3, m4
!
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   it = nt
      i  = icel(1,it)
      jmax  = jtmax(it)
      jmax1 = jmax - 1
!
      m3 = 2*nion + 1
      m4 = 2*nion + 2
!
      do 110 ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
      do 120 jw = 2, jmax1
      j = jcel(jw,it)
!
!--linearlized source terms   !  %%% 2002/11/13 (pxsorc)
!
!::q1a
      ss(m1a,m1a,jw) = ss(m1a,m1a,jw) - ssnv(j,i,ia)*hvol(j,i)
      ff(m1a,jw) = ff(m1a,jw)
     >  + ama(ia)*(ssnc(j,i,ia)+ssnv(j,i,ia)*vna(j,i,ia))*hvol(j,i)
!
!::q2a
      ss(m2a,m2a,jw) = ss(m2a,m2a,jw) - sspv(j,i,ia)*hvol(j,i)
      ff(m2a,jw) = ff(m2a,jw)
     >   + (sspc(j,i,ia)+sspv(j,i,ia)*q2a(j,i,ia))*hvol(j,i)
!
!::q3a
      ss(m3,m3,jw)  = ss(m3,m3,jw) - swiv(j,i)*hvol(j,i)*clpa(ia)
      ff(m3,jw) = ff(m3,jw)
     >  + (swic(j,i)+swiv(j,i)*q3(j,i))*hvol(j,i)*clpa(ia)
!
!::q4
      ss(m4,m4,jw)  = ss(m4,m4,jw)  - swev(j,i)*hvol(j,i)*clpa(ia)
      ff(m4,jw) = ff(m4,jw)
     >  + (swec(j,i)+swev(j,i)*q4(j,i))*hvol(j,i)*clpa(ia)
!
!::qcon ! %%% 2002/11/02
      qvl_sc(j,i,m1a) =
     >     ama(ia)*(ssnc(j,i,ia)+ssnv(j,i,ia)*vna(j,i,ia))*hvol(j,i)
      qvl_sc(j,i,m2a) =
     >        (sspc(j,i,ia)+sspv(j,i,ia)*q2a(j,i,ia))*hvol(j,i)
      qvl_sc(j,i,m3)  =
     >        (swic(j,i)+swiv(j,i)*q3(j,i))*hvol(j,i)*clpa(ia)
      qvl_sc(j,i,m4)  =
     >        (swec(j,i)+swev(j,i)*q4(j,i))*hvol(j,i)*clpa(ia)
!
 120  continue
 110  continue
!

      return
      end
