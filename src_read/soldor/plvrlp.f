!***********************************************************************
      subroutine plvrlp
!***********************************************************************
!
!        exclude swall & pwall   99/11/03
!
      use cplcom, only : dlpq1, dlpq2, dlpq3, dlpq4, dq1a, dq2a, dq3
     >    , dq4, elpq1, elpq2, elpq3, elpq4, ibcyl, nion, q1a, q2a, q3
     >    , q4, vva
      use cplmet, only : icel, itmax, itpve, jcel, jtmax, jtmin
      implicit none
!
      integer  itst, iten, it, ipbc, jwst, jwen, jw, j, i, ia
      real*8   sumq1, sumq2, sumq3, sumq4
      real*8   smdq1, smdq2, smdq3, smdq4
!
!::region
      itst = 2
      iten = itmax - 1
!xx      iten = itpve - 1   ! <=== /09/15
!
!::zero clear
      elpq1 = -1.0d20
      elpq2 = -1.0d20
      elpq3 = -1.0d20
      elpq4 = -1.0d20
      sumq1 = 0.0d0
      sumq2 = 0.0d0
      sumq3 = 0.0d0
      sumq4 = 0.0d0
      smdq1 = 0.0d0
      smdq2 = 0.0d0
      smdq3 = 0.0d0
      smdq4 = 0.0d0
!
      do 310 it = itst, iten
      if( it.eq.itpve ) goto 310
      ipbc = ibcyl(it)
      jwst = jtmin(it)
      jwen = jtmax(it)
      jwst = jwst + 1
      jwen = jwen - 1      ! Optchek  00/9/11
      if( ipbc.eq.1 ) then
      jwst = jtmin(it) + 1
      jwen = jtmax(it) - 1
      endif
!
      do 320 jw = jwst, jwen
      j = jcel(jw,it)
      i = icel(jw,it)
      do 325 ia = 1, nion
      elpq1 = dmax1( elpq1, dabs(dq1a(j,i,ia)/q1a(j,i,ia)) )
      if( dabs(vva(j,i,ia)).gt.1.0d2 )
     >elpq2 = dmax1( elpq2, dabs(dq2a(j,i,ia)/q2a(j,i,ia)) )
      sumq1 = sumq1 +  q1a(j,i,ia)**2
      sumq2 = sumq2 +  q2a(j,i,ia)**2
      smdq1 = smdq1 + dq1a(j,i,ia)**2
      smdq2 = smdq2 + dq2a(j,i,ia)**2
 325  continue
      elpq3 = dmax1( elpq3, dabs(dq3(j,i)/q3(j,i)) )
      elpq4 = dmax1( elpq4, dabs(dq4(j,i)/q4(j,i)) )
      sumq3 = sumq3 +  q3(j,i)**2
      sumq4 = sumq4 +  q4(j,i)**2
      smdq3 = smdq3 + dq3(j,i)**2
      smdq4 = smdq4 + dq4(j,i)**2
 320  continue
 310  continue
!
      dlpq1 = dsqrt(smdq1)/dsqrt(sumq1)
      dlpq2 = dsqrt(smdq2)/dsqrt(sumq2)
      dlpq3 = dsqrt(smdq3)/dsqrt(sumq3)
      dlpq4 = dsqrt(smdq4)/dsqrt(sumq4)
!
      return
      end
!
!***********************************************************************
      subroutine plwrlp(nf)
!***********************************************************************
      use cplcom, only : dlpq1, dlpq2, dlpq3, dlpq4,  dq1a, dq2a, dq3
     >     , dq4, elpq1, elpq2, elpq3, elpq4, q1a, q2a, q3, q4
      use cplmet, only : icmax, icmin, jcxp1, jtmax, jtmin
      use csize,  only : ndx
      use csonic, only : dtim, time
      use cunit,  only : lmspe, lmype
      implicit none
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer nf
      integer, intent(in) :: nf
!
!::local variables
      integer  ia, jst, jen, ist, ien, j, i, ih
      real*8   eps, zh
      character clin*(ndx), clev*(ndx)
!
      if( nf.le.0 ) return
      if( lmype.ne.lmspe ) return
!
!::species
      ia = 1
!
      write(nf,'(/2x,"***  plvrlp  ***   time,dtim =",1p2e12.3
     >  ,"  ia =",i3)') time, dtim, ia
!
!::region
      jst = jtmin(1)
      jen = jtmax(1)
      ist = icmin(jcxp1+1)
      ien = icmax(jcxp1+1)
!
!::clin
      clin = " "
      do 110 j = 1, jen
      clin(j:j) = "-"
      if( mod(j,10).eq.0 ) clin(j:j) = "I"
 110  continue
      clin(j:j) = "I"
!
!::Q1a
      write(nf,'(2x,"== Q1a ==  dlpq1,elpq1 =",1p2e10.2)') dlpq1,elpq1
      write(nf,'(7x,a)') clin(1:jen)
      do 210 i = ist, ien
      do 220 j = jst, jen
      eps = 0.0d0
      if( q1a(j,i,ia).ne.0.0d0 )
     >eps = dabs(dq1a(j,i,ia)/q1a(j,i,ia))
      zh  = eps/elpq1
      if( zh.le.0.001 ) then
        clev(j:j) = "."
      elseif( zh.le.0.01 ) then
        clev(j:j) = "-"
      elseif( zh.le.0.1 ) then
        clev(j:j) = "+"
      else
        zh = dmin1(zh*10.0d0+1.0d0,9.0d0)
        ih = int(zh)
        write(clev(j:j),'(i1)') ih
      endif
      if( eps.eq.0.0d0 ) clev(j:j) = "0"
 220  continue
      write(nf,'(2x,i3,2x,a,2x,i3)') i,clev(1:jen),i
 210  continue
!
!::Q2a
      write(nf,'(2x,"== Q2a ==  dlpq2,elpq2 =",1p2e10.2)') dlpq2,elpq2
      write(nf,'(7x,a)') clin(1:jen)
      do 310 i = ist, ien
      do 320 j = jst, jen
      eps = 0.0d0
      if( q2a(j,i,ia).ne.0.0d0 )
     >eps = dabs(dq2a(j,i,ia)/q2a(j,i,ia))
      zh  = eps/elpq2
      if( zh.le.0.001 ) then
        clev(j:j) = "."
      elseif( zh.le.0.01 ) then
        clev(j:j) = "-"
      elseif( zh.le.0.1 ) then
        clev(j:j) = "+"
      else
        zh = dmin1(zh*10.0d0+1.0d0,9.0d0)
        ih = int(zh)
        write(clev(j:j),'(i1)') ih
      endif
      if( eps.eq.0.0d0 ) clev(j:j) = "0"
 320  continue
      write(nf,'(2x,i3,2x,a,2x,i3)') i,clev(1:jen),i
 310  continue
!
!::Q3
      write(nf,'(2x,"== Q3 ==  dlpq3,elpq3 =",1p2e10.2)') dlpq3,elpq3
      write(nf,'(7x,a)') clin(1:jen)
      do 410 i = ist, ien
      do 420 j = jst, jen
      eps = 0.0d0
      if( q3(j,i).ne.0.0d0 )
     >eps = dabs(dq3(j,i)/q3(j,i))
      zh  = eps/elpq3
      if( zh.le.0.001 ) then
        clev(j:j) = "."
      elseif( zh.le.0.01 ) then
        clev(j:j) = "-"
      elseif( zh.le.0.1 ) then
        clev(j:j) = "+"
      else
        zh = dmin1(zh*10.0d0+1.0d0,9.0d0)
        ih = int(zh)
        write(clev(j:j),'(i1)') ih
      endif
      if( eps.eq.0.0d0 ) clev(j:j) = "0"
 420  continue
      write(nf,'(2x,i3,2x,a,2x,i3)') i,clev(1:jen),i
 410  continue
!
!::Q4
      write(nf,'(2x,"== Q4 ==  dlpq4,elpq4 =",1p2e10.2)') dlpq4,elpq4
      write(nf,'(7x,a)') clin(1:jen)
      do 510 i = ist, ien
      do 520 j = jst, jen
      eps = 0.0d0
      if( q4(j,i).ne.0.0d0 )
     >eps = dabs(dq4(j,i)/q4(j,i))
      zh  = eps/elpq4
      if( zh.le.0.001 ) then
        clev(j:j) = "."
      elseif( zh.le.0.01 ) then
        clev(j:j) = "-"
      elseif( zh.le.0.1 ) then
        clev(j:j) = "+"
      else
        zh = dmin1(zh*10.0d0+1.0d0,9.0d0)
        ih = int(zh)
        write(clev(j:j),'(i1)') ih
      endif
      if( eps.eq.0.0d0 ) clev(j:j) = "0"
 520  continue
      write(nf,'(2x,i3,2x,a,2x,i3)') i,clev(1:jen),i
 510  continue
!
      return
      end
