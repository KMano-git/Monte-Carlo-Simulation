!***********************************************************************
      subroutine plvrdt
!***********************************************************************
      use cplcom, only : ddtq1, ddtq2, ddtq3, ddtq4, edtq1, edtq2, edtq3
     >    , edtq4, ibcyl, nion, q1a, q2a, q3, q4, relmt, vte, wq1a, wq2a
     >    , wq3, wq4
      use cplmet, only : icel, itmax, itpve, jcel, jtmax, jtmin
      use csize,  only : ndeq
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  itst, iten, ndim, it, ipbck jwst, jwen, jw, ipbc, jwst
      integer  itst, iten, it, jwen, jw, ipbc, jwst
      integer  j, i, ia
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8   zmax2, emax2, emax3, emax4, sumq1, sumq2, sumq3, sumq4
      real*8   sumq1, sumq2, sumq3, sumq4
      real*8   smdq1, smdq2, smdq3, smdq4, eps1, eps2, eps3, eps4
      integer  m, m1a, m2a, m3, m4
      real*8   qmax(ndeq), dnm1, dnm2, dnm3, dnm4
      real*8   telmt
!
!::region
      itst = 2
      iten = itmax - 1
      telmt = 1.0d0
!
!::zero clear
      edtq1 = -1.0d20
      edtq2 = -1.0d20
      edtq3 = -1.0d20
      edtq4 = -1.0d20
! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik   emax2 = -1.0d20
!ik   emax3 = -1.0d20
!ik   emax4 = -1.0d20
      sumq1 = 0.0d0
      sumq2 = 0.0d0
      sumq3 = 0.0d0
      sumq4 = 0.0d0
      smdq1 = 0.0d0
      smdq2 = 0.0d0
      smdq3 = 0.0d0
      smdq4 = 0.0d0
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   ndim  = ndx*ndy
!
!::loop it
      do it = itst, iten
      if( it.eq.itpve )   cycle
      if( it.eq.itpve-1 ) cycle
!xx      if( chkoff(it).eq.1 ) cycle  ! KSFUJI
      ipbc = ibcyl(it)
      jwst = jtmin(it) + 1
      jwen = jtmax(it) - 1
      if( ipbc.eq.1 ) then
      jwst = jtmin(it) + 1
      jwen = jtmax(it) - 1
      endif
!
!::max value in flux tube
      do m = 1, ndeq
      qmax(m) = 0.0d0
      enddo
!
      do jw = jwst, jwen
      j = jcel(jw,it)
      i = icel(jw,it)
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      qmax(m1a) = dmax1( qmax(m1a), dabs(q1a(j,i,ia)) )
      qmax(m2a) = dmax1( qmax(m2a), dabs(q2a(j,i,ia)) )
      enddo
      m3 = 2*nion + 1
      m4 = 2*nion + 2
      qmax(m3) = dmax1( qmax(m3), dabs(q3(j,i)) )
      qmax(m4) = dmax1( qmax(m4), dabs(q4(j,i)) )
      enddo
!
!::variation
      do jw = jwst, jwen
      j = jcel(jw,it)
      i = icel(jw,it)
      if( vte(j,i).le.telmt ) cycle
      do ia = 1, nion
      m1a  = 2*ia - 1
      dnm1 = dabs(q1a(j,i,ia)) + relmt*qmax(m1a)
      eps1 = dabs((q1a(j,i,ia)-wq1a(j,i,ia)))/dnm1
      edtq1= dmax1(edtq1,eps1)
      m2a  = 2*ia
      dnm2 = dabs(q2a(j,i,ia)) + relmt*qmax(m2a)
      eps2 = dabs((q2a(j,i,ia)-wq2a(j,i,ia)))/dnm2
      edtq2= dmax1(edtq2,eps2)
      sumq1= sumq1 +  q1a(j,i,ia)**2
      sumq2= sumq2 +  q2a(j,i,ia)**2
      smdq1= smdq1 + (q1a(j,i,ia)-wq1a(j,i,ia))**2
      smdq2= smdq2 + (q2a(j,i,ia)-wq2a(j,i,ia))**2
      enddo
      m3   = 2*nion + 1
      dnm3 = dabs(q3(j,i)) + relmt*qmax(m3)
      eps3 = dabs((q3(j,i)-wq3(j,i)))/dnm3
      edtq3= dmax1(edtq3,eps3)
      m4   = 2*nion + 2
      dnm4 = dabs(q4(j,i)) + relmt*qmax(m4)
      eps4 = dabs((q4(j,i)-wq4(j,i)))/dnm4
      edtq4= dmax1(edtq4,eps4)
      sumq3= sumq3 +  q3(j,i)**2
      sumq4= sumq4 +  q4(j,i)**2
      smdq3= smdq3 + (q3(j,i)-wq3(j,i))**2
      smdq4= smdq4 + (q4(j,i)-wq4(j,i))**2
      enddo
      enddo
      ddtq1= dsqrt(smdq1)/dsqrt(sumq1)
      ddtq2= dsqrt(smdq2)/dsqrt(sumq2)
      ddtq3= dsqrt(smdq3)/dsqrt(sumq3)
      ddtq4= dsqrt(smdq4)/dsqrt(sumq4)
!
      return
      end
!
!***********************************************************************
      subroutine plvrdt_out(nf)
!***********************************************************************
      use cplcom, only : ddtq1, ddtq2, ddtq3, ddtq4, edtq1, edtq2, edtq3
     >    , edtq4, q1a, q2a, q3, q4, wq1a, wq2a, wq3, wq4
      use cplmet, only : icmax, icmin, jcxp1, jtmax, jtmin
      use csize,  only : ndx
      use csonic, only : dtim, time
      use cunit,  only : lmspe, lmype
      implicit none
!
      integer  nf
      integer  ia, jst, jen, ist, ien, j, i, ih
      real*8   eps, zh
!
      character clin*(ndx), clev*(ndx)
!
      if( nf.le.0 ) return
      if( lmype.ne.lmspe ) return
!
!::species
      ia = 1
!
      write(nf,'(/2x,"***  plvrdt  ***   time,dtim =",1p2e12.3
     >  ,"  ia =",i3)') time, dtim, ia
!
!::zero divide
      if( edtq1.le.0.0 .or. edtq2.le.0.0 .or. edtq3.le.0.0 .or.
     >    edtq4.le.0.0 ) then
      write(nf,'(2x,"return because of edtq = 0.0 ",1p4e12.3)')
     >  edtq1,edtq2,edtq3,edtq4
      return
      endif
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
      write(nf,'(2x,"== Q1a ==  ddtq1,edtq1 =",1p2e10.2)') ddtq1,edtq1
      write(nf,'(7x,a)') clin(1:jen)
      do 210 i = ist, ien
      do 220 j = jst, jen
      eps = 0.0d0
      if( q1a(j,i,ia).ne.0.0d0 )
     >eps = dabs((q1a(j,i,ia)-wq1a(j,i,ia))/q1a(j,i,ia))
      zh  = eps/edtq1
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
      write(nf,'(2x,"== Q2a ==  ddtq2,edtq2 =",1p2e10.2)') ddtq2,edtq2
      write(nf,'(7x,a)') clin(1:jen)
      do 310 i = ist, ien
      do 320 j = jst, jen
      eps = 0.0d0
      if( q2a(j,i,ia).ne.0.0d0 )
     >eps = dabs((q2a(j,i,ia)-wq2a(j,i,ia))/q2a(j,i,ia))
      zh  = eps/edtq2
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
      write(nf,'(2x,"== Q3 ==  ddtq3,edtq3 =",1p2e10.2)') ddtq3,edtq3
      write(nf,'(7x,a)') clin(1:jen)
      do 410 i = ist, ien
      do 420 j = jst, jen
      eps = 0.0d0
      if( q3(j,i).ne.0.0d0 )
     >eps = dabs((q3(j,i)-wq3(j,i))/q3(j,i))
      zh  = eps/edtq3
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
      write(nf,'(2x,"== Q4 ==  ddtq4,edtq4 =",1p2e10.2)') ddtq4,edtq4
      write(nf,'(7x,a)') clin(1:jen)
      do 510 i = ist, ien
      do 520 j = jst, jen
      eps = 0.0d0
      if( q4(j,i).ne.0.0d0 )
     >eps = dabs((q4(j,i)-wq4(j,i))/q4(j,i))
      zh  = eps/edtq4
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
