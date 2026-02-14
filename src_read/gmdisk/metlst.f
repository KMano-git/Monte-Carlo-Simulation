!***********************************************************************
      subroutine metlst
!***********************************************************************
!           No use Gmesh-code in planet
!-----------------------------------------------------------------------
      use cplmet, only : gare, gdsv, gwtm, gwtp, hare, hdsp, hdsv, hvol
     >    , hvsb, icel, icmax, icmin, itmax, itmps, itpvs, itsle, jcel
     >    , jcxp1, jtmax, jtmin, kce, kcn, kcs, kcw
      use csize,  only : ndx
      use csonic, only : lfdbg
      use cunit,  only : n6
      implicit none
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jmax, jmax1, jt, j, i, jst, jen, ist, ien, ih
      integer  it, jmax, jt, j, i, jst, jen, ist, ien, ih
      real*8   zm, zh
      character clin*(ndx),clev*(ndx)
!
      write(n6,'(/2x,"***  metlst  ***  lfdbg(1) =",i2)')  lfdbg(1)
!
      if( lfdbg(1).eq.2 ) goto 1000
!
!::list
      do 100 it = 2, itmax-1
!
      if( it.eq.itsle .or. it.eq.itpvs .or. it.eq.itmps ) then
      else
        goto 100
      endif
!
      jmax  = jtmax(it)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   jmax1 = jmax-1
!
!::mxmtrc
      write(n6,'(/2x,"mxmtrc   it =",i3,"  i =",i3)') it,icel(1,it)
      write(n6,602)
     >  "jt","jc","hdsp_E","hdsp_W","hdsp_N","hdsp_S","hdsv_E","hdsv_W"
     >  ,"hdsv_N","hdsp_S","hvol","hvsb"
 602  format(3x,a,3x,a,3x,a,6x,a,6x,a,6x,a,6x,a,6x,a,6x,a,6x,a,6x,a,
     >    6x,a,6x,a,8x,a,8x,a)
      do 210 jt = 1, jmax
      j = jcel(jt,it)
      i = icel(jt,it)
      write(n6,'(2x,2i5,1p10e12.3)')
     >  jt,j,hdsp(j,i,1),hdsp(j,i,2),hdsp(j,i,3),hdsp(j,i,4)
     >   ,hdsv(j,i,1),hdsv(j,i,2),hdsv(j,i,3),hdsv(j,i,4)
     >   ,hvol(j,i),hvsb(j,i)
!x   >   ,hwtp(j,i),hwtm(j,i)
 210  continue
!
!::mymtrc
      write(n6,'(/2x,"mymtrc   it =",i3,"  i =",i3)') it,icel(1,it)
      write(n6,602)
     >   "jt","jc","gdsv_E","gdsv_W","gdsv_N","gdsv_S","hvol  ","hare  "
     >  ,"hvsb  ","gare  ","gwtp","gwtm"
      do 220 jt = 1, jmax
      j = jcel(jt,it)
      i = icel(jt,it)
      write(n6,'(2x,2i5,1p10e12.3)')
     >  jt,j,gdsv(j,i,1),gdsv(j,i,2),gdsv(j,i,3),gdsv(j,i,4)
     >   ,hvol(j,i),hare(j,i),hvsb(j,i),gare(j,i)
     >   ,gwtp(j,i),gwtm(j,i)
 220  continue
!
 100  continue
!
      call wexit("metlst","debug list")
!
 1000  continue
      jst = jtmin(1)
      jen = jtmax(1)
      ist = icmin(jcxp1+1)
      ien = icmax(jcxp1+1)
!
!::clin
      clin = " "
      do 305 j = 1, jen
      clin(j:j) = "-"
      if( mod(j,10).eq.0 ) clin(j:j) = "I"
 305  continue
      clin(j:j) = "I"
!
!::hdsp
      write(n6,'(2x,"== hdsp ==  hdsp(N,S)/hdsp(E,W)")')
      write(n6,'(7x,a)') clin(1:jen)
      do 310 i = ist, ien
      do 320 j = jst, jen
      zm = dmax1(dabs(hdsp(j,i,kce)),dabs(hdsp(j,i,kcw)))
      if( zm.eq.0.0 ) then
      clev(j:j) = "0"
      else
      zh = dmax1(dabs(hdsp(j,i,kcn)),dabs(hdsp(j,i,kcs)))/zm
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
      endif
 320  continue
      write(n6,'(2x,i3,2x,a)') i,clev(1:jen)
 310  continue
!
!::hdsv
      write(n6,'(2x,"== hdsv ==  hdsv(E,W)/hdsv(N,S)")')
      write(n6,'(7x,a)') clin(1:jen)
      do 330 i = ist, ien
      do 340 j = jst, jen
      zm = dmax1(dabs(hdsv(j,i,kcn)),dabs(hdsv(j,i,kcs)))
      if( zm.eq.0.0 ) then
      clev(j:j) = "0"
      else
      zh = dmax1(dabs(hdsv(j,i,kce)),dabs(hdsv(j,i,kcw)))/zm
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
      endif
 340  continue
      write(n6,'(2x,i3,2x,a)') i,clev(1:jen)
 330  continue
!
!::gdsv
      write(n6,'(2x,"== gdsv ==  gdsv(E,W)/gdsv(N,S)")')
      write(n6,'(7x,a)') clin(1:jen)
      do 350 i = ist, ien
      do 360 j = jst, jen
      zm = dmax1(dabs(gdsv(j,i,kcn)),dabs(gdsv(j,i,kcs)))
      if (zm.eq.0.0 ) then
        clev(j:j) = "0"
      else
      zh = dmax1(dabs(gdsv(j,i,kce)),dabs(gdsv(j,i,kcw)))/zm
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
      endif
 360  continue
      write(n6,'(2x,i3,2x,a)') i,clev(1:jen)
 350  continue
!
      call wexit("metlst","debug list")
!
      end
