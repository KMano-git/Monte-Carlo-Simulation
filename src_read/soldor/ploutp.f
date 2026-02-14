!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   subroutine ploutp(kk,kt)
      subroutine ploutp(kk,it)
!***********************************************************************
      use cplcom, only : flmet, flmxe, q1a, q2a, q3, q4, vcs, vna, vni
     >    , vte, vti, vva, vve
      use cplmet, only : icel, itpve, itpvs, itsle, jcel, jcmax, jtmax
     >    , jtmin
      use csize,  only : ndx
      use csonic, only : time
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  kk, kt
      integer, intent(in) :: kk, it
!
!::local variables
      real*8  alnx(ndx)
      integer itb(10)
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer ia, it, jtst, jten, jt, j, i, jmax, l, nimx
!ik   real*8  zmch, zpcn
      integer ia, jtst, jten, jt, j, i, jmax, l, nimx
      real*8  zmch
!
!-----------------------------------------------------------------------
!::conservative variables
!-----------------------------------------------------------------------
      if(kk.eq.1) then
      write(n6,'(/2x,"***  ploutp (q1a,q2a,q3,q4) ***   time =",
     >  1pe12.4)') time
      write(n6,'( 2x,"it",3x,"ic",3x,"jc",5x,"lp",9x,"q1a",8x,"q2a",8x
     >    ,"q3",9x,"q4")')
      ia = 1
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   it = kt
      jtst = jtmin(it)
      jten = jtmax(it)
      call plenc(it,alnx)
      do 120 jt = jtst, jten
      j = jcel(jt,it)
      i = icel(jt,it)
      write(n6,'(1x,3i4,2x,1p5e11.3)')
     >    it,i,j,alnx(jt),q1a(j,i,ia),q2a(j,i,ia),q3(j,i),q4(j,i)
 120  continue
      endif
!
!-----------------------------------------------------------------------
!::non-conservative variables
!-----------------------------------------------------------------------
      if(kk.eq.2) then
      write(n6,'(/2x,"***  ploutp (ni,vi,ti,te) ***   time =",
     >   1pe12.4)') time
      ia = 1
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   it = kt
      jtst = jtmin(it)
      jten = jtmax(it)
      call plenc(it,alnx)
      write(n6,'(2x,"it =",i3,"  ic =",i3)') it,icel(jtst,it)
      write(n6,'(2x,"jt",3x,"jc",5x,"lp",10x,"ni",10x
     >    ,"vp",10x,"ti",7x,"te",7x,"mch")')
      do 220 jt = jtst, jten
      j = jcel(jt,it)
      i = icel(jt,it)
      zmch = vva(j,i,ia)/vcs(j,i)
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   zpcn = vne(j,i)*vte(j,i)+vni(j,i)*vti(j,i)
!ik  >         +ama(ia)*vni(j,i)*vva(j,i,ia)**2/cev
      write(n6,'(1x,2i5,2x,0pf10.3,1p5e12.3)')
     >    jt,j,alnx(jt),vna(j,i,ia),vva(j,i,ia),vti(j,i),vte(j,i)
     >      ,zmch
 220  continue
 210  continue
      endif
!
!-----------------------------------------------------------------------
!::flux limit of eta & Xe
!-----------------------------------------------------------------------
      if( kk.eq.3 ) then
      write(n6,'(/2x,"*** ploutp (flux limit of eta) ***  time =",
     >  1pe12.3)') time
      itb(1) = 2
      itb(2) = 10
      itb(3) = itsle
      itb(4) = itpvs
      itb(5) = itpve-2
      itb(6) = itpve-1
      jmax = jcmax
      write(n6,'(2x,3x,2x,6(3x,i3,4x))') (itb(l),l=1,6)
      do 310 j = 1, jmax
      write(n6,'(2x,i3,2x,1p6e10.2)') j,(flmet(j,itb(l)),l=1,6)
 310  continue
!
      write(n6,'(/2x,"*** ploutp (flux limit of Xe) ***  time =",
     >  1pe12.3)') time
      itb(1) = 2
      itb(2) = 10
      itb(3) = itsle
      itb(4) = itpvs
      itb(5) = itpve-2
      itb(6) = itpve-1
      jmax = jcmax
      write(n6,'(2x,3x,2x,6(3x,i3,4x))') (itb(l),l=1,6)
      do 320 j = 1, jmax
      write(n6,'(2x,i3,2x,1p6e10.2)') j,(flmxe(j,itb(l)),l=1,6)
 320  continue
      endif
!
!-----------------------------------------------------------------------
!::multi ion
!-----------------------------------------------------------------------
      if(kk.eq.4) then
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   it = kt
      jtst = jtmin(it)
      jten = jtmax(it)
      write(n6,'(/2x,"***  ploutp (ni,vi,na,va,ti,te) ***   time =",
     >   1pe12.4,"  it =",i3,"  i =",i3)') time, it, icel(jtst,it)
      write(n6,'(4x,"jt",3x,"jc",5x,"ni",10x,"n1",10x,"n2",
     >   10x,"ve",10x,"v1",10x,"v2",10x,"Ti",10x,"Te")')
      nimx = 2
      ia = 1
      do jt = jtst, jten
      j = jcel(jt,it)
      i = icel(jt,it)
      zmch = vva(j,i,ia)/vcs(j,i)
      write(n6,'(1x,2i5,2x,1p16e12.3)')
     >    jt,j,vni(j,i),(vna(j,i,ia),ia=1,nimx),
     >    vve(j,i),(vva(j,i,ia),ia=1,nimx),
     >    vti(j,i),vte(j,i)
      enddo
      endif
!
      return
      end
!
!***********************************************************************
      subroutine ploutp2(nft)
!***********************************************************************
      use cplcom, only : vna, vte, vti, vva
      use cplmet, only : icel, itmpe, itmps, itsle, jcel, jtmax, jtmin
      use csonic, only : time
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nft
      integer, intent(in) :: nft
!
!::local variables
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  i6, jmax, jmns, jmne, ia, it, jtst, jten, jt
!ik   integer  jc1, ic1, jc2, ic2, jc3, ic3, jc
      integer  jmns, jmne, ia, it, jtst, jten, jt
      integer  jc1, ic1, ic2, ic3, jc
      real*8   una, uva
!
      if( nft.le.0 ) return
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   i6 = nft
!
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   jmax = jtmax(itsle)
      jmns = jcel(2,itmps)
      jmne = jcel(jtmax(itmps)-1,itmps)
      ia = 1
      una = 1.0d19
      uva = 1.0d4
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   write(i6,'(/2x,"***  ploutp (ni,vi,ti,te) ***   time =",
      write(nft,'(/2x,"***  ploutp (ni,vi,ti,te) ***   time =",
     >   1pe12.4)') time
!
      it = itsle
      jtst = jtmin(it)
      jten = jtmax(it)
!
      do jt = jtst, jten
      jc1 = jcel(jt,itsle)
      ic1 = icel(jt,itsle)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   jc2 = jcel(jt,itmps)
      ic2 = icel(jt,itmps)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   jc3 = jcel(jt,itmpe)
      ic3 = icel(jt,itmpe)
      jc  = jc1
!
      if( jc.lt.jmns .or. jc.gt.jmne ) then
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   write(i6,'(1x,i3,3(2x,2i4,4f8.2))') jt
      write(nft,'(1x,i3,3(2x,2i4,4f8.2))') jt
     > ,jc,ic1,vna(jc,ic1,ia)/una,vva(jc,ic1,ia)/uva
     > ,        vti(jc,ic1),vte(jc1,ic1)
     > ,jc,ic2, 0.0d0, 0.0d0, 0.0d0, 0.0d0
     > ,jc,ic3, 0.0d0, 0.0d0, 0.0d0, 0.0d0
      else
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   write(i6,'(1x,i3,3(2x,2i4,4f8.2))') jt
      write(nft,'(1x,i3,3(2x,2i4,4f8.2))') jt
     > ,jc,ic1,vna(jc,ic1,ia)/una,vva(jc,ic1,ia)/uva
     > ,       vti(jc,ic1),vte(jc,ic1)
     > ,jc,ic2,vna(jc,ic2,ia)/una,vva(jc,ic2,ia)/uva
     > ,       vti(jc,ic2),vte(jc,ic2)
     > ,jc,ic3,vna(jc,ic3,ia)/una,vva(jc,ic3,ia)/uva
     > ,       vti(jc,ic3),vte(jc,ic3)
      endif
!
      enddo
!
      close(nft)
!
      return
      end
