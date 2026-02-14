!***********************************************************************
      subroutine plsprf
!***********************************************************************
      use cphcns, only : cev
      use cplcom, only : ama, vna, vne, vni, vte, vti, vva
      use cplmet, only : icel, itsle, jcel, jtmax, jtmin
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nfw, it, ia, np, jt, jtst, jten, j, i
      integer  nfw, ia, np, jt, jtst, jten, j, i
      real*8   pt, pt2
!
      write(n6,'(/2x,"*** plsprf ***  file = cmpprf")')
      nfw = 21
      open(unit=nfw,file="cmpprf")
!
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   it = itsle
      ia = 1
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   jtst = jtmin(it)
!ik   jten = jtmax(it)
      jtst = jtmin(itsle)
      jten = jtmax(itsle)
!
      write(nfw,'(2x,"*** plasma parameter ***  it =",i3,"  ia =",i2)')
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik  >  it,ia
     >  itsle,ia
      write(nfw,'(5x,"np",3x,"j",5x,"Na",13x,"Va",13x,
     >   "Te",13x,"Ti",13x,"Ptot",11x,"Pt2")')
      np=0
      do jt = jtst, jten
      if( jt-jtst.le.30 .or. jten-jt.le.30 ) then
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   j = jcel(jt,it)
!ik   i = icel(jt,it)
      j = jcel(jt,itsle)
      i = icel(jt,itsle)
      pt2 = (vne(j,i)*vte(j,i)+vni(j,i)*vti(j,i))*cev
      pt  = pt2+vna(j,i,ia)*ama(ia)*vva(j,i,ia)**2
      np=np+1
      write(nfw,'(2x,2i5,1p10e15.5)') np,j,
     >  vna(j,i,ia),vva(j,i,ia),vte(j,i),vti(j,i),pt,pt2
      endif
      enddo
!
      close(nfw)
!
      return
      end
