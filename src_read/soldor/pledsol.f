!**********************************************************************
      subroutine pledsol
!**********************************************************************
      use cplcom, only : nion, vna, vte, vti
      use cplmet, only : icel, itsle, itsls, jcel, jtmax
      use csize,  only : ndy
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer :: ia, it, jw, j, i, ihh
      integer, dimension(ndy) :: iswt
!
!::input data
      real*8 :: hnisw = 4.0d19
      real*8 :: htisw = 10.0d0
!
      write(n6,'(/2x,"*** pledsol ***")')
      write(n6,'(2x,"hnisw =",1pe12.3,"  htisw =",1pe12.3)')
     >    hnisw, htisw
!
!::switch
      ia  = 1
      do it = itsls, itsle
      ihh = 0
      jw  = 2
      j = jcel(jw,it)
      i = icel(jw,it)
      if( vna(j,i,ia).lt.hnisw ) ihh = ihh + 1
      jw = jtmax(it) - 1
      j = jcel(jw,it)
      i = icel(jw,it)
      if( vna(j,i,ia).lt.hnisw ) ihh = ihh + 1
      iswt(it) = 0
      if( ihh.gt.0 ) iswt(it) = 1
      enddo
!
!::debug write
      ia = 1
      jw = 2
      write(n6,'(2x,"Nad =",1p10e11.3)')
     >     (vna(jcel(jw,it),icel(jw,it),ia),it=itsls,itsle)
      write(n6,'(2x,"Tid =",1p10e11.3)')
     >     (vti(jcel(jw,it),icel(jw,it)),it=itsls,itsle)
      jw = jtmax(it) - 1
      write(n6,'(2x,"Nad =",1p10e11.3)')
     >     (vna(jcel(jw,it),icel(jw,it),ia),it=itsls,itsle)
      write(n6,'(2x,"Tid =",1p10e11.3)')
     >     (vti(jcel(jw,it),icel(jw,it)),it=itsls,itsle)
      write(n6,'(2x,"swt =",10(3x,i2,6x))')
     >          (iswt(it),it=itsls,itsle)
!
!::new value
      do it = itsls, itsle
      if( iswt(it).eq.0 ) cycle
      do jw = 1, jtmax(it)
      j = jcel(jw,it)
      i = icel(jw,it)
      do ia = 1, nion
      vna(j,i,ia) = hnisw
      enddo
      vti(j,i) = htisw
      vte(j,i) = htisw
      enddo
      enddo
!
!::debug write
      ia = 1
      jw = 2
      write(n6,'(2x,"Nad =",1p10e11.3)')
     >     (vna(jcel(jw,it),icel(jw,it),ia),it=itsls,itsle)
      write(n6,'(2x,"Tid =",1p10e11.3)')
     >     (vti(jcel(jw,it),icel(jw,it)),it=itsls,itsle)
      jw = jtmax(it) - 1
      write(n6,'(2x,"Nad =",1p10e11.3)')
     >     (vna(jcel(jw,it),icel(jw,it),ia),it=itsls,itsle)
      write(n6,'(2x,"Tid =",1p10e11.3)')
     >     (vti(jcel(jw,it),icel(jw,it)),it=itsls,itsle)
      write(n6,'(2x,"swt =",10(3x,i2,6x))')
     >          (iswt(it),it=itsls,itsle)
!
      return
      end
