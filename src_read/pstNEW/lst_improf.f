!**********************************************************************
      subroutine lst_improf(nt)
!***********************************************************************
      use csize
      use cimcom
      use cimden
      use cntcom
      use cntpls
      use cplcom
      use cplmet
      use cplwrd
      use csonic
      use cunit
      implicit none

!
!::argument
      integer  nt
!
!::local variables
      real*8  zero
      integer izmx, it, i6, jtst, jten, np, jt, j, i, ic, ix, iy, iz
      integer imox, imoy
      real*8  zne, znz, zav, zef, zwrds, zwrdd
      character  dsn*20
      real*8  alnx(ndx), zlmx 
!
!::init value
      zero = 1.0d-30
      izmx = min0( 10, ismax )
      izmx = ismax
      it = nt
!
!::check
      if( it.eq.itsle ) then
      write(n6,'(/2x,"--- lst_improf ---")')
      write(n6,'(2x,"wime =",1p10e12.3)') (wime(ic),ic=1,ncmax,100)
      write(n6,'(2x,"twrd =",1p10e12.3)') (twrd(ic),ic=1,ncmax,100)
      endif
!
!::output file
      i6 = 21
      write(dsn,'(a,i2.2,a)') "improf_", it, ".txt"
      open(unit=i6,file=trim(dsn))
!
      write(i6,'(2x,"*** lst_improf ***  it =",i3,"  itim =",i5,
     >  "  time =",1pe12.3)') it, itim, time
!
      write(i6,'(3x,"np",2x,"ix",2x,"iy",2x,"pl",7x,"pl2",6x,
     >  "dene",6x,"teme",6x,"temi",6x,"wrdD",6x,"wrdS",6x,
     >  "Ne",8x,"Nz",8x,"Nz/Ni",5x,"<Z>",7x,"Zef",7x,30(3x,i2,5x))')
     >   (iz,iz=0,izmx)
!
      jtst = jtmin(it)
      jten = jtmax(it)
      if( it.ge.itmps ) then
      jtst = jtst + 1
      jten = jten - 1
      endif
      call plenc(it,alnx)
      zlmx = alnx(jten)
!
      np = 0
      do jt = jtst, jten
      j  = jcel(jt,it)
      i  = icel(jt,it)
      ic = mcel(j,i)
      if( ic.le.0 ) cycle
!
      np = np + 1
      ix = imox(ic)
      iy = imoy(ic)
!
      zne = deni(ic,1)
      znz = 0.0d0
      zav = 0.0d0
      zef = deni(ic,1)
      do iz = 0, ismax
      zne = zne + dfloat(iz)*tdnz(iz,ic)
      znz = znz + tdnz(iz,ic)
      zav = zav + dfloat(iz)*tdnz(iz,ic)
      zef = zef + dfloat(iz)**2*tdnz(iz,ic)
      enddo
      if( znz.ne.0.0d0 ) zav = zav/znz
      if( zne.ne.0.0d0 ) zef = zef/zne
!
      zwrdd = twrd(ic)  ! imdisk
      zwrds = wime(ic)  ! soldor
!
      write(i6,'(1x,3i4,0p2f9.5,1p30e10.3)')
     >  np, ix, iy, alnx(jt), zlmx-alnx(jt),
     >  dmax1(zero,dene(ic)),dmax1(zero,teme(ic)),dmax1(zero,temi(ic)),
     >  dmax1(zero,zwrds),dmax1(zero,zwrdd),
     >  dmax1(zero,zne), dmax1(zero,znz), dmax1(zero,znz/deni(ic,1)), 
     >  dmax1(zero,zav), dmax1(zero,zef),
     >  (dmax1(zero,tdnz(iz,ic)),iz=0,izmx)
      enddo
      close(i6)
!
      return
      end

