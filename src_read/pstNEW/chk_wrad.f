!***********************************************************************
      subroutine chk_wrad(wrdm)
!***********************************************************************
!
!      wrdm : Monte cal. for Ar impurity
!
!-----------------------------------------------------------------------
      use csize
      use cntcom
      use cntpls
      use cplmet
      use cimctl
      use csonic
      use cunit
      implicit none
!
!::argument
      real*8  wrdm(ndmc)
!
!::local variables
      integer  itb(5)
!
!::local radiation
      integer  i6
      integer  ic, ir, it, jt, jtst, jten, j, i, np, ix, iy
      integer  imox, imoy
      real*8   zero, totwrd, trgwrd(10)
      character  dsn*20
!
      write(n6,'(/2x,"*** chk_wrad ***")')
!

      i6 = 21
      zero = 1.0d-25
!
      totwrd = 0.0d0
      trgwrd(1:10) = 0.0d0
!
      do ic = 1, ncmax
      ir = mrgn(ic)
      if( ir.gt.0 ) then
      if( ir.le.6 ) totwrd = totwrd + wrdm(ic)*volm(ic)
      trgwrd(ir) = trgwrd(ir) + wrdm(ic)*volm(ic)
      endif
      enddo
!
      write(n6,'(2x,"cradMC",1pe14.6,2x,1p10e12.3)') 
     > totwrd, (trgwrd(ir),ir=1,7)
!
      do it = itsls, itsle
      if( it.eq.itsls ) cycle
      write(dsn,'(a,i2.2)') "Wrdm_", it      
      open(unit=i6,file=dsn)
      jtst = jtmin(it)
      jten = jtmax(it)
      if( it.ge.itmps ) then
      jtst = jtst + 1
      jten = jten - 1
      endif
!
      np = 0
      do jt = jtst, jten
      j = jcel(jt,it)
      i = icel(jt,it)
      ic = mcel(j,i)
      if( ic.le.0 ) cycle
      np = np + 1
      ix = imox(ic)
      iy = imoy(ic)
!
      write(i6,'(1x,i4,i5,i4,i7,1p4e12.3)')
     >  np, ix, iy, ic, dmax1(zero,dene(ic)), dmax1(zero,teme(ic)),
     >   dmax1(zero,temi(ic)), dmax1(zero,dabs(wrdm(ic)))
      enddo
      close(i6)
      enddo
!
      return
      end
