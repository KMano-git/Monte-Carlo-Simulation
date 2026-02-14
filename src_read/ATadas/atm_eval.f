!***********************************************************************
      subroutine atm_eval2(ixs,vx1,vx2,vy,ns,ndm,kis)
!***********************************************************************
!
!   argument
!     vx1(Ne) [m3],  vx2(Te) [eV],  vy(<sig*v>) [m3/s]
!
!   assumption in catcom
!
!      rank = 2
!      spacing   vx1(Ne) & vx2(Te) & vy(sig*v)   log10
!  ==> spacing   irregular interval
!      cgs unit
!
!-----------------------------------------------------------------------
      use catcom, only : pin, pxa, pxd, pxi, pxs, xdaty, xdtx1, xdtx2
      implicit none

!::argument
      integer, intent(in)  :: ixs, kis, ndm
      integer, intent(out) :: ns
      real(8), intent(in)  :: vx1, vx2
      real(8),dimension(ndm), intent(out) ::  vy
! kis : impurity number
!
!::local variables
      real*8  zex1, zix1, zdx1; integer ix1a, ix1b
      real*8  zex2, zix2, zdx2; integer ix2a, ix2b
      real*8  zw, zy
      integer iyaa, iyab, iyba, iybb, i, iy0, ist
      integer, parameter :: i6 = 92
!
!-----------------------------------------------------------------------
!
!::number of charge state
      ns = pin(4,ixs,kis)
!-----
!
!::Ne [MKS ==> cgs unit]
      zex1 = vx1
      zex1 = zex1 * 1.0d-6
      zex1 = dmax1( zex1, pxi(1,ixs,kis) )
      zex1 = dmin1( zex1, pxa(1,ixs,kis) )
      zex1 = dlog10(zex1)
!
!::Te [eV]
      zex2 = vx2
      zex2 = dmax1( zex2, pxi(2,ixs,kis) )
      zex2 = dmin1( zex2, pxa(2,ixs,kis) )
      zex2 = dlog10(zex2)
!
!::index-Ne
      zix1 = (zex1-pxs(1,ixs,kis))/pxd(1,ixs,kis) + 1.0d0
      ix1a = int(zix1)
      ix1a = min0( ix1a, pin(2,ixs,kis)-1 )
      ix1b = ix1a + 1
      if( zex1.ge.xdtx1(ix1a,ixs,kis) .and.
     >    zex1.lt.xdtx1(ix1b,ixs,kis) ) goto 120
      ist  = ix1a + 1
      if( zex1.lt.xdtx1(ix1a,ixs,kis) ) ist = max0(ix1a-1,1)
      do i = ist, pin(2,ixs,kis)
        ix1b = i
        if( zex1.le.xdtx1(i,ixs,kis) ) goto 120
      enddo
      goto 900
 120  continue
      ix1a = ix1b - 1
      zdx1 = (zex1-xdtx1(ix1a,ixs,kis))
     >       /(xdtx1(ix1b,ixs,kis)-xdtx1(ix1a,ixs,kis))
!
!::index-Te
      zix2 = (zex2-pxs(2,ixs,kis))/pxd(2,ixs,kis) + 1.0d0
      ix2a = int(zix2)
      ix2a = min0( ix2a, pin(3,ixs,kis)-1 )
      ix2b = ix2a + 1
      if( zex2.ge.xdtx2(ix2a,ixs,kis) .and.
     >    zex2.lt.xdtx2(ix2b,ixs,kis) ) goto 140
      ist = ix2a + 1
      if( zex2.lt.xdtx2(ix2a,ixs,kis) ) ist = max0(ix2a-1,1)
      do i = ist, pin(3,ixs,kis)
        ix2b = i
        if( zex2.le.xdtx2(i,ixs,kis) ) goto 140
      enddo
      goto 900
 140  continue
      ix2a = ix2b - 1
      zdx2 = (zex2-xdtx2(ix2a,ixs,kis))
     >       /(xdtx2(ix2b,ixs,kis)-xdtx2(ix2a,ixs,kis))
!
      if( zdx1.lt.0.0d0 .or. zdx1.gt.1.0d0 .or.
     >    zdx2.lt.0.0d0 .or. zdx2.gt.1.0d0 ) then
        write(i6,'(2x,"serious error  zdx1, zdx2 ",1p2e12.3)')
     >    zdx1, zdx2
        goto 900
      endif
!
      do i = 1, ns
        iy0  = pin(2,ixs,kis)*pin(3,ixs,kis)*(i-1) + pin(1,ixs,kis)
        iyaa = ix1a + pin(2,ixs,kis)*(ix2a-1) + iy0
        iyba = ix1b + pin(2,ixs,kis)*(ix2a-1) + iy0
        iybb = ix1b + pin(2,ixs,kis)*(ix2b-1) + iy0
        iyab = ix1a + pin(2,ixs,kis)*(ix2b-1) + iy0
!
        zw  = xdaty(iyaa,kis)*(1.0d0-zdx1)*(1.0d0-zdx2) +
     >      xdaty(iyba,kis)*zdx1*(1.0d0-zdx2) +
     >      xdaty(iybb,kis)*zdx1*zdx2 +
     >      xdaty(iyab,kis)*(1.0d0-zdx1)*zdx2
        zy  = 10.0d0**zw
!
!::sigv [cgs ==> MKS unit]
        vy(i) = zy * 1.0d-6
      enddo
!
      return
!
!::error
 900  continue
      write(i6,'(/2x,"*** serious error at sub. atm_eval kis = ",i4)')
     >     kis
      write(i6,'(2x,"vx1 =",1pe12.3,"  min,max =",1p2e12.3)')
     >   vx1,10.0d0**xdtx1(1,ixs,kis)
     >  ,10.0d0**xdtx1(pin(2,ixs,kis),ixs,kis)
      write(i6,'(2x,"vx2 =",1pe12.3,"  min,max =",1p2e12.3)')
     >   vx2,10.0d0**xdtx2(1,ixs,kis)
     >  ,10.0d0**xdtx2(pin(3,ixs,kis),ixs,kis)
      call wexit("atm_eval","no found index")
!
      end
