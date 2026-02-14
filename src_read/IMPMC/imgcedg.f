!***********************************************************************
      subroutine imgcedg(cspt)
!***********************************************************************
      use cimcom, only : arbz, csbz, flbz, flximpamp, gfbz, icbz, nbz
     >    , ndbz, npemt, nrf, prbz, rhoimpsrc, sflux, snbz, tfbz, xpbz
     >    , ypbz
      use cphcns, only : cpi
      use cunit,  only : n6
      implicit none
!
!::argument
      character, intent(in) :: cspt*(*) ! dummy
!
!::local variables
      integer np, i, ii
      real*8  dl, x0, zsum, zmax, zfl, zfn
!
      write(n6,'(2x,"*** imgcedg ***")')
!
!::input data
      np  = nrf+1
!
      write(n6,'(2x,"np =",i4,"  rhoImpSrc =",f12.4," flxImpAmp =",e16.8,
     >   "  npemt =",i6)') np, rhoImpSrc, flxImpAmp, npemt
!
      call posman(rhoImpSrc, np, xpbz, ypbz, icbz)

      do i = 1, np-1
      dl = sqrt( (xpbz(i+1)-xpbz(i))**2+(ypbz(i+1)-ypbz(i))**2 )
      x0 = 0.5d0*(xpbz(i)+xpbz(i+1))
      arbz(i) = 2.0d0*cpi*x0*dl
      csbz(i) = (xpbz(i+1)-xpbz(i))/dl
      snbz(i) = (ypbz(i+1)-ypbz(i))/dl
      enddo
      arbz(np) = 0.0d0
      csbz(np) = 0.0d0
      snbz(np) = 0.0d0
!
!::birth profile
      zsum = 0.0d0
      zmax = -1.0d20
      ii = 0
      do i = 1, np-1
      zfl = 1.0d0
      zfn = zfl*arbz(i)
      zsum = zsum + zfn
      zmax = dmax1( zmax, zfn )
      ii = ii + 1
      if( ii.gt.ndbz ) goto 910
      prbz(ii) = zsum
      flbz(ii) = zfn
      enddo
!
!::normalization
      nbz = ii
      tfbz = zsum ! total (flux * area)
      gfbz = zmax
      do i = 1, nbz
      prbz(i) = prbz(i)/zsum
      enddo
      prbz(nbz) = 1.0001d0
      sflux = flxImpAmp
!
      write(n6,'(2x,"nbz =",i3,"  sflux =",1pe12.3)') nbz, sflux
      write(n6,'(3(2x,i3,3f10.5))') 1, xpbz(1), ypbz(1), prbz(1),
     >  nbz,xpbz(nbz),ypbz(nbz),prbz(nbz),nbz+1,xpbz(nbz+1),ypbz(nbz+1)
!
      return

 910  continue
      call wexit("imgcedg","ii.gt.ndbz")
      end
