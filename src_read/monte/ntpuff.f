!**********************************************************************
      subroutine ntpuff
!**********************************************************************
!
!   set neutral flux density of gas puf
!
!  input data
!     npf         : number of puff type
!     ipfnp(ndpf) : number of plasma surface (1:sol/2:idp/3:prv/4:odp)
!     pfmag(ndpf) : puff rate
!     pfpx1(ndpf) : region xmin
!     pfpx2(ndpf) : region xmax
!     pfpy1(ndpf) : region ymin
!     pfpy2(ndpf) : region ymax
!
!  output data
!     ipfmx       : number of puff segment
!     ipfmp(ndwp) : number of puff type
!     ipfiw(ndwp) : number of plasma surface segment
!     pfflx(ndwp) : flux density
!
!  How to use common variables
!
!      do i = 1, ipfmx
!        m  = ipfmp(i)
!        iw = ipfiw(i)
!        flux = pfflx(iw)*pfmag(m)  ! local flux density
!      enddo
!
!  How to set inout data
!
!  < sol Z > 2.0 >
!      npf = 2,
!      ipfnp(1) = 1,   pfmag(1) = 0.5d22
!      pfpx1(1)=-99.0, pfpx2(1)=99.0, pfpy1(1)=2.0, pfpy2(1)=99.0
!
!  <private region>
!      ipfnp(2) = 3,   pfmag(2) = 0.1d22
!      pfpx1(2)=-99.0, pfpx2(2)=99.0, pfpy1(2)=-99.0, pfpy2(2)=99.0
!
!----------------------------------------------------------------------
      use cntcom, only : ipfiw, ipfmp, ipfmx, ipfnp, ipps, ndpf, npep
     >    , npf, npsp, pfflx, pfmag, pfpx1, pfpx2, pfpy1, pfpy2, xpnt
     >    , ypnt
      use cphcns, only : cpi
      use csize,  only : ndwp
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer  m, nsiz, ii, ist, ien, nw, iw, iws, iwe, i
      integer  m, ii, ist, ien, nw, iw, iws, iwe, i
      integer  ip1, ip2
      real*8   zsum, x0, y0, dlx, dly, dln, zar, zwt, zflx
      real*8   zpmg(ndpf)
!
      write(n6,'(/2x,"*** ntpuff ***  npf =",i3)') npf
      if( npf.le.0 ) return
!
      do m = 1, npf
      write(n6,'(2x,i2,1pe12.3,2x,0p4f10.3)')
     >   m,pfmag(m),pfpx1(m),pfpx2(m),pfpy1(m),pfpy2(m)
      enddo
!
!::flux weight
! modified 2/1 lines organize local variables and include files by kamata 2021/06/16
!ik   nsiz = ndwp
!ik   call setd( pfflx, nsiz, 0.0d0 )
      call setd( pfflx, ndwp, 0.0d0 )
!
      ii = 0
      do m = 1, npf
      zsum = 0.0d0
      ist  = ii + 1
!
      nw  = ipfnp(m)
      iws = npsp(nw)
      iwe = npep(nw)
!
      do iw = iws, iwe-1
      ip1 = ipps(iw)
      ip2 = ipps(iw+1)
      x0  = 0.5d0*(xpnt(ip1)+xpnt(ip2))
      y0  = 0.5d0*(ypnt(ip1)+ypnt(ip2))
      dlx = xpnt(ip2)-xpnt(ip1)
      dly = ypnt(ip2)-ypnt(ip1)
      dln = dsqrt(dlx*dlx+dly*dly)
      zar = 2.0d0*cpi*x0*dln
      zwt = 0.0d0
      if( x0.ge.pfpx1(m) .and. x0.le.pfpx2(m) .and.
     >    y0.ge.pfpy1(m) .and. y0.le.pfpy2(m) ) zwt = 1.0d0
!-----
      write(n6,'(2x,"nw =",i2,"  iw =",i5,"  x0,y0 =",2f9.4,
     >  "  zwt =",f9.4)') nw, iw, x0, y0, zwt
!-----
      if( zwt.le.0.0d0 ) cycle
      ii = ii + 1
      if( ii.gt.ndwp ) goto 900
      ipfmp(ii) = m
      ipfiw(ii) = iw
      pfflx(ii) = zwt*zar
      zsum = zsum + pfflx(ii)
      enddo   ! loop(iw)
      ien = ii
      write(n6,'(2x,"ien =",i3)') ien
!
!::normalization
      do i = ist, ien
      pfflx(i) = pfflx(i)/zsum
      enddo
!
      enddo  ! loop(m)
      ipfmx = ien
!
!::debug write
! modified 2/1 lines organize local variables and include files by kamata 2021/06/16
!ik   nsiz = ndpf
!ik   call setd( zpmg, nsiz, 0.0d0 )
      call setd( zpmg, ndpf, 0.0d0 )
      do i = 1, ipfmx
      m  = ipfmp(i)
      iw = ipfiw(i)
      zflx = pfflx(i)*pfmag(m)
      zpmg(m) = zpmg(m) + zflx
      ip1 = ipps(iw)
      ip2 = ipps(iw+1)
      x0  = 0.5d0*(xpnt(ip1)+xpnt(ip2))
      y0  = 0.5d0*(ypnt(ip1)+ypnt(ip2))
      write(n6,'(2x,i4,i3,i5,2f9.3,1p2e12.3)')
     >   i, m, iw, x0, y0, zflx, zpmg(m)
      enddo
!
      write(n6,'(2x,"total puff rate ")')
      do m = 1, npf
      write(n6,'(2x,i3,1p2e12.3)') m, zpmg(m), pfmag(m)
      enddo
!
      return
!
!::dimension error
 900  continue
      call wexit("ntpuff","ii.gt.ndwp")
      end
