!**********************************************************************
      subroutine set_puff
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
      use cimpuf, only : ndpfz, pf_are, pf_deg, pf_flx, pf_iiw, pf_imx
     >    , pf_ity, pf_mag, pf_nty, pf_nwl, pf_px1, pf_px2, pf_py1
     >    , pf_py2, pf_rate
      use cntcom, only : ipps, ndpf, npep, npsp, xpnt, ypnt
      use cphcns, only : cpi
      use cpmpls, only : xmd1, xmd2, ymd1, ymd2
      use csize,  only : ndwp
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  m, nsiz, ii, ist, imx, nw, iw, iws, iwe, i
      integer  m, ii, ist, imx, nw, iw, iws, iwe, i
      integer  ip1, ip2
      real*8   zsum, x0, y0, dlx, dly, dln, zar, zwt, zflx
      real*8   zpmg(ndpf)
      real*8   xpc, ypc, dpc, dth
!
      write(n6,'(/2x,"*** set_puff ***  npf/ndpfz =",2i3)')
     >    pf_nty, ndpfz
      if( pf_nty.le.0 ) return
      if( pf_nty.gt.ndpfz ) call wexit("set_puff","pf_nty.gt.ndpfz")
!
      xpc = 0.5d0*(xmd1+xmd2)
      ypc = 0.5d0*(ymd1+ymd2)
      write(n6,'(2x,"xmd1,ymd1 =",2f10.5,"  xmd2,ymd2 =",2f10.5,
     >    "  xpc,ypc =",2f10.5)') xmd1,ymd1,xmd2,ymd2,xpc,ypc
!
!
      do m = 1, pf_nty
      write(n6,'(2x,2i3,1pe12.3," Pam3/s",2x,1pe12.3," 1/s",2x,
     >   "position ",0p4f10.3)')
     >   m, pf_nwl(m), pf_rate(m), pf_mag(m),
     >   pf_px1(m), pf_px2(m), pf_py1(m), pf_py2(m)
      enddo
!
!::flux weight
! modified 2/1 lines organize local variables and include files by kamata 2021/06/28
!ik   nsiz = ndwp
!ik   call setd( pf_flx, nsiz, 0.0d0 )
      call setd( pf_flx, ndwp, 0.0d0 )
!
      ii = 0
      do m = 1, pf_nty
      zsum = 0.0d0
      ist  = ii + 1
!
      nw  = pf_nwl(m)
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
      if( x0.ge.pf_px1(m) .and. x0.le.pf_px2(m) .and.
     >    y0.ge.pf_py1(m) .and. y0.le.pf_py2(m) ) zwt = 1.0d0
      dpc = sqrt((x0-xpc)**2+(y0-ypc)**2)
      dth = datan2( (y0-ypc)/dpc, (x0-xpc)/dpc )
      if( dth.lt.0.0d0 ) dth = dth + 2.0d0*cpi
      dth = dth/cpi*180.0d0
!-----
!x    write(n6,'(2x,"nw =",i2,"  iw =",i5,"  x0,y0 =",2f9.4,
!x   >  "  zwt =",f9.4,"  ar =",1pe12.3,"  th =",1pe12.3)')
!x   >   nw, iw, x0, y0, zwt, zar, dth
!-----
      if( zwt.le.0.0d0 ) cycle
      ii = ii + 1
      if( ii.gt.ndwp ) goto 900
      pf_ity(ii) = m
      pf_iiw(ii) = iw
      pf_flx(ii) = zwt*zar
      pf_are(ii) = zar
      pf_deg(ii) = dth
      zsum = zsum + pf_flx(ii)
      enddo   ! loop(iw)
      imx = ii
      write(n6,'(2x,"nw =",i2,"  imx =",i3)') nw, imx
!
!::normalization
      do i = ist, imx
      pf_flx(i) = pf_flx(i)/zsum
      enddo
!
      enddo  ! loop(m)
      pf_imx = imx
!
!::debug write
      write(n6,'(4x,"i",4x,"pf",3x,"nw",3x,"iw",3x,"x0",7x,"y0",7x,
     >   "flx",9x,"pmg",9x,"are",9x,"deg")')
! modified 2/1 lines organize local variables and include files by kamata 2021/06/28
!ik   nsiz = ndpf
!ik   call setd( zpmg, nsiz, 0.0d0 )
      call setd( zpmg, ndpf, 0.0d0 )
      do i = 1, pf_imx
      m  = pf_ity(i)
      iw = pf_iiw(i)
      zflx = pf_flx(i)*pf_mag(m)
      zpmg(m) = zpmg(m) + zflx
      ip1 = ipps(iw)
      ip2 = ipps(iw+1)
      x0  = 0.5d0*(xpnt(ip1)+xpnt(ip2))
      y0  = 0.5d0*(ypnt(ip1)+ypnt(ip2))
      write(n6,'(2x,4i5,2f9.4,1p4e12.3)')
     >   i, m, nw, iw, x0, y0, zflx, zpmg(m), pf_are(i), pf_deg(i)
      enddo
!
      write(n6,'(2x,"total puff rate ")')
      do m = 1, pf_nty
      write(n6,'(2x,i3,1p2e12.3)') m, zpmg(m), pf_mag(m)
      enddo
!
      return
!
!::dimension error
 900  continue
      call wexit("ntpuff","ii.gt.ndwp")
      end
