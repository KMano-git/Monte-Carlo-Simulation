!***********************************************************************
      subroutine set_bkgt
!***********************************************************************
!
!     set_bkgt : set distribution of backflow from gate
!
!::input data
!   bk_nty        : number of backflow type = nogt
!   bk_npt        : number of sample particles in MC-cal
!   bk_nwl(ndbkz) : number of plasma surface = 3(prv)
!   bk_mag(ndbkz) : puff rate

!   bk_px1(ndbkz), bk_px2(ndbkz) : dummy  gate poistion
!   bk_py1(ndbkz), bk_py2(ndbkz) : dummy  gate position
!
!::output data
!   bk_ity(ndwp) :  gate number (1/2)
!   bk_iiw(ndwp) :  wall segment (iw)
!   bk_imx       :  number of wall segment (50<ndwp)
!   bk_flx(ndwp) :  flux density
!   bk_wem(ndwp) :  emitted weight at wall segment
!   bk_are(ndwp) :  area (m2)
!   bk_deg(ndwp) :  angle of wall segment (deg)
!
!::Note
!   wall   : npsw, npew, ipwl, xpnt(ipwl(iw)) : imgphys,imgchem,imgbkgt
!   plasma : npsp, npep, ipps, xpnt(ipps(iw)) : imgpuff (Why?)
!
!-----------------------------------------------------------------------
      use cimpuf, only : bk_are, bk_deg, bk_flx, bk_iiw, bk_imx, bk_ity
     >    , bk_mag, bk_nty, bk_nwl
      use cntcom, only : igtw, ipwl, nogt, npew, npsw, xpnt, ypnt
     > , ncmax, icwl
      use cphcns, only : cpi
      use cpmpls, only : xmd1, xmd2, ymd1, ymd2
      use csize,  only : ndgt, ndwp
      use cunit,  only : n6
      use mod_externalgrid
      implicit none
!
!::argument
!
!::local variables
      integer :: m, nwl, ii, igt, ist, nw, iws, ien, iwe, iw, i
      integer :: ip1, ip2, ic, iw_end
      real(8) :: zsum, zar, zwt, zflx
      real(8) :: x0, y0, dlx, dly, dln, xpc, ypc, dpc, dth
      real(8) :: zpmg(ndgt)
      real(8) :: x1, y1, x2, y2
!
!
      write(n6,'(/2x,"*** set_bkgt ***  bk_nty/nogt =",2i3,
     >  "  bk_mag =",1p2e12.4)') bk_nty, nogt, (bk_mag(i),i=1,bk_nty)
      if( bk_nty.le.0 ) return
      if( nogt.le.0 ) then
        call wexit("set_bkgt","nogt <= 0")
      endif
!
!::(xpc,ypc) to define angle
      xpc = 0.5d0*(xmd1+xmd2)
      ypc = 0.5d0*(ymd1+ymd2)
!
!::flux weight
      call setd( bk_flx, ndwp, 0.0d0 )
!
      bk_nty = nogt  ! number of gate
      nwl = 3        ! private wall
!
      ii = 0
      do m = 1, bk_nty
        bk_nwl(m) = nwl
        igt = m        ! gate number
        zsum = 0.0d0
        ist  = ii + 1
!
        nw  = bk_nwl(m)
        iws = npsw(nw)
        iwe = npew(nw)
        if(use_exdata) then
          iw_end = iwe
        else
          iw_end = iwe-1
        endif
        do iw = iws, iw_end
          ip1 = ipwl(iw)
          if(use_exdata) then
            ip2 = ipwl2(iw)
          else
            ip2 = ipwl(iw+1)
          endif
          ic = icwl(iw)
          if(ic.le.ncmax .or. .not.use_exdata) then
            x1 = xpnt(ip1)
            y1 = ypnt(ip1)
            x2 = xpnt(ip2)
            y2 = ypnt(ip2)
          elseif(ic .le. ncmax+vac_ele_size)then
            x1 = vac_grid_x(ip1)
            y1 = vac_grid_y(ip1)
            x2 = vac_grid_x(ip2)
            y2 = vac_grid_y(ip2)
          elseif(ic .le. ncmax+vac_ele_size+pri_ele_size)then
            x1 = pri_grid_x(ip1)
            y1 = pri_grid_y(ip1)
            x2 = pri_grid_x(ip2)
            y2 = pri_grid_y(ip2)
          else
            x1 = vgx_EX(ip1)
            y1 = vgy_EX(ip1)
            x2 = vgx_EX(ip2)
            y2 = vgy_EX(ip2)
          endif
          x0  = 0.5d0*(x1+x2)
          y0  = 0.5d0*(y1+y2)
          dlx = x2-x1
          dly = y2-y1
          dln = dsqrt(dlx*dlx+dly*dly)
          zar = 2.0d0*cpi*x0*dln
!-----
          if( igtw(iw).eq.igt ) then
            zwt = 1.0d0
          else
            zwt = 0.0d0
          endif
!-----
          dpc = sqrt((x0-xpc)**2+(y0-ypc)**2)
          dth = datan2( (y0-ypc)/dpc, (x0-xpc)/dpc )
          if( dth.lt.0.0d0 ) dth = dth + 2.0d0*cpi
          dth = dth/cpi*180.0d0
          if( zwt.le.0.0d0 ) cycle
          ii = ii + 1
          if( ii.gt.ndwp ) goto 900
          bk_ity(ii) = m
          bk_iiw(ii) = iw
          bk_flx(ii) = zwt*zar
          bk_are(ii) = zar
          bk_deg(ii) = dth
          zsum = zsum + bk_flx(ii)
        enddo   ! loop(iw)
        ien = ii
!
!::normalization
        do i = ist, ien
          bk_flx(i) = bk_flx(i)/zsum
        enddo
!
      enddo  ! loop(m)
      bk_imx = ien
!
!::debug write
      write(n6,'(4x,"i",4x,"gt",3x,"nw",3x,"iw",3x,"x0",7x,"y0",7x,
     >   "xp1",6x,"yp1",6x,"xp2",6x,"yp2",6x,
     >   "flx",9x,"pmg",9x,"are",9x,"deg")')
!
      call setd( zpmg, ndgt, 0.0d0 )
!
      do i = 1, bk_imx
        m    = bk_ity(i)
        iw   = bk_iiw(i)
        zflx = bk_flx(i)*bk_mag(m)
        zpmg(m) = zpmg(m) + zflx
        ip1 = ipwl(iw)
        if(use_exdata) then
          ip2 = ipwl2(iw)
        else
          ip2 = ipwl(iw+1)
        endif
        ic = icwl(iw)
        if(ic.le.ncmax .or. .not.use_exdata) then
          x1 = xpnt(ip1)
          y1 = ypnt(ip1)
          x2 = xpnt(ip2)
          y2 = ypnt(ip2)
        elseif(ic .le. ncmax+vac_ele_size)then
          x1 = vac_grid_x(ip1)
          y1 = vac_grid_y(ip1)
          x2 = vac_grid_x(ip2)
          y2 = vac_grid_y(ip2)
        elseif(ic .le. ncmax+vac_ele_size+pri_ele_size)then
          x1 = pri_grid_x(ip1)
          y1 = pri_grid_y(ip1)
          x2 = pri_grid_x(ip2)
          y2 = pri_grid_y(ip2)
        else
          x1 = vgx_EX(ip1)
          y1 = vgy_EX(ip1)
          x2 = vgx_EX(ip2)
          y2 = vgy_EX(ip2)
        endif
        x0  = 0.5d0*(x1+x2)
        y0  = 0.5d0*(y1+y2)
        write(n6,'(2x,4i5,6f9.4,1p4e12.3)')
     >   i, m, nw, iw, x0,y0,x1,y1,x2,y2,
     >   zflx, zpmg(m), bk_are(i), bk_deg(i)
      enddo
!
      write(n6,'(2x,"sub. set_bkgt  Sum_bk_flx =",1p10e12.4)')
     >   (zpmg(m), m = 1, bk_nty)
!
      return
!
!::dimension error
 900  continue
      call wexit("set_bkgt","ii.gt.ndwp")
      end
