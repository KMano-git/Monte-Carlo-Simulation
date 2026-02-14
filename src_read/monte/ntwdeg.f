!**********************************************************************
      subroutine ntwdeg
!**********************************************************************
      use cntcom, only : icwl, iplx, iply, ipwl, next, npew, npsw, xpnt
     >    , ypnt
      use cntwfl, only : iywl, jxwl, xare, xdeg
      use cphcns, only : cpi
      use cpmpls, only : xmd1, xmd2, ymd1, ymd2
      use csize,  only : ndwp
      use cunit,  only : n6
      use mod_externalgrid, only:use_exdata
      implicit none
!
!::local variables
      real*8  xpc, ypc, x1, y1, x2, y2, x0, y0, dl, dd, th
      integer nw, iw, iws, iwe, ip, ip1
      integer ik, ic, ii, jx, iy

      if(use_exdata) then
        call ntwdeg_ex
        return
      endif
!
      write(n6,'(2x,"*** ntwdeg ***")')
!
!::plasma center
      xpc = 0.5d0*(xmd1+xmd2)
      ypc = 0.5d0*(ymd1+ymd2)
      write(n6,'(2x,"plasma center  (xpc,ypc) = ",2f10.5)') xpc, ypc
!
!::clear
      call setd( xare, ndwp, 0.0d0 )
      call setd( xdeg, ndwp, 0.0d0 )
!
!::xdeg, xare
      do nw = 1, 4
      iws = npsw(nw)
      iwe = npew(nw)
!
      do iw = iws, iwe-1
      ip  = ipwl(iw)
      ip1 = ipwl(iw+1)
      x1  = xpnt(ip)
      y1  = ypnt(ip)
      x2  = xpnt(ip1)
      y2  = ypnt(ip1)
      x0  = 0.5d0*(x1+x2)
      y0  = 0.5d0*(y1+y2)
      dl  = sqrt((x2-x1)**2+(y2-y1)**2)
      dd  = sqrt((x0-xpc)**2+(y0-ypc)**2)
      th  = datan2( (y0-ypc)/dd, (x0-xpc)/dd )
      if( th.lt.0.0d0 ) th = th + 2.0d0*cpi
      th = th/cpi*180.0d0
      xdeg(iw) = th
      xare(iw) = 2.0d0*cpi*x0*dl
      enddo
      enddo
!
!::jxwl,iywl
      write(n6,'(2x,"define  jxwl, iywl")')
      do nw = 1, 4
      iws = npsw(nw)
      iwe = npew(nw)
      if( nw.eq.1 ) ik = 3
      if( nw.eq.2 ) ik = 4
      if( nw.eq.3 ) ik = 1
      if( nw.eq.4 ) ik = 2
      do iw = iws, iwe-1
      ic = icwl(iw)
      ii = 0
 100  continue
      ii = ii + 1
      if( ii.gt.5 ) goto 120
      jx = iplx(ic)
      iy = iply(ic)
      if( jx.eq.0 .or. iy.eq.0 ) then
      ic = next(ic,ik)
      goto 100
      endif
!
 120  continue
      jxwl(iw) = jx
      iywl(iw) = iy
      enddo
      enddo
!
      return
      end

!**********************************************************************
      subroutine ntwdeg_ex
!**********************************************************************
      use cntcom, only : xpnt, ypnt, ncmax, npsw, npew,icwl,ipwl
     > ,iplx, iply
      use cntwfl, only : iywl, jxwl, xare, xdeg
      use cphcns, only : cpi
      use cpmpls, only : xmd1, xmd2, ymd1, ymd2
      use csize,  only : ndwp
      use cunit,  only : n6
      use mod_externalgrid
      implicit none
!
!::local variables
      real*8  xpc, ypc, x1, y1, x2, y2, x0, y0, dl, dd, th
      integer nw, iw, iws, iwe, ip, ip1,ic
!
      write(n6,'(2x,"*** ntwdeg_ex ***")')
!
!::plasma center
      xpc = 0.5d0*(xmd1+xmd2)
      ypc = 0.5d0*(ymd1+ymd2)
      write(n6,'(2x,"plasma center  (xpc,ypc) = ",2f10.5)') xpc, ypc
!
!::clear
      call setd( xare, ndwp, 0.0d0 )
      call setd( xdeg, ndwp, 0.0d0 )
!
!::xdeg, xare
      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe
          ic = icwl(iw)
          ip  = ipwl(iw)
          ip1 = ipwl2(iw)
          if(ic.le.ncmax) then
            x1  = xpnt(ip)
            y1  = ypnt(ip)
            x2  = xpnt(ip1)
            y2  = ypnt(ip1)
          elseif(ic.le.ncmax+vac_ele_size)then
            x1  = vac_grid_x(ip)
            y1  = vac_grid_y(ip)
            x2  = vac_grid_x(ip1)
            y2  = vac_grid_y(ip1)
          elseif(ic.le.ncmax+vac_ele_size+pri_ele_size)then
            x1  = pri_grid_x(ip)
            y1  = pri_grid_y(ip)
            x2  = pri_grid_x(ip1)
            y2  = pri_grid_y(ip1)
          else
            x1  = vgx_EX(ip)
            y1  = vgy_EX(ip)
            x2  = vgx_EX(ip1)
            y2  = vgy_EX(ip1)
          endif
          x0  = 0.5d0*(x1+x2)
          y0  = 0.5d0*(y1+y2)
          dl  = sqrt((x2-x1)**2+(y2-y1)**2)
          dd  = sqrt((x0-xpc)**2+(y0-ypc)**2)
          th  = datan2( (y0-ypc)/dd, (x0-xpc)/dd )
          if( th.lt.0.0d0 ) th = th + 2.0d0*cpi
          th = th/cpi*180.0d0
          xdeg(iw) = th
          xare(iw) = 2.0d0*cpi*x0*dl
        enddo
      enddo

!::jxwl,iywl
      write(n6,'(2x,"define  jxwl, iywl")')
      do nw = 1, 4
        ! only for diverter region
        if( nw.ne.2 .and. nw.ne.4) cycle
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe
          ip  = ipwl(iw)
          ip1 = ipwl2(iw)
          ic = icwl(iw)
          if(ic.gt.ncmax) cycle
          jxwl(iw) = iplx(ic)
          iywl(iw) = iply(ic)
        enddo
      enddo

      return
      end