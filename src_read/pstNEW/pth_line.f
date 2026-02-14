!**********************************************************************
      subroutine pth_line(xst,yst,xen,yen, ndim,ncol,posx,posy,ipos)
!**********************************************************************
!
!       path line  (xst,yst) ==> (xen,yen)
!
!----------------------------------------------------------------------         
      use csize
      use cntcom
      use csonic
      use cunit
      use com_pth_line, only : xwlp, ywlp, iwlp, nwlp
      implicit none
!
!::argument
      real*8  xst, yst, xen, yen
      integer ndim
      real*8  posx(ndim),posy(ndim)
      integer ncol, ipos(ndim)

      integer ncl
      real*8  xcl(2), ycl(2), tcl(2), scl(2)
      integer inon(4)
      data inon/3,4,1,2/
!
      integer iw, lon, ic, ii, k, j1, j2, kmax, ip1, ip2, kn, i, in
      real*8  vex, vey, zpx, zpy, zdx, zdy, zqx, zqy
      real*8  zaa, zss, ztt, wtm, wpx, wpy
      integer kmn
      real*8  wtmn
!
      integer idbg
      parameter (idbg=0)

!::start point
      vex = xen - xst
      vey = yen - yst
      call wlcold(xst,yst,xen,yen,ncl,xcl,ycl,tcl,scl,
     >            nwlp,xwlp,ywlp)
!          input:xst,yst, xen,yen, nwlp,xwlp,ywlp,   output:ncl,xcl,ycl,tcl,scl
!
      if(idbg.eq.1) then
      write(n6,'(2x,"----- trace of pth_line ----- ncl =",i2,
     >  "  xst,yst =",2f10.5,"  xen,yen =",2f10.5,"  ncl =",i2)')
     >   ncl, xst, yst, xen, yen, ncl
      write(n6,'(2x,4f10.5)') (xcl(i),ycl(i),tcl(i),scl(i),i=1,ncl)
      endif
!
!::no collision
      if(ncl.eq.0) then
          ncol = 0
      return
      endif
!
!::start
      wtm  = tcl(1)
      wpx  = xcl(1)
      wpy  = ycl(1)
      in   = scl(1)     ! <== Note
      iw   = iwlp(in)   ! <== Note
      lon  = -iw
      ic   = icwl(iw)
      ixpo = migx(ic)
      iypo = migy(ic)
      ii = 1
      posx(ii) = wpx
      posy(ii) = wpy
      ipos(ii) = ic
!
  200 continue
          if(idbg.eq.1)
     >        write(6,622) ii, ixpo, iypo, wtm, wpx, wpy, lon
  622     format(2x,3i5,3f10.5,i5)
!
          if(ii.gt.ndim) then
              write(6,601) ii,ndim
  601         format(2x,'***  sub pth_line  ii.gt.ndim  ',2i5)
              stop
          endif
!
          kmax = mseg(ic)
          do k = 1, kmax
              if( mknd(ic,k).eq.lon ) cycle
              kn  = k
              ip1 = mgrd(ic,k)
              ip2 = mgrd(ic,k+1)
              zdx = xpnt(ip2) - xpnt(ip1)
              zdy = ypnt(ip2) - ypnt(ip1)
              zpx = xst - xpnt(ip1)
              zpy = yst - ypnt(ip1)
              zaa = zdx*vey - zdy*vex

              if(zaa.eq.0.0d0) goto 910

              zss = (zpx*vey-zpy*vex)/zaa
              ztt = (zpx*zdy-zpy*zdx)/zaa
              if(zss.gt.0.0d0 .and. zss.lt.1.0d0) then
                  if(ztt.ge.0.0d0) go to 350
              endif
          enddo
!
          write(n6,'(2x,"pth_line  no coll. point")')
          ncol = -1
          return
!
  350     continue
          wtm = ztt
          wpx = xpnt(ip1) + zdx*zss
          wpy = ypnt(ip1) + zdy*zss
          lon = mknd(ic,kn)
!
          if(lon.gt.0) then
              ic   = next(ic,kn)
              ixpo = migx(ic)
              iypo = migy(ic)
              lon  = inon(lon)  !  center-cell  mseg(ic) = 3
          else
              ic = 0
              ixpo = 0
              iypo = 0
          endif
!
          ii = ii + 1
          posx(ii) = wpx
          posy(ii) = wpy
          ipos(ii) = ic
          if(lon.lt.0) go to 400
!xx   if( lon.lt.0 .or. lon.gt.10 ) goto 400
      go to 200
!
  400 continue
      if(idbg.eq.1)
     >    write(6,622) ii, ixpo, iypo, wtm, wpx, wpy, lon
      ncol = ii
      return
!
 910  continue
      write(n6,'(2x,"*** pth_line  zaa.eq.0.0    ",1pe12.3)') zaa
      call gdmod(0)
      stop
      end
!
!**********************************************************************
      subroutine set_pwall
!**********************************************************************
      use csize
      use cntcom
      use csonic
      use cunit
      use com_pth_line, only : xwlp, ywlp, iwlp, nwlp
      use mod_externalgrid, only : use_exdata
      implicit none

!::local variables
      integer  ii, k, iws, iwe, ip, iw, ic, nw, i
      integer  imox, imoy
!
      write(n6,'(/2x,"*** set_pwall ***")')
!
      ii = 0
      do k = 1, 4
        nw = k
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe
          if(.not.use_exdata .and. nw.ne.4 .and. iw.eq.iwe) cycle
          ii = ii + 1
          ip = ipwl(iw)
          xwlp(ii) = xpnt(ip)
          ywlp(ii) = ypnt(ip)
          iwlp(ii) = iw
        enddo
      enddo
      nwlp = ii
!
!::debug write
!      do i = 1, nwlp
!      iw = iwlp(i)
!      ic = icwl(iw)
!      write(n6,'(2x,i3,2x,i3,2x,i7,2x,i5,2x,i2,2x,2f10.5)')
!     >  i, iw, ic, imox(ic), imoy(ic), xwlp(i), ywlp(i)
!      enddo
!
      return
      end
