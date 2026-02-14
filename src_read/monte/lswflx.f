!**********************************************************************
      subroutine lswflx(cdsn)
!**********************************************************************
      use cntcom, only : ipps, ipwl, npep, npew, npsp, npsw, xpnt, ypnt
     > , ncmax, icwl
      use cntwfl, only : xare, xcfl, xdeg, xeema, xeemm, xees, xehta
     >    , xehtm, xfhta, xfhti, xfhtm, xsno, xtea, xtees, xtem, xtfa
     >    , xtfi, xtfm
      use cphcns, only : cpi
      use csonic, only : itim, time
      use cunit,  only : lmspe, lmype, n6
      use mod_externalgrid
      implicit none
!
!::argument
      character  cdsn*(*)
!
!::local variables
      integer :: i6 = 21
      integer    isrc, nw, iw, iws, iwe, ip, ip1, ic, iw_p
      real*8     x1, y1, x2, y2, x0, y0, th, ar
!
      write(n6,'(/2x,"*** lswflx ***  cdsn =",a,2x,"itim =",i7,
     >  "  time =",1pe14.6,"  lmype =",i5)') cdsn, itim, time, lmype
!
      if(lmype.ne.lmspe ) return
      write(n6,'(2x,"execute lswflx  lmype =",i5)') lmype
!
      open( unit=i6, file=trim(cdsn), status='replace' )
      write(i6,'(2x,"   itim =",i7,"  time =",1pe14.6,"  lmype =",i4)')
     >    itim, time, lmype
      write(i6,'(2x,"nw",3x,"iw",3x,"Xc",9x,"Yc",9x,"th",7x,"are",8x,
     >  "tot",8x,9("i",a3,6x,"a",a3,6x,"m",a3,6x:))')
     > (xcfl(isrc),xcfl(isrc),xcfl(isrc),isrc=0,xsno)
!
      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
!
        do iw = iws, iwe
          ip  = ipwl(iw)
          if(use_exdata) then
            ip1 = ipwl2(iw)
          else
            if(iw == iwe) cycle
            ip1 = ipwl(iw+1)
          endif
          ic = icwl(iw)
          if(ic.le.ncmax .or. .not.use_exdata) then ! SOL region
            x1  = xpnt(ip)
            y1  = ypnt(ip)
            x2  = xpnt(ip1)
            y2  = ypnt(ip1)
          elseif(ic .le. ncmax+vac_ele_size) then
            x1 = vac_grid_x(ip)
            y1 = vac_grid_y(ip)
            x2 = vac_grid_x(ip1)
            y2 = vac_grid_y(ip1)
          elseif(ic .le. ncmax+vac_ele_size+pri_ele_size) then
            x1 = pri_grid_x(ip)
            y1 = pri_grid_y(ip)
            x2 = pri_grid_x(ip1)
            y2 = pri_grid_y(ip1)
          else
            x1 = vgx_EX(ip)
            y1 = vgy_EX(ip)
            x2 = vgx_EX(ip1)
            y2 = vgy_EX(ip1)
          endif
          x0  = 0.5d0*(x1+x2)
          y0  = 0.5d0*(y1+y2)
          th  = xdeg(iw)
!
      ! xthti is not 0 only for diverter plate
          write(i6,'(2x,i2,i5,2f11.6,f9.3,1p30e10.3)')
     >      nw, iw, x0, y0, th, xare(iw),
     >      xfhti(iw,0)+xfhta(iw,0)+xfhtm(iw,0),
     >     (xfhti(iw,isrc),xfhta(iw,isrc),xfhtm(iw,isrc),isrc=0,xsno)
        enddo
      enddo
!
      ! xtfi is not 0 only for diverter plate
      write(i6,'(2x,i2,i5,2f11.6,f9.3,1p30e10.3)')
     >  0, 0, 0.0, 0.0, 0.0, 0.0, xtfi(0)+xtfa(0)+xtfm(0),
     >  (xtfi(isrc),xtfa(isrc),xtfm(isrc),isrc=0,xsno)
      close(i6)
!
! energy flux 200601KH
      open( unit=i6, file="Wneflx.txt", status='replace' )
      write(i6,*)"# Power flux to/from wall"
      write(i6,'("# itim =",i7,"  time =",1pe14.6,"  lmype =",i4)')
     >    itim, time, lmype
      write(i6,'("# Atom(in)  = ",30e11.3)')(xtea(isrc,1),isrc=0,xsno)
      write(i6,'("# Atom(out) = ",30e11.3)')(xtea(isrc,2),isrc=0,xsno)
      write(i6,'("# Atom(net) = ",30e11.3)')
     >     (xtea(isrc,1)-xtea(isrc,2),isrc=0,xsno)
      write(i6,*)
      write(i6,'("# Mol(in)   = ",30e11.3)')(xtem(isrc,1),isrc=0,xsno)
      write(i6,'("# Mol(out ) = ",30e11.3)')(xtem(isrc,2),isrc=0,xsno)
      write(i6,'("# Mol(net)  = ",30e11.3)')
     >     (xtem(isrc,1)-xtem(isrc,2),isrc=0,xsno)
      write(i6,*)
      write(i6,'(2x,"nw",3x,"iw",3x,"Xc",9x,"Yc",9x,"th",7x,"are",8x,
     >  "atm-in",5x,"mol-in",5x,"atm-ref",
     >  4x,"mol-ref",4x,"atm-net",4x,"mol-net")')
!
      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
!
        do iw = iws, iwe
          ip  = ipwl(iw)
          if(use_exdata) then
            ip1 = ipwl2(iw)
          else
            if(iw == iwe) cycle
            ip1 = ipwl(iw+1)
          endif
          ic = icwl(iw)
          if(ic.le.ncmax .or. .not.use_exdata) then ! SOL region
            x1  = xpnt(ip)
            y1  = ypnt(ip)
            x2  = xpnt(ip1)
            y2  = ypnt(ip1)
          elseif(ic .le. ncmax+vac_ele_size) then
            x1 = vac_grid_x(ip)
            y1 = vac_grid_y(ip)
            x2 = vac_grid_x(ip1)
            y2 = vac_grid_y(ip1)
          elseif(ic .le. ncmax+vac_ele_size+pri_ele_size) then
            x1 = pri_grid_x(ip)
            y1 = pri_grid_y(ip)
            x2 = pri_grid_x(ip1)
            y2 = pri_grid_y(ip1)
          else
            x1 = vgx_EX(ip)
            y1 = vgy_EX(ip)
            x2 = vgx_EX(ip1)
            y2 = vgy_EX(ip1)
          endif
          x0  = 0.5d0*(x1+x2)
          y0  = 0.5d0*(y1+y2)
          th  = xdeg(iw)
          ar  = xare(iw)
!
          write(i6,'(2x,i2,i5,2f11.6,f9.3,1p30e11.3)')
     >      nw, iw, x0, y0, th, ar,
     >      xehta(iw,0), xehtm(iw,0),
     >      -xeema(iw,0), -xeemm(iw,0),
     >      xehta(iw,0)-xeema(iw,0),xehtm(iw,0)-xeemm(iw,0)
        enddo
      enddo
      close(i6)
!
!  energy flux from plasma boundary (sol(nw=1)/prv(nw=3))
      open( unit=i6, file="Wneflxs.txt", status='replace' )
      write(i6,'("# Power flux from plasma surface" )')
      write(i6,'("# itim =",i7,"  time =",1pe14.6,"  lmype =",i4)')
     >    itim, time, lmype
      write(i6,'("# Power flux from (sol/idp/prv/odp) = ",4e11.3)')
     >     (xtees(nw),nw=1,4)
      write(i6,*)
      write(i6,'(2x,"np",3x,"Xc",6x,"Yc",6x,"are",5x)')
!
      do nw = 1,4
        if( nw.ne.1 .and. nw.ne.3 ) cycle
        iws = npsp(nw)
        iwe = npep(nw)
        do iw_p = iws, iwe-1
          ip  = ipps(iw_p)
          ip1 = ipps(iw_p+1)
          x1  = xpnt(ip)
          y1  = ypnt(ip)
          x2  = xpnt(ip1)
          y2  = ypnt(ip1)
          x0  = 0.5d0*(x1+x2)
          y0  = 0.5d0*(y1+y2)
          ar  = sqrt((x1-x2)**2+(y1-y2)**2) * x0 * 2.0d0*cpi
          write(i6,'(i5,3f8.4,e16.9)')
     >       iw_p,x0,y0,ar,xees(iw_p,0)/ar
        enddo
      enddo

      close(i6)

      return
      end
