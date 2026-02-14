!***********************************************************************
      subroutine ntwpls
!***********************************************************************
!
!      plasma surface in neut2d-code  <np> = 1, 2, 3, 4, 5
!
!                rcypwl   +-iaxs--------------------+     <--
!               pwal      !     <----------         !   pwal (3)
!            /!########## !~impl~~~~~~<5>~~~~~~~~~~~! ##########!/
!           |/!%%ppsf%%%% !                         ! %%%<3>%%%%!/
!         o |/!%        ! !     Main Plasma         ! !        %!/  i
!         d V/!%--------! !-------------------------! !--------%!/  d
!         p  /!%        ! !                         ! !  IN    %!/  p
!         l  /!<4>      ! !      Scrape-Off         ! !   DIV <%>/  l
!            /!%opsf    ! !                         ! !    ipsf%!/A
!        (4) /!%%%%%%%%%% %%%%%%%%%%%%%<1>%%%%%spsf%% %%%%%%%%%%!/| (2)
!    rcydpl  /!         ! !      Void               ! !         !/| rcydpl
!            /########### ########################### ##########!/
!                                    (1)  -->            swal
!                rcydwl             rcyswl                rcydwl
!
!                      Note.  counter clock wise
!
!                                                       ik      3
!     nops  :  number of plasma surface 5 <= ndwl          *********
!     npps  :  number of P-surf point     <= ndwp         4*       *
!     cpsf(5)  : "Spsf", "Ipsf", "Ppsf", "Opsf", "Mpsf"    *       *2
!     npsp(5)  : start point of P-surf                     *********
!     npep(5)  : end point                                     1
!     nbxp(5)  : x-index of plasma parameter at boundary
!     nbyp(5)  : y-index of plasma parameter at boundary
!
!     ikps(ndwp) : type of cell boundary (1/2/3/4)
!     ipps(ndwp) : grid number in monte carlo
!     icps(ndwp) : cell number in monte carlo (not dummy cell)
!     isxp(ndwp) : x-cell number of boundary cell in soldor
!     isyp(ndwp) : y-cell number of boundary cell in soldor
!     csps(ndwp) : cosine of P-surf segment
!     snps(ndwp) : sine of P-surf segment
!     rcps(ndwp) : recycling coefficients
!
!     Note
!        index of flux in soldor
!             j=isxp(iw),i=isyp(iw)  qfx_df(j,i,m1a)
!
!        index of plasma parameter at boundary
!             j=isxp(iw)+nbxp(nw),i=isyp(iw)+nbyp(nw)  vne(j,i)
!
!         ibc = 1 :  sol wall
!             = 2 :  inner divertor plate
!             = 3 :  private wall
!             = 4 :  outer divertor plate
!             = 5 :  core edge
!
!      index(only for diverter wall)  imwl(ndwp) : wall ==> plasma surface
!                                     imps(ndwp) : plasma surface ==> wall
!
!-----------------------------------------------------------------------
      use cntcom, only : cpsf, csps, icps, icwl, ikps, imps, imwl, iplx
     >    , iply, ipps, isxp, isyp, mcel, mgrd, migx, migy, mseg
     >    , nbxp, nbyp, nops, npep, npew, npps, npsp, npsw, rcps
     >    , snps, xpnt, ypnt
      use cplmet, only : icel, itmpe, itpve, jcel, jtmax
      use csize,  only : ndwp
      use cunit,  only : n6
      use mod_externalgrid, only : use_exdata
      implicit none
!
!::local variables
      integer  ii, np, ik, it, jt, j, i, ig, ic, k
      integer  ifx, ify, jf, if
      integer  mx1, mx2, mx3, mx4, my1, my2, my3, my4
      integer  nw, iw, ip1, ip2, ip, iwp
      integer  i2, icp, icw, iww, ixp, iyp
      real*8   dlx, dly, dln, x0, y0
      character  cmsg*1000
      integer  idbg/2/; save idbg
!
      write(n6,'(/2x,"*** ntwpls ***")')
!
!-----------------------------------------------------------------------
!::total number
!-----------------------------------------------------------------------
      ii = 0
!
!-----------------------------------------------------------------------
!::sol P-surf
!-----------------------------------------------------------------------
      if( idbg.ge.2 ) write(n6,'(2x,"Spsf")')
      np  = 1
      cpsf(np) = "Spsf"
      it = 2; ik  = 1; ifx = 0; ify = -1
!
      npsp(np) = ii + 1
      do jt = 2, jtmax(it)
      j = jcel(jt,it)
      i = icel(jt,it)
      jf = j + ifx
      if = i + ify
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
      ii = ii + 1
      if( ii.gt.ndwp ) goto 910
      ic = mcel(j,i)
      if( ic.gt.0 ) ig = mgrd(ic,ik)
      if( ic.eq.0 ) ig = mgrd(icps(ii-1),ik+1)
      ikps(ii) = ik
      ipps(ii) = ig
      icps(ii) = ic
      isxp(ii) = jf
      isyp(ii) = if
      enddo
      npep(np) = ii
      nbxp(np) = 0
      nbyp(np) = 0
!
!-----------------------------------------------------------------------
!::inner divertor P-surf
!-----------------------------------------------------------------------
      if( idbg.ge.2 ) write(n6,'(2x,"Ipsf")')
      np  = 2
      cpsf(np) = "Ipsf"
      ik = 2; ifx = 0; ify = 0
!
      npsp(np) = ii + 1
      do it = 2, itpve
      jt = jtmax(it) - 1
      j  = jcel(jt,it)
      i  = icel(jt,it)
      jf = j + ifx
      if = i + ify
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
      ii = ii + 1
      if( ii.gt.ndwp ) goto 910
      ic = mcel(j,i)
      if( ic.gt.0 ) ig = mgrd(ic,ik)
      if( ic.eq.0 ) ig = mgrd(icps(ii-1),ik+1)
      ikps(ii) = ik
      ipps(ii) = ig
      icps(ii) = ic
      isxp(ii) = jf
      isyp(ii) = if
      enddo
      npep(np) = ii
      nbxp(np) = 1
      nbyp(np) = 0
!
!-----------------------------------------------------------------------
!::prv P-surf
!-----------------------------------------------------------------------
      if( idbg.ge.2 ) write(n6,'(2x,"Ppsf")')
      np  = 3
      cpsf(np) = "Ppsf"
      it = itpve-1; ik  = 3; ifx = 0; ify = 0
!
      npsp(np) = ii + 1
      do jt = jtmax(it)-1, 1, -1
      j = jcel(jt,it)
      i = icel(jt,it)
      jf = j + ifx
      if = i + ify
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
      ii = ii + 1
      if( ii.gt.ndwp ) goto 910
      ic = mcel(j,i)
      if( ic.gt.0 ) ig = mgrd(ic,ik)
      if( ic.eq.0 ) ig = mgrd(icps(ii-1),ik+1)
      ikps(ii) = ik
      ipps(ii) = ig
      icps(ii) = ic
      isxp(ii) = jf
      isyp(ii) = if
      enddo
      npep(np) = ii
      nbxp(np) = 0
      nbyp(np) = 1
!
!-----------------------------------------------------------------------
!::outer divertor P-surf
!-----------------------------------------------------------------------
      if( idbg.ge.2 ) write(n6,'(2x,"Opsf")')
      np  = 4
      cpsf(np) = "Opsf"
      ik = 4; ifx = -1; ify = 0
!
      npsp(np) = ii + 1
      do it = itpve-1, 1, -1
      jt = 2
      j  = jcel(jt,it)
      i  = icel(jt,it)
      jf = j + ifx
      if = i + ify
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
      ii = ii + 1
      if( ii.gt.ndwp ) goto 910
      ic = mcel(j,i)
      if( ic.gt.0 ) ig = mgrd(ic,ik)
      if( ic.eq.0 ) ig = mgrd(icps(ii-1),ik+1)
      ikps(ii) = ik
      ipps(ii) = ig
      icps(ii) = ic
      isxp(ii) = jf
      isyp(ii) = if
      enddo
      npep(np) = ii
      nbxp(np) = 0
      nbyp(np) = 0
!
!-----------------------------------------------------------------------
!::core edge
!-----------------------------------------------------------------------
      if( idbg.ge.2 ) write(n6,'(2x,"Mpsf")')
      np  = 5
      cpsf(np) = "Mpsf"
      it = itmpe; ik = 3; ifx = 0; ify = -1
!
      npsp(np) = ii + 1
      do jt = jtmax(it)-1, 1, -1
      j = jcel(jt,it)
      i = icel(jt,it)
      jf = j + ifx
      if = i + ify
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
      ii = ii + 1
      if( ii.gt.ndwp ) goto 910
      ic = mcel(j,i)
      if( ic.gt.0 ) ig = mgrd(ic,ik)
      if( ic.eq.0 ) ig = mgrd(icps(ii-1),ik+1)  ! <== not dummy cell
      ikps(ii) = ik
      ipps(ii) = ig
      icps(ii) = ic
      isxp(ii) = jf
      isyp(ii) = if
      enddo
      npep(np) = ii
      nbxp(np) = 0
      nbyp(np) = 1
!
!-----------------------------------------------------------------------
!::point number
!-----------------------------------------------------------------------
      nops = 5
      npps = ii
!
!-----------------------------------------------------------------------
!::cos & sin of segment
!-----------------------------------------------------------------------
      do nw = 1, nops
      do iw = npsp(nw), npep(nw)-1
      ip1 = ipps(iw)
      ip2 = ipps(iw+1)
      dlx = xpnt(ip2)-xpnt(ip1)
      dly = ypnt(ip2)-ypnt(ip1)
      dln = dsqrt(dlx*dlx+dly*dly)
      csps(iw) = dlx/dln
      snps(iw) = dly/dln
      rcps(iw) = 1.0d0
      enddo
      iw = npep(nw)
      csps(iw) = 0.0d0
      snps(iw) = 0.0d0
      rcps(iw) = 0.0d0
      enddo
!
!-----------------------------------------------------------------------
!::plasma surface <==> wall (for diverter)
!-----------------------------------------------------------------------
      do nw = 2, 4, 2 !(nw=2:inner diverter, nw=4:outer diverter)
        do iw = npsp(nw), npep(nw)-1
          ic = icps(iw)
          do i2 = npsw(nw), npew(nw)
            iww = i2
            if( icwl(iww).eq.ic ) goto 170
          enddo
          ! fail to find cell
          ip1 = ipps(iw)
          ip2 = ipps(iw+1)
          write(cmsg,'("no found icwl(iww): ",3i6,4f13.5)') 
     >       nw, iw, ic
     >      ,xpnt(ip1), ypnt(ip1), xpnt(ip2), ypnt(ip2)
          call wexit("ntwpls",trim(cmsg))
 170      continue
          imwl(iww) = iw
          imps(iw)  = iww
        enddo  ! loop(i1)
      enddo    ! loop(nw)
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      if( idbg.ge.1 ) then
      write(n6,'(2x,"*** ntwpls ***")')
      write(n6,'(3x,"nw",4x,"iw",1x,"ik",3x,"ip",3x,"xp",8x,"yp",9x,
     > "ic",1x,"icx",1x,"icy",3x,"xc",8x,"yc",8x,"cs",8x,"sn",8x,"rc")')
      do nw = 1, nops
      do iw = npsp(nw), npep(nw)
      ic = icps(iw)
      ip = ipps(iw)
!-----
      x0 = 0.0d0
      y0 = 0.0d0
      if( ic.gt.0 ) then
      do k = 1, mseg(ic)
      x0 = x0 + xpnt(mgrd(ic,k))
      y0 = y0 + ypnt(mgrd(ic,k))
      enddo
      x0 = x0/dfloat(mseg(ic))
      y0 = y0/dfloat(mseg(ic))
      endif
!-----
      write(n6,'(2x,i3,i6,i3,i5,2f10.5, i6,2i4,5f10.5,2i5)')
     >  nw, iw, ikps(iw), ip, xpnt(ip), ypnt(ip),
     >  ic, iplx(ic), iply(ic), x0, y0,
     >  csps(iw), snps(iw), rcps(iw),
     >  isxp(iw)+nbxp(nw), isyp(iw)+nbyp(nw)
      enddo
      enddo
!
      write(n6,'(/2x,"*** ntwpls ***  ( wall <==> psurf )")')
      do nw = 2, 4, 2 !(nw=2:inner diverter, nw=4:outer diverter)
      do iw = npsw(nw), npew(nw)
      if(.not.use_exdata .and. iw==npew(nw)) cycle
      iww = iw
      icw = icwl(iw)
      iwp = imwl(iw)
      if( iwp.gt.0 ) then
      icp = icps(iwp)
      ixp = migx(icp); iyp = migy(icp)
      else
      icp = 0
      ixp = 0; iyp = 0
      endif
      write(n6,'(2x,i3,2x,"wall",2i6,2i4,4x,"pls",2i6,2i4)')
     >  nw, iww,icw,migx(icw),migy(icw), iwp,icp,ixp,iyp
      enddo
      enddo
!
      endif
!
      return
!
 910  continue
      write(cmsg,'("dimension error  ii.gt.ndwp  ",2i6)') ii,ndwp
      call wexit("plntpsf",trim(cmsg))
      end
