!**********************************************************************
      subroutine ntwarp(xp,yp,vx,vy,ic,ln,wt,wx,wy,ker)
!**********************************************************************
!)
!)   move particle from plasma region to vacume region at the gate
!)
!)   This routine is extended to treat in the case of
!)         hitting gate (ic=0) and  shevron (ic>0)
!)
!)     gate
!)                    plasma side
!)      |__|___|___|___|___|___|  gate
!)      |    V   |       |     |        mknd(ic,k) < 0
!)          ic=0      vaccum side       next(ic,k) = 0
!)
!)     shevron/partition
!)
!)      |__|_ _|__|__|__|__|__|__| S1
!)      |  | V |  |  |  |  |  |  | S2   mknd(ic,k) < 0  a kind of wall
!)                                      next(ic,k) > 0
!)
!)---------------------------------------------------------------------
      use cntcom, only : icgt, icwl, igtw, ipgt, iptl, istp, iwgt
     >   , migx, migy, mknd, mrgn, mseg, ncmax2, negt, next, nogt, nsgt
     >   , velx, vely, weit, xpnt, xpos, ypnt, ypos, ncmax
      use cunit,  only : n6
      use mod_externalgrid
      implicit none
!
!::argument
      real(8), intent(in)    :: xp, yp, vx, vy
      real(8), intent(inout) :: wt
      real(8), intent(out)   :: wx, wy
      integer, intent(inout) :: ic, ln
      integer, intent(out)   :: ker
!
!::local variables
      integer icnew, icold, ks, iw, ir, icx, icy, kg, kso
      integer ngs, nge, ng, ip1, ip2, k, kmax, ic_tmp
      real*8  zdx, zdy, zpx, zpy, zaa, zss, ztt, x1, y1, x2, y2
!
!::icold          ! 2003/11/21
      icnew = ic
      iw    = -ln
      ic    = icwl(iw)
      icold = ic
      kg = igtw(iw)
!
      wx = xpos + wt*velx
      wy = ypos + wt*vely
      istp = istp + 1
!
!---------------------------------------------------------------------
!::shevron/partition
!---------------------------------------------------------------------
      if( kg <= 0 .or. kg > nogt ) then
        ks = 2 ! subdiverter region
        if( icnew.gt.0 ) then
          ker = 0
          goto 215   ! 2003/11/21
        endif
        write(n6,'(2x,"ntwarp  kg =",i2,"  icnew =",i7)')
     >    kg, icnew
        call wexit("ntwarp","kgt(not gate)  icnew = 0")
      endif
!
!---------------------------------------------------------------------
!::gate
!---------------------------------------------------------------------
!::kg = 1/2  gate number
!::ks = 1 (plasma side)  = 2 (vacume side)
      if((.not. use_exdata) .and. (ic.le.0 .or. ic.gt.ncmax2) ) then
        call wexit("ntwarp","ic wrong")
      endif
      ks = 1
      ir = mrgn(ic)
      if( ir.eq.8 ) ks = 2
!
!
!::new ks
      kso = ks
      ks  = 3 - kso
!
!-----------------------------------------------------------------------
!::new cell
      ngs = nsgt(kg,ks)
      nge = negt(kg,ks)
!
      do ng = ngs, nge
        if(use_exdata) then
          ip1 = ipgt(ng)
          ip2 = ipgt2(ng)
          if(ks == 1) then ! plasma to vacume
            x1 = pri_grid_x(ip1)
            y1 = pri_grid_y(ip1)
            x2 = pri_grid_x(ip2)
            y2 = pri_grid_y(ip2)
          else ! vacume to plasma
            x1 = vgx_EX(ip1)
            y1 = vgy_EX(ip1)
            x2 = vgx_EX(ip2)
            y2 = vgy_EX(ip2)
          endif
        else
          if(ng.eq.nge) cycle
          ip1 = ipgt(ng)
          ip2 = ipgt(ng+1)
          x1 = xpnt(ip1)
          x2 = xpnt(ip2)
          y1 = ypnt(ip1)
          y2 = ypnt(ip2)
        endif
        zdx = x2 - x1
        zdy = y2 - y1
        zpx = xp - x1
        zpy = yp - y1
        zaa = zdx*vy - zdy*vx
        if( zaa.eq.0.0d0 ) then
          write(n6,'(2x,"*** sptim  zaa.eq.0.0    ",1pe12.3)') zaa
          ker = 2
          return
        endif
        zss = (zpx*vy-zpy*vx)/zaa
        if( zss.gt.0.0d0 .and. zss.le.1.0d0 ) then
          ztt = (zpx*zdy-zpy*zdx)/zaa
          if( ztt.ge.0.0d0 ) then
            wt = ztt
            wx = x1 + zdx*zss
            wy = y1 + zdy*zss
            ker = 0
            ic  = icgt(ng)
            iw  = iwgt(ng)
            icnew = ic
            goto 215
          endif
        endif
      enddo
!
      ker = 1
      if(use_exdata) then
        write(n6,'(2x,"ntwarp  no found col-point",
     >   "  iptl =",i6,"  ic=",i6,1pe12.3)')
     >     iptl,ic,weit
      else
        write(n6,'(2x,"ntwarp  no found col-point",
     >    "  iptl =",i6,"  ic,icx,icy =",3i6,1pe12.3)')
     >    iptl,ic,migx(ic),migy(ic),weit
      endif
      return
!
!-----------------------------------------------------------------------
 215  continue
      ic  = icnew
      if((.not.use_exdata) .and. (ic.le.0 .or. ic.gt.ncmax2) ) then
        call wexit("ntwarp","ic wrong")
      endif
!
!::check next to go
      if(use_exdata) then
        ! set kmax
        if(ks == 1) then ! plasma to vacume
          ic_tmp = ic - (ncmax+vac_ele_size)
          kmax = mseg_pri(ic_tmp)
        else ! vacume to plasma
          ic_tmp = ic - (ncmax+vac_ele_size+pri_ele_size)
          kmax = mseg_subdiv(ic_tmp)
        endif
        ! check
        do k = 1, kmax
          ln = mknd(ic,k)
          if( ln.eq.-iw ) goto 220
          if( next(ic,k).eq.icold ) goto 220  ! added 2002/12/13
        enddo
        !::no found next cell
        write(n6,'(2x,"ntwarp  no found next cell",
     >    "  iptl =",i6,"  ic =",i6)')
     >    iptl, ic
        do k = 1, kmax
          ln = mknd(ic,k)
          write(n6,'(2x,"  k =",i2,"  ln =",i6,"  iw =",i6,
     >    "  next =",i6,"  icold =",i6)') k,ln,iw,next(ic,k),icold
        enddo
      else ! ordinary mesh data
        do k = 1, mseg(ic)
          ln = mknd(ic,k)
          if( ln.eq.-iw ) goto 220
          if( next(ic,k).eq.icold ) goto 220  ! added 2002/12/13
        enddo
        icx = migx(ic)
        icy = migy(ic)
        ir  = mrgn(ic)
        write(n6,'(2x,"ntwarp  no found next cell",
     >    "  iptl =",i6,"  ic,icx,icy =",3i6,"  ir =",i2)')
     >    iptl, ic, icx, icy, ir
        do k = 1, mseg(ic)
          ln = mknd(ic,k)
          write(n6,'(2x,"  k =",i2,"  ln =",i6,"  iw =",i6,
     >    "  next =",i6,"  icold =",i6)') k,ln,iw,next(ic,k),icold
        enddo
      endif
      call wexit("ntwarp","ic wrong")
!
!:: success
 220  continue
      istp = istp + 1
      if( ic.le.0 .or. ic.gt.ncmax2 ) then
        call wexit("ntwarp","ic wrong")
      endif
      end
