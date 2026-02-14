!***********************************************************************
      subroutine tmwgti(r0,z0,ic,ln,r1,z1,ker)
!***********************************************************************
!
!     intersection 4 points  not rare case    2009/11/26
!     the same coll-time   ip1,ip2            2009/11/26
!
!
!          gate                (r1,z1)
!      *---------*---------*-----X---*---------*---------*
!      |-----------------------------> dsgt(i)
!                    A
!                    |         (r0,z0)
!      *-------------------------X-----------------------*
!
!----------------------------------------------------------------------
      use cntcom, only : dsgt, icgt, icwl, igtw, ipgt, iwgt, mgrd, mknd
     >    , mrgn, mseg, negt, nsgt, xpnt, ypnt
      use cunit,  only : n6
      use mod_externalgrid, only : use_exdata
      implicit none
!
!::argument
      real(8), intent(in)    :: r0, z0
      integer, intent(inout) :: ic, ln
      real(8), intent(out)   :: r1, z1
      integer, intent(out)   :: ker
!
!::local variables
      integer  iw, ir, ks, kg, ks0, ngs, nge, ips
      integer  ig, ign, k, kn, ip1, ip2
      real*8   ds, hs, hs0
!
!::for external grid version
      if(use_exdata) then
        r1 = r0
        z1 = z0
        ker = 0
        call tmwgti_grid(r0,z0,ic,ln)
        return
      endif

!::gate kg, ks
      iw = -ln
      ic = icwl(iw)
      ir = mrgn(ic)
      ks0 = 1                ! plasma side
      if( ir.eq.8 ) ks0 = 2  ! vacume size
      kg  = igtw(iw)
      ks  = 3 - ks0
      ngs = nsgt(kg,ks)       ! start point of gate
      nge = negt(kg,ks)       ! end point of gate
!
!::distance of coll. point from start point of gate
      ips = ipgt(ngs)
      ds  = sqrt((r0-xpnt(ips))**2+(z0-ypnt(ips))**2)
!
      do ig = ngs, nge
        ign = ig
        if( ds.le.dsgt(ign) ) goto 110
      enddo
      write(n6,'(2x,"error at sub. tmwari  no found ign")')
      ker = 1
      return
!
 110  continue
      ig = ign - 1      ! Note
      iw = iwgt(ig)
      ic = icgt(ig)
      do k = 1, mseg(ic)
        kn = k
        ln = mknd(ic,k)
        if( ln.eq.(-iw) ) goto 120
      enddo
      write(n6,'(2x,"error at sub. tmwari  no found ln")')
      return

 120  continue
      ln  = mknd(ic,kn) ! = -iw
      ip1 = mgrd(ic,kn)
      ip2 = mgrd(ic,kn+1)
      hs0 = (ds-dsgt(ig))/(dsgt(ig+1)-dsgt(ig))
      hs  = dmax1(hs0,0.00001d0)
      hs  = dmin1(hs0,0.99999d0)
      r1  = xpnt(ip1) + hs*(xpnt(ip2)-xpnt(ip1)) ! r1 = r0 ???
      z1  = ypnt(ip1) + hs*(ypnt(ip2)-ypnt(ip1)) ! z1 = z0 ???
      ker = 0
      return
!
      end subroutine tmwgti

!***********************************************************************
      subroutine tmwgti_grid(x0,y0,ic,ln)
!***********************************************************************
      use cntcom, only : icgt, icwl, igtw, ipgt, iwgt, mknd
     >    , mrgn, negt, nsgt, ncmax
      use mod_externalgrid

!::argument
      real(8), intent(in)    :: x0, y0
      integer, intent(inout) :: ic, ln
!
!::local variables
      integer  iw, ir, ks, kg, ks0, ngs, nge
      integer  ig, k, ip1, ip2, kmax, ic_tmp
      real*8   xg_1, yg_1, xg_2, yg_2

!::gate kg, ks
      iw = -ln
      ic = icwl(iw)
      ir = mrgn(ic)
      ks0 = 1                ! plasma side
      if( ir.eq.8 ) ks0 = 2  ! vacume size
      kg  = igtw(iw)
!
      ks  = 3 - ks0
      ngs = nsgt(kg,ks)       ! start point of gate
      nge = negt(kg,ks)       ! end point of gate
!
      igsx = 0
      igsy = 0
      do ig = ngs, nge
        ip1 = ipgt(ig)
        ip2 = ipgt2(ig)
        ! xg,yg: grid point in the destination region
        if(ks == 1) then ! vacume to plasma
          xg_1 = pri_grid_x(ip1)
          yg_1 = pri_grid_y(ip1)
          xg_2 = pri_grid_x(ip2)
          yg_2 = pri_grid_y(ip2)
        else ! plasma to vacume 
          xg_1 = vgx_EX(ip1)
          yg_1 = vgy_EX(ip1)
          xg_2 = vgx_EX(ip2)
          yg_2 = vgy_EX(ip2)
        endif
        sx  = (x0-xg_1)/(xg_2-xg_1)
        sy  = (y0-yg_1)/(yg_2-yg_1)
        if( sx.ge.0.0d0 .and. sx.le.1.0d0 ) igsx = ig
        if( sy.ge.0.0d0 .and. sy.le.1.0d0 ) igsy = ig
      enddo
!
      if( igsx.le.0 .or. igsx.ne.igsy ) then
        call wexit("tmwgti_grid","fail to find gate cell to go")
      endif
!
      ig = igsx
      iw = iwgt(ig)
      ic = icgt(ig)
      if(ks == 1) then
        ic_tmp = ic-(ncmax+vac_ele_size)
        kmax = mseg_pri(ic_tmp)
      else
        ic_tmp = ic-(ncmax+vac_ele_size+pri_ele_size)
        kmax = mseg_subdiv(ic_tmp)
      endif
      do k = 1, kmax
        ln = mknd(ic,k)
        if( ln.eq.(-iw) ) return
      enddo
      call wexit("tmwgti_grid","no found ln")
      end subroutine tmwgti_grid