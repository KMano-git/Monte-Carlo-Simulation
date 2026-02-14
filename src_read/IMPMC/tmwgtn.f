!***********************************************************************
      subroutine tmwgtn(x0,y0,z0,vtx,vty,vtz,ic,ln,ts,tm,x1,y1,z1,ker)
!***********************************************************************
!
!     intersection 4 points  not rare case    2009/11/26
!     the same coll-time   ip1,ip2            2009/11/26
!
!       donot use vz, ir  They are common bariables.
!         (vx,vy,vz) => (vtx,vty,vtz)
!           ir       =>  irg
!
!    ip   ig   iw  igt     ic    ix   iy   x1           y1           z1
!           vtx          vty          vtz          wgt
!          r1           z1
!
!----------------------------------------------------------------------
      use cimcom, only : amz, ipcm, wght
      use cimi60, only : i60_tmwgtn
      use cimpuf, only : lbkstw
      use cntcom, only : icgt, icwl, igtw, ipgt, iwgt, mrgn, negt
     >    , nsgt, xpnt, ypnt, mseg, ncmax
      use cphcns, only : cev
      use mod_shexe, only : impmc_model
      use mod_externalgrid
      use cunit, only : mype
      implicit none
!
!::argument
      real(8), intent(in)    :: x0, y0, z0, vtx, vty, vtz, ts
      real(8), intent(out)   :: x1, y1, z1, tm
      integer, intent(inout) :: ic, ln
      integer, intent(out)   :: ker
!
!::local variables
      integer ip
      integer :: i60
      integer ip1, ip2
      real*8  xs, ys, zs, rs, r0, ts0, eps, r1, e0
      real*8  sx, sy
      integer :: iw, irg, ks, kg, kso, ngs, nge
      integer :: ig, igsx, igsy, kcn, kmax, ic_tmp
      real(8) :: xg_1, xg_2, yg_1, yg_2
!
      real(8) :: dln, vl
      real(8) :: x2, y2, z2, r2
      integer :: lpbk = 0
      integer :: lpbkmx = 2000
!
! function
      integer    imox, imoy
!
      eps = -0.5d-8    ! <== cross with large angle  KSFUJI
      ts0 = ts + eps
!
!::gate
      ker = 0
      iw = -ln
      ic = icwl(iw)
      irg = mrgn(ic)
      ks = 1                  ! plasma side
      if( irg.eq.8 ) ks = 2   ! vacume size
      kg  = igtw(iw)
      kso = ks
      ks  = 3 - kso           ! ks=1:vacume to plasma, ks=2:plasma to vacume
      ngs = nsgt(kg,ks)       ! start point of gate
      nge = negt(kg,ks)       ! end point of gate
!
      xs = x0 + ts*vtx
      ys = y0 + ts*vty
      zs = z0 + ts*vtz
      rs = sqrt(xs**2+ys**2)
      r0 = sqrt(x0**2+y0**2) !debug
      ip1 = ipgt(ngs)
!
      igsx = 0
      igsy = 0
      do ig = ngs, nge
        ip1 = ipgt(ig)
        if(use_exdata) then
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
        else
          if(ig.eq.nge) cycle
          ip2 = ipgt(ig+1)
          xg_1 = xpnt(ip1)
          yg_1 = ypnt(ip1)
          xg_2 = xpnt(ip2)
          yg_2 = ypnt(ip2)
        endif
        sx  = (rs-xg_1)/(xg_2-xg_1)
        sy  = (zs-yg_1)/(yg_2-yg_1)
        if( sx.ge.0.0d0 .and. sx.le.1.0d0 ) igsx = ig
        if( sy.ge.0.0d0 .and. sy.le.1.0d0 ) igsy = ig
      enddo
!
!::next cell
      if( igsx.le.0 .or. igsx.ne.igsy ) then
        call go_through(xs,ys,zs,ts,vtx,vty,vtz,ic,tm,x1,y1,z1,ker)
        ln = 0
        if(ker .ne. 0) then
          x1 = xs
          y1 = ys
          z1 = zs
          tm = ts
          ic = icwl(iw)
        endif
        return
      endif
!
      ig = igsx
      ic = icgt(ig)
      ln = 0
      if(.not. use_exdata) then ! ordinary mesh
        call tmtrac(x0,y0,z0,vtx,vty,vtz,ic,ln,ts0,tm,x1,y1,z1,kcn)
      else ! the mesh made by Gmesh and external tool
        if(ic.le.ncmax) then ! SOL region
          kmax = mseg(ic)
          call tmtrac_ex(x0,y0,z0,vtx,vty,vtz,ic,ln
     >      ,ts0,tm,x1,y1,z1,kcn,kmax)
        elseif(ic .le. ncmax+vac_ele_size) then ! vacume region
          ic_tmp = ic-ncmax
          kmax = mseg_vacume(ic_tmp)
          call tmtrac_grid(x0,y0,z0,vtx,vty,vtz,ic,ln
     >      ,ts0,tm,x1,y1,z1,kcn
     >      ,vac_element,vac_ele_size,kmax
     >      ,vac_grid_size,vac_grid_x,vac_grid_y,ic_tmp)
        elseif(ic .le. ncmax+vac_ele_size+pri_ele_size) then ! private region
          ic_tmp = ic-(ncmax+vac_ele_size)
          kmax = mseg_pri(ic_tmp)
          call tmtrac_grid(x0,y0,z0,vtx,vty,vtz,ic,ln
     >      ,ts0,tm,x1,y1,z1,kcn
     >      ,pri_element,pri_ele_size,kmax
     >      ,pri_grid_size,pri_grid_x,pri_grid_y,ic_tmp)
        else ! sub-diverter
          ic_tmp = ic-(ncmax+vac_ele_size+pri_ele_size)
          kmax = mseg_subdiv(ic_tmp)
          call tmtrac_grid(x0,y0,z0,vtx,vty,vtz,ic,ln
     >      ,ts0,tm,x1,y1,z1,kcn
     >      ,subdiv_cell,cell_size,kmax
     >      ,grid_size,vgx_EX,vgy_EX,ic_tmp)
        endif
      endif
!
      ic = icgt(ig)  ! next cell ic = 0 in sub. tmtrac
      r1 = sqrt(x1**2+y1**2)
!
!::position and velocity at gate to plasma region
      if( lbkstw.eq.1 .and. ks.eq.1 ) then
      i60 = i60_tmwgtn
      if( i60.gt.0 .and. lpbk.le.lpbkmx ) then
      lpbk = lpbk + 1
      ip = ipcm
      dln = 0.05d0
      vl = sqrt(vtx**2+vty**2+vtz**2)
      x2 = x1 + dln*vtx/vl
      y2 = y1 + dln*vty/vl
      z2 = z1 + dln*vtz/vl
      r2 = sqrt(x2**2+y2**2)
      e0 = 0.5d0*amz*vl**2/cev
      write(i60,'(/2x,i8,3i5,i7,i6,i5,1p11e13.4)')
     >   ip, ig, iw, igtw(iw), ic, imox(ic), imoy(ic),
     >   x1, y1, z1, vtx, vty, vtz, wght(ip), r1, z1
      write(i60,'(2x,8x,15x,7x,6x,5x,78x,1p3e13.4)')  e0, r2, z2
      endif
      endif
!
      ker = kcn
!
!***********************************************************************
      contains
      subroutine go_through(xs,ys,zs,ts,vtx,vty,vtz,ic,tm
     > ,x1,y1,z1,ker )
!***********************************************************************
! retry for through gate.
! particle go forward 0.5mm(dist_delta) with its velocity.
! if there is no cell that the particle contains, error exit(e.g. go in the doom)
!
!::argument
      real*8, intent(in)  :: xs,ys,zs,ts,vtx,vty,vtz
      real*8, intent(out) :: tm,x1,y1,z1
      integer, intent(out) :: ic, ker
!
!::local variables
      real*8, parameter :: dist_delta = 0.0005d0 ! you can change this value
      real*8  t_delta, r_line
      integer ic_loop
!
      t_delta = dist_delta/sqrt(vtx**2+vty**2+vtz**2)
      tm = ts + t_delta
      x1 = xs + vtx*t_delta
      y1 = ys + vty*t_delta
      z1 = zs + vtz*t_delta
      r_line = sqrt(x1**2+y1**2)
!
      do ic_loop = 1, ncmax+vac_ele_size+pri_ele_size+cell_size
        call mchkin(r_line,z1,ic_loop,ker)
        if(ker.eq.0) then
          ic = ic_loop
          return
        endif
      enddo
      end subroutine go_through
!
      end subroutine tmwgtn