!***********************************************************************
      subroutine tmtrac_ex(x0,y0,z0,vlx,vly,vlz,ic,ln
     > ,ts,tm,x1,y1,z1,ker,kmax)
!***********************************************************************
      use cntcom, only:mgrd,xpnt,ypnt
      implicit none
!
!::argument
      real(8), intent(in)    :: x0, y0, z0, vlx, vly, vlz, ts
      real(8), intent(inout) :: x1, y1, z1, tm
      integer, intent(inout) :: ic, ln
      integer, intent(inout) :: ker
      integer, intent(in)    :: kmax
!
!::local variables
      integer  k,ip1,ip2
      real*8 x1_list(kmax),y1_list(kmax),x2_list(kmax),y2_list(kmax)

      do k = 1,kmax
          ip1 = mgrd(ic,k)
          ip2 = mgrd(ic,k+1)
          x1_list(k) = xpnt(ip1)
          y1_list(k) = ypnt(ip1)
          x2_list(k) = xpnt(ip2)
          y2_list(k) = ypnt(ip2)
      enddo
      call tmtrac_main(x0,y0,z0,vlx,vly,vlz,ic,ln,ts,tm
     >       ,x1,y1,z1,ker,x1_list,x2_list,y1_list,y2_list,kmax)
      end

!***********************************************************************
      subroutine tmtrac_grid(x0,y0,z0,vlx,vly,vlz,ic,ln
     > ,ts,tm,x1,y1,z1,ker
     > ,element,ele_size,kmax
     > ,grid_size,grid_x,grid_y,ic_in)
!***********************************************************************
      use csize, only:ndms
      use cunit, only:n6
      implicit none
!
!::argument
      real(8), intent(in)    :: x0, y0, z0, vlx, vly, vlz, ts
      real(8), intent(inout) :: x1, y1, z1, tm
      integer, intent(inout) :: ic, ln
      integer, intent(inout) :: ker
      integer, intent(in)    :: ele_size,grid_size,kmax,ic_in
      integer, intent(in)    :: element(ele_size,ndms)
      real*8,  intent(in)    :: grid_x(grid_size),grid_y(grid_size)
!
!::local variables
      integer  k, k1
      real*8 x1_list(kmax),y1_list(kmax),x2_list(kmax),y2_list(kmax)
!
      do k = 1,kmax
          k1 = k+1
          if( k.eq.kmax ) k1 = 1
          x1_list(k) = grid_x(element(ic_in,k))
          y1_list(k) = grid_y(element(ic_in,k))
          x2_list(k) = grid_x(element(ic_in,k1))
          y2_list(k) = grid_y(element(ic_in,k1))
      enddo
      call tmtrac_main(x0,y0,z0,vlx,vly,vlz,ic,ln,ts,tm
     >       ,x1,y1,z1,ker,x1_list,x2_list,y1_list,y2_list,kmax)
      end
!***********************************************************************
      subroutine tmtrac_main(x0,y0,z0,vlx,vly,vlz,ic,ln,ts,tm
     >       ,x1,y1,z1,ker,x1_list,x2_list,y1_list,y2_list,kmax)
!***********************************************************************
!
!     intersection 4 points  not rare case    2009/11/26
!     the same coll-time   ip1,ip2            2009/11/26
!
!    (vx,vy,vz) ==> (vlx,vly,vlz)  because cimcom contains vz
!                                             2010/05/13
!
!----------------------------------------------------------------------
      use cimcom, only : i6_trac, ipcm, lpcm
      use cntcom, only : mknd, next, ncmax
      use cunit,  only : n6, mype
      use mod_externalgrid
      implicit none
!
!::argument
      real*8, intent(in)     :: x0, y0, z0, vlx, vly, vlz, ts 
      real*8, intent(out)    :: x1, y1, z1, tm
      integer, intent(inout) :: ic, ln
      integer, intent(out)   :: ker
      integer, intent(in)    :: kmax
      real*8, intent(in)     :: x1_list(kmax),y1_list(kmax)
     > ,x2_list(kmax),y2_list(kmax)

!::local variables
      integer ii, k, i, kln, ln0, jj, nct, j, ic_new,k1
      integer klsv(20), onsv(20), inum(20), ipnt(20)
      real*8  a, b, c, a2, b2, c2, za, zb, zc, zd2, zd
      real*8  ztmin, zt(2), t, tmsv(20), prsv(20), pzsv(20), dltm(20)
      real*8  x, y, z, r
      real*8  xs, ys, zs, rs, r0, ts0, eps, r1
      integer  lmsg
      data   lmsg/0/
      save   lmsg
      real*8 r_1,r_2,z_1,z_2,r_line,z_line

! debug
      real*8 za_save(kmax),zd2_save(kmax)
! function
      integer    imox, imoy
!
      za_save = 777.0
      zd2_save =777.0
!
      if( lpcm.eq.1 ) lmsg = 0
!
      eps = 0.5d-12
      ln0 = ln
      if( ln0.eq.0 ) eps = 0.0d0
      ts0 = ts + eps
!
      ii  = 0
      kln = 0
      ztmin = 1.0d20
!
      do k = 1, kmax
!----- points of the cell edge
        r_1 = x1_list(k) ! r: poloidal radias
        r_2 = x2_list(k)
        z_1 = y1_list(k) ! z: tridal axis
        z_2 = y2_list(k)
!-----
        a   = z_2-z_1
        b   = -(r_2-r_1)
        c   = r_2*z_1-r_1*z_2

!----------------------------------------------------------------------
!::z_1=z_2 case. Assume the disk with a hole, instead of a cone.
!----------------------------------------------------------------------
        if( a.eq.0.0d0 ) then
          if( b.eq.0.0d0 ) cycle ! 2 points are at the same point
          if( vlz.eq.0.0d0 ) then
            cycle ! never cross
          else
            t = (z_1-z0)/vlz
          endif
          x = x0 + vlx*t
          y = y0 + vly*t
          z = z0 + vlz*t
          r = sqrt(x**2+y**2)
          ii = ii + 1
          tmsv(ii) = t
          klsv(ii) = k
          prsv(ii) = r
          pzsv(ii) = z
          onsv(ii) = 0
          if( (r-r_1)*(r-r_2).gt.0.0d0 ) cycle  ! <== KS
          onsv(ii) = 1
          if( t.gt.ts0 .and. t.lt.ztmin ) then       ! <== KS
            ztmin = t
            kln = k
          endif
          cycle
        endif

!----------------------------------------------------------------------------------------------------------------------
!:: solve quadratic equation for find the time(t) at which the trajectory and the cone intersect (including r1=r2 case).
!:: Then we get 2 real solution t+,t- (=zt)
!:: t+ = {-zb+sqrt(zb^2-za*zc)}/za, t- = {-zb-sqrt(zb^2-za*zc)}/za
!----------------------------------------------------------------------------------------------------------------------
        a2 = a*a
        b2 = b*b
        c2 = c*c
        za = a2*(vlx**2+vly**2)-b2*vlz**2
        zb = a2*(x0*vlx+y0*vly)-b2*z0*vlz-b*c*vlz
        zc = a2*(x0**2+y0**2)-b2*z0**2-2.0d0*b*c*z0-c2
        zd2 = zb*zb-za*zc
        za_save(k) = za
        zd2_save(k)= zd2
!
        if( zd2.lt.0.0d0 ) cycle ! no real solution
        if( za .eq.0.0d0 ) cycle ! only 1 real solution
!
        zd = sqrt(zd2)
        zt(1) = (-zb+zd)/za 
        zt(2) = (-zb-zd)/za
!
        do i = 1, 2 ! two real solution
          t = zt(i)
          x = x0 + vlx*t
          y = y0 + vly*t
          z = z0 + vlz*t
          r = sqrt(x**2+y**2)
          ii = ii + 1
          tmsv(ii) = zt(i)
          klsv(ii) = k
          prsv(ii) = r
          pzsv(ii) = z
          onsv(ii) = 0
!
          if( r_1.eq.r_2 ) then
            if( (z-z_1)*(z-z_2).gt.0.0d0 ) cycle
          else
            if( (r-r_1)*(r-r_2).gt.0.0d0 .or.
     >          (z-z_1)*(z-z_2).gt.0.0d0 ) cycle
          endif
!
          onsv(ii) = 1
          if( t.gt.ts0 .and. t.lt.ztmin ) then
            ztmin = t
            kln = k
          endif
        enddo  ! loop(i)
      enddo  ! loop(k)
!
!-----------------------------------------------------------------
      if( kln.eq.0 )  ztmin = 0.0d0
 100  continue
      tm = ztmin
      x1 = x0 + vlx*tm
      y1 = y0 + vly*tm
      z1 = z0 + vlz*tm
!-----------------------------------------------------------------
!
!::found coll-point
      if( kln.eq.0 ) then
!
!::no found coll-point (on line)
        if( ln0.ne.0 ) then
          jj = 0
          do i = 1, ii
            if( onsv(i).eq.0 ) cycle
            jj = jj + 1
            dltm(jj) = dabs(tmsv(i) - ts)
            inum(jj) = i
          enddo
          nct = jj
          if( nct.le.1 ) then
            write(n6,'(2x,"*** tmtrac ***  ln0 =",i4,"   ERR  nct =",
     >        i2,2x,"lpcm =",i10,"  ipcm =",i8)') ln0, nct, lpcm, ipcm
            goto 910
          endif
          call tsort(nct,dltm,ipnt)
          j = ipnt(2)
          i = inum(j)
          ztmin = tmsv(i)
          kln = klsv(i)
          lmsg = lmsg + 1
          if( lmsg.gt.100 ) then
            write(999,'(2x,"*** tmtrac ***  ln0 =",i4,"   ERR  nct =",
     >       i2,"  lpcm =",i10,"  ipcm =",i8)') ln0, nct, lpcm, ipcm
            ker = 1
            return
          endif
          goto 100
!
!::no found coll-point (inner point)
        else
          jj = 0
          do i = 1, ii
            if( onsv(i).eq.0 ) cycle
            if( tmsv(i).le.ts ) cycle
            jj = jj + 1
            dltm(jj) = tmsv(i) - ts
            inum(jj) = i
          enddo
          nct = jj
          if( nct.le.0 ) then
            write(n6,'(2x,"*** tmtrac ***  ln0 =",i4,"   ERR  nct =",
     >       i2,2x,i10,i8)') ln0, nct, lpcm, ipcm
            goto 910
          endif
          call tsort(nct,dltm,ipnt)
          j = ipnt(1)
          i = inum(j)
          ztmin = tmsv(i)
          kln = klsv(i)
          lmsg = lmsg + 1
          if( lmsg.gt.100 ) then
            write(999,'(2x,"tmtrac error  lmsg =",i5,"  lpcm =",i10,
     >       "  ipcm =",i8)') lmsg, lpcm, ipcm
            ker = 1
            return
          endif
          goto 100
        endif
      endif !not found cell
!
!::next cell
 200  continue
      ker = 0
      ln = mknd(ic,kln)
      ic_new = next(ic,kln)
!
      if(ic_new.eq.0 .and. ln.ge.0) then
        r_line = sqrt(x1**2+y1**2)
        z_line = z1
        call crossBorder(r_line,z_line,ic,ic_new,ker)
      endif
!
      ic = ic_new
      return
!
!::no found time
 910  continue
      ker = 1
      ! debug yamamoto 24/10/1
      write(40000+mype,'(2x,"*** tmtrac ***  ERR")')
      xs = x0 + vlx*ts
      ys = y0 + vly*ts
      zs = z0 + vlz*ts
      rs = sqrt(xs**2+ys**2)
      r0 = sqrt(x0**2+y0**2)
      write(40000+mype,'(2x,"ic =",i7,"  ix,iy =",2i5,"  ln =",i6)')
     >   ic, imox(ic), imoy(ic), ln
      write(40000+mype,'(2x,"birth",1p2e14.6)') r0, z0
      write(40000+mype,'(2x,"start",1p3e14.6)') rs, zs, ts
      write(40000+mype,'(2x,"veloc",1p3e14.6)') vlx, vly, vlz
      do k = 1, kmax
        write(40000+mype,'(2x,"grid ",1p4e14.6)') 
     >   x1_list(k),y1_list(k),za_save(k),zd2_save(k)
      enddo
      do i = 1, ii
        write(40000+mype,'(2x,"collp",1p4e14.6,1pe11.2,2i3)') 
     >   prsv(i),pzsv(i),tmsv(i),ts, tmsv(i)-ts,onsv(i), klsv(i)
      enddo
!
      return
      end