!***********************************************************************
      subroutine tmtrac(x0,y0,z0,vlx,vly,vlz,ic,ln,ts,tm,x1,y1,z1,ker)
!***********************************************************************
!
!     intersection 4 points  not rare case    2009/11/26
!     the same coll-time   ip1,ip2            2009/11/26
!
!    (vx,vy,vz) ==> (vlx,vly,vlz)  because cimcom contains vz
!                                             2010/05/13
!
!----------------------------------------------------------------------
! 3D axis(x,y,z) => xy:troidal direction, particle born at y=0(imstaw.f). z:toroidal axis
! 2D axis(r,z)   => poloidal surface defined in mesh(r:xpnt,z:ypnt).
!----------------------------------------------------------------------
      use cimcom, only : i6_trac, ipcm, isw2tr, iswtr, lpcm
      use cntcom, only : mgrd, mknd, mseg, next, xpnt, ypnt
      use cunit,  only : n6
      implicit none
!
!::argument
      real(8), intent(in)    :: x0, y0, z0, vlx, vly, vlz, ts
      real(8), intent(out)   :: x1, y1, z1, tm
      integer, intent(inout) :: ic, ln
      integer, intent(out)   :: ker

!::local variables
      integer ii, kmax, k, i, kln, ln0, ip, jj, nct, j
      integer ipa, ipb, ip1, ip2
      integer klsv(20), onsv(20), inum(20), ipnt(20)
      real*8  a, b, c, a2, b2, c2, za, zb, zc, zd2, zd
      real*8  ztmin, zt(2), t, tmsv(20), prsv(20), pzsv(20), dltm(20)
      real*8  x, y, z, r
      real*8  xs, ys, zs, rs, r0, ts0, eps, r_out
     > , r_1, r_2, z_1, z_2 
      integer  lmsg
      data   lmsg/0/
      save   lmsg
!
      integer  nxon(4); data  nxon/3,4,1,2/
      save     nxon
! function
      integer    imox, imoy
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
      kmax = mseg(ic)

!---------------------------------------------------------------------------------------------------------
!:: decide where to move from cell to cell. 
!---------------------------------------------------------------------------------------------------------
!:: Assume a cone with a toroidally rotated band of cell edges extended infinitely in the z-direction.
!:: Then, find the intersection of this cone and the particle trajectory.
!----------------------------------------------------------------------------------------------------------
      do k = 1, kmax
        ipa = mgrd(ic,k)
        ipb = mgrd(ic,k+1)
!-----
        ip1 = min0(ipa,ipb)
        ip2 = max0(ipa,ipb)
!----- points of the cell edge
        r_1 = xpnt(ip1) ! r: poloidal radias (xpnt)
        r_2 = xpnt(ip2)
        z_1 = ypnt(ip1) ! z: tridal axis (ypnt)
        z_2 = ypnt(ip2)
!-----
        a   = z_2-z_1
        b   = -(r_2-r_1)
        c   = r_2*z_1-r_1*z_2

!----------------------------------------------------------------------
!::z_1=z_2 case. Assume the disk with a hole, instead of a cone.
!----------------------------------------------------------------------
        if( a.eq.0.0d0 ) then
          if( b.eq.0.0d0 ) cycle ! ip1 and ip2 are at the same point
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
!
        if( zd2.lt.0.0d0 ) cycle
!
        zd = sqrt(zd2)
        zt(1) = (-zb+zd)/za
        zt(2) = (-zb-zd)/za
!
        do i = 1, 2
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
!
      enddo  ! loop(k)
!
      if( kln.eq.0 )  ztmin = 0.0d0
 100  continue
      tm = ztmin
      x1 = x0 + vlx*tm
      y1 = y0 + vly*tm
      z1 = z0 + vlz*tm
!
!-----------------------------------------------------------------
!::debug write  <==  KS
      if( iswtr.gt.0 .and. isw2tr.gt.0 ) then
        if( kln.ne.0 ) then
          write(i6_trac,'(2x,"*** tmtrac ***  O.K.  ",i10,i8)')
     >      lpcm, ipcm
        else
          write(i6_trac,'(2x,"*** tmtrac ***  ERR   ",i10,i8)')
     >      lpcm, ipcm
        endif
        ip = ipcm
!
        xs = x0 + vlx*ts
        ys = y0 + vly*ts
        zs = z0 + vlz*ts
        rs = sqrt(xs**2+ys**2)
        r0 = sqrt(x0**2+y0**2)
        r_out = sqrt(x1**2+y1**2)
!
        write(i6_trac,'(2x,"ic =",i7,"  ix,iy =",2i5,"  ln0 =",i6)')
     >   ic, imox(ic), imoy(ic), ln0
        write(i6_trac,'(2x,"Birth",1p3e14.6)') r0, z0, 0.0
        write(i6_trac,'(2x,"Start",1p3e14.6)') rs, zs, ts
        write(i6_trac,'(2x,"Hitp ",1p3e14.6,i3)') r_out, z1, tm, kln
        do k = 1, kmax+1
          write(i6_trac,'(2x,"grid ",1p2e14.6)')
     >     xpnt(mgrd(ic,k)),ypnt(mgrd(ic,k))
        enddo
        do i = 1, ii
          write(i6_trac,'(2x,"collp",1p4e14.6,1pe11.2,2i3)')
     >     prsv(i), pzsv(i), tmsv(i), ts, tmsv(i)-ts, onsv(i), klsv(i)
        enddo
      endif
!-----------------------------------------------------------------
!
!::found coll-point
      if( kln.ne.0 ) goto 200
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
     >     i2,2x,"lpcm =",i10,"  ipcm =",i8)') ln0, nct, lpcm, ipcm
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
     >     i2,"  lpcm =",i10,"  ipcm =",i8)') ln0, nct, lpcm, ipcm
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
     >     i2,2x,i10,i8)') ln0, nct, lpcm, ipcm
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
     >     "  ipcm =",i8)') lmsg, lpcm, ipcm
          ker = 1
          return
        endif
        goto 100
      endif
!
!::next cell
 200  continue
      ker = 0
      ln = mknd(ic,kln)
      ic = next(ic,kln)
      if( ln.ge.1 .and. ln.le.4 ) ln = nxon(ln)
      return
!
!::no found time
 910  continue
      ker = 1
      write(n6,'(2x,"*** tmtrac ***  ERR   ",i10,i8)') lpcm, ipcm
      xs = x0 + vlx*ts
      ys = y0 + vly*ts
      zs = z0 + vlz*ts
      rs = sqrt(xs**2+ys**2)
      r0 = sqrt(x0**2+y0**2)
      r_out = sqrt(x1**2+y1**2)
      write(n6,'(2x,"ic =",i7,"  ix,iy =",2i5,"  ln0 =",i6)')
     >   ic, imox(ic), imoy(ic), ln0
      write(n6,'(2x,"birth",1p3e14.6)') r0, z0, 0.0
      write(n6,'(2x,"start",1p3e14.6)') rs, zs, ts
      write(n6,'(2x,"hitp ",1p3e14.6,i3)') r_out, z1, tm, kln
      do k = 1, kmax+1
      ip = mgrd(ic,k)
      write(n6,'(2x,"grid ",1p2e14.6)') xpnt(ip), ypnt(ip)
      enddo
      do i = 1, ii
      write(n6,'(2x,"collp",1p4e14.6,1pe11.2,2i3)')
     >    prsv(i), pzsv(i), tmsv(i), ts, tmsv(i)-ts, onsv(i), klsv(i)
      enddo
!
      return
      end
