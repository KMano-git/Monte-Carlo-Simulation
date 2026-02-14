!**********************************************************************
      subroutine sptim(xp,yp,vx,vy,ic,ln,wt,wx,wy,ker)
!**********************************************************************
      use cntcom, only : iptl, mgrd, mknd, mseg, next, weit, xpnt, ypnt
      use cunit,  only : mype, n6
      implicit none
!
!::argument
      real(8), intent(in)    :: xp, yp, vx, vy
      real(8), intent(out)   :: wt, wx, wy
      integer, intent(inout) :: ic
      integer, intent(inout) :: ln
      integer, intent(out)   :: ker
!
!::local variables
      integer  kmax, k, ip1, ip2, kn
      real*8   zdx, zdy, zpx, zpy, zaa, zss, ztt
!
      integer  nxon(4)
      data nxon/3,4,1,2/
!
      kmax = mseg(ic)
!
      do k = 1, kmax
        if( mknd(ic,k).eq.ln ) cycle
        ip1 = mgrd(ic,k)
        ip2 = mgrd(ic,k+1)
        zdx = xpnt(ip2)-xpnt(ip1)
        zdy = ypnt(ip2)-ypnt(ip1)
!
        zpx = xp - xpnt(ip1)
        zpy = yp - ypnt(ip1)
        zaa = zdx*vy - zdy*vx
        if( zaa.eq.0.0d0 ) goto 910
        zss = (zpx*vy-zpy*vx)/zaa
        ztt = (zpx*zdy-zpy*zdx)/zaa
!
        if( zss.gt.0.0d0 .and. zss.le.1.0d0 ) then
          if( ztt.ge.0.0d0 ) then
            wt = ztt
            wx = xpnt(ip1) + zdx*zss
            wy = ypnt(ip1) + zdy*zss
            kn = k
            ln = mknd(ic,kn)
            goto 215
          endif
        endif
      enddo
!
      ker = 1
      write(n6,'(2x,"*** sptim  no found col-point  mype =",i5,
     > "  iptl =",i6,"  ic,ln =",2i5,1pe12.3)')
     >   mype,iptl,ic,ln,weit
      return
 215  continue
!
!::next cell
      ker  = 0
      ic = next(ic,kn)
      if( ln.ge.1 .and. ln.le.4 ) ln = nxon(ln)
!
      return
!
 910  continue
      write(n6,'(2x,"*** sptim  zaa.eq.0.0    ",1pe12.3)') zaa
      ker = 2
      return
      end
