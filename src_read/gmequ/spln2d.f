!
!   modified by K. Shimizu    wk(ndwk) worging area   2010/05/19
!
!*********************************************************************
      subroutine spln1( nx, ny, fx, fy, fz, ndx, ndy )
!*********************************************************************
!
!
!::initial set  derives
!
!     iwk          : length of the work area wk  6*max(nx,ny)
!     wk(iwk)      : work array
!     ux(ndx,ny)   : du/dx   at the mesh points
!     uy(ndx,ny)   : du/dy   at the mesh points
!     uxy(ndx,ndy) : du/dxdy at the mesh points
!
!::cell set  coeff
!      cbc(4,4)    : (4,4) array which contains the coef.
!
!
!::evaluation
!      xp  : fx(ip-1) <= xp <= fx(ip) ! Note
!      yp  : fy(jp-1) <= yp <= fy(jp) ! Note
!      s   : fz  at (xp,yp)
!      sx  : dfz/dx at (xp,yp)
!      sy  : dfz/dy at (xp,yp)
!      sxy : dfz/dxdy at (xp,yp)
!      sxx : dfz/dxdx at (xp,yp)
!      syy : dfz/dydy at (xp,yp)
!
!-----------------------------------------------------------------------
      use cunit,  only : n6
      use clocal, only : dbx, dby, ndspc, ndspx, ndspy, ntx, nty, tbcf
     >    , tbic, tbx, tby
      implicit none
!
!::argument
      integer, intent(in) :: nx, ny, ndx, ndy
      real(8), intent(in) :: fx(nx), fy(ny), fz(ndx,ndy)
!
!::temporary
      real*8   ux(ndx,ndy), uy(ndx,ndy), uxy(ndx,ndy)
      real*8   cbc(4,4)
!
!::local variables
      integer ip, jp, i, j, ii, kk, jc ,ic
!
      write(n6,'(/2x,"*** spln1 ***")')
!
!::dimension check
      write(n6,'(2x,"nx =",i4,"  ndx =",i4,"  ndspx =",i4)')
     >  nx, ndx, ndspx
      write(n6,'(2x,"ny =",i4,"  ndy =",i4,"  ndspy =",i4)')
     >  ny, ndy, ndspy
      if( nx.gt.ndx .or. nx.gt.ndspx ) then
      write(n6,'(2x,"Stop at sub. spln1 due to x-dimension error")')
      call wexit("spln1","x-dimension error")
      endif
      if( ny.gt.ndy .or. ny.gt.ndspy ) then
      write(n6,'(2x,"Stop at sub. spln1 due to y-dimension error")')
      call wexit("spln1","y-dimension error")
      endif
!
!
      call derivs( nx, ny, fx, fy, ndx, fz,
     >             ux, uy, uxy )
!
!::set common variables
      ntx = nx
      nty = ny
      do i = 1, ntx
      tbx(i) = fx(i)
      enddo
      do j = 1, nty
      tby(j) = fy(j)
      enddo
      dbx = (tbx(ntx)-tbx(1))/dfloat(ntx-1)
      dby = (tby(nty)-tby(1))/dfloat(nty-1)
!
      write(n6,'(2x,"tbx =",2(i3,1pe14.6,2x),"  dbx =",1pe14.6)')
     >   1, tbx(1), ntx, tbx(ntx), dbx
      write(n6,'(2x,"tby =",2(i3,1pe14.6,2x),"  dby =",1pe14.6)')
     >   1, tby(1), nty, tby(nty), dby
!
!::set coef for spline
      ii = 0
      do j = 2, ny
      do i = 2, nx
      ip = i
      jp = j
      call coeff( nx, ny, fx, fy, ndx, fz, ux, uy, uxy,
     >   ip, jp, cbc )
!
      ii = ii + 1
      tbic(i,j) = ii
!
      kk = 0
      do jc = 1, 4
      do ic = 1, 4
      kk = kk + 1
      tbcf(ii,kk) = cbc(ic,jc)
      enddo
      enddo
!
      enddo
      enddo
!
      write(n6,'(2x,"icmax =",i6,"  ndspc =",i6)') ii, ndspc
      if( ii.gt.ndspc ) then
      write(n6,'(2x,"Stop at sub. spln1 due to dimension error")')
      call wexit("spln1","spln1 error")
      endif
!
      return
      end
!
!*********************************************************************
      subroutine spln2( xp, yp, s, sx, sy, sxy, sxx, syy )
!*********************************************************************
      use clocal, only : dbx, dby, ntx, nty, tbcf, tbic, tbx, tby, xim1
     >    , yjm1
      use cunit,  only : n6
      implicit none
      real(8), intent(in)  :: xp, yp
      real(8), intent(out) :: s, sx, sy, sxy, sxx, syy
!
!::local variables
      integer ip, jp, ic, kk, i, j
      real*8  cbc(4,4)
!
 100  continue
      ip = int((xp-tbx(1))/dbx + 2.0d0)  ! Note  tbx(ip-1)<xp<tbx(ip)
      jp = int((yp-tby(1))/dby + 2.0d0)  ! Note  tby(jp-1)<yp<tby(jp)
!
!::check
      if( ip.le.1 ) then
        write(n6,'(2x,"warning spln2  xp.le.tbx(1)   ",1p2e12.3,i5)')
     >    xp, tbx(1), ip
        goto 910
      elseif( ip.gt.ntx ) then
        write(n6,'(2x,"warning spnl2  xp.ge.tbx(ntx) ",1p2e12.3,i5)')
     >    xp, tbx(ntx), ip
        goto 910
      elseif( jp.le.1 ) then
        write(n6,'(2x,"warning spln2  yp.le.tby(1)   ",1p2e12.3,i5)')
     >    yp, tby(1), jp
        goto 910
      elseif( jp.gt.nty ) then
        write(n6,'(2x,"warning spln2  yp.ge.tby(nty) ",1p2e12.3,i5)')
     >    yp, tby(nty), jp
        goto 910
      endif
!
      xim1 = tbx(ip-1)
      yjm1 = tby(jp-1)
      ic = tbic(ip,jp)
!
      kk = 0
      do j = 1, 4
      do i = 1, 4
      kk = kk + 1
      cbc(i,j) = tbcf(ic,kk)
      enddo
      enddo
!
      call bicube( cbc, xp, yp, s, sx, sy, sxy, sxx, syy )
!
      return
!
 910  continue
      call wexit("spln2","out of point")
      end
