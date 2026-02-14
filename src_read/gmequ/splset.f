!***********************************************************************
      subroutine splset
!***********************************************************************
      use com_eqdat, only : nr, nz, psi, rg, zg
      use cunit,     only : n6
      implicit none
!
!::local variables
! modified 2/2 lines organize local variables and include files by kamata 2021/06/16
!ik   integer  iz, ir, ii, i, j
!ik   real*8   xs, ys, ps, sx, sy, sxy, sxx, syy
      integer  iz, ir, i, j
      real*8   sx, sy, sxy, sxx, syy
      real*8   xa, ya, xb, yb, xc, yc, xd, yd, xp, yp
      real*8   pa, pb, pc, pd, pp, p2
! function
      real*8   fpsi
      fpsi(i,j) = psi(nr*(j-1) + i)
!
!::spline
      call spln1( nr, nz, rg, zg, psi, nr, nz )
!
      write(n6,'(/2x,"spline check  cell-center")')
      iz = nz/2
      do ir = 2, nr-1, 10
!
!::A (ir,iz)
      xa = rg(ir)
      ya = zg(iz)
      pa = fpsi(ir,iz)
!::B (ir+1,iz)
      xb = rg(ir+1)
      yb = zg(iz)
      pb = fpsi(ir+1,iz)
!::C (ir+1,iz+1)
      xc = rg(ir+1)
      yc = zg(iz+1)
      pc = fpsi(ir+1,iz+1)
!::D (ir,iz+1)
      xd = rg(ir)
      yd = zg(iz+1)
      pd = fpsi(ir,iz+1)
!::P
      xp = (xa+xb+xc+xd)/4.0d0
      yp = (ya+yb+yc+yd)/4.0d0
      pp = (pa+pb+pc+pd)/4.0d0
!::spline
      call spln2( xp, yp, p2, sx, sy, sxy, sxx, syy )
!
      write(n6,'(2x,2i4,4f11.5)')
     >  ir, iz, xp, yp, pp, p2
      enddo
!
      return
      end
