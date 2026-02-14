!**********************************************************************
      subroutine mcpnt(jc,ic,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
!**********************************************************************
!
!    point name
!                   ---t---         * 1,2,3,4
!            ---4--/      / --3---     r1 = grdx(mx1,my1)
!              /             /         z1 = grdy(mx1,my1)
!      ic     p     o       q       * s,t
!            /             /           rs = gcvx(mxs,mys)
!           /  ---s----   /            zs = gcvy(mxs,mys)
!       ---1--/       /--2---       * p,q
!                jc                    rp = 1/2*(r1+r4)
!                                      zp = 1/2*(z1+z4)
!                                   * o
!    Note  mxs = mx1, mys = my1        ro = 1/2*(rs+rt)
!          mxt = mx4, myt = my4        zo = 1/2*(zs+zt)
!
!----------------------------------------------------------------------
      use cplmet, only : kgdx, kgdy
      implicit none
!
      integer, intent(in)  :: jc, ic
      integer, intent(out) :: mx1, mx2, mx3, mx4, my1, my2, my3, my4
!
      if( jc.le.0 .or. ic.le.0 ) then
      write(6,'(2x,"mcpnt  jc, ic =",2i5)') jc, ic
      call wexit("mcpnt","invarid jc,ic")
      endif
!
      mx1 = kgdx(jc,ic,1)
      mx2 = kgdx(jc,ic,2)
      mx3 = kgdx(jc,ic,3)
      mx4 = kgdx(jc,ic,4)
!
      my1 = kgdy(jc,ic,1)
      my2 = kgdy(jc,ic,2)
      my3 = kgdy(jc,ic,3)
      my4 = kgdy(jc,ic,4)
!
      return
      end
!
!**********************************************************************
      subroutine mcpntp(jc,ic,gx1,gx2,gx3,gx4,gy1,gy2,gy3,gy4)
!**********************************************************************
      use cgdcom, only : grdx, grdy
      use cplmet, only : kgdx, kgdy
      implicit none
!
      integer, intent(in)  :: jc, ic
      real(8), intent(out) :: gx1, gx2, gx3, gx4, gy1, gy2, gy3, gy4
!
!::local vriables
      integer  mx1, mx2, mx3, mx4, my1, my2, my3, my4
!
      mx1 = kgdx(jc,ic,1)
      mx2 = kgdx(jc,ic,2)
      mx3 = kgdx(jc,ic,3)
      mx4 = kgdx(jc,ic,4)
!
      my1 = kgdy(jc,ic,1)
      my2 = kgdy(jc,ic,2)
      my3 = kgdy(jc,ic,3)
      my4 = kgdy(jc,ic,4)
!
      gx1 = grdx(mx1,my1)
      gx2 = grdx(mx2,my2)
      gx3 = grdx(mx3,my3)
      gx4 = grdx(mx4,my4)
!
      gy1 = grdy(mx1,my1)
      gy2 = grdy(mx2,my2)
      gy3 = grdy(mx3,my3)
      gy4 = grdy(mx4,my4)
!
      return
      end
