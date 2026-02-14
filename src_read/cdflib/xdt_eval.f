!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   function xdt_eval2(ixs,vx1,vx2)
      real(8) function xdt_eval2(ixs,vx1,vx2)
!***********************************************************************
      use cxdcom, only : pin, plg, pxa, pxd, pxi, pxs, xdaty
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/06/16
!ik   real*8 xdt_eval2, vx1, vx2
!ik   integer ixs
      real(8), intent(in) :: vx1, vx2
      integer, intent(in) :: ixs
!
!::local variables
      real*8  vy
      real*8  zex1, zix1, zdx1; integer ix1a, ix1b
      real*8  zex2, zix2, zdx2; integer ix2a, ix2b
      integer iyaa, iyab, iyba, iybb
!
!-----------------------------------------------------------------------
!
      zex1 = vx1
      zex1 = dmax1( zex1, pxi(1,ixs) )
      zex1 = dmin1( zex1, pxa(1,ixs) )
      if( plg(2,ixs).eq.1 ) zex1 = log(zex1)
!
      zex2 = vx2
      zex2 = dmax1( zex2, pxi(2,ixs) )
      zex2 = dmin1( zex2, pxa(2,ixs) )
      if( plg(3,ixs).eq.1 ) zex2 = log(zex2)
!
      zix1 = (zex1-pxs(1,ixs))/pxd(1,ixs) + 1.0d0
      ix1a = int(zix1)
      ix1a = min0( ix1a, pin(2,ixs)-1 )
      ix1b = ix1a + 1
      zdx1 = zix1 - dble( ix1a )
!
      zix2 = (zex2-pxs(2,ixs))/pxd(2,ixs) + 1.0d0
      ix2a = int(zix2)
      ix2a = min0( ix2a, pin(3,ixs)-1 )
      ix2b = ix2a + 1
      zdx2 = zix2 - dble( ix2a )
!
      iyaa = ix1a + pin(2,ixs)*(ix2a-1) + pin(1,ixs)
      iyba = ix1b + pin(2,ixs)*(ix2a-1) + pin(1,ixs)
      iybb = ix1b + pin(2,ixs)*(ix2b-1) + pin(1,ixs)
      iyab = ix1a + pin(2,ixs)*(ix2b-1) + pin(1,ixs)
!
      vy  = xdaty(iyaa)*(1.0d0-zdx1)*(1.0d0-zdx2) +
     >       xdaty(iyba)*zdx1*(1.0d0-zdx2) +
     >       xdaty(iybb)*zdx1*zdx2 +
     >       xdaty(iyab)*(1.0d0-zdx1)*zdx2
!
      if( plg(1,ixs).eq.1 ) vy = exp(vy)
!
      xdt_eval2 = vy
!
      return
      end
!
!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   function xdt_eval1(ixs,vx1)
      real(8) function xdt_eval1(ixs,vx1)
!***********************************************************************
      use cxdcom, only : pin, plg, pxa, pxd, pxi, pxs, xdaty
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/06/16
!ik   real*8 xdt_eval1, vx1
!ik   integer ixs
      real(8), intent(in) :: vx1
      integer, intent(in) :: ixs

!::local variables
      real*8  zex1, zix1, zdx1, vy
      integer ix1a, ix1b, iyaa, iyba
!
      zex1 = vx1
      zex1 = dmax1( zex1, pxi(1,ixs) )
      zex1 = dmin1( zex1, pxa(1,ixs) )
      if( plg(2,ixs).eq.1 ) zex1 = log(zex1)
!
      zix1 = (zex1-pxs(1,ixs))/pxd(1,ixs) + 1.0d0
      ix1a = int(zix1)
      ix1a = min0( ix1a, pin(2,ixs)-1 )
      ix1b = ix1a + 1
      zdx1 = zix1 - dble( ix1a )
!
      iyaa = ix1a + pin(1,ixs)
      iyba = ix1b + pin(1,ixs)
!
      vy  = xdaty(iyaa)*(1.0d0-zdx1) +
     >      xdaty(iyba)*zdx1
!
      if( plg(1,ixs).eq.1 ) vy = exp(vy)
!
      xdt_eval1 = vy
!
      return
      end
!
!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   function xdt_eval0(ixs)
      real(8) function xdt_eval0(ixs)
!***********************************************************************
      use cxdcom, only : pin, plg, xdaty
      implicit none
!
!::argument
! modified 2/1 lines organize local variables and include files by kamata 2021/06/16
!ik   real*8  xdt_eval0
!ik   integer ixs
      integer, intent(in) :: ixs

!::local variables
      integer  iyaa;  real*8 vy
!
      iyaa = 1 + pin(1,ixs)
      vy   = xdaty(iyaa)
!
      if( plg(1,ixs).eq.1 ) vy = exp(vy)
!
      xdt_eval0 = vy
!
      return
      end
