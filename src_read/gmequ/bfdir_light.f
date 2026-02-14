!***********************************************************************
      subroutine bfdir_light(bdir)
!***********************************************************************
      use cgdcom, only : grdx, grdy
      use cplmet, only : icel, itsle, jcel
      implicit none
!
!::argument
      real(8), intent(out) :: bdir
!
!::local variables
      integer  jt, j, i
      integer  mx1, mx2, mx3, mx4, my1, my2, my3, my4
      real*8   x1, y1, x2, y2, zln, esx, esy, zlb, ebx, eby
      real*8   bx, by, bz
!
      jt = 2
      j  = jcel(jt,itsle)
      i  = icel(jt,itsle)
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
!
      x1 = grdx(mx1,my1)
      y1 = grdy(mx1,my1)
      x2 = grdx(mx2,my2)
      y2 = grdy(mx2,my2)
      zln = dsqrt( (x2-x1)**2+(y2-y1)**2 )
      esx = (x2-x1)/zln
      esy = (y2-y1)/zln
!
      call mfdir(x1,y1,bx,by,bz)
      zlb = dsqrt( bx**2+by**2 )
      ebx = bx/zlb
      eby = by/zlb
!
      bdir = esx*ebx+esy*eby
!
      return
      end