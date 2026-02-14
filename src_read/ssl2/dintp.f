!***********************************************************************
      real(8) function dintp(x,np,xp,yp)
!***********************************************************************
      implicit none
      integer, intent(in) :: np
      real(8), intent(in) :: x, xp(np), yp(np)
!
      integer j, k
      real*8  y, fc
!
      if( x.le.xp(1) ) then
        y = yp(1)
        goto 100
      elseif( x.ge.xp(np) ) then
        y = yp(np)
        goto 100
      endif

      do j=2,np
        k=j
        if( (x-xp(j))*(x-xp(j-1)).le.0.0 ) go to 125
      enddo
      k=np
 125  continue
      fc = (x-xp(k-1))/(xp(k)-xp(k-1))
      y  = yp(k-1)+fc*(yp(k)-yp(k-1))
!
 100  continue
      dintp = y
      return
      end
