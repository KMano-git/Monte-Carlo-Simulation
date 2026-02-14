!***********************************************************************
      subroutine setd(a,n,c)
!***********************************************************************
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: c
      real*8, intent(out) :: a(n)
      a = c
      return
      end
