!***********************************************************************
      subroutine seti(ia,n,ic)
!***********************************************************************
      implicit none
! modified 1/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer n, ia(n), ic, i
      integer, intent(in)  :: n, ic
      integer, intent(out) :: ia(n)

! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   do 10 i = 1, n
!ik   ia(i) = ic
!ik 10   continue
      ia(1:n) = ic

      return
      end
