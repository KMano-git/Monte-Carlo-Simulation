!**********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!k    function lenx(clin)
      integer function lenx(clin)
!**********************************************************************
! modified 1/4 lines organize local variables and include files by kamata 2021/05/31
!ik   character clin*(*)
      implicit none
      character, intent(in) :: clin*(*)
! local variables
      integer    i, il, nmax
!
      nmax = len(clin)
      do 110 i = nmax, 1, -1
      il = i
      if( clin(i:i).ne." " ) goto 120
 110  continue
      il = 0
 120  continue
      lenx = il
!
      return
      end
