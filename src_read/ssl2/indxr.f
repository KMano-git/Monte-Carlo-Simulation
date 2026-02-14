!**********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   function indxr(clin,char)
      integer function indxr(clin,char)
!**********************************************************************
! modified 1/4 lines organize local variables and include files by kamata 2021/05/31
!ik   character clin*(*),char*1
      implicit none
      character, intent(in) :: clin*(*), char*1
! local variables
      integer    i, il, nmax
!
      nmax = len(clin)
      do 110 i = nmax, 1, -1
      il = i
      if( clin(i:i).eq.char ) goto 120
 110  continue
      il = 1
 120  continue
      indxr = il
!
      return
      end
