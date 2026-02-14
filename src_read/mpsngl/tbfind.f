!**********************************************************************
      subroutine tbfind(n,chtb,ino,item)
!**********************************************************************
      implicit none
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31, yamamoto 22/10/07
!ik   integer, intent(inout) :: n
!ik   character(*), dimension(n) :: chtb
      integer, intent(in) :: n
      character(*), dimension(n), intent(in) :: chtb
      character(*), intent(in) :: item
      integer, intent(out) :: ino

      integer :: i

      ino = -1
      do i = 1, n
        if( trim(item) == trim(chtb(i)) ) then
          ino = i
          exit
        endif
      enddo

! added 2 lines display of error occurrence by kamata 2021/08/18, yamamoto 22/10/07
      if( ino < 0 ) print '(a)', 'Error : derived data ' //
     >    trim( item ) // ' not found.'
 
      return
      end subroutine tbfind
