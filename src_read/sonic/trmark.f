!***********************************************************************
      subroutine  trmark(csub,cmsg)
!***********************************************************************
      use cunit, only : n6
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   character*(*) :: csub, cmsg
      character*(*), intent(in) :: csub, cmsg

      return

      if( csub(1:1) == "/" ) then
        write(n6,'(2x)')
        write(n6,'(a,a,a,a)') "*** ",trim(csub(2:)),"  ",trim(cmsg)
      else
        write(n6,'(a,a,a,a)') "*** ",trim(csub),"  ",trim(cmsg)
      endif
      call flush(n6)

      return
      end
