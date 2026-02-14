!**********************************************************************
      subroutine wf_exit(csub,cmsg)
!**********************************************************************

      implicit none
      character(*), intent(in) :: csub, cmsg
      integer :: errcode, ierr
      integer :: nf

      nf = 6

      write(nf,'(/2x,"*** wf_exit ***")')

      write(nf,'(5x,"stop at sub. ",a)') trim(csub)
      write(nf,'(5x,"because of   ",a)') trim(cmsg)
      write(nf,'(2x)')
      call flush(nf)

      stop
      end

