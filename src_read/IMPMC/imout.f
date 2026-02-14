!***********************************************************************
      subroutine imout_ipos(ip,ptmx,ptim)
!***********************************************************************
      use cimcom, only : i6_ipos
      use cimntl, only : stb
      use cunit,  only : lmspe, lmype
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/06/28
!ik   integer    ip
!ik   real*8     ptmx, ptim
      integer, intent(in) :: ip
      real(8), intent(in) :: ptmx, ptim
!
!::local variables
!
      if( lmype.ne.lmspe ) return
      if( i6_ipos.le.0 )   return
!
      call impdt(i6_ipos,ip,"nionS",ptmx,ptim,stb(ip))
!
      return
      end
