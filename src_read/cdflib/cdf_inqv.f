!**********************************************************************
      subroutine cdf_inqv( cvar, ivar )
!**********************************************************************
      use cdfcom, only : mxvr, nxvr, xvar
      implicit none
!
!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/16
!ik   character cvar*(*); integer ivar
      character, intent(in)  :: cvar*(*)
      integer,   intent(out) :: ivar

!::local variables
! modified 1/3 lines organize local variables and include files by kamata 2021/06/16
!ik   integer mj, i, lenx
      integer mj, i
! function
      integer    lenx
!
      mj = lenx(cvar)
!
      ivar = 0
      do i = 1, nxvr
      if( xvar(i)(1:mxvr(i)).eq.cvar(1:mj) ) then
        ivar = i
        exit
      endif
      enddo
!
!::error massage
      if( ivar.le.0 ) then
      write(6,'(/4x,"cdf_inqv  no found cvar ",a)') cvar(1:lenx(cvar))
      endif
!
      return
      end
