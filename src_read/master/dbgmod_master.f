!**********************************************************************
      subroutine figdat( cpath, ctime )
!**********************************************************************
! this is just dummy function for build
      implicit none
! arguments
      character, intent(in) :: cpath*(*), ctime*(*)
!
      return
      end

!**********************************************************************
      subroutine pst_qdpl( cpath, ctime )
!**********************************************************************
! this is just dummy function for build
      implicit none
! arguments
      character, intent(in) :: cpath*(*), ctime*(*)
!
      return
      end


!**********************************************************************
      subroutine timcntl( timtbl, dtmtbl, ntim, time, timb, judgement )
!**********************************************************************
! this is just dummy function for build
      implicit none
! arguments
      integer, intent(in)  :: ntim
      integer, intent(in) :: judgement
      real(8), intent(in)  :: dtmtbl(ntim), timb, time, timtbl(ntim)
!
      return
      end


!**********************************************************************
      subroutine rn2a( rnum, ndec, lntxt, atxt )
!**********************************************************************
! this is just dummy function for build
      implicit none
! arguments
      real(8),   intent(in)  :: rnum
      integer,   intent(in)  :: ndec
      integer,   intent(in)  :: lntxt
      character, intent(in) :: atxt*(lntxt)
!
      return
      end