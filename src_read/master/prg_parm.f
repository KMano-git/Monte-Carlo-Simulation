!**********************************************************************
      subroutine prg_parm(kitr,kcnd)
!**********************************************************************
      use mod_keylist, only : knorm
      implicit none
!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/22
!ik   integer :: kitr, kcnd
      integer, intent(in)  :: kitr ! dummy
      integer, intent(out) :: kcnd

!::local variables
! deleted 2 lines organize local variables and include files by kamata 2021/06/22
!ik   character(12) :: chnum
!ik   character(80) :: cmsg

      select case(kitr)
!--------------------------------------------------------------------
      case(1)
!--------------------------------------------------------------------

!--------------------------------------------------------------------
      case default
!--------------------------------------------------------------------
!::Note  master  multi-execution
      end select

      kcnd = knorm
      return
      end
