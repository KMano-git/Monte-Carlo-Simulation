!**********************************************************************
      subroutine prg_parm(kitr,kcnd)
!**********************************************************************
      use csonic,      only : lrand
      use mod_keylist, only : knorm
      implicit none

!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/16
!ik   integer :: kcnd, kitr
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   character(12) :: chnum
      character(80) :: cmsg

      select case(kitr)
!--------------------------------------------------------------------
      case(1)
!--------------------------------------------------------------------
      call ranini(lrand)
      call ntl_ini

!--------------------------------------------------------------------
      case default
!--------------------------------------------------------------------
      write(cmsg,'("kitr =",i2," > 1")') kitr
      call wexit("prg_init(monte)",trim(cmsg))
      end select

      kcnd = knorm
      return
      end
