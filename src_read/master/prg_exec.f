!***********************************************************************
      subroutine prg_exec(kitr,kcnd)
!***********************************************************************
      use csonic,      only : itend, itim, lstop
      use mod_keylist, only : klast, knorm, kstop
      implicit none

!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/22
!ik   integer :: kitr, kcnd
      integer, intent(in)  :: kitr ! dummy
      integer, intent(out) :: kcnd

!::local variables
! deleted 4 lines organize local variables and include files by kamata 2021/06/22
!ik   integer :: ichk, kcon, ifn
!ik   logical :: lex
!ik   real(8) :: tfac, cpmm, fclp, fcdt
!ik   character :: chnum*12, cmsg*80


      select case(kitr)
!--------------------------------------------------------------------
      case(1)
!--------------------------------------------------------------------

!--------------------------------------------------------------------
      case default
      end select

      if( lstop.eq.1 ) goto 200
      if( itim.eq.itend ) goto 200

!::loop
      kcnd = knorm
      return

 200  continue
      kcnd = klast
      if( lstop /= 0 ) kcnd = kstop
      return
      end
