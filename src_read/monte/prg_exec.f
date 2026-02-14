!***********************************************************************
      subroutine prg_exec(kitr,kcnd)
!***********************************************************************
      use csonic,      only : itend, itim, lstop
      use cunit,       only : n6
      use dbg_mod,     only : dbgcntl
      use mod_keylist, only : klast, knorm, kstop
      implicit none

!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/16
!ik   integer :: kitr, kcnd
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
! modified 4/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer :: ichk, kcon, ifn
!ik   logical :: lex
!ik   real(8) :: tfac, cpmm, fclp, fcdt
!ik   character :: chnum*12, cmsg*80
      character :: cmsg*80

      select case(kitr)
!--------------------------------------------------------------------
      case(1)
!--------------------------------------------------------------------
      call ntl_recd
      call ntl_pls
      call ntl_cal

! added 1 line periodic invocation of IMPMC/figdat by kamata 2023/08/04
      call dbgcntl( 2 )

      if( lstop.eq.1 ) goto 200
      if( mod(itim,100).eq.0 ) call flush(n6)
      if( itim.eq.itend ) goto 200  ! <== need to check

!::loop
      kcnd = knorm
      return

!::itim
 200  continue
      kcnd = klast
      if( lstop /= 0 ) kcnd = kstop
      return

!--------------------------------------------------------------------
      case default
!--------------------------------------------------------------------
      write(cmsg,'("kitr =",i2," > 1")') kitr
      call wexit("prg_exec(monte)",trim(cmsg))
      end select
      end
