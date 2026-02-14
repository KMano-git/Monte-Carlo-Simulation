!**********************************************************************
      subroutine prg_init(kitr,kcnd)
!**********************************************************************
      use cntmnt,      only : i6_src
      use cputim,      only : eltm00
      use csonic,      only : lpost
      use cunit,       only : n6
      use mod_keylist, only : knorm
      use mpi!,         only : mpi_wtime
      implicit none

!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/06/22
!ik   integer :: kitr, kcnd
!ik   character(12) :: chnum
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
! deleted 1 line organize local variables and include files by kamata 2021/06/22
!ik   character(80) :: cmsg

      select case(kitr)
!--------------------------------------------------------------------
      case(1)
!--------------------------------------------------------------------
      call clr_wftime
      eltm00 = MPI_WTIME()
      lpost = 0             ! <== simulation
      i6_src = 0
      call git_mesg(n6)

!--------------------------------------------------------------------
      case default
!--------------------------------------------------------------------
!::Note  master  multi-execution
      end select

      kcnd = knorm
      return
      end
