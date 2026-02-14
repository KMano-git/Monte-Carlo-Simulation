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
! modified 1/2 lines organize local variables and include files by kamata 2021/06/28
!ik   integer :: kitr, kcnd
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   character(12) :: chnum
      character(80) :: cmsg

      select case(kitr)
!--------------------------------------------------------------------
      case(1)
!--------------------------------------------------------------------
      eltm00 = MPI_WTIME()

      lpost = 0             ! <== simulation
      i6_src = 0

!::output files
!  (master:0,soldor:1,neut2d:2,impmc:3)
      call git_mesg(n6)

!::Read option flag & open fort.xx (Master PE)
      call stjob
      call opinpt
      call opnfrt

!::set physical const in all PE
      call phcnst

!::glb_tim/loc_tim
      call clr_wftime

      kcnd = knorm
      return

!--------------------------------------------------------------------
      case default
!--------------------------------------------------------------------
      write(cmsg,'("kitr =",i2," > 1")') kitr
      call wexit("prg_init(soldor)",trim(cmsg))
      end select
      end
