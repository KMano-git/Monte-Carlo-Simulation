!***********************************************************************
      subroutine prg_term(kitr,kcnd)
!***********************************************************************
      use cputim,      only : cputot, eltm0, eltm5
      use cunit,       only : n6
      use mod_keylist, only : knend
      use mod_loc_tim, only : msts
      use mpi!,         only : mpi_wtime
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

! deleted 1 line organize local variables and include files by kamata 2021/06/22
!ik   integer, parameter :: maxitr = 1

!xxx      if( kitr > maxitr ) then
!xxx        write(cmsg,'("kitr > maitr ",2i3)') kitr, maxitr
!xxx        call wf_exit("prg_term(master)",trim(cmsg))
!xxx      endif

      eltm5 = MPI_WTIME()
      cputot = eltm5-eltm0
!
 250  continue
      write(n6,'(/2x," total cpu = ", f10.3, " sec.",f6.2," hour")')
     >  cputot,  cputot/3600.0

      kcnd = knend  ! normal end
      msts = kcnd

      return
      end
