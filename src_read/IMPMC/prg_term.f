!***********************************************************************
      subroutine prg_term(kitr,kcnd)
!***********************************************************************
      use cputim,      only : cputot, eltm0, eltm5
      use cunit,       only : n6
      use dbg_mod,     only : cpathi
      use mod_keylist, only : knend
      use mod_loc_tim, only : msts
      use mod_shexe,   only : impmc_model
      use mpi!,         only : mpi_wtime
      implicit none

!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/07/04
!ik   integer :: kitr, kcnd
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
! deleted 1 line organize local variables and include files by kamata 2021/07/04
!ik   character(12) :: chnum
      character(80) :: cmsg

      integer, parameter :: maxitr = 1

      if( kitr > maxitr ) then
        write(cmsg,'("kitr > maitr ",2i3)') kitr, maxitr
        call wf_exit("prg_term(IMPMC)",trim(cmsg))
      endif

      eltm5 = MPI_WTIME()
      cputot = eltm5-eltm0
! deleted 1 line integration of ST and TD versionis of IMPMC by kamata 2022/05/14
!ik   write(n6,*) 'eltm5, eltm0', eltm5, eltm0
!
! added 3 lines integration of ST and TD versionis of IMPMC by kamata 2022/05/14
      if( impmc_model == 0 ) then
        call imp_end
      else
        call imp_wrd
      endif
      call figdat( cpathi, ' ' ) !DBG

      call pst_wimp
!
 250  continue
! added 1 line integration of ST and TD versionis of IMPMC by kamata 2022/05/14
      if( impmc_model == 1 ) then
        eltm5 = MPI_WTIME()
        cputot = eltm5-eltm0
! added 1 line integration of ST and TD versionis of IMPMC by kamata 2022/05/14
      endif
      write(n6,'(/2x," total cpu = ", f10.3, " sec.",f6.2," hour")')
     >  cputot,  cputot/3600.0

      kcnd = knend
      msts = kcnd

      return
      end
