!***********************************************************************
      subroutine prg_term(kitr,kcnd)
!***********************************************************************
      use cputim,      only : cputot, eltm0, eltm5
      use cunit,       only : n6
      use dbg_mod,     only : cpaths
      use mod_keylist, only : knend
      use mod_loc_tim, only : msts
      use mpi!,         only : mpi_wtime
      implicit none

!::argument
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
      character(80) :: cmsg

      integer, parameter :: maxitr = 1

      if( kitr > maxitr ) then
        write(cmsg,'("kitr > maitr ",2i3)') kitr, maxitr
        call wf_exit("prg_term(soldor)",trim(cmsg))
      endif

      eltm5 = MPI_WTIME()
      cputot = eltm5-eltm0
!
      call pls_end
      call figdat( cpaths, ' ' )
!
      write(n6,'(/5("=")," normal end prg. soldor")')
      write(n6,'(/2x," total cpu = ", f10.3, " sec.",f6.2," hour")')
     >  cputot,  cputot/3600.0
      write(n6,*)

      kcnd = knend
      msts = kcnd

      return
      end
