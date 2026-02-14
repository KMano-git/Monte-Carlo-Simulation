!***********************************************************************
      subroutine prg_term(kitr,kcnd)
!***********************************************************************
      use cputim,      only : cputot, eltm0, eltm5
      use cunit,       only : n6
      use dbg_mod,     only : cpathn
      use mod_keylist, only : knend
      use mod_loc_tim, only : msts
      use mpi!,         only : mpi_wtime
      use ntpfctl,     only : out_pfctl
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

      integer, parameter :: maxitr = 1

      if( kitr > maxitr ) then
        write(cmsg,'("kitr > maitr ",2i3)') kitr, maxitr
        call wf_exit("prg_term(neut2d)",trim(cmsg))
      endif

      eltm5 = MPI_WTIME()
      cputot = eltm5-eltm0
!
      call ntl_end
      call figdat( cpathn, ' ' )

      call out_pfctl
!
 250  continue
      write(n6,'(/2x," total cpu = ", f10.3, " sec. ",f6.2," hour")')
     >  cputot,  cputot/3600.0

      kcnd = knend
      msts = kcnd

      return
      end
