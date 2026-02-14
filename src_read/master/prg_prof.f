!**********************************************************************
      subroutine prg_prof(kitr,kcnd)
!**********************************************************************
!)
!)      sequential type    2015/07/23
!)
!)   Mnkntl_pls (vna,vti), (flps)  => NEUT2D
!)   Mnkimp_pls (vna,vti), (flps), (q1a), (N0, T0) => IMPMC
!)
!)
!)    order  code     IN-dtyp      => OUT-dtyp
!)   ----------------------------------------------------------
!)      1    soldor                => PLS_2B, PLS_2
!)      2    neut2d   PLS_2B       => NTL_2
!)      3    IMPMC    PLS_2        => IMP_2
!)      4    soldor   NTL_2, IMP_2 => PLS_2B, PLS_2
!)
!::SOLDOR
!)    PLS_2, PLS_2B contain flps calculated in sub. plpflx
!)    plmrgn Not depend on PLS_2(B)
!)    plwtfy depend on q1a
!)    plpflx depend on plwtfy and send to neut2d IMPMC
!)    plmflx tflne tflqi tflqe for output
!)
!----------------------------------------------------------------------
      use cputim,      only : eltm0, eltm00
      use cunit,       only : n6
      use mod_keylist, only : knorm
      use mpi!,         only : mpi_wtime
      implicit none

!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/22
!ik   integer :: kitr, kcnd
      integer, intent(in)  :: kitr ! dummy
      integer, intent(out) :: kcnd

!::local variables
! deleted 4 lines organize local variables and include files by kamata 2021/06/22
!ik   character(12) :: chnum
!ik   character(80) :: cmsg
!ik   real(8) :: tfac
!ik   integer :: ierr, k

      select case(kitr)
!--------------------------------------------------------------------
      case(1)
!--------------------------------------------------------------------

!--------------------------------------------------------------------
      case default
      end select

      eltm0 = MPI_WTIME()
      write(n6,'(2x," kitr = ",i1,"  initial cpu ",f10.3,"  sec.")')
     >   kitr, eltm0-eltm00
      kcnd = knorm
      return
      end
