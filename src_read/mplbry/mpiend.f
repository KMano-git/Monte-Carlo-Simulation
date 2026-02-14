!**********************************************************************
      subroutine mpiend
!**********************************************************************
      use cunit, only : mype, n6
      use mpi!,   only : mpi_finalize
      implicit none

      integer ierr

      write(n6,'(//2x,"=== mpiend ===   normal end  mype =",i5)') mype

      call MPI_Finalize( ierr )

      return
      end
