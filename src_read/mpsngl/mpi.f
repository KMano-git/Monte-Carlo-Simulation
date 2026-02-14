      module mpi

      integer    MPI_INTEGER, MPI_BYTE
      integer    MPI_CHARACTER
      real(8)    MPI_REAL8, MPI_SUM!, MPI_WTIME

      integer, parameter :: MPI_Address_Kind = 4
      real(4) MPI_REAL4

      integer    MPI_STATUS_SIZE
      parameter (MPI_STATUS_SIZE=1)

      contains
      
      real*8 function MPI_WTIME()
      MPI_WTIME = 0.0d0
      return
      end


      end module mpi
