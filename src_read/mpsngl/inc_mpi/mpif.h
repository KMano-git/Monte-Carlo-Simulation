!----------------------------------------------------------------------
!::dummy    mpsngl/inc_mpi/mpif.h
!----------------------------------------------------------------------
!::MPI type definitions
      integer  MPI_BYTE, MPI_CHARACTER, MPI_INTEGER, MPI_REAL8
      integer  MPI_SUM, MPI_STATUS_SIZE
      parameter ( MPI_STATUS_SIZE = 8 )
      real*8       MPI_WTIME
!
      integer  MPI_MODE_CREATE, MPI_MODE_RDWR, MPI_INFO_NULL
      integer  MPI_SEEK_SET
      integer  MPI_COMM_WORLD
