! added dynamic allocation of arrays by kamata 2022/05/29
! separated from strct_comblk  'CMQNOW'
! bcast variables in cplcom.f
!**********************************************************************
      subroutine cplcom_lsr
!**********************************************************************
      use cplcom, only : q1a, q2a, q3, q4, qtim
      use csize,  only : ndsp, ndx, ndy
      use cunit,  only : lmspe, mywld
      use mpi
      implicit none
! local variables
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 405, ns

! 'CMQNOW'
! ../../soldor/inc/cplcom  /cmqnow/
        ns = ndx * ndy * ndsp
        call mpi_bcast( q1a,   ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( q2a,   ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        ns = ndx * ndy
        call mpi_bcast( q3,    ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( q4,    ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( qtim,  1,    MPI_REAL8,     lmspe
     >                , mywld, merr )

      return
      end
