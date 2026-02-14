! added dynamic allocation of arrays by kamata 2022/05/29
! moved from strct_comblk 'CPLMDPRF'
! bcast variables in cpmpls.f
!**********************************************************************
      subroutine cpmpls_lsr
!**********************************************************************
      use cpmpls, only : mdp_ni, mdp_rh, mdp_ro, mdp_te, mdp_ti
      use csize,  only : ndy
      use cunit,  only : lmspe, mywld
      use mpi
      implicit none
! local variables
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 409

! 'CPLMDPRF'
! ../../soldor/inc/cpmpls  /com_plmdprf/
      call mpi_bcast( mdp_ni,  ndy,  MPI_REAL8,     lmspe
     >              , mywld, merr )
      call mpi_bcast( mdp_rh,  ndy,  MPI_REAL8,     lmspe
     >              , mywld, merr )
      call mpi_bcast( mdp_ro,  ndy,  MPI_REAL8,     lmspe
     >              , mywld, merr )
      call mpi_bcast( mdp_te,  ndy,  MPI_REAL8,     lmspe
     >              , mywld, merr )
      call mpi_bcast( mdp_ti,  ndy,  MPI_REAL8,     lmspe
     >              , mywld, merr )

      return
      end
