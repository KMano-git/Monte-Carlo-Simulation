! added dynamic allocation of arrays by kamata 2022/05/29
! separated from def_IMP_2  'IMPx_2'
! sending and receiving variables in cimden.f
      subroutine cimden_sr( kk, peno )
      use cimden,      only : tdnz, tfrz, tionZ, tradiZ, tradliZ, tradrZ
     >    , trecZ, tthz, tvlz, twci, twne, twrd
      use csize,       only : ndmc, ndmis
      use mod_mpicomm, only : nwld_cmm, nwld_grp
      use mpi
      implicit none
! arguments
      integer, intent(in) :: kk, peno 
! kk   : processing flag
!        = 1 : send data to all PEs in the group ( bcast )
!        = 2 : send data to PE = peno
!        = 3 : receive data from PE = peno
! peno : send or receive PE number

! local variables
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 412, ns

      select case( kk )
      case( 1 )
! ../../IMPMC/inc/cimden  /cimden/
        call mpi_bcast( twci,    ndmc, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( twne,    ndmc, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( twrd,    ndmc, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        ns = ( ndmis + 1 ) * ndmc
        call mpi_bcast( tdnz,    ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
! ../../IMPMC/inc/cimden  /cimden_2/
        call mpi_bcast( tfrz,    ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( tionZ,   ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( trecZ,   ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( tthz,    ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( tvlz,    ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
! ../../IMPMC/inc/cimden  /cimden_3/
        call mpi_bcast( tradiZ,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( tradliZ, ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( tradrZ,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
      case( 2 )
! ../../IMPMC/inc/cimden  /cimden/
        call mpi_send( twci,    ndmc, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( twne,    ndmc, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( twrd,    ndmc, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        ns = ( ndmis + 1 ) * ndmc
        call mpi_send( tdnz,    ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
! ../../IMPMC/inc/cimden  /cimden_2/
        call mpi_send( tfrz,    ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( tionZ,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( trecZ,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( tthz,    ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( tvlz,    ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
! ../../IMPMC/inc/cimden  /cimden_3/
        call mpi_send( tradiZ,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( tradliZ, ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( tradrZ,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
      case( 3 )
! ../../IMPMC/inc/cimden  /cimden/
        call mpi_recv( twci,    ndmc, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( twne,    ndmc, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( twrd,    ndmc, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        ns = ( ndmis + 1 ) * ndmc
        call mpi_recv( tdnz,    ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
! ../../IMPMC/inc/cimden  /cimden_2/
        call mpi_recv( tfrz,    ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( tionZ,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( trecZ,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( tthz,    ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( tvlz,    ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
! ../../IMPMC/inc/cimden  /cimden_3/
        call mpi_recv( tradiZ,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( tradliZ, ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( tradrZ,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
      end select

      return
      end
