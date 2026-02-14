! added dynamic allocation of arrays by kamata 2022/05/29
! separated from def_PLS_1 'PLS_1' and strct_comblk 'CPMAIN' 
! sending and receiving variables in cpmpls.f
!**********************************************************************
      subroutine cpmpls_sr( kk, peno )
!**********************************************************************
      use cpmpls,      only : arhmp, armp, drmp, dvmp, rohmp, romp, sfmp
     >    , vlmp, wdmp, wfmp, wlmp
      use csize,       only : ndy
      use mod_mpicomm, only : nwld_cmm, nwld_grp
      use mpi
      implicit none
! arguments
      integer,   intent(in) :: kk, peno 
! kk   : processing flag
!        = 1 : send data to all PEs in the group ( bcast )
!        = 2 : send data to PE = peno
!        = 3 : receive data from PE = peno
! peno : send or receive PE number

! local variables
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 408

! def_PLS_1, strct_comblk
      select case( kk )
      case( 1 )
! ../../soldor/inc/cpmpls  /cpmain/
        call mpi_bcast( arhmp, ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( armp,  ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( drmp,  ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( dvmp,  ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( rohmp, ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( romp,  ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( sfmp,  ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( vlmp,  ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( wdmp,  ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( wfmp,  ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( wlmp,  ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )

      case( 2 )
! ../../soldor/inc/cpmpls  /cpmain/
        call mpi_send( arhmp, ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( armp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( drmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( dvmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( rohmp, ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( romp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( sfmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( vlmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( wdmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( wfmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( wlmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )

      case( 3 )
! ../../soldor/inc/cpmpls  /cpmain/
        call mpi_recv( arhmp, ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( armp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( drmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( dvmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( rohmp, ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( romp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( sfmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( vlmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( wdmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( wfmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( wlmp,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
      end select

      return
      end
