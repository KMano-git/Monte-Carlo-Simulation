! added dynamic allocation of arrays by kamata 2022/05/29
! separated from def_PLS_1 'PLS_1' and strct_comblk 'CGMSH' 
! sending and receiving variables in cgdcom.f
!**********************************************************************
      subroutine cgdcom_sr( kk, peno )
!**********************************************************************
      use cgdcom,      only : grdx, grdy, hbr, hbt, hbz, psman, psprv
     >    , pssol
      use csize,       only : ndq => ndx, ndp => ndy
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
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 414, ns

! def_PLS_1, strct_comblk
! ../../sonic/inc/size_XX/cgdcom  /cgmsh/
      select case( kk )
      case( 1 )
        ns = ndq * ndp
        call mpi_bcast( grdx,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( grdy,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hbr,   ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hbt,   ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hbz,   ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( psman, ndp,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( psprv, ndp,  MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( pssol, ndp,  MPI_REAL8,     peno
     >                , nwld_grp, merr )

      case( 2 )
        ns = ndq * ndp
        call mpi_send( grdx,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( grdy,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hbr,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hbt,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hbz,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( psman, ndp,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( psprv, ndp,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( pssol, ndp,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )

      case( 3 )
        ns = ndq * ndp
        call mpi_recv( grdx,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( grdy,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hbr,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hbt,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hbz,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( psman, ndp,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( psprv, ndp,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( pssol, ndp,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
      end select

      return
      end
