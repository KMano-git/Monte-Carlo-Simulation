! added dynamic allocation of arrays by kamata 2022/05/29
! separated from def_PLS_1 'PLS_1' and strct_comblk 'CMETRC' 
! sending and receiving variables in cplmet.f
!**********************************************************************
      subroutine cplmet_sr( kk, peno )
!**********************************************************************
      use cplmet,      only : gare, gdsv, gwtm, gwtp, hare, hdsp, hdsv
     >    , hdxm, hdxp, hgdx, hgdy, hpit, hvol, hvsb, hwtm, hwtp, icel
     >    , icmax, icmin, jcel, jnxm, jnxp, jtmax, jtmin, kgdx, kgdy
     >    , kreg, romn, vlmn
      use csize,       only : ndx, ndy
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
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 406, ns

! def_PLS_1, strct_comblk
! ../../soldor/inc/cplcom  /cmetrc/
      select case( kk )
      case( 1 )
        ns = ndx * ndy
        call mpi_bcast( gare,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( gwtm,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( gwtp,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hare,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hdxm,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hdxp,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hgdx,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hgdy,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hpit,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hvol,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hvsb,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hwtm,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hwtp,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( icel,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( jcel,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( jnxm,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( jnxp,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( kreg,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )

        ns = ndx * ndy * 4
        call mpi_bcast( gdsv,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hdsp,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( hdsv,  ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( kgdx,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( kgdy,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )

        call mpi_bcast( icmax, ndx,  MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( icmin, ndx,  MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( vlmn,  ndx,  MPI_REAL8,     peno
     >                , nwld_grp, merr )

        call mpi_bcast( jtmax, ndy,  MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( jtmin, ndy,  MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( romn,  ndy,  MPI_REAL8,     peno
     >                , nwld_grp, merr )

      case( 2 )
        ns = ndx * ndy
        call mpi_send( gare,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( gwtm,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( gwtp,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hare,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hdxm,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hdxp,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hgdx,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hgdy,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hpit,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hvol,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hvsb,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hwtm,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hwtp,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( icel,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( jcel,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( jnxm,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( jnxp,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( kreg,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )

        ns = ndx * ndy * 4
        call mpi_send( gdsv,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hdsp,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( hdsv,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( kgdx,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( kgdy,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )

        call mpi_send( icmax, ndx,  MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( icmin, ndx,  MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( vlmn,  ndx,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )

        call mpi_send( jtmax, ndy,  MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( jtmin, ndy,  MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( romn,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )

      case( 3 )
        ns = ndx * ndy
        call mpi_recv( gare,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( gwtm,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( gwtp,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hare,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hdxm,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hdxp,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hgdx,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hgdy,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hpit,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hvol,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hvsb,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hwtm,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hwtp,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( icel,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( jcel,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( jnxm,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( jnxp,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( kreg,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )

        ns = ndx * ndy * 4
        call mpi_recv( gdsv,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hdsp,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( hdsv,  ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( kgdx,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( kgdy,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )

        call mpi_recv( icmax, ndx,  MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( icmin, ndx,  MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( vlmn,  ndx,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )

        call mpi_recv( jtmax, ndy,  MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( jtmin, ndy,  MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( romn,  ndy,  MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
      end select

      return
      end
