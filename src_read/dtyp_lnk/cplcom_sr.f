! added dynamic allocation of arrays by kamata 2022/05/29
! separated from def_PLS_1 'PLS_1', def_PLS_2 'PLS_2', def_PLS_2B 'PLS_2B'
! sending and receiving variables in cplcom.f
!**********************************************************************
      subroutine cplcom_sr( dtyp, kk, peno )
!**********************************************************************
      use cplcom,      only : anamp, anapv, anasl, anemp, anepv, anesl
     >    , animp, anipv, anisl, atemp, atepv, atesl, atimp, atipv
     >    , atisl, flmxe, flmxi, flps, flvl
     >    , heat_flux_para_by_kappa_para, heat_flux_para_elec, q1a, q2a
     >    , q3, q4, tN0, tNg, tT0, tTg, tV0, vcs, vea, vna, vne, vnezef
     >    , vni, vte, vti, vva, vve, vzf, wdnz, wrad
      use csize,       only : ndgs, ndmc, ndsp, ndwp, ndx, ndy
      use mod_mpicomm, only : nwld_cmm, nwld_grp
      use mpi
      implicit none
! arguments
      character, intent(in) :: dtyp*(*)
      integer,   intent(in) :: kk, peno 
! dtyp : data group name
! kk   : processing flag
!        = 1 : send data to all PEs in the group ( bcast )
!        = 2 : send data to PE = peno
!        = 3 : receive data from PE = peno
! peno : send or receive PE number

! local variables
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 404, ns

      select case( dtyp )
! def_PLS_1
      case( 'PLS_1' )
        select case( kk )
        case( 1 )
! ../../soldor/inc/cplcom  /cmcrad/
          ns = ndx * ndy
          call mpi_bcast( wdnz,  ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( wrad,  ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
        case( 2 )
! ../../soldor/inc/cplcom  /cmcrad/
          ns = ndx * ndy
          call mpi_send( wdnz,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( wrad,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
        case( 3 )
! ../../soldor/inc/cplcom  /cmcrad/
          ns = ndx * ndy
          call mpi_recv( wdnz,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( wrad,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
        end select

! def_PLS_2
      case( 'PLS_2' )
        select case( kk )
        case( 1 )
! ../../soldor/inc/cplcom  /cmauxv/
          ns = ndy * ndsp
          call mpi_bcast( anamp,  ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anapv,  ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anasl,  ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anemp,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anepv,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anesl,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( animp,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anipv,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anisl,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atemp,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atepv,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atesl,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atimp,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atipv,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atisl,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          ns = ndx * ndy
          call mpi_bcast( vcs,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vne,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vnezef, ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vni,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vte,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vti,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vve,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vzf,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          ns = ndx * ndy * ndsp
          call mpi_bcast( vea,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vna,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vva,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )

! ../../soldor/inc/cplcom  /cmflm2/
          ns = ndx * ndy
          call mpi_bcast( flmxe,  ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( flmxi,  ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )

! ../../soldor/inc/cplcom  /cmpflx/
          ns = ndwp * ndsp
          call mpi_bcast( flps,   ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          ns = ndx * ndy * ndsp
          call mpi_bcast( flvl,   ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )

! ../../soldor/inc/cplcom  /cmfth/
          ns = ndx * ndy
          call mpi_bcast( heat_flux_para_by_kappa_para,  ns,   MPI_REAL8
     >                  , peno, nwld_grp, merr )
          call mpi_bcast( heat_flux_para_elec,           ns,   MPI_REAL8
     >                  , peno, nwld_grp, merr )

! ../../soldor/inc/cplcom  /cmqnow/
          ns = ndx * ndy * ndsp
          call mpi_bcast( q1a,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( q2a,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          ns = ndx * ndy
          call mpi_bcast( q3,     ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( q4,     ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )

! ../../soldor/inc/cplcom  /cmntot/
          ns = ( ndmc + 1 ) * ndgs
          call mpi_bcast( tN0,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( tT0,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( tV0,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          ns = ( ndmc + 1 ) * 2
          call mpi_bcast( tNg,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( tTg,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )

        case( 2 )
! ../../soldor/inc/cplcom  /cmauxv/
          ns = ndy * ndsp
          call mpi_send( anamp,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anapv,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anasl,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anemp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anepv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anesl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( animp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anipv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anisl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atemp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atepv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atesl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atimp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atipv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atisl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          ns = ndx * ndy
          call mpi_send( vcs,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vne,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vnezef, ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vni,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vte,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vti,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vve,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vzf,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          ns = ndx * ndy * ndsp
          call mpi_send( vea,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vna,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vva,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )

! ../../soldor/inc/cplcom  /cmflm2/
          ns = ndx * ndy
          call mpi_send( flmxe,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( flmxi,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )

! ../../soldor/inc/cplcom  /cmpflx/
          ns = ndwp * ndsp
          call mpi_send( flps,   ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          ns = ndx * ndy * ndsp
          call mpi_send( flvl,   ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )

! ../../soldor/inc/cplcom  /cmfth/
          ns = ndx * ndy
          call mpi_send( heat_flux_para_by_kappa_para,  ns,   MPI_REAL8
     >                 , peno, mtag, nwld_cmm, merr )
          call mpi_send( heat_flux_para_elec,           ns,   MPI_REAL8
     >                 , peno, mtag, nwld_cmm, merr )

! ../../soldor/inc/cplcom  /cmqnow/
          ns = ndx * ndy * ndsp
          call mpi_send( q1a,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( q2a,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          ns = ndx * ndy
          call mpi_send( q3,     ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( q4,     ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )

! ../../soldor/inc/cplcom  /cmntot/
          ns = ( ndmc + 1 ) * ndgs
          call mpi_send( tN0,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( tT0,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( tV0,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          ns = ( ndmc + 1 ) * 2
          call mpi_send( tNg,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( tTg,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )

        case( 3 )
! ../../soldor/inc/cplcom  /cmauxv/
          ns = ndy * ndsp
          call mpi_recv( anamp,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anapv,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anasl,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anemp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anepv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anesl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( animp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anipv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anisl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atemp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atepv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atesl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atimp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atipv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atisl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          ns = ndx * ndy
          call mpi_recv( vcs,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vne,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vnezef, ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vni,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vte,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vti,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vve,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vzf,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          ns = ndx * ndy * ndsp
          call mpi_recv( vea,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vna,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vva,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )

! ../../soldor/inc/cplcom  /cmflm2/
          ns = ndx * ndy
          call mpi_recv( flmxe,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( flmxi,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )

! ../../soldor/inc/cplcom  /cmpflx/
          ns = ndwp * ndsp
          call mpi_recv( flps,   ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          ns = ndx * ndy * ndsp
          call mpi_recv( flvl,   ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )

! ../../soldor/inc/cplcom  /cmfth/
          ns = ndx * ndy
          call mpi_recv( heat_flux_para_by_kappa_para,  ns,   MPI_REAL8
     >                 , peno, mtag, nwld_cmm, mstat, merr )
          call mpi_recv( heat_flux_para_elec,           ns,   MPI_REAL8
     >                 , peno, mtag, nwld_cmm, mstat, merr )

! ../../soldor/inc/cplcom  /cmqnow/
          ns = ndx * ndy * ndsp
          call mpi_recv( q1a,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( q2a,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          ns = ndx * ndy
          call mpi_recv( q3,     ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( q4,     ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )

! ../../soldor/inc/cplcom  /cmntot/
          ns = ( ndmc + 1 ) * ndgs
          call mpi_recv( tN0,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( tT0,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( tV0,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          ns = ( ndmc + 1 ) * 2
          call mpi_recv( tNg,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( tTg,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
        end select

! def_PLS_2B
      case( 'PLS_2B' )
        select case( kk )
        case( 1 )
! ../../soldor/inc/cplcom  /cmauxv/
          ns = ndy * ndsp
          call mpi_bcast( anamp,  ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anapv,  ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anasl,  ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anemp,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anepv,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anesl,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( animp,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anipv,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( anisl,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atemp,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atepv,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atesl,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atimp,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atipv,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( atisl,  ndy,  MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          ns = ndx * ndy
          call mpi_bcast( vcs,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vne,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vnezef, ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vni,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vte,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vti,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vve,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vzf,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          ns = ndx * ndy * ndsp
          call mpi_bcast( vea,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vna,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( vva,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )

! ../../soldor/inc/cplcom  /cmpflx/
          ns = ndwp * ndsp
          call mpi_bcast( flps,   ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          ns = ndx * ndy * ndsp
          call mpi_bcast( flvl,   ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )

! ../../soldor/inc/cplcom  /cmqnow/
          ns = ndx * ndy * ndsp
          call mpi_bcast( q1a,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( q2a,    ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          ns = ndx * ndy
          call mpi_bcast( q3,     ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )
          call mpi_bcast( q4,     ns,   MPI_REAL8,     peno
     >                  , nwld_grp, merr )

        case( 2 )
! ../../soldor/inc/cplcom  /cmauxv/
          ns = ndy * ndsp
          call mpi_send( anamp,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anapv,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anasl,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anemp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anepv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anesl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( animp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anipv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( anisl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atemp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atepv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atesl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atimp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atipv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( atisl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          ns = ndx * ndy
          call mpi_send( vcs,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vne,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vnezef, ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vni,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vte,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vti,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vve,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vzf,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          ns = ndx * ndy * ndsp
          call mpi_send( vea,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vna,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( vva,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )

! ../../soldor/inc/cplcom  /cmpflx/
          ns = ndwp * ndsp
          call mpi_send( flps,   ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          ns = ndx * ndy * ndsp
          call mpi_send( flvl,   ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )

! ../../soldor/inc/cplcom  /cmqnow/
          ns = ndx * ndy * ndsp
          call mpi_send( q1a,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( q2a,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          ns = ndx * ndy
          call mpi_send( q3,     ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )
          call mpi_send( q4,     ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, merr )

        case( 3 )
! ../../soldor/inc/cplcom  /cmauxv/
          ns = ndy * ndsp
          call mpi_recv( anamp,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anapv,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anasl,  ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anemp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anepv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anesl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( animp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anipv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( anisl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atemp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atepv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atesl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atimp,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atipv,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( atisl,  ndy,  MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          ns = ndx * ndy
          call mpi_recv( vcs,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vne,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vnezef, ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vni,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vte,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vti,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vve,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vzf,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          ns = ndx * ndy * ndsp
          call mpi_recv( vea,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vna,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( vva,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )

! ../../soldor/inc/cplcom  /cmpflx/
          ns = ndwp * ndsp
          call mpi_recv( flps,   ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          ns = ndx * ndy * ndsp
          call mpi_recv( flvl,   ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )

! ../../soldor/inc/cplcom  /cmqnow/
          ns = ndx * ndy * ndsp
          call mpi_recv( q1a,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( q2a,    ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          ns = ndx * ndy
          call mpi_recv( q3,     ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
          call mpi_recv( q4,     ns,   MPI_REAL8,     peno
     >                 , mtag, nwld_cmm, mstat, merr )
        end select
      end select

      return
      end
