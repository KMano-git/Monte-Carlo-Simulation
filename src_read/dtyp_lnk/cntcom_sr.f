! added dynamic allocation of arrays by kamata 2022/05/29
! separated from def_NTL_1  'NTL_1'
! sending and receiving variables in cntcom.f
      subroutine cntcom_sr( kk, peno )
      use cntcom,      only : bvx, bvy, bvz, csps, cswl, e0wl, icps
     >    , icwl, igtw, ihwl, ikps, ikwl, imps, imwl, inwl, ipfiw, ipfmp
     >    , iplx, iply, ipmp, ipps, ipwl, isxp, isxw, isyp, isyw
     >    , mcel, mgrd, migx, migy, mknd, mrgn, mseg, next, pfflx, rcps
     >    , rcwl, snps, snwl, tywl, volm, xpnt, ypnt
     >    , dsgt, iwgt, ipgt, icgt
      use csize,       only : ndmc, ndmg, ndms, ndwp, ndx, ndy, ndgtp
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
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 402, ns

      select case( kk )
      case( 1 )
! ../../monte/inc/cntcom  /cntgrd/
        ns = ndmc + 1
        call mpi_bcast( bvx,   ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( bvy,   ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( bvz,   ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( iplx,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( iply,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        ns = ( ndx + 1 ) * ( ndy + 1 ) 
        call mpi_bcast( mcel,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        ns = ndmc * ( ndms + 1 )
        call mpi_bcast( mgrd,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        ns = ndmc * ndms
        call mpi_bcast( mknd,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( next,  ns,   MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( migx,  ndmc, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( migy,  ndmc, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( mrgn,  ndmc, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( mseg,  ndmc, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( volm,  ndmc, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( xpnt,  ndmg, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( ypnt,  ndmg, MPI_REAL8,     peno
     >                , nwld_grp, merr )

! ../../monte/inc/cntcom  /cntpsf/
        call mpi_bcast( csps,  ndwp, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( icps,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ikps,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( imps,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( imwl,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ipps,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( isxp,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( isyp,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( rcps,  ndwp, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( snps,  ndwp, MPI_REAL8,     peno
     >                , nwld_grp, merr )

! ../../monte/inc/cntcom  /cntpuf/
        call mpi_bcast( ipfiw, ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ipfmp, ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( pfflx, ndwp, MPI_REAL8,     peno
     >                , nwld_grp, merr )

! ../../monte/inc/cntcom  /cntwal/
        call mpi_bcast( cswl,  ndwp, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( e0wl,  ndwp, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( icwl,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( igtw,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ihwl,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ikwl,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( inwl,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ipmp,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ipwl,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( isxw,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( isyw,  ndwp, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( rcwl,  ndwp, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( snwl,  ndwp, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        ns = ndwp * 4
        call mpi_bcast( tywl,  ns,   MPI_CHARACTER, peno
     >                , nwld_grp, merr )

        call mpi_bcast( dsgt,  ndgtp, MPI_REAL8,    peno
     >                , nwld_grp, merr )
        call mpi_bcast( iwgt,  ndgtp, MPI_INTEGER,  peno
     >                , nwld_grp, merr )
        call mpi_bcast( ipgt,  ndgtp, MPI_INTEGER,  peno
     >                , nwld_grp, merr )
        call mpi_bcast( icgt,  ndgtp, MPI_INTEGER,  peno
     >                , nwld_grp, merr )
      case( 2 )
! ../../monte/inc/cntcom  /cntgrd/
        ns = ndmc + 1
        call mpi_send( bvx,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( bvy,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( bvz,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( iplx,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( iply,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        ns = ( ndx + 1 ) * ( ndy + 1 ) 
        call mpi_send( mcel,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        ns = ndmc * ( ndms + 1 )
        call mpi_send( mgrd,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        ns = ndmc * ndms
        call mpi_send( mknd,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( next,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( migx,  ndmc, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( migy,  ndmc, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( mrgn,  ndmc, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( mseg,  ndmc, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( volm,  ndmc, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xpnt,  ndmg, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ypnt,  ndmg, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )

! ../../monte/inc/cntcom  /cntpsf/
        call mpi_send( csps,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( icps,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ikps,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( imps,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( imwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ipps,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( isxp,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( isyp,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( rcps,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( snps,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )

! ../../monte/inc/cntcom  /cntpuf/
        call mpi_send( ipfiw, ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ipfmp, ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( pfflx, ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )

! ../../monte/inc/cntcom  /cntwal/
        call mpi_send( cswl,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( e0wl,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( icwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( igtw,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ihwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ikwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( inwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ipmp,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ipwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( isxw,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( isyw,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( rcwl,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( snwl,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        ns = ndwp * 4
        call mpi_send( tywl,  ns,   MPI_CHARACTER, peno
     >               , mtag, nwld_cmm, merr )

        call mpi_send( dsgt, ndgtp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( iwgt, ndgtp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ipgt, ndgtp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( icgt, ndgtp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
      case( 3 )
! ../../monte/inc/cntcom  /cntgrd/
        ns = ndmc + 1
        call mpi_recv( bvx,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( bvy,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( bvz,   ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( iplx,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( iply,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        ns = ( ndx + 1 ) * ( ndy + 1 ) 
        call mpi_recv( mcel,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        ns = ndmc * ( ndms + 1 )
        call mpi_recv( mgrd,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        ns = ndmc * ndms
        call mpi_recv( mknd,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( next,  ns,   MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( migx,  ndmc, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( migy,  ndmc, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( mrgn,  ndmc, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( mseg,  ndmc, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( volm,  ndmc, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xpnt,  ndmg, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ypnt,  ndmg, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )

! ../../monte/inc/cntcom  /cntpsf/
        call mpi_recv( csps,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( icps,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ikps,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( imps,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( imwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ipps,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( isxp,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( isyp,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( rcps,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( snps,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )

! ../../monte/inc/cntcom  /cntpuf/
        call mpi_recv( ipfiw, ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ipfmp, ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( pfflx, ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )

! ../../monte/inc/cntcom  /cntwal/
        call mpi_recv( cswl,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( e0wl,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( icwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( igtw,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ihwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ikwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( inwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ipmp,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ipwl,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( isxw,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( isyw,  ndwp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( rcwl,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( snwl,  ndwp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        ns = ndwp * 4
        call mpi_recv( tywl,  ns,   MPI_CHARACTER, peno
     >                , mtag, nwld_cmm, mstat, merr )

        call mpi_recv( dsgt,  ndgtp, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( iwgt,  ndgtp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ipgt,  ndgtp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( icgt,  ndgtp, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, mstat, merr )
      end select

      return
      end
