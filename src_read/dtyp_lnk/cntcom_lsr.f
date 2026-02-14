! added dynamic allocation of arrays by kamata 2022/05/29
! separated from strct_comblk  'CNTGRD', 'CNTWAL'
! bcast variables in cntcom.f
!**********************************************************************
      subroutine cntcom_lsr( kd )
!**********************************************************************
      use cntcom, only : bvx, bvy, bvz, cswl, e0wl, icwl, igtw, ihwl
     >    , ikwl, inwl, iplx, iply, ipmp, ipwl, isxw, isyw, mcel
     >    , mgrd, migx, migy, mknd, mrgn, mseg, next, rcwl, snwl, tywl
     >    , volm, xpnt, ypnt, mclx, mcly
     >    , dsgt, iwgt, ipgt, icgt, nocl
      use csize, only : ndmc, ndmg, ndms, ndwp, ndx, ndy, ndgtp
      use cunit, only : lmspe, mywld
      use mpi
      implicit none
! arguments
      integer, intent(in) :: kd
! kd   : processing flag
!        = 1 : 'CNTGRD', = 2 : 'CNTWAL'

! local variables
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 403, ns

      select case( kd )
      case( 1 ) ! 'CNTGRD'
! ../../monte/inc/cntcom  /cntgrd/
        ns = ndmc + 1
        call mpi_bcast( bvx,   ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( bvy,   ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( bvz,   ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( iplx,  ns,   MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( iply,  ns,   MPI_INTEGER,   lmspe
     >                , mywld, merr )
        ns = ( ndx + 1 ) * ( ndy + 1 )
        call mpi_bcast( mcel,  ns,   MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( nocl,  ns,   MPI_INTEGER,   lmspe
     >                , mywld, merr )
        ns = ndmc * ( ndms + 1 )
        call mpi_bcast( mgrd,  ns,   MPI_INTEGER,   lmspe
     >                , mywld, merr )
        ns = ndmc * ndms
        call mpi_bcast( mknd,  ns,   MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( next,  ns,   MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( migx,  ndmc, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( migy,  ndmc, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( mrgn,  ndmc, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( mseg,  ndmc, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( mclx,  ndx,  MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( mcly,  ndy,  MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( volm,  ndmc, MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( xpnt,  ndmg, MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( ypnt,  ndmg, MPI_REAL8,     lmspe
     >                , mywld, merr )
      case( 2 ) ! 'CNTWAL'
! ../../monte/inc/cntcom  /cntwal/
        call mpi_bcast( cswl,  ndwp, MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( e0wl,  ndwp, MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( icwl,  ndwp, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( igtw,  ndwp, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( ihwl,  ndwp, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( ikwl,  ndwp, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( inwl,  ndwp, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( ipmp,  ndwp, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( ipwl,  ndwp, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( isxw,  ndwp, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( isyw,  ndwp, MPI_INTEGER,   lmspe
     >                , mywld, merr )
        call mpi_bcast( rcwl,  ndwp, MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( snwl,  ndwp, MPI_REAL8,     lmspe
     >                , mywld, merr )
        ns = ndwp * 4
        call mpi_bcast( tywl,  ns,   MPI_CHARACTER, lmspe
     >                , mywld, merr )

        call mpi_bcast( dsgt,  ndgtp, MPI_REAL8,    lmspe
     >                , mywld, merr )
        call mpi_bcast( iwgt,  ndgtp, MPI_INTEGER,  lmspe
     >                , mywld, merr )
        call mpi_bcast( ipgt,  ndgtp, MPI_INTEGER,  lmspe
     >                , mywld, merr )
        call mpi_bcast( icgt,  ndgtp, MPI_INTEGER,  lmspe
     >                , mywld, merr )
      end select

      return
      end
