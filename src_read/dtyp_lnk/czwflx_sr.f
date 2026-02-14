! added treat 4 or more impurities with IMPMC by kamata 2022/04/21
! base def_IMP_2B.f, def_IMP3_2B.f
! sending and receiving variables in czwflx.f
      subroutine czwflx_sr( kk, peno, kis, kir )
      use czwflx, only : zw_catmz, zw_eflx, zw_ismax, zw_pflx
      use csize,  only : ndmis, ndwp
      use mod_mpicomm, only : nwld_cmm, nwld_grp
      use mpi
      implicit none
! arguments
      integer, intent(in) :: kir, kis, kk, peno 
! kir  : impurity number on the receiving side
! kis  : impurity number on the sending side
! kk   : processing flag
!        = 1 : send data to all PEs in the group ( bcast ), kis = kir
!        = 2 : send data to PE = peno, use kis
!        = 3 : receive data from PE = peno, use kir
! peno : send or receive PE number

! local variables
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 401, ns

      select case( kk )
      case( 1 )
! ../../IMPMC_com/inc/czwflx  /czwflx/
        call mpi_bcast( zw_catmz(kis),    2,  MPI_CHARACTER, peno
     >                , nwld_grp, merr )
        ns = ( ndmis + 1 ) * ndwp
        call mpi_bcast( zw_eflx(0,1,kis), ns, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( zw_pflx(0,1,kis), ns, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( zw_ismax(kis),    1,  MPI_INTEGER,   peno
     >                , nwld_grp, merr )
      case( 2 )
! ../../IMPMC_com/inc/czwflx  /czwflx/
        call mpi_send( zw_catmz(kis),    2,  MPI_CHARACTER, peno
     >               , mtag, nwld_cmm, merr )
        ns = ( ndmis + 1 ) * ndwp
        call mpi_send( zw_eflx(0,1,kis), ns, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( zw_pflx(0,1,kis), ns, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( zw_ismax(kis),    1,  MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
      case( 3 )
! ../../IMPMC_com/inc/czwflx  /czwflx/
        call mpi_recv( zw_catmz(kir),    2,  MPI_CHARACTER, peno
     >                , mtag, nwld_cmm, mstat, merr )
        ns = ( ndmis + 1 ) * ndwp
        call mpi_recv( zw_eflx(0,1,kir), ns, MPI_REAL8,     peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( zw_pflx(0,1,kir), ns, MPI_REAL8,     peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( zw_ismax(kir),    1,  MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
      end select

      return
      end
