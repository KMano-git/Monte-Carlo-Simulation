! added treat 4 or more impurities with IMPMC by kamata 2022/04/15
! base def_IMP_1.f, def_IMP2_1.f, def_IMP3_1.f
! sending and receiving variables in catcom.f
      subroutine catcom_sr( kk, peno, kis, kir )
      use catcom, only : aionz, catmz, ip_cxr, ip_ion, ip_plt
     >    , ip_prb, ip_prc, ip_rad, ip_rec, natmz, ndxdt, ndxkn, ndxrw
     >    , ndxx1, ndxx2, nxs, nyp
     >    , xdaty, xdatz, xdevl, xdnam, xdrnk, xdsn, xdtx1, xdtx2, xjend
     >    , xjnum, xjsta, xwmax, xwspc, xwunt, xwvar, xwmin, xwmlt
     >    , xwnum
      use csize,  only : ndmis
      use mod_mpicomm, only : nwld_cmm, nwld_grp
      use mpi
      implicit none
! arguments
      integer, intent(in) :: kir, kis, kk, peno 
! kir  : impurity number on the receiving side
! kis  : impurity number on the sending side
! kk   : processing flag
!        = 1 : send data to all PEs in the group ( bcast ), kis = kir = 1
!        = 2 : send data to PE = peno, use kis
!        = 3 : receive data from PE = peno, use kir
! peno : send or receive PE number

! local variables
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 400, ns

      select case( kk )
      case( 1 )
! ../../ATadas/inc/catcom  /catcom/
        call mpi_bcast( xdaty(1,kis),   ndxdt, MPI_REAL8,     peno
     >                , nwld_grp, merr )
        ns = ndxx1 * ndxkn
        call mpi_bcast( xdtx1(1,1,kis), ns,    MPI_REAL8,     peno
     >                , nwld_grp, merr )
        ns = ndxx2 * ndxkn
        call mpi_bcast( xdtx2(1,1,kis), ns,    MPI_REAL8,     peno
     >                , nwld_grp, merr )
        ns = ndmis * ndxkn
        call mpi_bcast( xdatz(1,1,kis), ns,    MPI_REAL8,     peno
     >                , nwld_grp, merr )
        ns = ndxrw * ndxkn
        call mpi_bcast( xwmin(1,1,kis), ns,    MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( xwmax(1,1,kis), ns,    MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( xwmlt(1,1,kis), ns,    MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( ip_ion(kis),    1,     MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ip_rec(kis),    1,     MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ip_cxr(kis),    1,     MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ip_rad(kis),    1,     MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ip_plt(kis),    1,     MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ip_prb(kis),    1,     MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( ip_prc(kis),    1,     MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( nxs(kis),       1,     MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( nyp(kis),       1,     MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( xdrnk(1,kis),   ndxkn, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        ns = ndxrw * ndxkn
        call mpi_bcast( xwnum(1,1,kis), ns,    MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( xjsta(1,kis),   ndxkn, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( xjend(1,kis),   ndxkn, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( xjnum(1,kis),   ndxkn, MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        ns = ndxkn * 80
        call mpi_bcast( xdsn(1,kis),    ns,    MPI_CHARACTER, peno
     >                , nwld_grp, merr )
        ns = ndxkn * 20
        call mpi_bcast( xdnam(1,kis),   ns,    MPI_CHARACTER, peno
     >                , nwld_grp, merr )
        call mpi_bcast( xdevl(1,kis),   ns,    MPI_CHARACTER, peno
     >                , nwld_grp, merr )
        ns = ndxrw * ndxkn * 12
        call mpi_bcast( xwvar(1,1,kis), ns,    MPI_CHARACTER, peno
     >                , nwld_grp, merr )
        call mpi_bcast( xwunt(1,1,kis), ns,    MPI_CHARACTER, peno
     >                , nwld_grp, merr )
        call mpi_bcast( xwspc(1,1,kis), ns,    MPI_CHARACTER, peno
     >                , nwld_grp, merr )
! ../../ATadas/inc/catcom  /catcom3/
        call mpi_bcast( aionz(kis),     1,     MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( natmz(kis),     1,     MPI_INTEGER,   peno
     >                , nwld_grp, merr )
        call mpi_bcast( catmz(kis),     2,     MPI_CHARACTER, peno
     >                , nwld_grp, merr )
      case( 2 )
! ../../ATadas/inc/catcom  /catcom/
        call mpi_send( xdaty(1,kis),   ndxdt, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        ns = ndxx1 * ndxkn
        call mpi_send( xdtx1(1,1,kis), ns,    MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        ns = ndxx2 * ndxkn
        call mpi_send( xdtx2(1,1,kis), ns,    MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        ns = ndmis * ndxkn
        call mpi_send( xdatz(1,1,kis), ns,    MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        ns = ndxrw * ndxkn
        call mpi_send( xwmin(1,1,kis), ns,    MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xwmax(1,1,kis), ns,    MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xwmlt(1,1,kis), ns,    MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ip_ion(kis),    1,     MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ip_rec(kis),    1,     MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ip_cxr(kis),    1,     MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ip_rad(kis),    1,     MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ip_plt(kis),    1,     MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ip_prb(kis),    1,     MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( ip_prc(kis),    1,     MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( nxs(kis),       1,     MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( nyp(kis),       1,     MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xdrnk(1,kis),   ndxkn, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        ns = ndxrw * ndxkn
        call mpi_send( xwnum(1,1,kis), ns,    MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xjsta(1,kis),   ndxkn, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xjend(1,kis),   ndxkn, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xjnum(1,kis),   ndxkn, MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        ns = ndxkn * 80
        call mpi_send( xdsn(1,kis),    ns,    MPI_CHARACTER, peno
     >               , mtag, nwld_cmm, merr )
        ns = ndxkn * 20
        call mpi_send( xdnam(1,kis),   ns,    MPI_CHARACTER, peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xdevl(1,kis),   ns,    MPI_CHARACTER, peno
     >               , mtag, nwld_cmm, merr )
        ns = ndxrw * ndxkn * 12
        call mpi_send( xwvar(1,1,kis), ns,    MPI_CHARACTER, peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xwunt(1,1,kis), ns,    MPI_CHARACTER, peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xwspc(1,1,kis), ns,    MPI_CHARACTER, peno
     >               , mtag, nwld_cmm, merr )
! ../../ATadas/inc/catcom  /catcom3/
        call mpi_send( aionz(kis),     1,     MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( natmz(kis),     1,     MPI_INTEGER,   peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( catmz(kis),     2,     MPI_CHARACTER, peno
     >               , mtag, nwld_cmm, merr )
      case( 3 )
! ../../ATadas/inc/catcom  /catcom/
        call mpi_recv( xdaty(1,kir),   ndxdt, MPI_REAL8,     peno
     >                , mtag, nwld_cmm, mstat, merr )
        ns = ndxx1 * ndxkn
        call mpi_recv( xdtx1(1,1,kir), ns,    MPI_REAL8,     peno
     >                , mtag, nwld_cmm, mstat, merr )
        ns = ndxx2 * ndxkn
        call mpi_recv( xdtx2(1,1,kir), ns,    MPI_REAL8,     peno
     >                , mtag, nwld_cmm, mstat, merr )
        ns = ndmis * ndxkn
        call mpi_recv( xdatz(1,1,kir), ns,    MPI_REAL8,     peno
     >                , mtag, nwld_cmm, mstat, merr )
        ns = ndxrw * ndxkn
        call mpi_recv( xwmin(1,1,kir), ns,    MPI_REAL8,     peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xwmax(1,1,kir), ns,    MPI_REAL8,     peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xwmlt(1,1,kir), ns,    MPI_REAL8,     peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ip_ion(kir),    1,     MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ip_rec(kir),    1,     MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ip_cxr(kir),    1,     MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ip_rad(kir),    1,     MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ip_plt(kir),    1,     MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ip_prb(kir),    1,     MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( ip_prc(kir),    1,     MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( nxs(kir),       1,     MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( nyp(kir),       1,     MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xdrnk(1,kir),   ndxkn, MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        ns = ndxrw * ndxkn
        call mpi_recv( xwnum(1,1,kir), ns,    MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xjsta(1,kir),   ndxkn, MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xjend(1,kir),   ndxkn, MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xjnum(1,kir),   ndxkn, MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        ns = ndxkn * 80
        call mpi_recv( xdsn(1,kir),    ns,    MPI_CHARACTER, peno
     >                , mtag, nwld_cmm, mstat, merr )
        ns = ndxkn * 20
        call mpi_recv( xdnam(1,kir),   ns,    MPI_CHARACTER, peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xdevl(1,kir),   ns,    MPI_CHARACTER, peno
     >                , mtag, nwld_cmm, mstat, merr )
        ns = ndxrw * ndxkn * 12
        call mpi_recv( xwvar(1,1,kir), ns,    MPI_CHARACTER, peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xwunt(1,1,kir), ns,    MPI_CHARACTER, peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xwspc(1,1,kir), ns,    MPI_CHARACTER, peno
     >                , mtag, nwld_cmm, mstat, merr )
! ../../ATadas/inc/catcom  /catcom3/
        call mpi_recv( aionz(kir),     1,     MPI_REAL8,     peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( natmz(kir),     1,     MPI_INTEGER,   peno
     >                , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( catmz(kir),     2,     MPI_CHARACTER, peno
     >                , mtag, nwld_cmm, mstat, merr )
      end select

      return
      end
