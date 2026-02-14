! added dynamic allocation of arrays by kamata 2022/05/29
! separated from def_PLS_1 'PLS_2'
! sending and receiving variables in cntwfl.f
!**********************************************************************
      subroutine cntwfl_sr( kk, peno )
!**********************************************************************
      use cntwfl,      only : xfhta, xfhti, xfhtm
      use csize,       only : ndmsr, ndwp
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
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 413, ns

! def_PLS_1, strct_comblk
      select case( kk )
      case( 1 )
! ../../soldor/inc/cpmpls  /cpmain/
        ns = ndwp * ( ndmsr + 1 )
        call mpi_bcast( xfhta, ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( xfhti, ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
        call mpi_bcast( xfhtm, ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )

      case( 2 )
! ../../soldor/inc/cpmpls  /cpmain/
        ns = ndwp * ( ndmsr + 1 )
        call mpi_send( xfhta, ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xfhti, ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
        call mpi_send( xfhtm, ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )

      case( 3 )
! ../../soldor/inc/cpmpls  /cpmain/
        ns = ndwp * ( ndmsr + 1 )
        call mpi_recv( xfhta, ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xfhti, ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
        call mpi_recv( xfhtm, ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
      end select

      return
      end
