! added dynamic allocation of arrays by kamata 2022/05/29
! separated from def_NTL_2  'NTL_2'
! sending and receiving variables in cntmnt.f
      subroutine cntmnt_sr( kk, peno )
      use cntmnt,      only : vmwork
      use csize,       only : ndmsr, nwkmp_nt
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
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 407, ns

      select case( kk )
      case( 1 )
! ../../monte/inc/cntmnt  /cntmnt/
        ns = nwkmp_nt * ndmsr
        call mpi_bcast( vmwork, ns,   MPI_REAL8,     peno
     >                , nwld_grp, merr )
      case( 2 )
! ../../monte/inc/cntmnt  /cntmnt/
        ns = nwkmp_nt * ndmsr
        call mpi_send( vmwork, ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
      case( 3 )
! ../../monte/inc/cntmnt  /cntmnt/
        ns = nwkmp_nt * ndmsr
        call mpi_recv( vmwork, ns,   MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
      end select

      return
      end
