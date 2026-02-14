! added dynamic allocation of arrays by kamata 2022/05/29
! separated from def_IMP_1  'IMPx_1'
! sending and receiving variables in cimcom.f
      subroutine cimcom_sr( kk, peno )
      use cimcom,      only : dfcf
      use csize,       only : ndmc
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
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 410

      select case( kk )
      case( 1 )
! ../../IMPMC/inc/cimcom  /cimcom_8/
        call mpi_bcast( dfcf,  ndmc, MPI_REAL8,     peno
     >                , nwld_grp, merr )
      case( 2 )
! ../../IMPMC/inc/cimcom  /cimcom_8/
        call mpi_send( dfcf,  ndmc, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, merr )
      case( 3 )
! ../../IMPMC/inc/cimcom  /cimcom_8/
        call mpi_recv( dfcf,  ndmc, MPI_REAL8,     peno
     >               , mtag, nwld_cmm, mstat, merr )
      end select

      return
      end
