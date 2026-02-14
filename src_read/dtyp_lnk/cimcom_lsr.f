! added dynamic allocation of arrays by kamata 2022/05/29
! separated from strct_comblk  'CIMCOM_11', 'CIMCOMT11', 'CIMCOM_11Y', 'CIMCOMT11Y'
! bcast variables in cimcom.f
!***********************************************************************
      subroutine cimcom_lsr( kd )
!***********************************************************************
      use cimcom, only : denZ, friZ, ionZ, recZ, temZ, thfZ, vzpZ, wsct
      use csize,  only : ndmc, ndmis
      use cunit,  only : lmspe, mywld
      use mod_shexe, only : impmc_model
      use mpi
      implicit none
! arguments
      integer, intent(in) :: kd
! kd   : processing flag
!        = 1 : 'CIMCOM_11', = 2 : 'CIMCPM_11Y'

! local variables
      integer :: merr, mstat(MPI_STATUS_SIZE), mtag = 411, ns

      select case( kd )
      case( 1 ) ! 'CIMCOM_11', 'CIMCOMT11'
! ../../IMPMC/inc/cimcom  /cimcom_11/, /cimcom_11y/
        ns = ( ndmis + 1 ) * ( ndmc + 1 )
        call mpi_bcast( denZ,  ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( temZ,  ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        if( impmc_model == 1 )
     >    call mpi_bcast( vzpZ,  ns,   MPI_REAL8,     lmspe
     >                  , mywld, merr )
        ns = ndmc + 1
        call mpi_bcast( wsct,  ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )

      case( 2 ) ! 'CIMCOM_11Y', 'CIMCOMT11Y'
! ../../IMPMC/inc/cimcom  /cimcom_11y/
        ns = ( ndmis + 1 ) * ( ndmc + 1 )
        call mpi_bcast( friZ,  ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( thfZ,  ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        if( impmc_model == 0 )
     >     call mpi_bcast( vzpZ,  ns,   MPI_REAL8,     lmspe
     >                  , mywld, merr )
        call mpi_bcast( ionZ,  ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
        call mpi_bcast( recZ,  ns,   MPI_REAL8,     lmspe
     >                , mywld, merr )
      end select

      return
      end
