c**********************************************************************
      subroutine ntl_pls
c**********************************************************************
!
!      call lnkntl_pls  1 => 2                       ! lmspe
!      call ntl_pls
!         call  mcplas  vne(SOLDOR) => dene(MONTE)   ! lmspe
!         call  MPI_Bcast  vna                       ! grp member
!         call  MPI_Bcast  tfps                      ! grp member
!         call  MPI_Bcast  dene                      ! grp member
!         call  ntplas                               ! grp member
!
!::New structure
!      call lnkntl_pls  1 => 2
!      call lnkntl_pls  2 => 2 for all grp member
!      call ntl_pls
!         call  mcplas
!         call  ntplas
!                                        2015/04/22/  K.Shimizu
!
!----------------------------------------------------------------------
! delete 1 line for use plimp
!      use cputim, only : cntpls
      use cunit,  only : cdgr, lmype, mygrp, n6
      use mpi!,    only : mpi_wtime
      implicit none
!
!::local variables
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   integer  ia, it
! delete 1 line for use plimp
!      real(8) :: cptm0, cptm1
!
!::mpi variables
!-----
!x      real(8) :: MPI_WTIME
!x      integer, parameter ::  MPI_STATUS_SIZE = 1
!-----
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   integer istat(MPI_STATUS_SIZE), nls, nle, nbf, ierr
!
!
      call trmark("ntl_plsN","start")

! delete 1 line for use plimp
!      cptm0 = MPI_WTIME()
      if( mygrp.ne.cdgr(2) ) goto 900

      write(n6,'(2x,"call mcplas  mygrp, lmype =",i3,i6)') mygrp,lmype

      call mcplas
      call ntplas
      call trmark("ntl_plsN","return")
!
 900  continue
! delete 1 line for use plimp
!      cptm1 = MPI_WTIME()
!      cntpls = cptm1-cptm0
!
      return
      end
