!**********************************************************************
      subroutine ntl_pls
!**********************************************************************
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
      use cputim, only : cntpls
      use cunit,  only : cdgr, lmype, mygrp, n6
      use mpi!,    only : mpi_wtime
      implicit none
!
!::local variables
      real(8) :: cptm0, cptm1
!-----
!
      call trmark("ntl_plsN","start")

      cptm0 = MPI_WTIME()
      if( mygrp.ne.cdgr(2) ) goto 900

      write(n6,'(2x,"call mcplas  mygrp, lmype =",i3,i6)') mygrp,lmype

      call mcplas
      call ntplas
      call trmark("ntl_plsN","return")
!
 900  continue
      cptm1 = MPI_WTIME()
      cntpls = cptm1-cptm0
!
      return
      end
