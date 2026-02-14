! added replace all include files with module files by kamata 2021/08/18
      module cunit
      implicit none
!-----------------------------------------------------------------------
!    ngrp, mygrp 
!    conflict with module mpicomm   
!
!    How to resolve   2015/02/06  
!    use mpicomm, ngrpM => ngrp, mygrpM => mygrp     
!
!-----------------------------------------------------------------------

!::number of group
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   integer, parameter :: ndgr = 6
      integer, parameter :: ndgr = 103   ! IMPMC = 3 --> 100 same as mod_sizedef

!::code name & group name
      character :: cdnm(ndgr)*6 = ' ', grnm(ndgr)*6 = ' '
! added 2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
      character(10) :: cprg = ' '
! cprg : code name

!::variables of mpi
      integer :: ncode = 0, ngrp = 0, cdgr(ndgr) = 0
      integer :: nope = 0, mspe = 0, mype = 0, mygrp = 0, mywld = 0
     >    , mpeout = 0
      integer :: gnope(ndgr) = 0, gmspe(ndgr) = 0, gmpbg(ndgr) = 0
     >    , gmped(ndgr) = 0, glwld = 0
      integer :: lnope = 0, lmspe = 0, lexpe = 0, lmype = 0
      integer :: lmceq(ndgr) = 0
      character :: cmype*8 = ' '
      integer :: mjpe = 0

!::unit of I/O and file name
      character :: ufl5*80 = ' ', ufl6*80 = ' ', ufl7*80 = ' '
      character :: stcpy*8 = ' '

!::unit at the present code
      integer :: n5 = 0, n6 = 0, n7 = 0

!::information about this run
      character :: cdrout*80 = ' '
      character :: wh_day*9 = ' ', wh_tim*8 = ' ', wh_lod*80 = ' '
     >    , wh_dir*80 = ' ', wh_job*10 = ' '

!::[MPI_Bcast in lnkall_ini]  cunit/comstjb/ (cdrou#t,emrk)  10/04/21

!::outfile (pls/ntl/imp) unit number
      integer :: ftgrp(3) = 0

      end module cunit
