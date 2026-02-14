!**********************************************************************
      subroutine spmd_init
!**********************************************************************
!
!    ngrp(mod_mpicomm) = 4 (master, soldor, neut2d, impmc)
!    gnope, gmspe, gmpbg, gmped   !  bg:begin  ed:end
!
!    conflict variables
!       ngrp, mygrp, n5, n6, n7
!
!     mod_mpicomm : ngrpM = 4              0      1      2      3
!                   grpprg(0:ngrpM-1) = "mst", "sol", "mon", "imp"
!
!     cunit : ngrp = 3
!             cdnm(1:ngrp) = "soldor", "monte", "impmc"
!
!  *** mpista ***  mype = 3 mygrp = 2  mspe  =  0  nope = 32  glwld = 0
!  group name: [GRP1]   gnope =  1  gmspe =  0  gmpbg =  0 gmped =  0
!  group name: [GRP2]   gnope = 31  gmspe =  1  gmpbg =  1 gmped = 31
!
!  code = 3               sol, mon, imp
!    (A) ngrp = 2 case    sol, mon/imp
!    (B) ngrp = 3 case    sol, mon, imp
!
!----------------------------------------------------------------------
      use cunit,       only : cdgr, cdnm, cmype, gmpbg, gmped, gmspe
     >    , gnope, glwld, grnm, lexpe, lmceq, lmspe, lmype, lnope, mjpe
     >    , mpeout, mspe, mygrp, mype, mywld, n5, n6, n7, ncode, ngrp
     >    , nope, stcpy, ufl5, ufl6, ufl7
      use mod_mpicomm, only : grpprg, m6, mygrpM => mygrp, n5M => n5
     >    , n6M => n6, nfnam, ngrpM => ngrp, nmst_cgr, nmst_grp
     >    , nope_cgr, nope_grp, nrnk_cgr, nrnk_grp, nwld_cgr, nwld_grp
     >    , pegrp
      implicit none

!::local variables
      integer :: ii, i, j, igrp, ipmn, ipmx, ipe
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   integer :: igmx, ipno
      character(5) :: ctmp

!:ufl5, ufl6, ugl7, stcpy, n5, n6, n7
      ufl5 = "inppls"
      ufl6 = nfnam
      ufl7 = "erstop"
      stcpy = "delete"
      n5 = n5M
      n6 = n6M
      n7 = m6

!::ngrp, ncode
      ii = 0
      do i = 1, ngrpM-1
        ii = ii + 1
        cdnm(ii) = grpprg(i)
        write(grnm(ii),'(a,i1)')  "GRP", ii
      enddo
      ngrp  = ii
      ncode = ii

!::cdgr
      do i = 1, ngrp
        cdgr(i) = i
      enddo

!::global world
! modified 3/3 lines generalization of MPI communicator by kamata 2020/10/25
!ik   nope  = nope_all
!ik   mspe  = nmst_all
!ik   mype  = nrnk_all
      nope  = nope_cgr
      mspe  = nmst_cgr
      mype  = nrnk_cgr
      mygrp = mygrpM
      mpeout = 3
! modified 1/1 lines generalization of MPI communicator by kamata 2020/10/25
!ik   glwld = nwld_all
      glwld = nwld_cgr

!::gnope, gmspe, gmpbg, gmped in glwld
      do j = 1, ngrp
        igrp = j
        ipmn = +10000
        ipmx = -10000
! modified 1/1 lines generalization of MPI communicator by kamata 2020/10/25
!ik     do i = 1, nope_all
        do i = 1, nope_cgr
          ipe = i-1
          if( pegrp(ipe) /= igrp ) cycle
          ipmn = min0(ipmn,ipe)
          ipmx = max0(ipmx,ipe)
        enddo
        gmspe(igrp) = ipmn
        gmpbg(igrp) = ipmn
        gmped(igrp) = ipmx
        gnope(igrp) = ipmx-ipmn+1
      enddo

!::lnope, lmspe, lmype
      lnope = nope_grp
      lmspe = nmst_grp
      lmype = nrnk_grp
      mywld = nwld_grp

!::no use
      lexpe = 0
      lmceq(1:ngrp) = 0

!::cmype, mjpe
      write(cmype,'(i8.5)') mype
      do i = 1, 8
      mjpe = i
      if( index("123456789",cmype(i:i)).gt.0 ) exit
      enddo

!::wh_day etc

!::debug write
      write(n6,'(/2x,"*** spmd_init ***")')
      if( mygrp == 0 ) then
        ctmp = "master"
      else
        ctmp = cdnm(mygrp)
      endif
      write(n6,'(2x,"ngrp =",i2,"  mygrp =",i2,"  cdnm =",a,
     >  "  mype/lmype =",2i5)') ngrp, mygrp, ctmp, mype, lmype

      return
      end
