!%%mkspt +
!***********************************************************************
      subroutine opnfrt
!***********************************************************************
      use cunit, only : cdgr, cdrout, lmspe, lmype, mygrp, n6
      implicit none
!
!::local variables
! modified 1/3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer   mj, lenx
      integer   mj
! function
      integer   lenx
!
      if( lmype.ne.lmspe )  return
      write(n6,'(2x,"*** opnfrt ***  mygrp =",i2,"  lmype =",i4,
     >   "  n6 =",i4)') mygrp, lmype, n6
!
!::drout
      call getenv( "DROUT", cdrout )
      mj = lenx(cdrout)
      if( mj.eq.0 ) cdrout = "./"
      if( mj.gt.0 ) cdrout = cdrout(1:mj)//"/"
      mj = lenx(cdrout)
!
!-----------------------------------------------------------------------
!::SOLDOR
!-----------------------------------------------------------------------
      if( mygrp.eq.cdgr(1) ) then
!KH      open( unit=60, file=cdrout(1:mj)//"zcpu_S" )
      open( unit=61, file=cdrout(1:mj)//"zpld" )
      open( unit=62, file=cdrout(1:mj)//"zflx" )
!xx      open( unit=64, file=cdrout(1:mj)//"zvrf" )
!xx      open( unit=65, file=cdrout(1:mj)//"zvrc" )
!
      write(n6,'(/2x,"*** opnfrt ***  in PLS")')
!KH      write(n6,602) 60,"zcpu_S","cpu-time"
      write(n6,602) 61,"zpld","monitering data"
      write(n6,602) 62,"zflx","flux"
!KH      write(n6,602) 64,"zvrf","variation of final result"
!KH      write(n6,602) 65,"zvrc","variation of non converged result"
      endif
!
!-----------------------------------------------------------------------
!::NEUT2D
!-----------------------------------------------------------------------
      if( mygrp.eq.cdgr(2) ) then
!KH      open( unit=60, file=cdrout(1:mj)//"zcpu_N" )
      open( unit=68, file=cdrout(1:mj)//"zcnd" )
!
      write(n6,'(/2x,"*** opnfrt ***  in NTL")')
!KH      write(n6,602) 60,"zcpu_N","cpu-time"
      write(n6,602) 68,"zcnd","conductance"
      endif
!
 602  format(4x,i3,2x,a," : ",a)
!
      return
      end
