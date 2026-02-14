!***********************************************************************
      subroutine pls_cal
!***********************************************************************
!
!          Do not change itend !                              May/09
!
!-----------------------------------------------------------------------
      use cntmnt, only : pfl_ion, psm_abs, psm_man, tnpuf
      use csonic, only : itim, lcstg, lntl, lstop
      use cunit,  only : n6
!     220824, for soldor time Series
      use mod_soldorTimeSeries
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  ia, itcnv, itend0, lcstgb
      integer  lcstgb
      real*8   zfinp, zfexh
      data lcstgb/1/
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   save lcstgb, itcnv
      save lcstgb
!
!xx   catmz:undef  ismax:O.K.  azmas:undef    2012/06/20
!xx      if( itim.le.2 ) then
!xx        write(n6,'(/2x,"impurity : ",a2,"  ismax =",i3,"  azmas =",
!xx     >  f10.3)') catmz, ismax, azmas
!xx      endif
!
      call out_sol2
      call plstep

!     220824, for soldor time Series
      interNow = interNow + 1
      if(interNow.eq.interNum) then
            call calcLatestVal
            call saveLatest
            interNow = 0
      end if
!
      if( lstop.ne.0 ) then
        return
      endif
      call plauxv
!
!::judegment of conversion
      if( lntl.eq.0 )   return
      if( itim.le.200 ) return
!
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   ia    = 1
      zfinp = -pfl_ion(5)+tnpuf
      zfexh = psm_abs+psm_man
      if( lcstg.eq.1 .and. zfexh.gt.zfinp*0.70d0 ) lcstg = 2
      if( lcstg.eq.2 .and. zfexh.gt.zfinp*0.90d0 ) lcstg = 3
      if( lcstg.eq.3 .and. zfexh.gt.zfinp*0.98d0 ) lcstg = 4
!
      if( lcstg.ne.lcstgb ) then
      write(n6,'(2x,"=== change lcstg =",i3," ==>",i3,"  flux ="
     > ,1p2e12.3)') lcstgb,lcstg,zfinp,zfexh
      endif
      lcstgb = lcstg
!
      return
      end