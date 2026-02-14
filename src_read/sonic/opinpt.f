!***********************************************************************
      subroutine opinpt
!***********************************************************************
      use csonic,    only : limp, lntl, lpls, lqick, lrand, lstep, mxcpu
      use ctopics,   only : ktop_cr, ktop_di, ktop_ese, ktop_esi
     >    , ktop_fni, ktop_fqe, ktop_fqi, ktop_hp, ktop_ne, ktop_ni
     >    , ktop_pp, ktop_ps, ktop_sn, ktop_te, ktop_ti, ktop_we
     >    , ktop_wi, ktop_xe, ktop_xi
      use cunit,     only : gmspe, mspe, mype, n6
      use mod_mpicomm, only : grpprg, mygrp
      use mod_shexe, only : use_namelist, INP_PLS, set_shexe_variables
      implicit none
!
      character cmsg(10)*80,cmsg2(0:5)*16
      character cdsn*80,clin*120
      logical lex
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer nls, nle, nbf, ierr, nft
! modified 1/1 lines addition and reconsider of passed data by kamata 2022/02/23
!ik   integer nft
      integer    ios, nft
      integer :: mhdpe
!
      namelist /usonic/ lstep,lpls,lntl,limp,lqick,lrand,mxcpu
! added 4 lines addition and reconsider of passed data by kamata 2022/02/23
      namelist /utopics/ktop_cr, ktop_di, ktop_ese, ktop_esi, ktop_fni
     >    , ktop_fqe, ktop_fqi, ktop_hp, ktop_ne, ktop_ni, ktop_pp
     >    , ktop_ps, ktop_sn, ktop_te, ktop_ti, ktop_we, ktop_wi
     >    , ktop_xe, ktop_xi
!
!-----------------------------------------------------------------------
!::process to get steady state
!-----------------------------------------------------------------------
      cmsg(1) = "soldor + analytic neutral"
      cmsg(2) = "soldor + neut2d"
      cmsg(3) = "impmc to get refl-coefficients"
      cmsg(4) = "impmc to get initial distribution"
      cmsg(5) = "soldor + neut2d + impmc to increase gradually rad-loss"
      cmsg(6) = "soldor + neut2d + impmc to get steady state"
      cmsg(7) = "check flux"
!
!-----------------------------------------------------------------------
!::meaning of lpls,lntl,limp
!-----------------------------------------------------------------------
      cmsg2(0) = "no effect       "
      cmsg2(1) = "read disk file  "
      cmsg2(2) = "analytical model"
      cmsg2(3) = "run a full model"
!
      call trmark("opinpt","start")
      mhdpe = gmspe(1)   ! master pe in soldor_group

      write(n6,'(4x,"mype,mspe =",2i5,"  mhdpe =",i5)') mype,mspe,mhdpe
!
!-----------------------------------------------------------------------
!::input /usonic/
!-----------------------------------------------------------------------
      if( mype.eq.mhdpe ) then
        write(n6,'(2x,"open file  option of inppls  mype =",i5)') mype
        lqick = 0
        nft = 21
        call set_shexe_variables
        if(use_namelist) then
           call nopen(nft,INP_PLS,"text",cdsn,lex)
        else
           call nopen(nft,"INP_PLS","text",cdsn,lex)
        endif
        read(nft,usonic)
        call check_IMPMC_for_limp
        if( grpprg(mygrp)(1:6) == 'soldor' ) then
          read(nft,utopics,iostat=ios)
          if( ios /= 0 ) then
            write(n6,*) 'namelist data utopics not defined.'
          endif
        endif
        close(nft)
        if( mxcpu.le.0 ) mxcpu = 30
      endif
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      write(n6,'(/2x,"*** S O N I C - code ***   lstep =",i2,$)') lstep
      if( lstep.lt.1 .or. lstep.gt.7 ) goto 110
      write(n6,'(4x,a)') cmsg(lstep)
      write(n6,602) "soldor   lpls =",lpls,cmsg2(lpls)
      write(n6,602) "neut2d   lntl =",lntl,cmsg2(lntl)
      write(n6,602) "impmc    limp =",limp,cmsg2(limp)
      write(n6,603) "random   lrand=",lrand
      write(n6,603) "lqick    lqick=",lqick
      write(n6,604) "cputime  mxcpu=",mxcpu
! added 20 lines addition and reconsider of passed data by kamata 2022/02/23
      if( grpprg(mygrp)(1:6) == 'soldor' ) then
! check if the flag value used by plheat is 0 or 1
        if( ktop_cr  /= 1 ) ktop_cr  = 0
        if( ktop_ese /= 1 ) ktop_ese = 0
        if( ktop_esi /= 1 ) ktop_esi = 0
        if( ktop_ps  /= 1 ) ktop_ps  = 0
        if( ktop_sn  /= 1 ) ktop_sn  = 0
        if( ktop_we  /= 1 ) ktop_we  = 0
        if( ktop_wi  /= 1 ) ktop_wi  = 0
        write(n6,900) 'ktop_ne  = ', ktop_ne,  'ktop_ni  = ', ktop_ni
     >              , 'ktop_te  = ', ktop_te,  'ktop_ti  = ', ktop_ti
        write(n6,900) 'ktop_fni = ', ktop_fni, 'ktop_fqe = ', ktop_fqe
     >              , 'ktop_fqi = ', ktop_fqi, 'ktop_di  = ', ktop_di
        write(n6,900) 'ktop_xe  = ', ktop_xe,  'ktop_xi  = ', ktop_xi
     >              , 'ktop_pp  = ', ktop_pp,  'ktop_hp  = ', ktop_hp
        write(n6,900) 'ktop_sn  = ', ktop_sn,  'ktop_ese = ', ktop_ese
     >              , 'ktop_esi = ', ktop_esi, 'ktop_cr  = ', ktop_cr
        write(n6,900) 'ktop_ps  = ', ktop_ps,  'ktop_we  = ', ktop_we
     >             ,  'ktop_wi  = ', ktop_wi
      endif
 602  format(5x,a,i2,"  : ",a,2x,a,i3,2x,a,1pe11.2)
 603  format(5x,a,i2)
 604  format(5x,a,i5)
      return
! added 1 line addition and reconsider of passed data by kamata 2022/02/23
  900 format( 4( 2x, a, i5 ) )
!
!-----------------------------------------------------------------------
!::error
!-----------------------------------------------------------------------
 110  continue
      write(clin,'("wrong lstep =",i3)') lstep
      call wexit( "aasonic",clin )
      end subroutine opinpt

!***********************************************************************
      subroutine check_IMPMC_for_limp
!***********************************************************************
      use mod_mpicomm, only : grpprg
      use mod_sizedef, only : ndgrp
      use csonic,    only : limp
      implicit none

!:: local variable
      integer i

      do i = 0, ndgrp
        if(grpprg(i)(1:5) == 'IMPMC') then
          limp = 3
          return
        endif
      enddo

      limp = 0
      end subroutine check_IMPMC_for_limp