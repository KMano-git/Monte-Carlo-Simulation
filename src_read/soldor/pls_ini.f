!***********************************************************************
!
!    S O L D O R - code
!
!   << solve 2d fluid equation  >>
!
!     d(ma*na)/dt + div(ma*na*vpa*b) + div(ma*na*vda) - ma*Sna = 0
!
!     d(ma*na*vpa)/dt + div((ma*na*vpa**2+Pa)*b)
!                     + div(ma*na*vpa*Vda)
!                     - div(4/3*eta*grad(vpa))
!                     + (za*na/ne)*div(Pe*b) - (pa+za*na/ne*Pe)*div(b)
!                     + sum_b(mab*na*nuab*(vpa-vpb))
!                     - (za/zeff-1)*za*na*ci*(div(Ti*b)-Ti*div(b))
!                     - (za/zeff-1)*za*na*ce*(div(Te*b)-Te*div(b))
!                     - Spa = 0
!
!     d(sum_a(Ea))/dt + div(sum_a(Ea+Pa)*vpa*b)
!                     + div(sum_a(Ea+Pa)*Vda)
!                     - div(4/3*eta*vpa*grad(vpa))
!                     - div(ni*Xicl*grad(Ti)) - div(ni*Xian*grad(Ti))
!                     + div(Pe*vpe*b) - Pe*div(vpe*b)
!                     + div(Pe*Vde)   - Pe*div(Vde)
!                     + ne*nueq*(Ti-Te)
!                     - Sei = 0
!
!     d(Ee)/dt        + div((Ee+Pe)*vpe*b)
!                     + div((Ee+Pe)*Vde)
!                     - div(ne*Xecl*grad(Te)) - div(ne*Xean*grad(Te))
!                     - div(Pe*vpe*b) + Pe*div(vpe*b)
!                     - div(Pe*Vde)   + Pe*div(Vde)
!                     - ne*nueq*(Ti-Te)
!                     - See = 0
!
!       q1a=m*na,  q2a=m*na*vpa,  q3=sum_a(Ea),  q4=Ee
!
!       Ea = 3/2*na*Ti + 1/2*ma*na*vpa*vpa,  Pa = na*Ta
!       Ee = 3/2*ne*Te,  Pe = ne*Te
!       Vda = -Da/na*grad(na)  Vde = sum_a(za*Vda)
!
!
!           explicit maccromack   95/4/28 - 95/5/09
!           1-ion 1-d tvd scheme  95/5/10 - 95/5/25
!           multi-ion 2-d tvd scheme   98/10/03 -
!           monte carlo neut2d         99/10/05 -
!
!                               coded by k.shimizu
!
!::pls_ini  : read input data / load restart file
!::pls_cal  : set plasma parameter
!                           and calculation of neut2d or impmc
!::pls_src  : calculation source terms every loop (dummy)
!::pls_out  : output
!::pls_end  : epilog
!
!::initial profile
!xx      call plrprf(lscal)    !  1d initial profile
!xx
!xx      call plinit                !  2d initial profile
!xx      call pldisk(nftr,2,cftr)   !  plasma profile from disk data
!xx      call plauxv                !  q => na,va,te,ti
!xx      call plmprf                !  main plasma profile
!xx
!xx      call pledprf               !  modify profile
!xx      call out_mprf(1)           !  Ni, Ti, Xi in main plasma
!xx
!xx      close (n5, status=stcpy)
!xx
!***********************************************************************
      subroutine pls_ini_t
!***********************************************************************
      use cplcom,    only : cftr, nftr, nlp
      use cplmet,    only : itsle
      use csonic,    only : itim, lcstg, lfdbg, lstop
      use cunit,     only : cdgr, cdrout, lmspe, lmype, mygrp, n5, n6
     >    , stcpy
      use mod_shexe, only : use_namelist, MESH_MTR, MESH_NTL, INP_PLS
!::   soldor time series
      use mod_soldorTimeSeries, only:timeNum,interNum,
     > dataNum,icNum,timeSeries_o,timeSeries_i
      implicit none
!
!::local variables
      character cdsn*80
      integer   nft, nf, lscal, i6
      logical   lex
! function
      integer    lenx
!
      i6 = 7000 + 20
!
      if( mygrp.ne.cdgr(1) ) return
      write(n6,'(/2x,"*** pls_ini_t.f ***")')
!
!::clear variables
      call plcler
!
!::numerical constant
      call plcset
!
!::metric data
      nft = 21
      if( use_namelist ) then ! use namelist
        call nopen( nft, MESH_MTR, 'binary', cdsn, lex )
      else ! use environment variable
        call nopen(nft,"MESH_MTR","binary",cdsn,lex)
      endif
      call mtrdsk(nft,"read")
      close (nft)
      call plmtrc     !  %%% 2002/12/20
!
!::KSFUJI  need data of mcel for plauxv
      if(use_namelist) then
!..   use variable in inppls
        call nopen(nft,MESH_NTL,"binary",cdsn,lex)
      else ! use environment variable
        call nopen(nft,"MESH_NTL","binary",cdsn,lex)
      endif
      call ntgdsk(nft,"read")
      close (nft)
!
!::mesh in core
      call rcvftop( 1 ) ! recieve nro from TOPICS
      call set_top_sr( 2 )
      call sndttop( 1 ) ! send mesh data to TOPICS
!
!::input data
      nft = 21
      if( use_namelist ) then ! use namelist
        call ncopy( nft, INP_PLS, 'text', cdsn, lex )
      else ! use environment variable
        call ncopy(nft,"INP_PLS","text",cdsn,lex)
      endif
      open( unit=n5, file=cdsn )
      call pldfin           !  default value
      call plinpt(n5)       !  input  ( refer itmax )

!::   save time series
      allocate( timeSeries_o(timeNum,dataNum,icNum) )
      allocate( timeSeries_i(timeNum,dataNum,icNum) )

!
c::loop control of soldor
      itim  = 0             !  loop number
      nlp   = 0             !  execution flag of pwstep
      lcstg = 1             !  conversion flag
      lstop = 0             !  KSFUJI
!
!::set some variables
      call plmpvl           !  volume in main plasma   move
      call plupwd           !  up-wind flag for y-diffusion
      call pldfan           !  anomalous diffusion
      call plfset           !  option flag
!
!::preparation for plrprf and plrecm and ntsrcm
      write(n6,'(4x)')
      write(n6,'(4x,"call atomhy : ionization & recombination",
     >  " data of Hydrogen  degas2  00/08/29")')
      call atomhy
      call molehy
!
!::disk initialize
      call plhist(1)
!
!::initial set
      call plpset           !  mass, charge of plsma ion
      call plmidp           !  distance in mid-plane
!
!
!::new type
      open(unit=i6,file="chk_manprf")
      write(i6,'(/2x,"-----------exd_tokrcv")')
      call plrprf_t(lscal)          ! exd_tokrcv  atimp
      call plmanprf(1)              ! set man_ti
      write(i6,'(/2x,"-----------plrprf")')
      call plmanprf(3)              ! out vti, man_ti
      call plinit_t                 ! atimp ==> vti
      write(i6,'(/2x,"-----------plinit")')
      call plmanprf(3)

      call pldisk(nftr,2,cftr)      ! q3
      call plauxv                   ! vti
      write(i6,'(/2x,"-----------pldisk")')
      call plmanprf(3)
!
      call plmanprf(2)              ! atimp  correction
      call plconsvq                 ! q3
      call plauxv                   ! vti
      write(i6,'(/2x,"-----------plconsbq")')
      call plmanprf(3)              ! out vti, man_ti
      close (n5, status=stcpy)
!::check
      if( lmype.eq.lmspe ) then
      nf = 21
      open( unit=nf, file=cdrout(1:lenx(cdrout))//"outsta" )
      call pllist(nf)       !  initial profile
      close (nf)
      endif
      if( lfdbg(1).gt.0 ) call metlst
      if( lfdbg(2).gt.0 ) call ploutp(2,itsle)
!
      write(n6,'(2x,"=== pls_ini_t === END")')
!
      return
      end
