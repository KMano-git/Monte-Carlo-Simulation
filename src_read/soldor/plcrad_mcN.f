!**********************************************************************
      subroutine plcrad_mc
!**********************************************************************
!
!    energy loss due to interactions with impurity
!
!         wtime, witim, wnlp, wfac
!         wmci0, wmce0    IMPMC
!         wmci,  wmce     soldor   normalization
!                wcre     corona   wcre = W + Ar rad
!         wimi,  wime     total energy loss due to impurity
!
!  sign of radiation
!     IMPMC  +5.0 MW/m3 : twrd
!     soldor -5.0 MW/m3 : wcre, wmce, wmci, wime, wimi
!
!
!::Note
!     IMPMC/calwrd.f
!         twrd(ic)  iz + il + rc
!         twci(ic)  collision with ions
!
!   Note  mc0 : IMPMC
!         mc  : soldor normalization
!         im  : impurity
!         nt  : neutral
!         tt  : total
!
!         wrg_mci0(10) : twci  in impmc
!         wrg_mce0(10) : twrd  in impmc
!
!         wrg_mce(10) : radiation in soldor
!
!         wrg_cre(10) : corona
!
!         wrg_imi(10) : wimi in soldor   fac*twci
!         wrg_ime(10) : wime in soldor   corona + fac*twrd
!
!         wrg_nti(10) : swic+swiv*q3  with neutral particles
!         wrg_nte(10) : swec+swev*q4  with neutral partciles
!
!         wrg_tti(10) : swic+swiv*q3  with neutral and impurity
!         wrg_tte(10) : swec+swev*q4  with neutral and impurity
!
!          ir = mrgn(ic) ==> ir2 = mrgnp(ic)
!         plasma parameter    radiation
!
!
!    mdl_wrd = 0   no control
!    mdl_wrd = 1   control of total radiation
!    mdl_wrd = 2   control of radiation in SOL/Edge
!    mdl_wrd = 3   control of radiation in Div/Prv
!
!----------------------------------------------------------------------
      use catcom,  only : catmz
      use cntcom,  only : iplx, iply, mrgn, ncmax, volm
      use cntmnt,  only : dotn, i6_src, sdotn2, totmwe, totwe
      use cntsrc,  only : trgwe
      use cphcns,  only : cev
      use cplcom,  only : anamp, animp, atemp, atimp, crad_totwe, nlp
     >    , nlpmn, q3, q4, mdl_wrd, swec, swev, swic, swiv, trad
     >    , trad_max, trad_min, tradrg, wfac_lv 
      use cplimp,  only : nzmxl, wmc_nty
      use cplmet,  only : hvol, icel, itmax, itpve, jcel, jtmax, kreg
      use cplqcn,  only : mrgnp
      use cplwrd,  only : wcre, wfac, wime, wimi, witim, wmce, wmci
     >    , wnlp, wtime
      use cplwrd2, only : wmc_zwe, wmc_zwi
      use csize,   only : ndmc, nzsmx
      use csonic,  only : itim, time
      use cunit,   only : cdgr, mygrp, n6
      use cpmpls,  only : flxqi, flxqe
     > , nty_ctl, rad_const_step, rad_const_fac
      implicit none
!
      integer :: lsmwrd
!
!::local common
      real*8,dimension(10)   :: wrg_cre
      real*8,dimension(10)   :: wrg_mce
      real*8,dimension(10)   :: wrg_ime
      real*8,dimension(10)   :: wrg_nti, wrg_nte  ! ntl
      real*8,dimension(10)   :: wrg_tti, wrg_tte  ! ntl + imp
      real*8,dimension(10,nzsmx)   :: wrg_mci0, wrg_mce0
!
      real*8 :: smwrd(0:nzsmx) ! 0=CR model, 1~=IMPMC
      real*8 :: totwrd, fac(0:nzsmx), trgwrd, totwrd2
      real*8 :: tot_cre, tot_mce, tot_ime
      real*8 :: tot_nti, tot_nte, tot_tti, tot_tte
      real*8 ::  tot_mci0(nzsmx), tot_mce0(nzsmx)
      real*8 :: zwe, zwec, zwev, zwi, zwic, zwiv
!
!::local variables
      integer  :: ic, j, i, it, jmax, jt, ift, ldbg
      integer  :: ir, ir2, nty
      real*8   :: hne, hni, hte, hti,hq3, hq4
!
!::KSFUJI 0519
      if( nzmxL(1).le.0 ) then
        write(n6,'(2x,"check plrad_mc  return  itim =",i7)') itim
        return
      endif
!
!::zero clear
      wcre(1:ndmc) = 0.0d0
      wmci(1:ndmc) = 0.0d0
      wmce(1:ndmc) = 0.0d0
      wimi(1:ndmc) = 0.0d0
      wime(1:ndmc) = 0.0d0
!
      wrg_cre(1:10) = 0.0d0
      wrg_mce(1:10) = 0.0d0
      wrg_ime(1:10) = 0.0d0
      wrg_nti(1:10) = 0.0d0
      wrg_nte(1:10) = 0.0d0
      wrg_tti(1:10) = 0.0d0
      wrg_tte(1:10) = 0.0d0
      wrg_mci0(1:10,:) = 0.0d0
      wrg_mce0(1:10,:) = 0.0d0
!
!::total radiation in IMPMC
!::energy loss due to neutral particles
      do ic = 1, ncmax    ! KSFUJI
        j  = iplx(ic)
        i  = iply(ic)
        if( j.le.0 .or. i.le.0 ) cycle
        ir2 = mrgnp(ic)
        wrg_nti(ir2)=wrg_nti(ir2)
     >   +(swic(j,i)+swiv(j,i)*q3(j,i))*hvol(j,i)
        wrg_nte(ir2)=wrg_nte(ir2)
     >   +(swec(j,i)+swev(j,i)*q4(j,i))*hvol(j,i)
      enddo
!
!
!::corona-model for background impurity (mcr)
      call plwrdcrG(wcre)   ! <===== 2015/06/11
!-----
!xx      wcr1(1:ndmc) = wcr_wrd(1:ndmc,1)  ! W
!xx      wcr2(1:ndmc) = wcr_wrd(1:ndmc,2)  ! Ar
!xx      wcre(1:ndmc) = wcr1(1:ndmc) + wcr2(1:ndmc)
!-----
!
!::IMPMC results
      call plwrdMC ! nty=1 to wmc_nty wmc_zwi, wmc_zwe
!
      fac(0:nzsmx) = 1.0d0
!
!******************************************************************************************
! During the first few steps, force the limit to impurity radiation power.
! Assume that impurity of nty = 1 is controlled by trad_min/max
!******************************************************************************************
      if( itim <= rad_const_step ) then
!::evaluate radiation
        smwrd   = 0.0d0
        totwrd2  = 0.0d0
        trgwrd  = 0.0d0
!
        do ic = 1, ncmax  ! KSFUJI
          j = iplx(ic)
          i = iply(ic)
          if( j.le.0 .or. i.le.0 ) cycle
          ir2 = mrgnp(ic)
!
!::total radiation in the system (except the hot core)
          if( ir2.le.6 ) then
            smwrd(0) = smwrd(0) + wcre(ic)*volm(ic)
            do nty = 1, wmc_nty
              if(mdl_wrd.ge.1 .and. nty == nty_ctl) cycle
              smwrd(nty) = smwrd(nty) + wmc_zwe(ic,nty)*volm(ic)
            enddo
          endif
        enddo ! ic
!
!:: limit of impurity radiation
        smwrd = -smwrd
!
        do nty = 0, wmc_nty ! 0:CR model
            if(nty == nty_ctl) cycle
            totwrd2 = totwrd2 + smwrd(nty)
        enddo

        ! force trgwrd = flxqi+flxqe
        trgwrd = (flxqi+flxqe)*rad_const_fac
        if(totwrd2.ge.trgwrd)then
          fac = trgwrd/totwrd2
          if(mdl_wrd.ge.1) fac(nty_ctl) = 1.0d0
        endif
      endif
!******************************************************************************************

!:: mdl_wrd 1,2,3
      if( mdl_wrd.ge.1) then
        !::evaluate radiation
        smwrd   = 0.0d0
        totwrd  = 0.0d0
        totwrd2 = 0.0d0
        !
        do ic = 1, ncmax  ! KSFUJI
          j = iplx(ic)
          i = iply(ic)
          if( j.le.0 .or. i.le.0 ) cycle
          ir2 = mrgnp(ic)
          !
          !::total radiation in the system (except the hot core)
          lsmwrd = 0
          if(mdl_wrd.eq.1 ) then
            if( ir2.le.6 ) lsmwrd = 1
          elseif(mdl_wrd.eq.2) then
            if( ir2.eq.2 .or. ir2.eq.6  ) lsmwrd = 1
          elseif(mdl_wrd.eq.3) then
            if( ir2.le.5 .and. ir2.ne.2 ) lsmwrd = 1 ! ir2 = 1, 3, 4, 5
          endif
          !
          if( lsmwrd.eq.1)then
            smwrd(0) = smwrd(0) + wcre(ic)*volm(ic)*fac(0)
            do nty = 1, wmc_nty
              smwrd(nty) = smwrd(nty)
     >          + wmc_zwe(ic,nty)*volm(ic)*fac(nty)
            enddo
          endif
        enddo ! ic
!
!:: limit of impurity radiation
        smwrd = -smwrd

!:: sum of smwrd
        do nty = 0, wmc_nty ! 0:CR model
          totwrd = totwrd + smwrd(nty)
          if(nty.ne.nty_ctl) totwrd2 = totwrd2 + smwrd(nty)
        enddo
!
        if( totwrd.le.trad_min ) totwrd = trad_min
        if( totwrd.ge.trad_max ) totwrd = trad_max
        if( smwrd(nty_ctl).le.0.0d0 .or. totwrd.le.totwrd2 ) then
          fac(nty_ctl) = 0.0d0
        else
          fac(nty_ctl) = (totwrd-totwrd2)/smwrd(nty_ctl)
        endif
!
!KH130725 add limit of fac variation per internal iteration
        if(itim.gt.1 .and. wfac_lv.gt.0.0d0)then
          fac(nty_ctl) = min(fac(nty_ctl), (1.0d0+wfac_lv)*wfac)
          fac(nty_ctl) = max(fac(nty_ctl), (1.0d0-wfac_lv)*wfac)
        endif
      endif
!:: end mdl_wrd 1,2,3
!
!
!::set energy loss with normalization
      do nty = 1, wmc_nty
        wmci(1:ncmax) = wmci(1:ncmax) + wmc_zwi(1:ncmax,nty)*fac(nty)
        wmce(1:ncmax) = wmce(1:ncmax) + wmc_zwe(1:ncmax,nty)*fac(nty)
      enddo
      wimi(1:ncmax) = wmci(1:ncmax)   ! not taken into accoutn, see below
      wime(1:ncmax) = wcre(1:ncmax) + wmce(1:ncmax)
!
!::setup source terms & integ energy loss
 100  continue
      do ic = 1, ncmax
        j = iplx(ic)
        i = iply(ic)
        ir  = mrgn(ic)
        ir2 = mrgnp(ic)
        if( j.le.0 .or. i.le.0 ) cycle
!
        if( ir.eq.7 ) then
          hni = animp(i)
          hne = anamp(i,1)
          hti = atimp(i)
          hte = atemp(i)
          hq3 = 1.5d0*hni*hti*cev
          hq4 = 1.5d0*hne*hte*cev
        else
          hq3 = q3(j,i)
          hq4 = q4(j,i)
        endif
!
!::radiation
        zwe  = wime(ic)
        zwec = 0.0d0
        zwev = zwe/hq4
        swec(j,i)  = swec(j,i) + zwec
        swev(j,i)  = swev(j,i) + zwev
!
!::exchange enrgy due to scattering process with impurity ions
!xx     zwi  = wimi(ic)
        zwi  = 0.0d0     ! heating ion in almost region  2009/01/12
        zwic = 0.0d0
        zwiv = zwi/hq3
        swic(j,i) = swic(j,i) + zwic
        swiv(j,i) = swiv(j,i) + zwiv
!
!::region
        wrg_cre(ir2) = wrg_cre(ir2) + wcre(ic)*volm(ic)*fac(0)
        do nty = 1, wmc_nty
          wrg_mci0(ir2,nty) = wrg_mci0(ir2,nty)
     >     + wmc_zwi(ic,nty)*volm(ic)*fac(nty)
          wrg_mce0(ir2,nty) = wrg_mce0(ir2,nty)
     >     + wmc_zwe(ic,nty)*volm(ic)*fac(nty)
        enddo
        wrg_mce(ir2) = wrg_mce(ir2) + wmce(ic)*volm(ic)
        wrg_ime(ir2) = wrg_ime(ir2) + wime(ic)*volm(ic)
!
!::energy loss due to neutral and impurity
        wrg_tti(ir2)=wrg_tti(ir2)
     >   +(swic(j,i)+swiv(j,i)*q3(j,i))*hvol(j,i)
        wrg_tte(ir2)=wrg_tte(ir2)
     >   +(swec(j,i)+swev(j,i)*q4(j,i))*hvol(j,i)
      enddo
!
!::total radiation in the system
      tot_mci0 = 0.0d0
      tot_mce0 = 0.0d0
      tot_cre  = 0.0d0
      tot_mce  = 0.0d0
      tot_ime  = 0.0d0
      tot_nti  = 0.0d0
      tot_nte  = 0.0d0
      tot_tti  = 0.0d0
      tot_tte  = 0.0d0
!
      do ir2 = 1, 6  !  6:core edge 7:hot core
        do nty = 1, wmc_nty
          tot_mci0(nty) = tot_mci0(nty) + wrg_mci0(ir2,nty)
          tot_mce0(nty) = tot_mce0(nty) + wrg_mce0(ir2,nty)
        enddo
        tot_cre  = tot_cre  + wrg_cre(ir2)
        tot_mce  = tot_mce  + wrg_mce(ir2)
        tot_ime  = tot_ime  + wrg_ime(ir2)
        tot_nti  = tot_nti  + wrg_nti(ir2)
        tot_nte  = tot_nte  + wrg_nte(ir2)
        tot_tti  = tot_tti  + wrg_tti(ir2)
        tot_tte  = tot_tte  + wrg_tte(ir2)
      enddo
!
!::trad (common variables)
      trad = tot_ime
!
!::conv-check
      totwe(1:10) = 0.0d0
      if( mygrp.eq.cdgr(2) ) then
        do ir = 1, 7
          totwe(ir) = trgwe(ir)
        enddo
      endif
!
      do it = 2, itmax-1
        if( it.eq.itpve ) cycle
        jmax = jtmax(it)
        do jt = 2, jmax-1
          j  = jcel(jt,it)
          i  = icel(jt,it)
          ir = kreg(j,i)
          totwe(ir) = totwe(ir) 
     >              + (swec(j,i)+swev(j,i)*q4(j,i))*hvol(j,i)
        enddo
      enddo
!
!                neutral   radiation
!xx   totwe(7) = totmwe + tradrg(7)
      totwe(7) = totmwe + wrg_ime(7)  ! region_impurity_electron
!
      totwe(10) = 0.0d0
      do ir = 1, 7
        totwe(10) = totwe(10) + totwe(ir)
      enddo
!
!::total energy
      do ir = 1, 10
        crad_totwe(ir) = totwe(ir)
      enddo
!
!::total radiation (corona+mc)
      do ir = 1, 7
        tradrg(ir) = wrg_ime(ir)
      enddo
!
!::variables of cmwrad
      wtime = time
      witim = itim
      wnlp  = nlp
      wfac  = fac(nty_ctl)
!
      ift = i6_src
      ift = n6
      if( ift.eq.0 ) return
!
      ldbg = 0
      if(itim.le.1) ldbg = 1 !SYamoto
      if( mod(itim,10).eq.0 .and. nlp.ge.nlpmn ) ldbg = 1
!
      if( ldbg.eq.1 ) then
      write(ift,602) "NEUT","Swi","--",dotn,1.0d0,tot_nti,
     >                                  (wrg_nti(ir2),ir2=1,7)
      write(ift,602) "NEUT","Swe","--",dotn,1.0d0,tot_nte,
     >                                  (wrg_nte(ir2),ir2=1,7)
      do nty = 1, wmc_nty
      write(ift,602) "cradMC","Swi",catmz(nty), dotn,fac(nty),
     >               tot_mci0(nty),(wrg_mci0(ir2,nty),ir2=1,7)
      write(ift,602) "cradMC","Swe",catmz(nty), dotn,fac(nty),
     >               tot_mce0(nty),(wrg_mce0(ir2,nty),ir2=1,7)
      enddo
      write(ift,602) "cradCR","Swe","--",dotn,1.0d0,tot_cre,
     >                                  (wrg_cre(ir2),ir2=1,7)
      write(ift,602) "cradMC","Swe","--",dotn,fac(1),tot_mce,
     >                                  (wrg_mce(ir2),ir2=1,7)
      write(ift,602) "cradTT","Swe","--",dotn,fac(1),tot_ime,
     >                                  (wrg_ime(ir2),ir2=1,7)
      write(ift,602) "NTLIMP","Swi","--",dotn,1.0d0,tot_tti,
     >                                  (wrg_tti(ir2),ir2=1,7)
      write(ift,602) "NTLIMP","Swe","--",dotn,1.0d0,tot_tte,
     >                                  (wrg_tte(ir2),ir2=1,7)
!
      write(ift,602) "total","SWe","--",sdotn2,1.0d0,
     >   crad_totwe(10),(crad_totwe(ir),ir=1,7)
 602  format(2x,a6,2x,a3,2x,a2,x,1pe11.3,1pe11.3,1pe14.6,1x,1p12e11.3)
      ! output for time dependent analysis 
      if(nlp.eq.nlpmn+1) then
      write(2023,'(2x,1pe15.6,1x,1pe11.3,1pe11.3,1pe14.6,1x,1p12e11.3)')
     >   time,dotn,fac(1),tot_ime,(wrg_ime(ir2),ir2=1,7)
      endif
      endif
!
      return
      end