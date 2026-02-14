!***********************************************************************
      subroutine iminpt_x(nft,kis)
!***********************************************************************
!
!     ldist :  ON(1)/OFF(0)  imdist
!     lpcut :  ON(1) : suppression
!
!     lwrfl : No use (wall reflection)
!     lmrfl : main reflection
!     lgstk : stick at gate ==> lwstk(ih) = 0/1
!     fwstp : = 0.10  cal. stop at the condition of (Wcal/Wemt < fwstp)
!
!     he0w(ndh) : reflection energy of Impurity neutral  imwrfn.f
!     Erfl = Ein  reflection energy of Impurity ion      imwrfi.f
!                    Under development
!
!                                               2012/07/02  K.Shimizu
!     edg_flx : impurity flux at the core edge  2014/09/17  K.Shimizu
!-----------------------------------------------------------------------
      use catcom, only : catmz
      use cimcom, only : aimas, ami, amz, an0, ane, at0, ate, ati, avf
     >    , azmas, cdif, cdifrg, cftrz, cftwz, dlmtn, dtcut, dtemt
     >    , dtimz, dtout, dtscr, dttgt, e0inz, e0bkz, e0tmz, flximpamp
     >    , fwstp, he0w, hrfw, hv0w, irmax, ishez, ishte, ishti, ishvf
     >    , ismax, jrmax, latom, lcutz, ldifs, ldist, lfgti
     >    , lfgtzf, lgstk, lmrfl, lorbt, lpcut, lpdif, lpmax, lscat
     >    , lsctmx, lslf, ltdst, ltest, lthst, ltmax, ltout, lufrm
     >    , lwstk, mdl_felec, mdl_fth_q_base, mdl_fthe_qe_base
     >    , mdl_fthi_wo_limiter, nctst, ndhr, nhr
     >    , nhsiz, npcut, npemt, npmax, nrf, rf1_ro, rf2_ro, rflidp
     >    , rflodp, rflpmp, rflves, rflwal, rhoimpsrc, sgmtn, tmcut
     >    , tynsmp, typspt, wrfmin, ychem, yfomt, yfphy, yfslf
     >    , use_neut_refl
     >    , cdifrg_imp, cdif_imp, mrgn_max, match_solimpdiff
     >    , lerp_points_cdif,lerp_r_cdif,lerp_diff_cdif
     >    , mdl_roth, rothfit_flx, rothfit_eps
      use cimctl, only : dtimp, ftav, lstdy, nimp, nprg
      use cimden, only : csput, npsput, nsput
      use cimfls, only : fls_roh, lfls
      use cimpuf, only : bk_mag, bk_npt, bk_nty, ndpfz, pf_mag, pf_nty
     >    , pf_nwl, pf_px1, pf_px2, pf_py1, pf_py2, pf_rate
      use cntcom, only : chwl, iransd, ncmax, ncmax2, nhwl, temwl
     >    , void_ne, void_ni, void_te, void_ti
      use cntpls, only : dene, teme, temi, vflw
      use cntsrc, only : tden0, teng0
      use cphcns, only : cev, cmp
      use csize,  only : ndwh
      use cunit,  only : grnm, lmype, lnope, mygrp, n6, wh_dir
      use mod_shexe, only : impmc_model
      use cplimp, only : aimasL,azmasL,amiL,amzL,ismaxL
      use dbg_mod,   only : dtmdbgimp, kimpsy, timdbgimp
      implicit none
!
!::argument
      integer, intent(in) :: nft, kis

!
!::local variables
      integer  ic, nran, i, ii, ih
      integer  nsmp, ihvs
      real*8     zrf, ze0, ztw, zv0, pfrt0
      character*4   zwl
      real*8   xran
      integer   k, kmx, mjs(30), mje(30)
      character(80) ::  cline, clin1, clin2
      character(80) ::  cdatm
      character(10) ::  cslf
      integer ido1, ido2, ido3
      integer nsmp_total
      real(8), parameter :: pump_ref_def = 0.98d0

!:: dummy
      integer lcalZ, ndskZ, nftrZ, nftwZ

! function
      integer    lenx
      real(8)    random
!
      namelist /uplzsp/ cdatm
      namelist /uiminp/
     >  cdatm,
     >  lcalZ, lpmax, ltmax, ltout, dtimZ, lstdy,
     >  dtemt, dtscr, dttgt,
     >  ltdst, lthst,
     >  nimp,  dtimp,
     >  aimas, azmas, ismax, npmax, npemt,
     >  typspt, tynsmp,
     >  ltest, lufrm, ane, ate, ati, avf, an0, at0,
     >  lscat, lorbt, ldifs, latom, lpdif, lsctmx,
     >  ishez, ishte, ishti, ishvf, lfgti, lfgtzf,
     >  mdl_fth_q_base, mdl_fthi_wo_limiter,
     >  mdl_fthe_qe_base, mdl_felec,
     >  rf1_ro, rf2_ro, nrf,
     >  ychem, yfphy, yfslf, yfomt,
     >  cdif,  cdifrg, e0inZ,
     >  nran,
     >  nftrZ, cftrZ, nftwZ, cftwZ, ndskZ,
     >  lpcut, ldist, sgmtn, dlmtn,
     >  lmrfl, lgstk, fwstp,
     >  rflodp, rflidp, rflwal, rflves,
     >  rflpmp, wrfmin,
     >  lslf, ftav, nsmp_total,
     >  pf_nty, pf_nwl, pf_rate, pf_mag, pf_px1,pf_px2,pf_py1,pf_py2,
     >  bk_nty, bk_mag, bk_npt,
     >  rhoImpSrc, flxImpAmp,
     >  lfls, fls_roh, iransd,
     >  lcutz, nhr, jrmax, npcut, nhsiz, tmcut, dtcut, nctst
     >  ,void_ne, void_ni, void_te, void_ti
     >  , dtmdbgimp, kimpsy, timdbgimp, use_neut_refl
     >  , mdl_roth, rothfit_flx, rothfit_eps
      namelist /uiminp2/
     >  cdatm,
     >  lcalZ, lpmax, ltmax, ltout, dtimZ, lstdy,
     >  dtemt, dtscr, dttgt,
     >  ltdst, lthst,
     >  nimp,  dtimp,
     >  aimas, azmas, ismax, npmax, npemt,
     >  typspt, tynsmp,
     >  ltest, lufrm, ane, ate, ati, avf, an0, at0,
     >  lscat, lorbt, ldifs, latom, lpdif, lsctmx,
     >  ishez, ishte, ishti, ishvf, lfgti, lfgtzf,
     >  mdl_fth_q_base, mdl_fthi_wo_limiter,
     >  mdl_fthe_qe_base, mdl_felec,
     >  rf1_ro, rf2_ro, nrf,
     >  ychem, yfphy, yfslf, yfomt,
     >  cdif,  cdifrg, e0inZ,
     >  nran,
     >  nftrZ, cftrZ, nftwZ, cftwZ, ndskZ,
     >  lpcut, ldist, sgmtn, dlmtn,
     >  lmrfl, lgstk, fwstp,
     >  rflodp, rflidp, rflwal, rflves,
     >  rflpmp, wrfmin,
     >  lslf, ftav, nsmp_total,
     >  pf_nty, pf_nwl, pf_rate, pf_mag, pf_px1,pf_px2,pf_py1,pf_py2,
     >  bk_nty, bk_mag, bk_npt,
     >  rhoImpSrc, flxImpAmp,
     >  lfls, fls_roh, iransd,
     >  lcutz, nhr, jrmax, npcut, nhsiz, tmcut, dtcut, nctst
     >  ,void_ne, void_ni, void_te, void_ti
     >  , dtmdbgimp, kimpsy, timdbgimp, use_neut_refl
     >  , mdl_roth, rothfit_flx, rothfit_eps
      namelist /uiminp3/
     >  cdatm,
     >  lcalZ, lpmax, ltmax, ltout, dtimZ, lstdy,
     >  dtemt, dtscr, dttgt,
     >  ltdst, lthst,
     >  nimp,  dtimp,
     >  aimas, azmas, ismax, npmax, npemt,
     >  typspt, tynsmp,
     >  ltest, lufrm, ane, ate, ati, avf, an0, at0,
     >  lscat, lorbt, ldifs, latom, lpdif, lsctmx,
     >  ishez, ishte, ishti, ishvf, lfgti, lfgtzf,
     >  mdl_fth_q_base, mdl_fthi_wo_limiter,
     >  mdl_fthe_qe_base, mdl_felec,
     >  rf1_ro, rf2_ro, nrf,
     >  ychem, yfphy, yfslf, yfomt,
     >  cdif,  cdifrg, e0inZ,
     >  nran,
     >  nftrZ, cftrZ, nftwZ, cftwZ, ndskZ,
     >  lpcut, ldist, sgmtn, dlmtn,
     >  lmrfl, lgstk, fwstp,
     >  rflodp, rflidp, rflwal, rflves,
     >  rflpmp, wrfmin,
     >  lslf, ftav, nsmp_total,
     >  pf_nty, pf_nwl, pf_rate, pf_mag, pf_px1,pf_px2,pf_py1,pf_py2,
     >  bk_nty, bk_mag, bk_npt,
     >  rhoImpSrc, flxImpAmp,
     >  lfls, fls_roh, iransd,
     >  lcutz, nhr, jrmax, npcut, nhsiz, tmcut, dtcut, nctst
     >  ,void_ne, void_ni, void_te, void_ti
     >  , dtmdbgimp, kimpsy, timdbgimp, use_neut_refl
     >  , mdl_roth, rothfit_flx, rothfit_eps
      namelist /uiminp4/
     >  cdatm,
     >  lcalZ, lpmax, ltmax, ltout, dtimZ, lstdy,
     >  dtemt, dtscr, dttgt,
     >  ltdst, lthst,
     >  nimp,  dtimp,
     >  aimas, azmas, ismax, npmax, npemt,
     >  typspt, tynsmp,
     >  ltest, lufrm, ane, ate, ati, avf, an0, at0,
     >  lscat, lorbt, ldifs, latom, lpdif, lsctmx,
     >  ishez, ishte, ishti, ishvf, lfgti, lfgtzf,
     >  mdl_fth_q_base, mdl_fthi_wo_limiter,
     >  mdl_fthe_qe_base, mdl_felec,
     >  rf1_ro, rf2_ro, nrf,
     >  ychem, yfphy, yfslf, yfomt,
     >  cdif,  cdifrg, e0inZ,
     >  nran,
     >  nftrZ, cftrZ, nftwZ, cftwZ, ndskZ,
     >  lpcut, ldist, sgmtn, dlmtn,
     >  lmrfl, lgstk, fwstp,
     >  rflodp, rflidp, rflwal, rflves,
     >  rflpmp, wrfmin,
     >  lslf, ftav, nsmp_total,
     >  pf_nty, pf_nwl, pf_rate, pf_mag, pf_px1,pf_px2,pf_py1,pf_py2,
     >  bk_nty, bk_mag, bk_npt,
     >  rhoImpSrc, flxImpAmp,
     >  lfls, fls_roh, iransd,
     >  lcutz, nhr, jrmax, npcut, nhsiz, tmcut, dtcut, nctst
     >  ,void_ne, void_ni, void_te, void_ti
     >  , dtmdbgimp, kimpsy, timdbgimp, use_neut_refl
     >  , mdl_roth, rothfit_flx, rothfit_eps
      namelist /uiminp5/
     >  cdatm,
     >  lcalZ, lpmax, ltmax, ltout, dtimZ, lstdy,
     >  dtemt, dtscr, dttgt,
     >  ltdst, lthst,
     >  nimp,  dtimp,
     >  aimas, azmas, ismax, npmax, npemt,
     >  typspt, tynsmp,
     >  ltest, lufrm, ane, ate, ati, avf, an0, at0,
     >  lscat, lorbt, ldifs, latom, lpdif, lsctmx,
     >  ishez, ishte, ishti, ishvf, lfgti, lfgtzf,
     >  mdl_fth_q_base, mdl_fthi_wo_limiter,
     >  mdl_fthe_qe_base, mdl_felec,
     >  rf1_ro, rf2_ro, nrf,
     >  ychem, yfphy, yfslf, yfomt,
     >  cdif,  cdifrg, e0inZ,
     >  nran,
     >  nftrZ, cftrZ, nftwZ, cftwZ, ndskZ,
     >  lpcut, ldist, sgmtn, dlmtn,
     >  lmrfl, lgstk, fwstp,
     >  rflodp, rflidp, rflwal, rflves,
     >  rflpmp, wrfmin,
     >  lslf, ftav, nsmp_total,
     >  pf_nty, pf_nwl, pf_rate, pf_mag, pf_px1,pf_px2,pf_py1,pf_py2,
     >  bk_nty, bk_mag, bk_npt,
     >  rhoImpSrc, flxImpAmp,
     >  lfls, fls_roh, iransd,
     >  lcutz, nhr, jrmax, npcut, nhsiz, tmcut, dtcut, nctst
     >  ,void_ne, void_ni, void_te, void_ti
     >  , dtmdbgimp, kimpsy, timdbgimp, use_neut_refl
     >  , mdl_roth, rothfit_flx, rothfit_eps
      namelist /uiminp6/
     >  cdatm,
     >  lcalZ, lpmax, ltmax, ltout, dtimZ, lstdy,
     >  dtemt, dtscr, dttgt,
     >  ltdst, lthst,
     >  nimp,  dtimp,
     >  aimas, azmas, ismax, npmax, npemt,
     >  typspt, tynsmp,
     >  ltest, lufrm, ane, ate, ati, avf, an0, at0,
     >  lscat, lorbt, ldifs, latom, lpdif, lsctmx,
     >  ishez, ishte, ishti, ishvf, lfgti, lfgtzf,
     >  mdl_fth_q_base, mdl_fthi_wo_limiter,
     >  mdl_fthe_qe_base, mdl_felec,
     >  rf1_ro, rf2_ro, nrf,
     >  ychem, yfphy, yfslf, yfomt,
     >  cdif,  cdifrg, e0inZ,
     >  nran,
     >  nftrZ, cftrZ, nftwZ, cftwZ, ndskZ,
     >  lpcut, ldist, sgmtn, dlmtn,
     >  lmrfl, lgstk, fwstp,
     >  rflodp, rflidp, rflwal, rflves,
     >  rflpmp, wrfmin,
     >  lslf, ftav, nsmp_total,
     >  pf_nty, pf_nwl, pf_rate, pf_mag, pf_px1,pf_px2,pf_py1,pf_py2,
     >  bk_nty, bk_mag, bk_npt,
     >  rhoImpSrc, flxImpAmp,
     >  lfls, fls_roh, iransd,
     >  lcutz, nhr, jrmax, npcut, nhsiz, tmcut, dtcut, nctst
     >  ,void_ne, void_ni, void_te, void_ti
     >  , dtmdbgimp, kimpsy, timdbgimp, use_neut_refl
     >  , mdl_roth, rothfit_flx, rothfit_eps
      namelist /uiminp7/
     >  cdatm,
     >  lcalZ, lpmax, ltmax, ltout, dtimZ, lstdy,
     >  dtemt, dtscr, dttgt,
     >  ltdst, lthst,
     >  nimp,  dtimp,
     >  aimas, azmas, ismax, npmax, npemt,
     >  typspt, tynsmp,
     >  ltest, lufrm, ane, ate, ati, avf, an0, at0,
     >  lscat, lorbt, ldifs, latom, lpdif, lsctmx,
     >  ishez, ishte, ishti, ishvf, lfgti, lfgtzf,
     >  mdl_fth_q_base, mdl_fthi_wo_limiter,
     >  mdl_fthe_qe_base, mdl_felec,
     >  rf1_ro, rf2_ro, nrf,
     >  ychem, yfphy, yfslf, yfomt,
     >  cdif,  cdifrg, e0inZ,
     >  nran,
     >  nftrZ, cftrZ, nftwZ, cftwZ, ndskZ,
     >  lpcut, ldist, sgmtn, dlmtn,
     >  lmrfl, lgstk, fwstp,
     >  rflodp, rflidp, rflwal, rflves,
     >  rflpmp, wrfmin,
     >  lslf, ftav, nsmp_total,
     >  pf_nty, pf_nwl, pf_rate, pf_mag, pf_px1,pf_px2,pf_py1,pf_py2,
     >  bk_nty, bk_mag, bk_npt,
     >  rhoImpSrc, flxImpAmp,
     >  lfls, fls_roh, iransd,
     >  lcutz, nhr, jrmax, npcut, nhsiz, tmcut, dtcut, nctst
     >  ,void_ne, void_ni, void_te, void_ti
     >  , dtmdbgimp, kimpsy, timdbgimp, use_neut_refl
     >  , mdl_roth, rothfit_flx, rothfit_eps
      namelist /uiminp8/
     >  cdatm,
     >  lcalZ, lpmax, ltmax, ltout, dtimZ, lstdy,
     >  dtemt, dtscr, dttgt,
     >  ltdst, lthst,
     >  nimp,  dtimp,
     >  aimas, azmas, ismax, npmax, npemt,
     >  typspt, tynsmp,
     >  ltest, lufrm, ane, ate, ati, avf, an0, at0,
     >  lscat, lorbt, ldifs, latom, lpdif, lsctmx,
     >  ishez, ishte, ishti, ishvf, lfgti, lfgtzf,
     >  mdl_fth_q_base, mdl_fthi_wo_limiter,
     >  mdl_fthe_qe_base, mdl_felec,
     >  rf1_ro, rf2_ro, nrf,
     >  ychem, yfphy, yfslf, yfomt,
     >  cdif,  cdifrg, e0inZ,
     >  nran,
     >  nftrZ, cftrZ, nftwZ, cftwZ, ndskZ,
     >  lpcut, ldist, sgmtn, dlmtn,
     >  lmrfl, lgstk, fwstp,
     >  rflodp, rflidp, rflwal, rflves,
     >  rflpmp, wrfmin,
     >  lslf, ftav, nsmp_total,
     >  pf_nty, pf_nwl, pf_rate, pf_mag, pf_px1,pf_px2,pf_py1,pf_py2,
     >  bk_nty, bk_mag, bk_npt,
     >  rhoImpSrc, flxImpAmp,
     >  lfls, fls_roh, iransd,
     >  lcutz, nhr, jrmax, npcut, nhsiz, tmcut, dtcut, nctst
     >  ,void_ne, void_ni, void_te, void_ti
     >  , dtmdbgimp, kimpsy, timdbgimp, use_neut_refl
     >  , mdl_roth, rothfit_flx, rothfit_eps
      namelist /uiminp9/
     >  cdatm,
     >  lcalZ, lpmax, ltmax, ltout, dtimZ, lstdy,
     >  dtemt, dtscr, dttgt,
     >  ltdst, lthst,
     >  nimp,  dtimp,
     >  aimas, azmas, ismax, npmax, npemt,
     >  typspt, tynsmp,
     >  ltest, lufrm, ane, ate, ati, avf, an0, at0,
     >  lscat, lorbt, ldifs, latom, lpdif, lsctmx,
     >  ishez, ishte, ishti, ishvf, lfgti, lfgtzf,
     >  mdl_fth_q_base, mdl_fthi_wo_limiter,
     >  mdl_fthe_qe_base, mdl_felec,
     >  rf1_ro, rf2_ro, nrf,
     >  ychem, yfphy, yfslf, yfomt,
     >  cdif,  cdifrg, e0inZ,
     >  nran,
     >  nftrZ, cftrZ, nftwZ, cftwZ, ndskZ,
     >  lpcut, ldist, sgmtn, dlmtn,
     >  lmrfl, lgstk, fwstp,
     >  rflodp, rflidp, rflwal, rflves,
     >  rflpmp, wrfmin,
     >  lslf, ftav, nsmp_total,
     >  pf_nty, pf_nwl, pf_rate, pf_mag, pf_px1,pf_px2,pf_py1,pf_py2,
     >  bk_nty, bk_mag, bk_npt,
     >  rhoImpSrc, flxImpAmp,
     >  lfls, fls_roh, iransd,
     >  lcutz, nhr, jrmax, npcut, nhsiz, tmcut, dtcut, nctst
     >  ,void_ne, void_ni, void_te, void_ti
     >  , dtmdbgimp, kimpsy, timdbgimp, use_neut_refl
     >  , mdl_roth, rothfit_flx, rothfit_eps
      namelist /uiminp10/
     >  cdatm,
     >  lcalZ, lpmax, ltmax, ltout, dtimZ, lstdy,
     >  dtemt, dtscr, dttgt,
     >  ltdst, lthst,
     >  nimp,  dtimp,
     >  aimas, azmas, ismax, npmax, npemt,
     >  typspt, tynsmp,
     >  ltest, lufrm, ane, ate, ati, avf, an0, at0,
     >  lscat, lorbt, ldifs, latom, lpdif, lsctmx,
     >  ishez, ishte, ishti, ishvf, lfgti, lfgtzf,
     >  mdl_fth_q_base, mdl_fthi_wo_limiter,
     >  mdl_fthe_qe_base, mdl_felec,
     >  rf1_ro, rf2_ro, nrf,
     >  ychem, yfphy, yfslf, yfomt,
     >  cdif,  cdifrg, e0inZ,
     >  nran,
     >  nftrZ, cftrZ, nftwZ, cftwZ, ndskZ,
     >  lpcut, ldist, sgmtn, dlmtn,
     >  lmrfl, lgstk, fwstp,
     >  rflodp, rflidp, rflwal, rflves,
     >  rflpmp, wrfmin,
     >  lslf, ftav, nsmp_total,
     >  pf_nty, pf_nwl, pf_rate, pf_mag, pf_px1,pf_px2,pf_py1,pf_py2,
     >  bk_nty, bk_mag, bk_npt,
     >  rhoImpSrc, flxImpAmp,
     >  lfls, fls_roh, iransd,
     >  lcutz, nhr, jrmax, npcut, nhsiz, tmcut, dtcut, nctst
     >  ,void_ne, void_ni, void_te, void_ti
     >  , dtmdbgimp, kimpsy, timdbgimp, use_neut_refl
     >  , mdl_roth, rothfit_flx, rothfit_eps

! not define gaspuff_IMP_time as imgpuff.f is not used in plimp
! not define upldif as imdifus.f is not used in plimp

      write(n6,'(/2x,"*** iminpt ***")')

!::default value
      iransd = 98815  !  initial random number
      ychem = 0.05d0
      yfphy = 2.00d0
      yfslf = 1.50d0
      yfomt = 1.0d0/100.0d0
      if( impmc_model == 0 ) then
        npemt = 100
        dtimp = 1.0d-6       ! No use
      else
        npemt = 1       !SY reduced 100 => 1 because 100 is too much
        dtemt = 1.0d-5
        dtscr = 1.0d-5
        dttgt = 2.0d0*dtscr
        dtimp = dtscr ! dtimp is used in prg_call.f
      endif
      lpcut = 0            ! No use
      lmrfl = 0            ! reflection at surface in core edge
      lgstk = 0            ! stcik at gate
      fwstp = 0.0d0        ! continue to cal. untill fwcal < fwstp
      rflodp = 0.60d0
      rflidp = 0.60d0
      rflwal = 0.60d0
      rflpmp(1:4) = pump_ref_def
      wrfmin = 0.10d0
      ltdst(1:50) = 0
      lthst(1:50) = 0
      tynsmp = "  "
      typspt = "  "
      lslf = 0 ! 20121107 KH
      csput(1:5) = "  "
      npsput(1:5) = 0
      lfgtzf = 0
      pf_rate(1:ndpfz) = 0.0d0
      pf_mag (1:ndpfz) = 0.0d0
      rhoImpSrc = 0.95d0
      flxImpAmp = 1.06d21
      ftav = 0.0d0
      lfls = 0
      fls_roh = (/0.998d0, 0.98d0, 0.96d0, 0.95d0, 0.94d0/)
      nsmp_total = 1 ! 0: for each PE, 1:for all PE
      cdatm = "---"
      mdl_fth_q_base = 0  ! if 1 --> use thermal force formule
!                           based on b.g. heat flux (Homma JCP2013).
      mdl_fthi_wo_limiter = 0 ! if 1 AND mdl_fth_q_base=1,
!                           --> Kinetic thermal force from ion but not flux limited.
      mdl_fthe_qe_base = 1 ! if 1 AND mdl_fth_q_base=1,
!                         --> Thermal force from electron (kinetic): ON
      mdl_felec = 1   ! if 1 AND mdl_fth_q_base=1,
!                         --> Force from electro-static field: ON
      void_ne = 0.0d0
      void_ni = 0.0d0
      void_te = 0.0d0
      void_ti = 0.0d0

      if( impmc_model == 1 ) then
!::imp_cutbk
        lcutz = 1
        nhr   = 50
        jrmax = 48
        npcut = 50
        nhsiz = 10
        tmcut = 0.001d0
        dtcut = 0.001d0
        nctst = 5000
      endif
!
!::input
      select case( kis )
      case(  2 )
        read(nft,uiminp2)
      case(  3 )
        read(nft,uiminp3)
      case(  4 )
        read(nft,uiminp4)
      case(  5 )
        read(nft,uiminp5)
      case(  6 )
        read(nft,uiminp6)
      case(  7 )
        read(nft,uiminp7)
      case(  8 )
        read(nft,uiminp8)
      case(  9 )
        read(nft,uiminp9)
      case( 10 )
        read(nft,uiminp10)
      case default
        read(nft,uiminp)
      end select

      if(cdatm=="---") then
        read(nft,uplzsp)
        if(cdatm.ne."---")then
          write(n6,*)"cdatm should be moved from uplzsp to uiminp"
          write(*,*) "cdatm should be moved from uplzsp to uiminp"
        else
          call wexit("iminpt","cdatm is not defined")
        endif
      endif
!
      call set_atomz_x(cdatm,kis)
      write(n6,'(/2x,"impurity : ",a2,"  ismax =",i3,"  azmas =",
     >   f10.3)') catmz(kis), ismax, azmas

      ami = aimas*cmp
      amz = azmas*cmp
      irmax = ncmax
      lpmax = ltmax*ltout
      dtout = dtimZ*dfloat(ltout)
      aimasL(kis) = aimas
      azmasL(kis) = azmas
      amiL(kis)   = ami
      amzL(kis)   = amz
      ismaxL(kis) = ismax
!xxx      npmax = 0     ! <=== 2007.11.22
!
!::check value
      if( lsctmx.le.0 ) lsctmx = 20
!
!::type of impurity generation  typspt => csput
      if( len_trim(typspt).le.0 ) then
        call wexit("iminpt","typspt No data")
      endif
      clin1 = typspt
      cline = typspt
      call linsep( cline, " ,", kmx, mjs, mje, 30 )
      ii = 0
      typspt = "  "
      do k = 1, kmx
      if( cline(mjs(k):mjs(k)).eq."-" ) cycle
      ii = ii + 1
      typspt(lenx(typspt)+2:) = cline(mjs(k):mje(k))
      csput(ii) = cline(mjs(k):mje(k))
      enddo
      nsput = ii
!
!::number of sample partcile  tynsmp => npsput
      if( len_trim(tynsmp).le.0 ) then
        call wexit("iminpt","tynsmp No data")
      endif
      clin2 = tynsmp
      cline = tynsmp
      call linsep( cline, " ,", kmx, mjs, mje, 30 )
      ii = 0
      tynsmp = "  "
      do k = 1, kmx
        if( cline(mjs(k):mjs(k)).eq."-" ) cycle
        ii = ii + 1
        read(cline(mjs(k):mje(k)),*) nsmp
        if( impmc_model == 0 ) then
          if( nsmp_total .eq. 1 ) then
            if( lnope.gt.1 ) then
              ido1 = lmype+1; ido2 = nsmp; ido3 = lnope
              nsmp=0
              do i = ido1, ido2, ido3
                nsmp = nsmp + 1
              enddo
              if(nsmp.eq.0)nsmp=1
            endif
          endif
        endif
        npsput(ii) = nsmp
      enddo
!
!::check lgstk when Arbk is on
!SY Backflow calc. particle # is defined as the same way as steady state impmc
      do k=1, nsput
      if( csput(k).eq."Bflw" .or. csput(k).eq."Arbk" ) then
        if(csput(1).ne."Arpf" .and. csput(1).ne."Cedg")then
          write(n6,*)
     &   "#1 of csput should be Arbk or Cedg when backflow mode is used"
          call wexit("iminpt","wrong csput with backflow model")
        endif
        if(lgstk.eq.0) then
          write(n6,'(2x,"Arbk is set. lgstk is forced to be 1.")')
          lgstk=1
        endif
        if( nsmp_total .eq. 1 ) then
          if(bk_npt.gt.0)then
          if( lnope.gt.1 ) then
            ido1 = lmype+1; ido2 = bk_npt; ido3 = lnope
            nsmp=0
            do i = ido1, ido2, ido3
              nsmp = nsmp + 1
            enddo
            if(nsmp.eq.0)nsmp=1
            bk_npt = nsmp
          endif
          endif
        endif
      endif
      enddo
!
!::lthst
      if( lthst(1).gt.ltmax ) then
        lthst(1) = ltmax
        lthst(2) = 0
      endif
      if( ltdst(1).gt.ltmax ) then
        ltdst(1) = ltmax
        ltdst(2) = 0
      endif
!
!
      write(n6,'(2x,"ncmax =",i6,"  ncmax2 =",i6,"  irmax =",i6)')
     >   ncmax, ncmax2, irmax
      write(n6,'(2x,"typspt (inp) = ",a)') trim(clin1)
      write(n6,'(2x,"tynsmp (inp) = ",a)') trim(clin2)
      if( ii.ne.nsput ) call wexit("iminpt","wrong  typspt & tynsmp")
!
      write(n6,'(2x,"csput  (cal) =",5(1x,a4,2x))')
     >                                 (trim(csput(i)),i=1,nsput)
      write(n6,'(2x,"npsput (cal) =",5i7)') (npsput(i),i=1,nsput)
!
      write(n6,'(2x,"lstdy =",i2)') lstdy
      write(n6,'(2x,"ltmax =",i7,"  ltout =",i8,"  lpmax =",i15,
     >   "  dtout =",1pe10.2)') ltmax, ltout, lpmax, dtout
      write(n6,'(2x,"lthst =",5i6)') lthst(1:5)
      write(n6,'(2x,"nimp  =",i5,"  dtimp =",1pe10.2)') nimp, dtimp
      write(n6,'(2x,"aimas =",f8.4,"  azmas =",f8.4,"  ismax =",i3)')
     >   aimas, azmas, ismax
      write(n6,'(2x,"ami =",1pe10.2,"  amz =",1pe10.2)') ami, amz
      write(n6,'(2x,"lpdif =",i2,"  lsctmx =",i7)')
     >     lpdif, lsctmx
      write(n6,'(2x,"lscat =",i2,"  lorbt =",i2,"  ldifs =",i2,
     >  "  latom =",i2)') lscat,lorbt,ldifs,latom
      write(n6,'(2x,"ishez =",i2,"  ishte =",i2,"  ishti =",i2,
     >  "  ishvf =",i2)') ishez, ishte, ishti, ishvf
      write(n6,'(2x,"lfgti =",i2,"  0 (fluid)/1 (kinetic type)",
     >   2x,"lfgtzf =",i2,"  0(Neuhauser)/1(B2-classic)")')
     >    lfgti, lfgtzf
      write(n6,'(2x,"cdif  =",f8.4)') cdif
      write(n6,'(2x,"cdifrg(if<0, it is ignored)  =",8f8.4)') 
     > cdifrg(1:8)
!
      write(n6,'(2x,"lmrfl =",i2,"  rf1_ro =",f8.4,"  rf2_ro =",f8.4,
     >  "  nrf =",i4)') lmrfl, rf1_ro, rf2_ro, nrf
      write(n6,'(2x,"lgstk =",i3,"  fwstp =",f8.3)') lgstk, fwstp
!
      write(n6,'(2x,"lpcut =",i3,"  ldist =",i3)') lpcut, ldist
      write(n6,'(2x,"rfl_(odp,idp,wal) =",3f8.4)')
     >    rflodp, rflidp, rflwal
      write(n6,'(2x,"rfl_pmp =",4f8.4)') rflpmp
      write(n6,'(2x,"wrfmin =",f8.4)')   wrfmin
      write(n6,'(2x,"rfl_(odp,idp,wal,pmp) =",4f8.4,"  wrfmin =",
     >   f8.4)') rflodp, rflidp, rflwal, rflpmp, wrfmin
!
      write(n6,'(2x,"ychem =",f8.4,"  yfphy =",f8.4,"  yfslf =",f8.4,
     >   "  yfomt =",f8.4)') ychem, yfphy, yfslf, yfomt
      write(n6,'(2x,"sgmtn=",f8.4,"  dlmtn =",f8.4)') sgmtn, dlmtn
      write(n6,'(2x,"lslf =", i3)')lslf
!
      write(n6,'(2x,"BackFlow model works: bk_npt = ",i6 )') bk_npt
!
      if(ftav.eq.0.0d0) then
      elseif(ftav.gt.0.0d0 .and. ftav.lt.1.0d0) then
        write(n6,'(2x,"time average model works : ftav = ",f6.2)')
     >      ftav    !KH130409
      else
        write(n6,'(2x,"ftav is wrong value: ",f6.2)') ftav
        write(n6,'(2x,"ftav is forced to be 0.0d0")')
        ftav = 0.0d0
      endif

      write(n6,'(/2x,"lfls =",i2,"  fls_roh =",5f9.4)') lfls, fls_roh
      write(n6,'(/2x,"mdl_fth_q_base (1-->ON) =",i2)') mdl_fth_q_base
      write(n6,'(/2x,"mdl_fthi_wo_limiter (1-->ON) =",i2)')
     >     mdl_fthi_wo_limiter
      write(n6,'(/2x,"mdl_fthe_qe_base (1-->ON) =",i2)')
     >     mdl_fthe_qe_base
      write(n6,'(/2x,"mdl_felec (1-->ON) =",i2)') mdl_felec

      write(n6,'(/2x,"ne/ni/te/ti in void =",4e12.3)')
     >    void_ne, void_ni, void_te,void_ti

      if( irmax <= 0 ) call wexit("iminpt","irmax <= 0")
      if( impmc_model == 1 .and. dtemt .gt. dtscr )
     >   call wexit("iminpt","dtemt>dtscr")
!
!::random
      do i = 1, nran
        xran = random(0)
      enddo
      write(n6,'(2x,"nran =",i12,"  xran =",f12.8)') nran, xran

!::uniform plasma
      write(n6,'(2x,"lufrm =",i2)') lufrm
      if( lufrm.eq.1 ) then
        write(n6,'(2x,"Ne =",1pe10.2,"  N0 =",1pe10.2,"  Te =",1pe10.2,
     >    "  Ti =",1pe10.2,"  T0 =",1pe10.2,"  Vf =",1pe10.2)')
     >    ane, an0, ate, ati, at0, avf
!
        irmax = 1
        ic = 1
        dene(ic) = ane
        teme(ic) = ate
        temi(ic) = ati
        vflw(ic,1) = avf
        tden0(ic,1) = an0
        teng0(ic,1) = at0
      endif
!
!::reflection at wall
      do ih = 1, nhwl
      zrf = -1.0d0
      zwl = chwl(ih)
      if( zwl(1:3).eq."W2d" ) zrf = rflidp
      if( zwl(1:3).eq."W4d" ) zrf = rflodp
      if( zwl(3:3).eq."g" )   zrf = 0.0d0
      if( zwl(3:4).eq."a1" )  zrf = rflpmp(1)
      if( zwl(3:4).eq."a2" )  zrf = rflpmp(2)
      if( zwl(3:4).eq."a3" )  zrf = rflpmp(3)
      if( zwl(3:4).eq."a4" )  zrf = rflpmp(4)
      if( zwl(1:1).eq."g" )   zrf = 0.0d0
      if( zwl(1:2).eq."a1" )  zrf = rflpmp(1)
      if( zwl(1:2).eq."a2" )  zrf = rflpmp(2)
      if( zwl(1:2).eq."a3" )  zrf = rflpmp(3)
      if( zwl(1:2).eq."a4" )  zrf = rflpmp(4)
!-----
      if( zwl(1:1).eq."p" )   zrf = rflves
      if( zwl(1:1).eq."s" )   zrf = 0.5d0
!-----
      if( zwl(3:3).eq." " )   then
        if( zwl(1:1).eq."W" ) zrf = rflwal
        if( zwl(1:1).eq."V" ) zrf = rflves
      endif
      if( zrf.lt.0.0d0 ) then
         write(n6,'(2x,"ih =",i3,"  zwl =",a,"  zrf =",f9.4)')
     >     ih, trim(zwl), zrf
         goto 910
      endif
      hrfw(ih) = zrf
      enddo
!
!::reflection energy at wall
      do ih = 1, nhwl
      ze0 = (temwl(ih)+273.0d0)/11600.0d0
      ze0 = 1.5d0*ze0
      if( ze0.lt.0.0d0 ) then
         write(n6,'(2x,"ih =",i3,"  zwl =",a,"  ze0 =",f9.4)')
     >     ih, trim(zwl), ze0
         goto 910
      endif
      he0w(ih) = ze0
      enddo
!
!::stick at gate
      lwstk(1:ndwh) = 0
      if( lgstk.eq.1 ) then
      do ih = 1, nhwl
      if( chwl(ih)(3:3).eq."g" ) lwstk(ih) = 1
      if( chwl(ih)(1:1).eq."g" ) lwstk(ih) = 1
      enddo
      endif
!
!::e0inZ, e0bkZ, e0tmZ
      ihvs = 0
      do ih = 1, nhwl
      if( chwl(ih)(1:4).eq."V1  " ) then
        ihvs = ih
        exit
      endif
      enddo
      if( ihvs.eq.0 ) then
        call wexit("iminpt","no found W3g1")
      endif
!
      e0bkZ = he0w(ihvs)  ! back flow
      e0tmZ = 0.0d0       ! Thompson model
!
      write(n6,'(2x,"cev =",1pe12.3,"  amz =",1pe12.3)') cev, amz
      write(n6,'(2x,"e0inZ =",1pe12.3,"  e0bkZ =",1pe12.3,
     > "  e0tmZ =",1pe12.3)') e0inZ, e0bkZ, e0tmZ
!
      write(n6,'(/2x,"hrfw : reflection for impurity")')
      write(n6,'(2x,"ih",2x,"chwl",2x,"wstk",5x,"hrfw",5x,"tmwl",
     >  3x,"e0w",9x,"v0w")')
      do ih = 1, nhwl
      ztw = he0w(ih)/1.5d0*11600.0d0 - 273.0d0
      zv0 = sqrt(2.0d0*he0w(ih)*cev/amz)
      hv0w(ih) = zv0
      write(n6,'(2x,i2,2x,a4,2x,i4,f9.1,f9.1,1p2e12.3)')
     >  ih, chwl(ih), lwstk(ih), hrfw(ih), ztw, he0w(ih), zv0
      enddo
!
!::pf_rate [Pam3/s] ==> pf_mag [1/s]
      pfrt0 = 0.24d21 ! this coefficient comes from the ideal gas equation
      write(n6,'(2x,"aimas =",1pe11.3,"  azmas =",1pe11.3,"  pfrt0 =",
     >  1pe11.3)') aimas, azmas, pfrt0
      do i = 1, pf_nty
        if(pf_rate(i).ne.0.0d0)then  !KH130419
          pf_mag(i) = pfrt0*pf_rate(i)
        endif
      enddo
!
      write(n6,'(2x,"pf_nty  =",i2)') pf_nty
      if( pf_nty.gt.0 ) then
      write(n6,'(2x,"pf_rate =",1p5e11.3)') pf_rate(1:pf_nty)
      write(n6,'(2x,"pf_mag  =",1p5e11.3)') pf_mag(1:pf_nty)
      endif
!
!xx   self-sputtering of C
      if(lslf.eq.1) then
        cslf=csput(1)
        if(cslf(1:1).ne."C") then
          write(n6,'(2x,"lslf is turend off because of csput = ",a4)')
     >           cslf
          lslf=0
        else
          write(n6,'(2x,"lslf is turend on for C")')
        endif
      endif
!
      return
!
!::error
 910  continue
      call wexit("iminpt","under development  //zwal")
      end
