!**********************************************************************
      subroutine dtypfile(ctyp,dsn,kact,is_imprst)
!**********************************************************************
!)
!)      dtyp-data => file  for restart run
!)
!)
!) @soldor_0: RSTF  mygrp = 1  ../wxdr_57/RSTF-PLS_2
!)                    mlp/itim =    24    21  time =  3.192762E-01
!) @soldor_0: RSTF  PLS_2   write  ../wxdr_57/RSTF-PLS_2
!)                                  24  1  vti(5,5)=  7.72406172E+01
!)
!)    currentry PLS_2 and PLS_2B for restart data
!)
!)  Note.
!)    when you change dtyp-data, change dtypfile.f at the same time
!)
!---------------------------------------------------------------------
      use catcom,      only : ip_cxr, ip_ion, ip_plt, ip_prb, ip_prc
     >    , ip_rad, ip_rec, ndxdt, ndxkn, ndxrw, ndxx1, ndxx2, nxs, nyp
     >    , xdaty, xdatz, xdevl, xdnam, xdrnk, xdsn, xdtx1, xdtx2, xjend
     >    , xjnum, xjsta, xwmax, xwmin, xwmlt, xwnum, xwspc, xwunt
     >    , xwvar
      use cgdcom,      only : grdx, grdy, hbr, hbt, hbz, mpax, mpman
     >    , mpprv, mpsol, mpsp, mpw1, mpw2, mqd1, mqd2, mqh1, mqh2
     >    , mqs1, mqs2, mqx1, mqx2, npmx, nqmx, psman, psprv, pssol
      use cimcom,      only : aimas, ami, amz, azmas, cdif, dfcf, ismax
     >    , lwrad
      use cimctl,      only : cdirz, dtimp, ftav, icalz, kimp, lstdy
     >    , ndirz, nimp, nsvtm, svdsn, svitm, tmimp
      use cimden,      only : bkspt, csput, eipot, fsput, ncmx, npsput
     >    , nsput, nsrc_spt, nwspt, nzmx, tdnz, twci, twne, twrd, wtsput
      use cntcom,      only : bvx, bvy, bvz, chwl, cpsf, csps, cswl
     >    , cwal, dsgt, e0in, e0pfa, e0pfm, e0wl, flxin, grax, grxp
     >    , gzax, gzxp, iaxs, icgt, icps, icwl, igtw, ihgw, ihwl, ikps
     >    , ikwl, impl, imps, imwl, inwl, ipfiw, ipfmp, ipfmx, ipfnp
     >    , ipgt, iplx, iply, ipmp, ipps, ipwl, irfw, ispx, isxp, isxw
     >    , isyp, isyw, ivb1, ivb2, ivs1, ivs2, iwgt, iwl1, iwl2
     >    , jdp1, jdp2, jsl1, jsl2, jvb1, jvb2, jxp1, jxp2, ldpmd, lelmd
     >    , lemdp, lempf, lemwl, lntmd, lrcmd, lsmax, lsrmd, lwlmd, mcel
     >    , mgrd, migx, migy, mknd, mrgn, mseg
     >    , myit, nbxp, nbyp, ncmax, ncmax2, negt, next, ngas
     >    , ngmax, ngmax2, nhwl, nocl, nogd, nogdv, nogt, noia
     >    , nops, nowl, nowl2, npep, npew, npf, npmp, npps, npsp, npsw
     >    , nptldp, nptlvl, nptlwl, npwl, npwl2, nsgt, pfflx, pfmag
     >    , pfpx1, pfpx2, pfpy1, pfpy2, rcps, rcwl, rcywl, rmas, snps
     >    , snwl, temin_ion, temin_rec, temwl, tlmt_el, tywl, v0in, volm
     >    , wtmin, xpnt, ypnt
      use cntctl,      only : dtntl, imxnt, itmnt, itmntb, kntl
     >    , mdl_ntsm, nclnt, nhfl, nhunt, nini_cal, nntl, tmntl
      use cntmnt,      only : mcflx, mfmax, mnknd, mnsmp, mnvwk, vcsrc
     >    , vflux, visrc, vitim, vitnt, vkflx, vksty, vmwork, vnsmp
     >    , vnsrc, vsmno, vsmty, vtime
      use cntpls,      only : dene, teme
      use cplcom,      only : aion, ama, anamp, anapv, anasl, anemp
     >    , anepv, anesl, animp, anipv, anisl, atemp, atepv, atesl
     >    , atimp, atipv, atisl, aza, cimp, cimprg, cimprg2, flps, flvl
     >    , fvpmp, gfvl, mdl_wrd, nion, q1a, q2a, q3, q4, qtim, rcydt
     >    , rcysp, rcytm, snflx, tn0, tng, tt0, ttg, tv0, tfps, tfvl
     >    , trad, trad_fac, trad_max, trad_min, tradrg, vcs, vea, vna
     >    , vne, vnezef, vni, vva, vve, vte, vti, vzf, wdnz, wrad
      use cplmet,      only : gare, gdsv, gwtm, gwtp, hare, hdsp, hdsv
     >    , hdxm, hdxp, hgdx, hgdy, hpit, hvol, hvsb, hwtm, hwtp, icaxs
     >    , icel, icmax, icmin, icmpe, icspx, icmps, icwl1, icwl2, itmax
     >    , itmpe, itmps, itpve, itpvs, itsle, itsls, jcdp1, jcdp2, jcel
     >    , jcmax, jcxp1, jcxp2, jnxm, jnxp, jtmax, jtmin, kce, kcn, kcs
     >    , kcw, kgdx, kgdy, kreg, nompl, noprv, nosol, romn, vlmn
      use cplwrd,      only : wfac
      use cpmpls,      only : arhmp, armp, drmp, dvmp, fflna, flxni
     >    , flxqe, flxqi, fprna, imd1, imd2, jmd1, jmd2, lnedg, prfna
     >    , prfni, prfte, prfti, r0mp, ramp, rohmp, romp, sfmp, vlmp
     >    , wdmp, wfmp, wlmp, xmd1, xmd2, ymd1, ymd2
      use csize,       only : ndmis, ndwp, ndgt, ndmc, nvxd, nvyd
      use csonic,      only : dtim, itend, itim, kdsk, khst, kpcn, kpls
     >    , lcode, lcstg, lfdbg, lfopt, limp, lmstd, lntl, lpchk, lpcn
     >    , lpls, lpost, lpst_bprf, lqick, lrand, lstep, lstop, mstop
     >    , mxcpu, mxdsk, mxhst, mxprf, ndsk, nfdbg, nfopt, nhsav, nhst
     >    , npcn, npls, nprf, tend, time
      use cunit,       only : cdrout, n6, wh_day, wh_dir, wh_job, wh_lod
     >    , wh_tim
      use czwflx,      only : zw_catmz, zw_eflx, zw_ismax, zw_pflx
      use mod_dtflw,   only : mlp, ntyp, typ_dnam
      use mod_mpicomm, only : wf_exit
      implicit none
!::argument
      character(*), intent(in) :: ctyp, dsn
      integer,      intent(in) :: kact
      logical,      intent(out) :: is_imprst

!::local variables
      integer :: ityp, ift = 281
      integer :: maxz, i, j, ic, iz
      integer         nc
      character(2)    np
! nc  : code number
! np  : processing number
      real(8) :: zne, zte

!::temporary variables
      integer, parameter :: f_size = 5
      integer :: nsput_d, nzmx_d, ncmx_d, npsput_d(f_size)
      integer :: nsrc_spt_d
      real(8) :: fsput_d(f_size), wtsput_d(f_size)
      real(8) :: bkspt_d(ndgt)
      real(8) :: eipot_d(0:ndmis)
      real(8) :: twrd_d(ndmc), tdnz_d(0:ndmis,ndmc)
     > , twci_d(ndmc), twne_d(ndmc)
      character :: csput_d(f_size)*4, nwspt_d*4
!
!::dummy
      integer noclv(0:nvxd,0:nvyd)
!
      noclv = 0
      is_imprst = .true.
!
      call tbfind( ntyp, typ_dnam, ityp, ctyp )
      if( ityp <= 0 ) then
      call wf_exit("dtypfile","no found ctyp in typ_dnam "//trim(ctyp))
      endif

!::send  => write
      if( kact == 1 ) then

      open(unit=ift,file=dsn,action="write",form="unformatted")

      select case(trim(ctyp))

!::MST_1 : ../dtyp_lnk/def_MST_1.f
      case("MST_1")
       write(ift)
     >  lstep,  lcode,  ndsk,  nhst,  nhsav,  nprf,  mxcpu,  mxdsk,
     >  mxhst,  mxprf,  lfopt,  lfdbg,  lpst_bprf,  nfopt,  nfdbg,
     >  lpls,  npls,  kpls,  lpcn,  npcn,  kpcn,  lpchk,  kdsk,
     >  khst,  lntl,  limp,  lmstd,  lqick,  lpost,  lrand,
     >  time, tend,  dtim,  itim,  itend,  mstop,  lstop, lcstg,
     >  cdrout, wh_day,  wh_tim,  wh_lod,  wh_dir,  wh_job
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,i6)') "RSTF",
     > "MST_1 ","write",trim(dsn),mlp,"itim =",itim

!::MST_2 : ../dtyp_lnk/def_MST_2.f
      case("MST_2")
       write(ift)
     >  time,  tend,  dtim,  itim,  itend,  mstop,  lstop,  lcstg,
     >  kntl,  kimp,  kpcn,  kdsk,  khst,  kpls
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     > "MST_2 ","write",trim(dsn),mlp,"time =",time

!::PLS_1 : ../dtyp_lnk/def_PLS_1.f
      case("PLS_1")
       write(ift)
     >  aion,  aza,  ama,  nion,  grdx,  grdy,  hbr,  hbz,  hbt,
     >  pssol,  psprv,  psman,  mqd1,  mqx1,  mqs1,  mqs2,  mqh1,
     >  mqh2,  mqx2,  mqd2,  nqmx,  mpw1,  mpsp,  mpw2,  mpax,
     >  npmx,  mpsol,  mpprv,  mpman,  vlmn,  romn,  hvol,  hgdx,
     >  hgdy,  hdxm,  hdxp,  hdsp,  hdsv,  hvsb,  hare,  hpit,
     >  hwtm,  hwtp,  gdsv,  gare,  gwtm,  gwtp,  nosol,  noprv,
     >  nompl,  jcdp1,  jcdp2,  jcxp1,  jcxp2,  jcmax,  icwl1,
     >  icspx,  icwl2,  icmps,  icmpe,  icaxs,  itsls,  itsle,
     >  itpvs,  itpve,  itmps,  itmpe,  itmax,  icmin,  icmax,
     >  jtmin,  jtmax,  kgdx,  kgdy,  kreg,  jcel,  icel,  jnxp,
     >  jnxm,  kce,  kcw,  kcn,  kcs,  flxni,  flxqi,  flxqe,
     >  fflna,  prfna,  prfni,  prfti,  prfte,  fprna,  r0mp,
     >  ramp,  armp,  arhmp,  romp,  rohmp,  vlmp,  sfmp,  dvmp,
     >  drmp,  wdmp,  wfmp,  wlmp,  xmd1,  ymd1,  xmd2,  ymd2,
     >  jmd1,  imd1,  jmd2,  imd2,  lnedg,  wrad,  wdnz,  tradrg,
     >  trad,  trad_min,  trad_max,  trad_fac,  cimp,  cimprg,
     >  cimprg2,  rcydt,  rcysp,  rcytm,  mdl_wrd
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "PLS_1 ","write",trim(dsn),mlp,"flxni =",flxni

!::PLS_1B : ../dtyp_lnk/def_PLS_1B.f
      case("PLS_1B")
       write(ift)
     >   mdl_wrd
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,i3)') "RSTF",
     >  "PLS_1B","write",trim(dsn),mlp,"mdl_wrd=",mdl_wrd

!::PLS_2 : ../dtyp_lnk/def_PLS_2.f
      case("PLS_2")
       write(ift)
     >  vna,  vne,  vni,  vnezef,  vzf,  vva,  vve,  vti,  vte,
     >  vcs,  vea,  anamp,  animp,  anemp,  atimp,  atemp,  fvpmp,
     >  anasl,  anisl,  anesl,  atisl,  atesl,  anapv,  anipv,
     >  anepv,  atipv,  atepv,  tfps,  flps,  tfvl,  gfvl,  flvl,
     >  q1a,  q2a,  q3,  q4,  qtim,  tN0,  tT0,  tV0,  tNg,  tTg,
     >  snflx
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "PLS_2 ","write",trim(dsn),mlp,"vti(5,5)=",vti(5,5)

!::PLS_2B : ../dtyp_lnk/def_PLS_2B.f
      case("PLS_2B")
       write(ift)
     >  vna,  vne,  vni,  vnezef,  vzf,  vva,  vve,  vti,  vte,
     >  vcs,  vea,  anamp,  animp,  anemp,  atimp,  atemp,  fvpmp,
     >  anasl,  anisl,  anesl,  atisl,  atesl,  anapv,  anipv,
     >  anepv,  atipv,  atepv,  tfps,  flps,  tfvl,  gfvl,  flvl
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "PLS_2B","write",trim(dsn),mlp,"vti(5,5)=",vti(5,5)

!::PLS_3 :
      case("PLS_3")
       write(ift)
     >  wfac
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "PLS_3","write",trim(dsn),mlp,"wfac=",wfac

!::NTL_1 : ../dtyp_lnk/def_NTL_1.f
      case("NTL_1")
       write(ift)
     >  xpnt,  ypnt,  volm,  bvx,  bvy,  bvz,  grax,  gzax,  grxp,
     >  gzxp,  myit,  nogd,
     >  nocl,  mcel,  nogdv,  noclv,  mseg,  mgrd,  mknd,  mrgn,
     >  migx,  migy,  next,  iplx,  iply,  ngmax,  ngmax2,  ncmax,
     >  ncmax2,  jdp1,  jxp1,  jsl1,  jsl2,  jxp2,  jdp2,  iwl1,
     >  ispx,  iwl2,  impl,  iaxs,  ivs1,  ivs2,  jvb1,  jvb2,
     >  ivb1,  ivb2,  rcywl,  temwl,  rcwl,  e0wl,  cswl,  snwl,
     >  cwal,  chwl,  tywl,  npsw,  npew,  irfw,  ihgw,  inwl,
     >  ikwl,  ipwl,  icwl,  isxw,  isyw,  igtw,  ihwl,
     >  ipmp,  nowl,  npwl,  npmp,  nowl2,  npwl2,  nhwl,
     >  dsgt,  nsgt,  negt,  iwgt,  ipgt,
     >  icgt,  nogt,  rcps,  csps,  snps,  cpsf,  npsp,  npep,
     >  nbxp,  nbyp,  ikps,  ipps,  icps,  isxp,  isyp,  imps,
     >  imwl,  nops,  npps,  pfflx,  ipfmp,  ipfiw,  ipfmx,  tmntl,
     >  dtntl,  nntl,  imxnt,  nhunt,  nhfl,  nini_cal,  mdl_ntsm,
     >  kntl,  nclnt,  itmnt,  itmntb,  mfmax,  mnknd,  mnsmp,
     >  mnvwk,  mcflx,  rmas,  e0in,  v0in,  flxin,  wtmin,  pfmag,
     >  pfpx1,  pfpx2,  pfpy1,  pfpy2,  e0pfa,  e0pfm,  temin_ion,
     >  temin_rec,  tlmt_el,  ngas,  noia,  nptldp,  nptlwl,  nptlvl,
     >  lsmax,  lemwl,  lemdp,  lempf,  ipfnp,  npf,  lntmd,  ldpmd,
     >  lwlmd,  lrcmd,  lelmd,  lsrmd
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "NTL_1 ","write",trim(dsn),mlp,"xpnt(15)=",xpnt(15)


!::NTL_2 : ../dtyp_lnk/def_NTL_2.f
      case("NTL_2")
       write(ift)
     >  vtime,  vflux,  vmwork,  vitim,  vitnt,  vnsrc,  vsmty,
     >  vsmno,  visrc,  vkflx,  vksty,  vnsmp,  vcsrc
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "NTL_2 ","write",trim(dsn),mlp,"vflux(2)=",vflux(2)
!::
      case default
        if(  ctyp(1:3) == 'IMP' ) then
! get code number and processing number
          call getncnp( ctyp, nc, np )

! execution
          if( nc > 0 .and. np /= ' ' ) then
            select case( np )
! IMP(nc)_1
            case( '1 ' )
              write(ift)
     >          xdaty(1:ndxdt,nc), xdtx1(1:ndxx1,1:ndxkn,nc)
     >        , xdtx2(1:ndxx2,1:ndxkn,nc), xdatz(1:ndmis,1:ndxkn,nc)
     >        , xwmin(1:ndxrw,1:ndxkn,nc), xwmax(1:ndxrw,1:ndxkn,nc)
     >        , xwmlt(1:ndxrw,1:ndxkn,nc), ip_ion(nc), ip_rec(nc)
     >        , ip_cxr(nc), ip_rad(nc), ip_plt(nc), ip_prb(nc)
     >        , ip_prc(nc), nxs(nc), nyp(nc)
     >        , xdrnk(1:ndxkn,nc), xwnum(1:ndxrw,1:ndxkn,nc)
     >        , xjsta(1:ndxkn,nc), xjend(1:ndxkn,nc), xjnum(1:ndxkn,nc)
     >        , xdsn(1:ndxkn,nc), xdnam(1:ndxkn,nc), xdevl(1:ndxkn,nc)
     >        , xwvar(1:ndxrw,1:ndxkn,nc), xwunt(1:ndxrw,1:ndxkn,nc)
     >        , xwspc(1:ndxrw,1:ndxkn,nc)
     >        , tmimp, dtimp, ftav, nimp, kimp, icalZ, lstdy, svitm
     >        , nsvtm, ndirZ, svdsn, cdirZ, aimas, azmas, ami, amz
     >        , ismax, cdif, dfcf, lwrad
              write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)')
     >          "RSTF ", trim( ctyp ), "write ", trim( dsn ), mlp
     >        , "tmimp =",tmimp
! IMP(nc)_2
            case( '2 ' )
              write(ift)
     >          nsput, nzmx, ncmx, npsput, nsrc_spt, fsput, wtsput
     >        , bkspt, eipot, twrd, tdnz, twci, twne, csput, nwspt

              write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)')
     >          "RSTF ", trim( ctyp ), "write ", trim( dsn ), mlp
     >        , "fsput(1) =", fsput(1)

              maxz = min0( nzmx, 3 )
              do ic = 100, 300, 50
                zne = 0.0d0
                zte = 0.0d0
                j = iplx(ic)
                i = iply(ic)
                if( j.gt.0 .and. i.gt.0 ) then
                  zne = vne(j,i)
                  zte = vte(j,i)
                endif

                write(n6,'(2x,i6,2i5,"  Ne =",1p2e12.3,"  Te =",1p2e12.3
     >          , " twne,twrd,tdnz =",1p2e12.3,1p5e12.3)')
     >            ic, j, i, zne, dene(ic), zte, teme(ic), twne(ic)
     >          ,  twrd(ic), (tdnz(iz,ic),iz=0,maxz)
              enddo
! IMP(nc)_2B
            case( '2B' )
              write(ift)
     >          zw_eflx(0:ndmis,1:ndwp,nc), zw_pflx(0:ndmis,1:ndwp,nc)
     >        , zw_ismax(nc), zw_catmz(nc)
              write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,a)')
     >          "RSTF ", trim( ctyp ), "write ", trim( dsn ), mlp
     >        , "zw_catmz =", zw_catmz(nc)
            end select
          endif
        else 
           call wf_exit("dtypfile","no found MPI_Bcast("//trim(ctyp))
        endif
      end select
      close(ift)

      endif        !  write

!::recv  => read
      if( kact == 2 ) then
      open(unit=ift,file=dsn,action="read",form="unformatted")
      select case(trim(ctyp))

!::MST_1 : ../dtyp_lnk/def_MST_1.f
      case("MST_1")
       read(ift)
     >  lstep,  lcode,  ndsk,  nhst,  nhsav,  nprf,  mxcpu,  mxdsk,
     >  mxhst,  mxprf,  lfopt,  lfdbg,  lpst_bprf,  nfopt,  nfdbg,
     >  lpls,  npls,  kpls,  lpcn,  npcn,  kpcn,  lpchk,  kdsk,
     >  khst,  lntl,  limp,  lmstd,  lqick,  lpost,  lrand,
     >  time, tend,  dtim,  itim,  itend,  mstop,  lstop, lcstg,
     >  cdrout, wh_day,  wh_tim,  wh_lod,  wh_dir,  wh_job
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,i6)') "RSTF",
     > "MST_1 ","read ",trim(dsn),mlp,"itim =",itim

!::MST_2 : ../dtyp_lnk/def_MST_2.f
      case("MST_2")
       read(ift)
     >  time,  tend,  dtim,  itim,  itend,  mstop,  lstop,  lcstg,
     >  kntl,  kimp,  kpcn,  kdsk,  khst,  kpls
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     > "MST_2 ","read ",trim(dsn),mlp,"time =",time

!::PLS_1 : ../dtyp_lnk/def_PLS_1.f
      case("PLS_1")
       read(ift)
     >  aion,  aza,  ama,  nion,  grdx,  grdy,  hbr,  hbz,  hbt,
     >  pssol,  psprv,  psman,  mqd1,  mqx1,  mqs1,  mqs2,  mqh1,
     >  mqh2,  mqx2,  mqd2,  nqmx,  mpw1,  mpsp,  mpw2,  mpax,
     >  npmx,  mpsol,  mpprv,  mpman,  vlmn,  romn,  hvol,  hgdx,
     >  hgdy,  hdxm,  hdxp,  hdsp,  hdsv,  hvsb,  hare,  hpit,
     >  hwtm,  hwtp,  gdsv,  gare,  gwtm,  gwtp,  nosol,  noprv,
     >  nompl,  jcdp1,  jcdp2,  jcxp1,  jcxp2,  jcmax,  icwl1,
     >  icspx,  icwl2,  icmps,  icmpe,  icaxs,  itsls,  itsle,
     >  itpvs,  itpve,  itmps,  itmpe,  itmax,  icmin,  icmax,
     >  jtmin,  jtmax,  kgdx,  kgdy,  kreg,  jcel,  icel,  jnxp,
     >  jnxm,  kce,  kcw,  kcn,  kcs,  flxni,  flxqi,  flxqe,
     >  fflna,  prfna,  prfni,  prfti,  prfte,  fprna,  r0mp,
     >  ramp,  armp,  arhmp,  romp,  rohmp,  vlmp,  sfmp,  dvmp,
     >  drmp,  wdmp,  wfmp,  wlmp,  xmd1,  ymd1,  xmd2,  ymd2,
     >  jmd1,  imd1,  jmd2,  imd2,  lnedg,  wrad,  wdnz,  tradrg,
     >  trad,  trad_min,  trad_max,  trad_fac,  cimp,  cimprg,
     >  cimprg2,  rcydt,  rcysp,  rcytm,  mdl_wrd
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "PLS_1 ","read ",trim(dsn),mlp,"flxni =",flxni

!::PLS_1B : ../dtyp_lnk/def_PLS_1B.f
      case("PLS_1B")
       read(ift)
     >   mdl_wrd
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,i3)') "RSTF",
     >  "PLS_1B","read ",trim(dsn),mlp,"mdl_wrd=",mdl_wrd

!::PLS_2 : ../dtyp_lnk/def_PLS_2.f
      case("PLS_2")
       read(ift)
     >  vna,  vne,  vni,  vnezef,  vzf,  vva,  vve,  vti,  vte,
     >  vcs,  vea,  anamp,  animp,  anemp,  atimp,  atemp,  fvpmp,
     >  anasl,  anisl,  anesl,  atisl,  atesl,  anapv,  anipv,
     >  anepv,  atipv,  atepv,  tfps,  flps,  tfvl,  gfvl,  flvl,
     >  q1a,  q2a,  q3,  q4,  qtim,  tN0,  tT0,  tV0,  tNg,  tTg,
     >  snflx
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "PLS_2 ","read ",trim(dsn),mlp,"vti(5,5)=",vti(5,5)

!::PLS_2B : ../dtyp_lnk/def_PLS_2B.f
      case("PLS_2B")
       read(ift)
     >  vna,  vne,  vni,  vnezef,  vzf,  vva,  vve,  vti,  vte,
     >  vcs,  vea,  anamp,  animp,  anemp,  atimp,  atemp,  fvpmp,
     >  anasl,  anisl,  anesl,  atisl,  atesl,  anapv,  anipv,
     >  anepv,  atipv,  atepv,  tfps,  flps,  tfvl,  gfvl,  flvl
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "PLS_2B","read",trim(dsn),mlp,"vti(5,5)=",vti(5,5)

!::PLS_3 :
      case("PLS_3")
       read(ift)
     >  wfac
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "PLS_3","read",trim(dsn),mlp,"wfac=",wfac

!::NTL_1 : ../dtyp_lnk/def_NTL_1.f
      case("NTL_1")
       read(ift)
     >  xpnt,  ypnt,  volm,  bvx,  bvy,  bvz,  grax,  gzax,  grxp,
     >  gzxp,  myit,  nogd,
     >  nocl,  mcel,  nogdv,  noclv,  mseg,  mgrd,  mknd,  mrgn,
     >  migx,  migy,  next,  iplx,  iply,  ngmax,  ngmax2,  ncmax,
     >  ncmax2,  jdp1,  jxp1,  jsl1,  jsl2,  jxp2,  jdp2,  iwl1,
     >  ispx,  iwl2,  impl,  iaxs,  ivs1,  ivs2,  jvb1,  jvb2,
     >  ivb1,  ivb2,  rcywl,  temwl,  rcwl,  e0wl,  cswl,  snwl,
     >  cwal,  chwl,  tywl,  npsw,  npew,  irfw,  ihgw,  inwl,
     >  ikwl,  ipwl,  icwl,  isxw,  isyw,  igtw,  ihwl,
     >  ipmp,  nowl,  npwl,  npmp,  nowl2,  npwl2,  nhwl,
     >  dsgt,  nsgt,  negt,  iwgt,  ipgt,
     >  icgt,  nogt,  rcps,  csps,  snps,  cpsf,  npsp,  npep,
     >  nbxp,  nbyp,  ikps,  ipps,  icps,  isxp,  isyp,  imps,
     >  imwl,  nops,  npps,  pfflx,  ipfmp,  ipfiw,  ipfmx,  tmntl,
     >  dtntl,  nntl,  imxnt,  nhunt,  nhfl,  nini_cal,  mdl_ntsm,
     >  kntl,  nclnt,  itmnt,  itmntb,  mfmax,  mnknd,  mnsmp,
     >  mnvwk,  mcflx,  rmas,  e0in,  v0in,  flxin,  wtmin,  pfmag,
     >  pfpx1,  pfpx2,  pfpy1,  pfpy2,  e0pfa,  e0pfm,  temin_ion,
     >  temin_rec,  tlmt_el,  ngas,  noia,  nptldp,  nptlwl,  nptlvl,
     >  lsmax,  lemwl,  lemdp,  lempf,  ipfnp,  npf,  lntmd,  ldpmd,
     >  lwlmd,  lrcmd,  lelmd,  lsrmd
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "NTL_1 ","read ",trim(dsn),mlp,"xpnt(15)=",xpnt(15)


!::NTL_2 : ../dtyp_lnk/def_NTL_2.f
      case("NTL_2")
       read(ift)
     >  vtime,  vflux,  vmwork,  vitim,  vitnt,  vnsrc,  vsmty,
     >  vsmno,  visrc,  vkflx,  vksty,  vnsmp,  vcsrc
      write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)') "RSTF",
     >  "NTL_2 ","read ",trim(dsn),mlp,"vflux(2)=",vflux(2)
      case default
!::
        if(  ctyp(1:3) == 'IMP' ) then
! get code number and processing number
          call getncnp( ctyp, nc, np )

! execution
          if( nc > 0 .and. np /= ' ' ) then
            select case( np )
! IMP(nc)_1
            case( '1 ' )
              read(ift)
     >          xdaty(1:ndxdt,nc), xdtx1(1:ndxx1,1:ndxkn,nc)
     >        , xdtx2(1:ndxx2,1:ndxkn,nc), xdatz(1:ndmis,1:ndxkn,nc)
     >        , xwmin(1:ndxrw,1:ndxkn,nc), xwmax(1:ndxrw,1:ndxkn,nc)
     >        , xwmlt(1:ndxrw,1:ndxkn,nc), ip_ion(nc), ip_rec(nc)
     >        , ip_cxr(nc), ip_rad(nc), ip_plt(nc), ip_prb(nc)
     >        , ip_prc(nc), nxs(nc), nyp(nc)
     >        , xdrnk(1:ndxkn,nc), xwnum(1:ndxrw,1:ndxkn,nc)
     >        , xjsta(1:ndxkn,nc), xjend(1:ndxkn,nc), xjnum(1:ndxkn,nc)
     >        , xdsn(1:ndxkn,nc), xdnam(1:ndxkn,nc), xdevl(1:ndxkn,nc)
     >        , xwvar(1:ndxrw,1:ndxkn,nc), xwunt(1:ndxrw,1:ndxkn,nc)
     >        , xwspc(1:ndxrw,1:ndxkn,nc)
     >        , tmimp, dtimp, ftav, nimp, kimp, icalZ, lstdy, svitm
     >        , nsvtm, ndirZ, svdsn, cdirZ, aimas, azmas, ami, amz
     >        , ismax, cdif, dfcf, lwrad

              write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)')
     >          "RSTF ", trim( ctyp ), "read ", trim( dsn ), mlp
     >        , "tmimp =",tmimp
! IMP(nc)_2
            case( '2 ' )
              read(ift)
     >          nsput_d, nzmx_d, ncmx_d, npsput_d, nsrc_spt_d, fsput_d
     >        , wtsput_d
     >        , bkspt_d, eipot_d, twrd_d, tdnz_d, twci_d, twne_d
     >        , csput_d, nwspt_d
!
!:: check the impurity source is consistent or not
              if(nsput_d==nsput .and. nzmx_d==ismax) then ! at now nzmx = 0
                do i = 1, f_size
                  if(trim(csput_d(i)) /= trim(csput(i))) then
                    is_imprst = .false.
                    close(ift)
                    return
                  endif
                enddo
              else
                is_imprst = .false.
                close(ift)
                return
              endif
!
              rewind(ift)
              read(ift)
     >          nsput, nzmx, ncmx, npsput, nsrc_spt, fsput, wtsput
     >        , bkspt, eipot, twrd, tdnz, twci, twne, csput, nwspt

              write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,1pe16.8)')
     >          "RSTF ", trim( ctyp ), "read ", trim( dsn ), mlp
     >        , "fsput(1) =", fsput(1)
! IMP(nc)_2B
            case( '2B' )
              read(ift)
     >          zw_eflx(0:ndmis,1:ndwp,nc), zw_pflx(0:ndmis,1:ndwp,nc)
     >        , zw_ismax(nc), zw_catmz(nc)
              write(n6,'(2x,a,2x,a,2x,a,2x,a,2x,i5,2x,a,a)')
     >          "RSTF ", trim( ctyp ), "read ", trim( dsn ), mlp
     >        , "zw_catmz =", zw_catmz(nc)
            end select
          endif
        else 
          call wf_exit("dtypfile","no found MPI_Bcast("//trim(ctyp))
        endif
      end select

      close(ift)
      endif        !  read

      call flush(n6)
      return
      end
