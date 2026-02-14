!foc::impurity
      module cimcom
      use csize, only : ndgt, ndis => ndmis, ndwh
      implicit none

! for IMPMC
      real(8) :: aimas = 0.0_8, ami = 0.0_8, amz = 0.0_8, azmas = 0.0_8
      integer :: ismax = 0
      integer, parameter :: impmc_max = 10

!::[MPI_Send in lnkimp_iniSnd]   cimcom/cimcom_1/ (ismax)  10/04/21
!::[MPI_Recv in lnkimp_iniSnd]   cimcom/cimcom_1/ (ismax)  10/04/21

!::uniform plasma
      integer :: irmax = 0, ltest = 0, lufrm = 0
      real(8) :: an0 = 0.0_8, ane = 0.0_8, at0 = 0.0_8, ate = 0.0_8
     >    , ati = 0.0_8, avf = 0.0_8

      integer :: ishez = 0, ishte = 0, ishti = 0, ishvf = 0, lfgti = 0
     >    , lfgtzf = 0
      integer :: mdl_fth_q_base = 0, mdl_fthi_wo_limiter = 0
      integer :: mdl_felec = 0, mdl_fthe_qe_base = 0

!::particle information
      integer, parameter :: ndmp = 500000  ! NEW TimeDep
!      integer, parameter :: ndmp = 40000  ! NEW IAEA 
!      integer, parameter :: ndmp = 20000  ! OLD PSI
!      Change  ndmp2 in cimntl

      integer :: jpmax = 0, npemt = 0, npmax = 0
      integer :: il(ndmp) = 0, ir(ndmp) = 0, is(ndmp) = 0
     >    , ivd(ndmp) = 0, ivs(ndmp) = 0
      integer :: igi(ndmp) = 0, igs(ndmp) = 0
     >    , itg(ndmp) = 0
      integer :: ien(ndmp) = 0, igtem(ndmp) = 0, igtno(ndmp) = 0
     >    , iml(ndmp) = 0
      real(8) :: fi(ndmp) = 0.0_8, rr(ndmp) = 0.0_8, tt(ndmp) = 0.0_8
     >    , zz(ndmp) = 0.0_8
      real(8) :: vr(ndmp) = 0.0_8, vv(ndmp) = 0.0_8, vz(ndmp) = 0.0_8
     >    , wght(ndmp) = 0.0_8
      real(8) :: v(ndmp) = 0.0_8, vvr(ndmp) = 0.0_8, vvz(ndmp) = 0.0_8
      real(8) :: fforc(ndmp) = 0.0_8, tforc(ndmp) = 0.0_8

!::force
      real(8), allocatable :: gdez(:), gdte(:), gdti(:), glte(:)
     >    , glti(:)

      real(8) :: cfez(0:ndis) = 0.0_8, cfgte(0:ndis) = 0.0_8
     >    , cfgti(0:ndis) = 0.0_8

      real(8), allocatable :: heat_flux_para_at_ic_center(:)
     >    ,  heat_flux_para_e_ic_center(:), flmxe_at_impmc(:)
     >    ,  flmxi_at_impmc(:)

!::cross section
      real(8), allocatable :: cdi(:,:), cdi0(:), cdr(:,:), cdLz_i(:,:)
     >    , cdLz_r(:,:), cdLz_cx(:,:)
      real(8) :: dlmtn = 0.0_8, sgmtn = 0.0_8

!::scattering
      integer, parameter :: ndfm = 1001
      integer :: ifm = 0
      real(8), allocatable :: slnv(:,:), slw0(:,:), slw1(:,:)
      real(8) :: fxm1(ndfm,2) = 0.0_8, fxm2(ndfm,2) = 0.0_8
     >    , fxm3(ndfm,2) = 0.0_8, fxvu(ndfm,2) = 0.0_8
      real(8) :: dfm = 0.0_8, fxam(2) = 0.0_8, fxxs(2) = 0.0_8

!::reflecion point in the main plasma
      integer, parameter :: ndrf = 501
      integer :: nrf = 0
      real(8) :: rf1_r(ndrf) = 0.0_8, rf1_z(ndrf) = 0.0_8
     >    , rf2_r(ndrf) = 0.0_8, rf2_z(ndrf) = 0.0_8
      real(8) :: frf = 0.0_8, rf1_ro = 0.0_8, rf2_ro = 0.0_8
      integer :: rf1_c(ndrf) = 0, rf2_c(ndrf) = 0

!::reflection at wall
      real(8) :: rflidp = 0.0_8, rflodp = 0.0_8, rflpmp(4) = 0.0_8
     >    , rflves = 0.0_8, rflwal = 0.0_8, wrfmin = 0.0_8
      real(8) :: he0w(ndwh) = 0.0_8, hrfw(ndwh) = 0.0_8
     >    , hv0w(ndwh) = 0.0_8
      integer :: lwstk(ndwh) = 0
      integer :: lslf = 0  !20121107 KH
      logical :: use_neut_refl = .true. ! true: overwrite rflpmp by rcy_al(neut2d refrection coefficient)

!::diffusion   [MPI_send lnkimp_iniSnd]  cimcom/cimcom_7b/ cdif
      integer, parameter :: mrgn_max = 10
      real(8) :: cdif = 0.0_8                      ! uniform diffusion coeffient input
     > , cdifrg(mrgn_max) = -99.9d0                ! diffusion coeffient input for each segment
     > , cdifrg_imp(mrgn_max,impmc_max) = -99.9d0  ! diffusion coeffient input for each segment for each IMPMC
     > , cdif_imp(impmc_max) = 0.0_8               ! uniform diffusion coeffient input for each IMPMC
      real(8), allocatable :: dfcf(:,:)            ! diffusion coeffient for each cell, index1:valence, index2:cell number
      logical :: match_solimpdiff = .true.

!************************************************!
! set diffusion coefficient by liner interporation along r_omd (in src/IMPMC/imp_pls.f)
!************************************************!
! number of data point. When lerp_points <= 0, it is ignored
! lerp_r    : position for interporation
! lerp_diff : diffusion coefficient value at lerp_r
      integer :: lerp_points_cdif(impmc_max) = 0
      real(8) :: lerp_r_cdif(100,impmc_max)=0.0_8
     > , lerp_diff_cdif(100,impmc_max)=0.0_8

! v_imp_conv:convection verocity in IMPMC. index1:valence, index2:cell number
! This is not given in SONIC at now, but it might given by IMPACT.
      real(8), allocatable :: v_imp_conv(:,:)

!::time
      integer :: itimZ = 0, ltmax = 0, lpmax = 0, ltout = 0
      integer :: lpsct = 0
      real(8) :: dtimZ = 0.0_8, dtout = 0.0_8, timeZ = 0.0_8
     >    , tout = 0.0_8

!::restart file
      character*80, parameter :: cftrZ_def = "DummyDefCftrName0120XX"
      character :: cftrZ*80 = cftrZ_def
      character :: cftwZ*80 = ' '

!::screing
!  mimic IMPV5/inc/cimcom
      real(8) :: sflux = 0.0_8, sptcl = 0.0_8, swtot = 0.0_8
      real(8) :: sdtmz = 0.0_8, sitmp = 0.0_8, stimp = 0.0_8
     >    , stimz = 0.0_8
      real(8) :: sdmy1 = 0.0_8, sdmy2 = 0.0_8, sdmy3 = 0.0_8
      real(8), allocatable :: denZ(:,:), temZ(:,:), wsct(:)
! variables for imsave
      real(8), allocatable :: sflux_save(:),swtot_save(:)
     > , sptcl_save(:),stimp_save(:), sdtmz_save(:)
     > , sdmy1_save(:), sdmy2_save(:), sdmy3_save(:),hcnpt_save(:,:)
     > , stimz_save(:), denZ_save(:,:,:), temZ_save(:,:,:)
     > , wsct_save(:,:), hrgpc_save(:,:)
     > , phwrd_save(:,:), phwci_save(:,:)
     > , phdnz_save(:,:,:), phvlz_save(:,:,:)
     > , hgrpp_save(:,:), hgrpw_save(:,:)
     > , hregp_save(:,:), hregw_save(:,:)
     > , friZ_save(:,:,:), thfZ_save(:,:,:)
     > , vzpZ_save(:,:,:), ionZ_save(:,:,:)
     > , recZ_save(:,:,:), sitmp_save(:)
     > , phfrz_save(:,:,:), phthz_save(:,:,:)
     > , phionz_save(:,:,:), phrecz_save(:,:,:)
     > , hcnwt_save(:,:), hrgpa_save(:,:)
      character(len=4), allocatable :: sptyc_save(:)
      integer :: nwkmp_im1 = 0

      real(8), allocatable :: friZ(:,:), thfZ(:,:), vzpZ(:,:), ionZ(:,:)
     > , recZ(:,:)
      integer :: nwkmp_im1y = 0

      integer :: sptyi = 0
      character :: sptyc*4 = ' '

!-----
!::whit(,,3) => whit(,,4)    ! 2012/04/19
      real(8) :: wpemt = 0.0_8, wtemt = 0.0_8
      real(8) :: wexh(0:ndis,30) = 0.0_8
     > , whit(0:ndis,30,4) = 0.0_8
      real(8), allocatable :: denFlxZtoWall(:,:), eneFlxZtoWall(:,:)
     >    , weflx(:,:,:), wpflx(:,:,:)

      integer, parameter :: nwkmp_im2
     > = 2 + (ndis+1)*30*4 + (ndis+1)*30

!-----
!::[MPI_Send in imhist]  cimcom_12 (hmype,hcnwt)        10/05/26
!::[MPI_Recv in imhist]  cimcom_12 (hmype,hcnwt)        10/05/26
      real(8) :: hcpum = 0.0_8, hmype = 0.0_8, hpmax = 0.0_8
     > , hstep = 0.0_8, htout = 0.0_8
     > , hrgpa(0:mrgn_max) = 0.0_8, hrgpc(0:mrgn_max) = 0.0_8! region of partciles
     > , hcnpt(0:mrgn_max) = 0.0_8, hcnwt(0:mrgn_max) = 0.0_8! condition of particle

      integer, parameter :: nwkmp_im3 = 5 + 15 * 3 + 10 * 3 ! IMPMC_TD value
      integer, parameter :: idemt = 1  ! emit
      integer, parameter :: idcal = 2  ! cal
      integer, parameter :: idexh = 3  ! exhaust
      integer, parameter :: idstk = 4  ! stick (NO trace)
      integer, parameter :: iderr = 5  ! error
      integer, parameter :: idetc = 6  ! et cetra

!-----
      integer, parameter :: ndext = ndgt + 5
      real(8) :: wvaci(0:ndgt,0:ndext) = 0.0_8
     >    , wvacn(0:ndgt,0:ndext) = 0.0_8
      integer, parameter :: nwkmp_im4 = (ndgt+1)*(ndext+1)*2
!-----

!::conductance in vac-region
      integer, parameter :: ndbkf = ndgt + 1
      real(8) :: bkflux(0:ndbkf) = 0.0_8, bkfrac(0:ndbkf) = 0.0_8

!::convergence of calculation 
      real(8) :: fwcal = 0.0_8, fwstp = 0.0_8
!-----
      real(8) :: e0emtZ = 0.0_8, v0emtZ = 0.0_8
      real(8) :: e0bkZ = 0.0_8, e0inZ = 0.0_8, e0tmZ = 0.0_8
!-----
!::impurity flux sputtered from wall
      integer, parameter :: ndbz = 501
      real(8) :: arbz(ndbz) = 0.0_8, flbz(ndbz) = 0.0_8
     >    , prbz(ndbz) = 0.0_8
      real(8) :: csbz(ndbz) = 0.0_8, snbz(ndbz) = 0.0_8
     >    , xpbz(ndbz) = 0.0_8, ypbz(ndbz) = 0.0_8
      real(8) :: gfbz = 0.0_8, tfbz = 0.0_8
      integer :: icbz(ndbz) = 0, ikbz(ndbz) = 0, iwbz(ndbz) = 0
      integer :: iebz(ndbz) = 0, isbz(ndbz) = 0, nbz = 0

      real(8) :: ychem = 0.0_8, yfomt = 0.0_8, yfphy = 0.0_8
     >    , yfslf = 0.0_8

!-----
      real(8) :: r1lim = 0.0_8, r2lim = 0.0_8, z1lim = 0.0_8
     >    , z2lim = 0.0_8

      integer :: ip1tr = 0, ip2tr = 0, iswtr = 0, isw2tr = 0, lp1tr = 0
     >    , lp2tr = 0, mpetr = 0
      integer :: ipcm = 0, lpcm = 0
      integer :: i6_epos = 0, i6_ipos = 0, i6_repo = 0, i6_trac = 0
     >    , i6_wrfi = 0, i6_wrfn = 0

!::impdt   trace of test particle
      integer, parameter :: ndclng = 256
      character :: cstrg*(ndclng) = ' '

!-----
!::option
      integer :: latom = 0, ldifs = 0, lorbt = 0, lpdif = 0
     >    , lscat = 0, lsctmx = 0
      integer :: ldist = 0, lgstk = 0, lmrfl = 0, lpcut = 0, lwrad = 0
      integer :: lhst = 0, lprf = 0, ltdst(50) = 0, lthst(50) = 0

!::calculation condition
      character :: tynsmp*80 = ' ', typspt*80 = ' '
      integer, parameter :: ndspt = 5
   
      character(4), dimension(ndspt) :: twsput = ' '
      integer, dimension(ndspt) :: tlgstk = 0, tnpmax = 0, tnskip = 0
      real(8), dimension(ndspt) :: te0emt = 0.0_8

!::He source on Core Edge
      real(8) :: flxImpAmp = 0.0_8, rhoImpSrc = 0.0_8

!::chemical sputtering model option
!:: Y = ychem/(1+(f/f_0)^e)
!:: (J. Roth et al NF (2004) 44 L21-25 eqs (1) and (2))
      logical :: mdl_roth = .false. ! true: use Roth model
      real*8 :: rothfit_flx = 6.0d21, rothfit_eps=0.54d0! parameter for Roth model(f_0,e)

!::only in IMPMC_TD::!!
!::uniform plasma
      integer :: npmaxB = 0, npmaxL(5) = 0, npmaxdummy = 0
      integer :: ntags(ndmp) = 0
      real(8) :: pemt(ndmp) = 0.0_8

!::particle information
      real(8) :: wstp(ndmp) = 0.0_8
!
!::time
      real(8) :: tmsim = 0.0_8, tmmcz = 0.0_8, tmref = 0.0_8
     >  , dtemt = 0.0_8, dtscr = 0.0_8, dttgt = 0.0_8
     >  , dtexe = 0.0_8, dtsiz = 0.0_8
      integer :: itmmc = 0

!  mimic IMPV5/inc/cimcom
      real(8) :: sitmz = 0.0_8, scpum = 0.0_8

!      real(8) :: phwrd(ndmc) = 0.0_8
!     >, phwci(ndmc) = 0.0_8
!     >, phdnz(0:ndis,ndmc) = 0.0_8
!     >, phvlz(0:ndis,ndmc) = 0.0_8
!     >, phfrz(0:ndis,ndmc) = 0.0_8
!     > ,phthz(0:ndis,ndmc) = 0.0_8
!     > ,phionz(0:ndis,ndmc) = 0.0_8
!     > ,phrecz(0:ndis,ndmc) = 0.0_8
      real(8), allocatable :: phdnz(:,:), phwci(:), phwrd(:), phvlz(:,:)
     >    , phfrz(:,:), phthz(:,:), phionz(:,:), phrecz(:,:)

!     keep the vals. identical as much as possible
!     reference: IMPV5/inc/cimcom
      real(8) :: emflx = 0.0_8, emwgt = 0.0_8, emptl = 0.0_8

!::sample number of 1PE   not NPE      
      real(8) :: myptl(15) = 0.0_8, mywgt(15) = 0.0_8
     > , mypcn(15) = 0.0_8

!::define igxxx  g: group classified by ien(ip)
      integer, parameter :: igemt = 1  ! emit
      integer, parameter :: igcal = 2  ! cal
      integer, parameter :: igsgt = 3  ! sticl at gate
      integer, parameter :: igsce = 4  ! stick at core edge
      integer, parameter :: igmin = 5  ! wght < wmin
      integer, parameter :: igerr = 6  ! error
      integer, parameter :: igexh = 7  ! exhaust at wall
      integer, parameter :: igpmp = 8  ! pump
      integer, parameter :: ignul = 9  ! null  ???  no use
      integer, parameter :: igrun = 10 ! run

      integer, parameter :: igtst = 11  ! tst: at start time
      integer, parameter :: igten = 12  ! ten: at end time
      integer, parameter :: igtav = 13  ! tav: Intg(<Nz>*dv)
      integer, parameter :: igtbl = 14  ! tbl: balance check = Wten

!-----
!::[MPI_Send in imhist]  cimcom_12 (hmype,hcnwt)        10/05/26
!::[MPI_Recv in imhist]  cimcom_12 (hmype,hcnwt)        10/05/26
      real(8) :: hgrpp(15) = 0.0_8, hgrpw(15) = 0.0_8
     > , hgrpn(15) = 0.0_8 ! conditon
      real(8) :: hregp(mrgn_max) = 0.0_8, hregw(mrgn_max) = 0.0_8
     > , hregn(mrgn_max) = 0.0_8 ! region

!
!::convergence of calculation 
      real(8) :: twemt = 0.0_8, twcal = 0.0_8, twabs = 0.0_8
     > , twerr = 0.0_8
!
!-----
      integer :: kspcm = 0, ichcm = 0, ichmx = 0
     > , istcm = 0, istmx = 0
      integer :: ip_trak = 0, i6_trak = 0, i6_pcon = 0
     > , i6_scat = 0, i6_pcnT = 0, i6_flxb = 0
      integer :: i6_wrad = 0, i6_hflxb = 0, i6_mscat = 0

!::partcile reduction (to be implemented)
      integer,parameter :: ndhr = 50
      integer :: nhr = 0
      real(8) :: hdr = 0.0_8, hroh(ndhr) = 0.0_8
     > ,hstz(ndhr,0:ndis) = 0.0_8
      real(8) :: crho(ndmp) = 0.0_8

      integer, parameter :: ndmct = 100
      integer :: lcutz = 0, npcut = 0, jrmax = 0
     > , nccnd = 0, nhsiz = 0
      integer :: cutjr(ndmct) = 0, cutiz(ndmct) = 0

      integer :: nctst = 0
      real(8) :: tmcut = 0.0_8, dtcut = 0.0_8

      integer:: dumlt = 0
! end   used only in IMPMC_TD

! rename ( nwkmp_im* --> nwkmps_im* ) IMPMC values
      integer :: nwkmps_im1 = 0
      integer :: nwkmps_im1y = 0
      integer, parameter :: nwkmps_im3 = 5 + 11 * 4
! end   used only in IMPMC
      end module cimcom

! Note : difference between IMPMC and IMPMC_TD
! only IMPMC :
! hcnpt, hcnwt, hrgpa, hrgpc
! idemt, idcal, idexh, idstk, iderr, idetc

! only IMPMC_TD :
! npmaxB, npmaxdummy, npmaxL
! ntags, pemt
! dtemt, dtexe, dtscr, dtsiz, dttgt, tmmcz, tmref, tmsim
! itmmc
! phdnz, phwci, phwrd, phvlz, phfrz, phthz, phionz, phrecz
! scpum, sitmz
! emflx, emptl, emwgt
! mypcn, myptl, mywgt
! igemt, igcal, igsgt, igsce, igmin, igerr, igexh, igpmp, ignul, igrun,
! igtst, igten, igtav, igtbl
! hgrpn, hgrpp, hgrpw, hregn, hregp, hregw
! twabs, twcal, twemt, twerr
! ichcm, ichmx, istcm, istmx, kspcm
! i6_flxb, i6_pcnT, i6_pcon, i6_scat, i6_trak, ip_trak, i6_hflxb, i6_mscat, i6_wrad
! ndhr, nhr, hdr, hroh, hstz, crho, ndmct, jrmax, lcutz, nccnd, nhsiz, npcut,
! cutiz, cutjr, nctst, dtcut, tmcut, dumlt, ltmax

! deleted from IMPMC_TD
! ltmaxXX

! deleted from IMPMC
! wspt ( only initial value setting with imclear )

! rename IMPMC
! wspt --> wstp, ltmax --> ltmaxXX
! nwkmp_im1 = 10+(ndis+1)*(ndmc+1)*2+(ndmc+1) -->
! nwkmp_im1 = 12+(ndis+1)*(ndmc+1)*3+(ndmc+1)
! nwkmp_im1y = (ndis+1)*(ndmc+1)*5 -->
! nwkmp_im1y = (ndis+1)*(ndmc+1)*4
! nwkmp_im3 = 5 + 11*4 -->
! nwkmp_im3 = 5 + 15*3 +10*3
