! added replace all include files with module files by kamata 2021/08/18
!----------------------------------------------------------------------
!::grid information for neut2d     soldori/cntcom
!----------------------------------------------------------------------
      module cntcom
      use csize, only : ndgs, ndgt, ndgtp, ndpm, ndwh, ndwl
      implicit none

      real(8) :: grax = 0.0_8, gzax = 0.0_8, grxp = 0.0_8, gzxp = 0.0_8
      real(8), allocatable :: xpnt(:), ypnt(:), volm(:), bvx(:), bvy(:)
     >    , bvz(:)
      integer, allocatable :: mclx(:), mcly(:), mcly2(:), mxjw(:)
     >    , mxjw2(:), myit(:), nogd(:,:), nocl(:,:), mcel(:,:)
     >    , nogdv(:,:), mseg(:), mgrd(:,:), mknd(:,:)
     >    , mrgn(:), migx(:), migy(:), next(:,:), iplx(:), iply(:)
      integer :: ngmax = 0, ngmax2 = 0, ncmax = 0, ncmax2 = 0, jdp1 = 0
     >    , jxp1 = 0, jsl1 = 0, jsl2 = 0, jxp2 = 0, jdp2 = 0, iwl1 = 0
     >    , ispx = 0, iwl2 = 0, impl = 0, iaxs = 0, ivs1 = 0, ivs2 = 0
     >    , jvb1 = 0, jvb2 = 0, ivb1 = 0, ivb2 = 0

!::[MPI_Bcast in ntgdsk]   cntcom (xpnt,emrk)  10/04/21
!::[MPI_Send  in lnkntl_iniSnd]  cntcom/cntgrd/ (xpnt,emrk)  10/04/21
!::[MPI_Recv  in lnkntl_iniSnd]  cntcom/cntgrd/ (xpnt,emrk)  10/04/21

      real(8) :: rcywl(ndwh) = 0.0_8, temwl(ndwh) = 0.0_8
      real(8), allocatable :: rcwl(:), e0wl(:), cswl(:), snwl(:)
      character :: cwal(ndwl)*4 = ' ', chwl(ndwh)*4 = ' '
      character, allocatable :: tywl(:)*4
      integer :: npsw(ndwl) = 0, npew(ndwl) = 0, irfw(ndwh) = 0
     >    , ihgw(ndwh) = 0
      integer, allocatable :: inwl(:), ikwl(:), ipwl(:), icwl(:)
     >    , isxw(:), isyw(:), ixyw(:), igtw(:), ihwl(:), ipmp(:)
      integer :: nowl = 0, npwl = 0, npmp = 0, nowl2 = 0, npwl2 = 0
     >    , nhwl = 0

!::[MPI_Bcast in ntgdsk]  cntcom (rcywl,emrk)  10/04/21
!::[MPI_Send  in lnkntl_iniSnd]  cntcom/cntwal/ (rcywl,emrk)  10/04/21
!::[MPI_Recv  in lnkntl_iniSnd]  cntcom/cntwal/ (rcywl,emrk)  10/04/21

!::[MPI_Bcast in ntgdsk]   cntcom (xsgt,emrk)  10/04/21
      real(8) :: xsgt(ndgt,2) = 0.0_8, xegt(ndgt,2) = 0.0_8
     >    , ysgt(ndgt,2) = 0.0_8, yegt(ndgt,2) = 0.0_8
      real(8), allocatable :: dsgt(:)
      integer :: nsgt(ndgt,2) = 0, negt(ndgt,2) = 0
      integer, allocatable :: iwgt(:), ipgt(:), icgt(:)
      integer :: nogt = 0

!::conductance
      integer :: nopm = 0, nwpm(ndpm) = 0
      real(8) :: cflx(ndpm) = 0.0_8, cqfl(ndpm) = 0.0_8
     >    , cpre(ndpm) = 0.0_8, ccnd(ndpm) = 0.0_8, cvol(ndpm) = 0.0_8

!----------------------------------------------------------------------
!::Plasma surface
!----------------------------------------------------------------------
      real(8), allocatable :: rcps(:), csps(:), snps(:)
      character :: cpsf(ndwl)*6 = ' '
      integer :: npsp(ndwl) = 0, npep(ndwl) = 0, nbxp(ndwl) = 0
     >    , nbyp(ndwl) = 0
      integer, allocatable :: ikps(:), ipps(:), icps(:), isxp(:)
     >    , isyp(:), imps(:), imwl(:)
      integer :: nops = 0, npps = 0

!::[MPI_Send  in lnkntl_iniSnd]  cntcom/cntpsd/ (rcps,emrk)  10/04/21
!::[MPI_Recv  in lnkntl_iniSnd]  cntcom/cntpsf/ (rcps,emrk)  10/04/21

      real(8), allocatable :: eips(:), e0ps(:)
      real(8) :: e0em = 0.0_8, csem = 0.0_8, snem = 0.0_8

!::strike point
      integer :: ichti = 0, ichto = 0

!----------------------------------------------------------------------
!::neutral model
!----------------------------------------------------------------------
!::[MPI_Send  in lnkntl_iniSnd] cntcom/cntmdl/ (lntmd,emrk)  10/04/21
!::[MPI_Recv  in lnkntl_iniSnd] cntcom/cntmdl/ (lntmd,emrk)  10/04/21

      integer :: lntmd = 0, ldpmd = 0, lwlmd = 0, lrcmd = 0, lelmd = 0
     >    , lnnel = 0, lsrmd = 0

!----------------------------------------------------------------------
!::neutral parameter
!----------------------------------------------------------------------
      integer, parameter :: ndpf = 10
      real(8) :: rmas(ndgs) = 0.0_8, e0in = 0.0_8, v0in = 0.0_8
     >    , flxin = 0.0_8, wtmin = 0.0_8
      real(8) :: pfmag(ndpf) = 0.0_8, pfpx1(ndpf) = 0.0_8
     >    , pfpx2(ndpf) = 0.0_8, pfpy1(ndpf) = 0.0_8
     >    , pfpy2(ndpf) = 0.0_8, e0pfa = 0.0_8, e0pfm = 0.0_8
      real(8) :: temin_ion = 0.0_8, temin_rec = 0.0_8, tlmt_el = 0.0_8
      integer :: ngas = 0, noia(ndgs) = 0, nptldp = 0, nptlwl = 0
     >    , nptlvl = 0, lsmax = 0
      integer :: lemwl = 0, lemdp = 0, lemd2 = 0
      integer :: lempf = 1 ! Type of gas puff. =0:D0 3 ev inward, =1:D2
      integer :: ipfnp(ndpf) = 0, npf = 0

!::time depending puff rate setting
      logical :: use_gaspufftime = .false.
      integer, parameter :: pufftime_max = 20 ! max size of pufftime
      real*8 :: pufftime(pufftime_max) = 0.0_8 ! time when change puff rate
      real*8 :: pfmag_time(ndpf,pufftime_max) = 0.0_8 ! time depending puff rate
      integer :: pufftime_size = 0

!::input data
      real(8) :: cvrel(ndgs) = 0.0_8, cerel = 0.0_8
      real(8), allocatable :: vion(:,:)
      character :: cstyp*12 = ' '
      real(8), allocatable :: den0(:,:), eng0(:,:), vlp0(:,:)
     >    , deng(:,:), engg(:,:)
      real(8), allocatable :: nfl0x(:,:),nfl0y(:,:),nfl0z(:,:)

!----------------------------------------------------------------------
!::particle variables
!----------------------------------------------------------------------
      real(8) :: xpos = 0.0_8, ypos = 0.0_8, zpos = 0.0_8, velx = 0.0_8
     >    , vely = 0.0_8, velz = 0.0_8, velp = 0.0_8, vel = 0.0_8
     >    , vel2 = 0.0_8, evel = 0.0_8, tmas = 0.0_8, weit = 0.0_8
      integer :: icpo = 0, ixpo = 0, iypo = 0, izpo = 0, ilon = 0
     >    , ikon = 0, igas = 0, istp = 0, ievt = 0, ihit = 0, iptl = 0
      real(8) :: velxb = 0.0_8, velyb = 0.0_8, velzb = 0.0_8
     >    , velpb = 0.0_8, velb = 0.0_8, vel2b = 0.0_8, evelb = 0.0_8
     >    , tmasb = 0.0_8, weitb = 0.0_8
      integer :: igasb = 0
      integer :: is_atom_or_mole = 0 ! use in ntfolw, =0:atom. =1:molecular.
      integer :: ntmont_flg = 0 ! use in ntmont for ntfolw end flag

!----------------------------------------------------------------------
!::cross section table of ele-ionization
!----------------------------------------------------------------------
      real(8), allocatable :: tbsg_ei(:), tbel_ei(:)
     >    , tbsg_ds(:,:), tbel_ds(:,:), tbde_ds(:,:)

!----------------------------------------------------------------------
!::cross section
!----------------------------------------------------------------------
      integer, parameter :: ndrc = 7 * ndgs, ndrmx = 15

      character :: lbrc(ndrmx)*4 = ' ', ctrc(ndrmx)*80 = ' '
      real(8) :: dmfp = 0.0_8, trct = 0.0_8, trcx = 0.0_8, wrct = 0.0_8
     >    , frct(ndrc) = 0.0_8, srct(ndrc) = 0.0_8, elrc(ndrc) = 0.0_8
     >    , derc(ndrc) = 0.0_8
      integer :: igrc(ndrc) = 0, irct(ndrc) = 0, nrct = 0, mrct = 0
     >    , nrcmx = 0
      real(8) :: rcfc(ndrmx) = 0.0_8, scfc(ndrmx) = 0.0_8

!----------------------------------------------------------------------
!::gass puff
!----------------------------------------------------------------------
      real(8), allocatable :: pfflx(:)
      integer, allocatable :: ipfmp(:), ipfiw(:)
      integer :: ipfmx = 0

!::[MPI_Send  in lnkntl_iniSnd]  cntcom/cntpuf/ (pfflx,emrk)  10/04/21
!::[MPI_Recv  in lnkntl_iniSnd]  cntcom/cntpuf/ (pfflx,emrk)  10/04/21

!----------------------------------------------------------------------
!::birth profile
!----------------------------------------------------------------------
      integer :: kbr = 0, lbr = 0, nbr = 0
      integer, allocatable :: icbr(:), ikbr(:), isbr(:), iebr(:)
     >    , iwbr(:)
      real(8), allocatable :: csbr(:), snbr(:), eibr(:), e0br(:)
     >    , prbr(:), flbr(:)!, fmbr(:) was deleted
      real(8) :: tfbr = 0.0_8
      integer, allocatable :: iwpbr(:), iwpbr_wall(:)
      character :: cfbr*6 = ' '

!----------------------------------------------------------------------
!::wall information
!----------------------------------------------------------------------
      real*8 :: tflpm = 0.0_8, tflex = 0.0_8
      character :: wmate*2 = ' '

!----------------------------------------------------------------------
!::normal random   2 set (CX, EL)
!----------------------------------------------------------------------
      integer, parameter :: ndnrm = 4001, ndnrm2 = 4501
      real(8) :: tnor(ndnrm) = 0.0_8, fnor = 0.0_8
      real(8) :: tnor2(ndnrm2) = 0.0_8, fnor2 = 0.0_8
      integer :: iransd = 0

!----------------------------------------------------------------------
!::trigonometric function
!----------------------------------------------------------------------
      integer, parameter :: ndthe = 2001, ndfai = 4001
      real(8) :: fnth = 0.0_8, tcth(ndthe) = 0.0_8, tsth(ndthe) = 0.0_8
     >    , fnfi = 0.0_8, tcfi(ndfai) = 0.0_8, tsfi(ndfai) = 0.0_8
      integer :: nthe = 0, nfai = 0
      real(8) :: fvnth = 0.0_8, tvcth(ndthe) = 0.0_8
     >    , tvsth(ndthe) = 0.0_8

!----------------------------------------------------------------------
!::flag
!----------------------------------------------------------------------
      integer :: ntrc = 0, ltrc = 0, lscrg = 0
      character :: cftnr*80 = ' ', cftnw*80 = ' '
      integer :: nftnr = 0, nftnw = 0

!----------------------------------------------------------------------
!::analytical neutral transport code
!----------------------------------------------------------------------
      real(8) :: eioni = 0.0_8, eion = 0.0_8
      real(8) :: famp1 = 0.0_8, famp2 = 0.0_8, faref = 0.0_8
     >    , fawid = 0.0_8

!----------------------------------------------------------------------
!::density and temperature in void region
!----------------------------------------------------------------------
      real(8) :: void_ne = 0.0_8, void_ni = 0.0_8, void_te = 0.0_8
     >    , void_ti = 0.0_8

      end module cntcom
