! added replace cplcom from include file to module file by kamata 2021/08/18
!----------------------------------------------------------------------
!::numerical constant
!         phisical const (cpi,cev,cmp) move to cphcns
!----------------------------------------------------------------------
      module cplcom
      use csize, only : ndeq, ndmfl, ndsp, ndwl
      use cimcom, only : mrgn_max
      implicit none

      real(8) :: c00 = 0.0_8, c11 = 0.0_8, c12 = 0.0_8, c13 = 0.0_8
     >    , c14 = 0.0_8, c15 = 0.0_8, c16 = 0.0_8, c10 = 0.0_8
     >    , c21 = 0.0_8, c23 = 0.0_8, c25 = 0.0_8, c31 = 0.0_8
     >    , c32 = 0.0_8, c34 = 0.0_8, c43 = 0.0_8, c51 = 0.0_8
     >    , c52 = 0.0_8, c53 = 0.0_8, c54 = 0.0_8, c56 = 0.0_8
     >    , c58 = 0.0_8, c73 = 0.0_8, calf = 0.0_8
      real(8) :: clpa(ndsp) = 0.0_8, crod(ndsp,ndsp) = 0.0_8

!----------------------------------------------------------------------
!::ions species
!----------------------------------------------------------------------
      real(8) :: aion(ndsp) = 0.0_8, ama(ndsp) = 0.0_8
     >    , aza(ndsp) = 0.0_8
      integer :: nion = 0

!::[MPI_Bcast in ntl_ini]  cplcom/cmispc/ (aion,nion)  10/04/21
!::[MPI_Send  in lnkntl_ini]  cplcom/cmispc/ (aion,nion)  10/04/21
!::[MPI_Recv  in lnkntl_ini]  cplcom/cmispc/ (aion,nion)  10/04/21
!----------------------------------------------------------------------
!::initial profile
!----------------------------------------------------------------------
      real(8) :: exsx(5) = 0.0_8, exni(5) = 0.0_8, exvp(5) = 0.0_8
     >    , exte(5) = 0.0_8, exti(5) = 0.0_8
      real(8) :: fcna(ndsp) = 0.0_8
      integer :: nex = 0
!xxx      real*8    pni(ndx),pvp(ndx),pti(ndx),pte(ndx)

      real(8) :: nas0 = 0.0_8, tis0 = 0.0_8, tes0 = 0.0_8, nad0 = 0.0_8
     >    , tid0 = 0.0_8, ted0 = 0.0_8
!----------------------------------------------------------------------
!::decay length
!----------------------------------------------------------------------
      real(8), allocatable :: dmid(:), xmid(:), ymid(:)
      integer :: nmid = 0
!----------------------------------------------------------------------
!::aux. parameter  from q(N+1,l)
!----------------------------------------------------------------------
      real(8), allocatable :: vna(:,:,:), vne(:,:), vni(:,:)
     >    , vnezef(:,:)                   ! <== impurity effect
     >    , vcs(:,:), vea(:,:,:), vte(:,:), vti(:,:), vva(:,:,:)
     >    , vve(:,:), vzf(:,:) 
      real(8), allocatable :: anamp(:,:), anapv(:,:), anasl(:,:)
     >    , anemp(:), anepv(:), anesl(:), animp(:), anipv(:), anisl(:)
     >    , atemp(:), atepv(:), atesl(:), atimp(:), atipv(:), atisl(:)
     >    , fvpmp(:)
!::[MPI_Bcast in ntl_pls]  cplcom/cmauxv/ (vna,emrk)  10/04/21
!::[MPI_Send  in lnkimp_pls]  cplcom/cmauxv/ (vna,emrk)   10/04/21
!::[MPI_Recv  in lnkimp_pls]  cplcom/cmauxv/ (vna,emrk)   10/04/21
!::[MPI_Send  in lnkntl_pls]  cplcom/cmauxv/ (vna,emrk)  10/04/21
!::[MPI_Recv  in lnkntl_pls]  cplcom/cmauxv/ (vna,emrk)  10/04/21
!----------------------------------------------------------------------
!::aux. parameter at grid point (j+1/2,i+1/2)
!----------------------------------------------------------------------
      real(8), allocatable :: vnag(:,:,:), vteg(:,:), vtig(:,:)
     >    , vvag(:,:,:)

!**********************************************************************
!::move /cplpfl/ from monte/inc/cntmnt to solodr/inc/cplcom  2015/04/23
!-----------------------------------------------------------------------
!::plasma flux & volume source
!-----------------------------------------------------------------------
      real(8) :: tfps(ndwl,ndsp) = 0.0_8
      real(8) :: gfvl(ndsp) = 0.0_8, tfvl(ndsp) = 0.0_8
      real(8), allocatable :: flps(:,:), flvl(:,:,:)

!::[MPI_Bcast in ntl_pls]     cntmnt/cplpfl/ (tfps,emrk)  10/04/21
!::[MPI_Send  in lnkimp_pls]  cntmnt/cplpfl/ (tfps,emrk)  10/04/21
!::[MPI_Recv  in lnkimp_pls]  cntmnt/cplpfl/ (tfps,emrk)  10/04/21
!::[MPI_Send  in lnkntl_pls]  cntmnt/cplpfl/ (tfps,emrk)  10/04/21
!::[MPI_Recv  in lnkntl_pls]  cntmnt/cplpfl/ (tfps,emrk)  10/04/21
!**********************************************************************

!----------------------------------------------------------------------
!::conservative variables  q(N+1,l)
!----------------------------------------------------------------------
      real(8), allocatable :: q1a(:,:,:), q2a(:,:,:), q3(:,:), q4(:,:)
      real(8) :: qtim = 0.0_8, qtim_ini = 0.0_8

!::[MPI_Send  in lnkimp_pls]  cplcom/cmqnow/ (q1a,qtim)   10/04/21
!::[MPI_Recv  in lnkimp_pls]  cplcom/cmqnow/ (q1a,qtim)   10/04/21
!::[MPI_Bcast in lnkimp_pls]  cplcom/cmqnow/ (q1a,qtim)   10/04/21
!::[MPI_Bcast in pldisk]  cplcom/cmqnow/ (q1a,qtim)  10/04/21
!----------------------------------------------------------------------
!::conservative variables  q(N)
!----------------------------------------------------------------------
      real(8), allocatable :: wq1a(:,:,:), wq2a(:,:,:), wq3(:,:)
     >    , wq4(:,:)

!----------------------------------------------------------------------
!::variation of conv-variables   dlt-q = q(N+1,l+1) - q(N+1,l)
!----------------------------------------------------------------------
      real(8), allocatable :: dq1a(:,:,:), dq2a(:,:,:), dq3(:,:)
     >    , dq4(:,:)
!----------------------------------------------------------------------
!::collision frequency
!----------------------------------------------------------------------
      real(8) :: cfab(ndsp,ndsp) = 0.0_8, cfeb(ndsp) = 0.0_8
     >    , cfrab(ndsp,ndsp) = 0.0_8, ckab(ndsp,ndsp) = 0.0_8
      real(8) :: cthi = 0.0_8, cthe = 0.0_8
!----------------------------------------------------------------------
!::diffusion coefficient
!----------------------------------------------------------------------
      real(8) :: vlda(mrgn_max) = 0.0_8, vlet(6) = 0.0_8
     >    , vlxi(6) = 0.0_8, vlxe(6) = 0.0_8
      real(8) :: fcda(ndsp) = 0.0_8, fcet(ndsp) = 0.0_8
      real(8), allocatable :: vdda(:,:,:), vdet(:,:,:), vdxe(:,:)
     >    , vdxi(:,:)
      real(8) :: fdps = 0.0_8
      integer :: ldps = 0, ldluc = 0
      real(8) :: clzef = 0.0_8

      real(8), allocatable :: daan(:,:), etan(:,:), xian(:), xean(:)
     >    , dacl(:,:), etcl(:,:), xicl(:), xecl(:)
      real(8), allocatable :: dacfE(:,:), dacfW(:,:), ehcfE(:,:)
     >    , ehcfW(:,:), etcfE(:,:), etcfW(:,:), xecfE(:), xecfW(:)
     >    , xicfE(:), xicfW(:)
      real(8), allocatable :: dfl1a(:,:), dfl2a(:,:), dfl3(:), dfl4(:)
     >    , hea(:,:), hia(:,:), r1ae(:,:), r1ai(:,:), r1av(:,:)
     >    , r1aw(:,:), r2ai(:,:), r2av(:,:), r2aw(:,:), r3i(:), r4e(:)
     >    , xste(:,:), xsti(:,:), xsva(:,:) 
      real(8) :: temndf = 0.0_8, timndf = 0.0_8, timneta = 0.0_8
     >    , temn_dpsl = 0.0_8, temn_dprv = 0.0_8

      real(8), allocatable :: temndp(:)
      real(8) :: fdsl2 = 0.0_8
      integer :: itsl2 = 0

      integer :: mdl_bal = 0
      integer :: jst_bal = 0, jen_bal = 0
      integer :: mdl_bale = 0
      real(8) :: factor_bal = 0.0_8
      real(8) :: fdeg2 = 1.0 ! diffusion enhance factor for edge

      real(8), allocatable :: heat_flux_para_by_kappa_para(:,:)
     >    , heat_flux_para_elec(:,:), kappa_helander_test(:,:)
     >    , xicl_test(:,:)

      !************************************************!
      ! set diffusion coefficient by liner interporation along r_omd (in src/soldor/pldfan.f)
      !************************************************!
      ! number of data point. When lerp_points <= 0, upldif is ignored
      integer :: lerp_points_vldar = 0, lerp_points_vletr = 0
     >          ,lerp_points_vlxir = 0, lerp_points_vlxer = 0
      ! lerp_r    : position for interporation
      ! lerp_diff : diffusion coefficient value at lerp_r
      real(8) :: lerp_r_vldar(100)=0.0_8, lerp_diff_vldar(100)=0.0_8, 
     >           lerp_r_vletr(100)=0.0_8, lerp_diff_vletr(100)=0.0_8,
     >           lerp_r_vlxir(100)=0.0_8, lerp_diff_vlxir(100)=0.0_8,
     >           lerp_r_vlxer(100)=0.0_8, lerp_diff_vlxer(100)=0.0_8
      !************************************************!

! modified 1/1 lines dynamic allocation of arrays by kamata 2022/05/29
!!$omp threadprivate( daan, etan, xian, xean, dacl, etcl, etch, xicl
!$omp threadprivate( daan, etan, xian, xean, dacl, etcl, xicl
!$omp&             , xecl, dacfE, dacfW, etcfE, etcfW, ehcfE, ehcfW
!$omp&             , xicfE, xicfW, xecfE, xecfW, dfl1a, dfl2a, dfl3
!$omp&             , dfl4, hia, hea, xsva, xsti, xste, r1av, r2av, r1aw
! modified 1/1 lines dynamic allocation of arrays by kamata 2022/05/29
!!$omp&             , r2aw, r1ai, r2ai, r3i, r1ae, r2ae, r4e )
!$omp&             , r2aw, r1ai, r2ai, r3i, r1ae, r4e )
!----------------------------------------------------------------------
!::equi-partition
!----------------------------------------------------------------------
      real(8), allocatable :: alnr(:), feqp(:)

!$omp threadprivate( alnr, feqp )
!----------------------------------------------------------------------
!::flux limit
!----------------------------------------------------------------------
      real(8) :: flimi = 0.0_8, flime = 0.0_8, flimv = 0.0_8
      real(8), allocatable :: flmxe(:,:), flmxi(:,:), flmet(:,:)

!$omp threadprivate( flmxi, flmxe, flmet )
!----------------------------------------------------------------------
!::time
!----------------------------------------------------------------------
!
!    move (time,tend,dtim,lcstg) in cplcom-/cmtime/ ==> csonic-/cscntl/
!                                      ==> /cmdtim/
!
!        dtrateq(ndx,ndy,ndeq)
!                          mxx      mesh
!             (j,i,m1a)   density    (j,i)
!             (j,i,m2a)   momentum
!             (j,i,m3)    Ti
!             (j,i,m4)    Te
!
!         region   irg   (odiv,SOL,idiv,oprv,iprv,edge,core)
!
!        fine dtim just after IMPMC calculation   2013/01/07
!
!----------------------------------------------------------------------
      real(8) :: dtnxt = 0.0_8, dtmax = 0.0_8, dtmin = 0.0_8
     >    , dtimb = 0.0_8, dtlmt = 0.0_8
      real(8) :: dtrod = 0.0_8, dtrid = 0.0_8, dtrmp = 0.0_8
     >    , dtrsl = 0.0_8, dtrwl(10) = 0.0_8
      real(8) :: dtna = 0.0_8, dtva = 0.0_8, dtte = 0.0_8, dtti = 0.0_8
      real(8) :: dteq(ndeq) = 0.0_8
      real(8), allocatable :: dtrateq(:,:,:)
      real(8) :: dtprf = 0.0_8, dtdsk = 0.0_8, dtevl = 0.0_8
      real(8) :: tmtb(31) = 0.0_8, dttb(31) = 0.0_8, eltb(31) = 0.0_8
     >    , edtb(31) = 0.0_8
      real(8) :: rcydt = 0.0_8, rcysp = 0.0_8, rcytm = 0.0_8
      integer :: lttb(31) = 0, lttbr(31) = 0, nttb = 0, ittb = 0
     >    , itfix = 0
      integer :: itcn = 0

      real(8) :: dtfn(0:30) = 0.0_8
      integer :: ntfn = 0
!----------------------------------------------------------------------
!::time step control
!----------------------------------------------------------------------
      real(8) :: elpmx = 0.0_8, edtmx = 0.0_8, edtmn = 0.0_8
     >    , relmt = 0.0_8
      integer :: nlpmx = 0, nlpmn = 0, nlp = 0
      real(8) :: dlpq = 0.0_8, dlpq1 = 0.0_8, dlpq2 = 0.0_8
     >    , dlpq3 = 0.0_8, dlpq4 = 0.0_8, ddtq = 0.0_8, ddtq1 = 0.0_8
     >    , ddtq2 = 0.0_8, ddtq3 = 0.0_8, ddtq4 = 0.0_8, elpq = 0.0_8
     >    , elpq1 = 0.0_8, elpq2 = 0.0_8, elpq3 = 0.0_8, elpq4 = 0.0_8
     >    , edtq = 0.0_8, edtq1 = 0.0_8, edtq2 = 0.0_8, edtq3 = 0.0_8
     >    , edtq4 = 0.0_8
      real(8) :: rxflp(51) = 0.0_8

!----------------------------------------------------------------------
!::matrix
!----------------------------------------------------------------------
      real(8), allocatable :: cc(:,:,:), dd(:,:,:), ee(:,:,:), ff(:,:)
     >    , gg(:,:), ss(:,:,:)
!$omp threadprivate( cc, dd, ee, ss, ff, gg )
!----------------------------------------------------------------------
!::flux due to convection
!----------------------------------------------------------------------
      real(8), allocatable :: fl1a(:,:), fl2a(:,:), fl3(:), fl4(:)
     >    , fl4m(:), fl4p(:)
! modified 2/1 lines dynamic allocation of arrays by kamata 2022/05/29
!!$omp threadprivate( fl1a, fl2a, fl3, fl4, fl1ap, fl2ap, fl3p, fl4p
!!$omp&             , fl1am, fl2am, fl3m, fl4m )
!$omp threadprivate( fl1a, fl2a, fl3, fl4, fl4p, fl4m )
!----------------------------------------------------------------------
!::source term  (linearlization)
!----------------------------------------------------------------------
      real(8), allocatable :: ssnc(:,:,:), ssnv(:,:,:), sspc(:,:,:)
     >    , sspv(:,:,:), swec(:,:), swev(:,:), swic(:,:), swiv(:,:)

!----------------------------------------------------------------------
!::total neutral desnity (N0) send to IMPMC code
!----------------------------------------------------------------------
      real(8), allocatable :: tN0(:,:), tT0(:,:), tV0(:,:), tNg(:,:)
     >    , tTg(:,:)
      real(8) :: snflx(ndmfl) = 0.0_8

!----------------------------------------------------------------------
!::radiation of corona-model
!----------------------------------------------------------------------
!::[MPI_Send in lnkntl_ini]  cplcom/cmcrad/ /cmtrcy/ temp(30) 10/04/21
!::[MPI_Recv in lnkntl_ini]  cplcom/cmcrad/ /cmtrcy/ temp(30) 10/04/21
!::[MPI_Bcast in imp_ini]    cplcom/cmcrad/ /cmtrcy/ temp(30) 10/04/21
!::[MPI_Bcast in ntl_ini]    cplcom/cmcrad/ /cmtrcy/ temp(30) 10/04/21

      integer, parameter :: nwkmp_pl1 = 30  ! > 24  KSFUJI
      real(8), allocatable :: wdnz(:,:), wrad(:,:)
      real(8) :: tradrg(10) = 0.0_8
      real(8) :: trad = 0.0_8, trad_min = 0.0_8, trad_max = 0.0_8
     >    , trad_fac = 0.0_8
      real(8) :: cimp = 0.0_8, cimprg(10) = 0.0_8, cimprg2(10) = 0.0_8

!::ionization loss of impurities
      real(8) :: tlosrg(10) = 0.0_8, tlos = 0.0_8
      real(8) :: crad_totwe(10) = 0.0_8

      integer :: xplm_opt = 0, xplm_jc = 0, xplm_ic = 0
      real(8) :: xplm_ne = 0.0_8, xplm_te = 0.0_8
!----------------------------------------------------------------------
!::TVD scheme & up-wind scheme for y-diffusion
!----------------------------------------------------------------------
      real(8) :: fai = 0.0_8, faip = 0.0_8, faim = 0.0_8, beta = 0.0_8
      integer :: lordr = 0, lmuscl = 0
      real(8), allocatable :: xwtam(:,:,:), xwtap(:,:,:), xfdfa(:,:,:)
      integer, allocatable :: iupwd(:,:)
      integer :: lupwd = 0

!----------------------------------------------------------------------
!::Boundary condition
!----------------------------------------------------------------------
      real(8) :: cdli = 0.0_8, cdle = 0.0_8, gcsi_sl = 0.0_8
     >    , gcsi_pv = 0.0_8, gcse_sl = 0.0_8, gcse_pv = 0.0_8
      real(8) :: fcbcvl = 0.0_8
      real(8), allocatable :: gcse(:), gcsi(:)
      real(8) :: gwsni = 0.0_8, gwsti = 0.0_8, gwste = 0.0_8
     >    , bwsni = 0.0_8, bwsti = 0.0_8, bwste = 0.0_8, gwpni = 0.0_8
     >    , gwpti = 0.0_8, gwpte = 0.0_8, bwpni = 0.0_8, bwpti = 0.0_8
     >    , bwpte = 0.0_8, bmpni = 0.0_8, bmpti = 0.0_8, bmpte = 0.0_8
      real(8) :: bmpfn(ndsp) = 0.0_8, bwpfn(ndsp) = 0.0_8
     >    , bwsfn(ndsp) = 0.0_8
      integer :: lbcgd = 0, lbcsw = 0, lbcpw = 0
      integer, allocatable :: ibcyl(:)
      real(8), allocatable :: qb1a(:,:,:), qb2a(:,:,:), qb3(:,:)
     >    , qb4(:,:)
      real(8) :: qw1a(ndsp,4) = 0.0_8, qw2a(ndsp,4) = 0.0_8
     >    , qw3(4) = 0.0_8, qw4(4) = 0.0_8
!cKH
      integer :: lfbbc = 0
!----------------------------------------------------------------------
!::particle flux & energy flux from the plasma edge
!----------------------------------------------------------------------
      real(8) :: tflna(ndsp) = 0.0_8
      real(8) :: tflni = 0.0_8, tflne = 0.0_8, tflqi = 0.0_8
     >    , tflqe = 0.0_8
!----------------------------------------------------------------------
!::file 
!----------------------------------------------------------------------
      integer :: nftr = 0, nftw = 0, nfth = 0, nftp = 0, nfti = 0
      character :: cftr*80 = ' ', cftw*80 = ' ', cfth*80 = ' '
     >    , cftp*80 = ' ', cfti*80 = ' '
!----------------------------------------------------------------------
!::Flags
!----------------------------------------------------------------------
      character :: cfopt(100)*80 = ' ', cfdbg(100)*80 = ' '
      character :: caim*80 = ' ', cprm*80 = ' '
      integer :: mdl_wrd = 0, mdl_eqp = 0, mdl_edt = 0, mdl_srw = 0
     >    , mdl_srp = 0, mdl_vis = 0, mdl_ini = 0, mdl_hcv = 0
     >    , mdl_cgen = 0, mdl_fimp = 0
      real(8) :: temin_aux = 0.0_8, timin_aux = 0.0_8, timin_vis = 0.0_8
     >    , nimin_aux = 0.0_8
      real(8) :: chfpr = 0.0_8, chfcv = 0.0_8, cevpr = 0.0_8
      real(8) :: clmdnz = 0.0_8
      integer :: nsmp_wrad = 0
      real(8) :: wfac_lv = 0.0_8

      end module cplcom
