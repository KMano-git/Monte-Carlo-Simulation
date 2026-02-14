! added replace all include files with module files by kamata 2021/08/18
!-----------------------------------------------------------------------
!::cal results of monte-carlo
!::KSFUJI   vwork(nwksz,ndmsr) ==> vmwork(nwkmp_nt,ndmsr)  10/04/21
!-----------------------------------------------------------------------
      module cntmnt
      use csize, only : ndmfl, ndptl, ndsp
      implicit none

      real(8) :: vtime = 0.0_8, vflux(ndmfl) = 0.0_8
      real(8), allocatable :: vmwork(:,:)
      integer :: vitim = 0, vitnt = 0, vnsrc = 0, vsmty = 0, vsmno = 0
      integer :: visrc(ndmfl) = 0, vkflx(ndmfl) = 0, vksty(ndmfl) = 0
     >    , vnsmp(ndmfl) = 0
      character :: vcsrc(ndmfl)*6 = ' '
      character :: cntmnt_emrk*1 = ' '

!::[MPI_Bcast in ntdisk]  cntmnt /cntmnt/ (vtime,emrk)  10/04/21
!::[MPI_Recv  in lnkntl_cal]  cntmnt/cntmnt/ (vtime,emrk)  10/04/21
!::[MPI_Send  in lnkntl_cal]  cntmnt/cntmnt/ (vtime,emrk)  10/04/21

!-----------------------------------------------------------------------
!::cal. condition of monte-carlo  (initial set)
!-----------------------------------------------------------------------
      integer :: mfmax = 0, mnknd(ndmfl) = 0, mnsmp(ndmfl) = 0
     >    , mnvwk(ndmfl) = 0
      character :: mcflx(ndmfl)*6 = ' '
!**********************************************************************

!-----------------------------------------------------------------------
!::particle conservation
!-----------------------------------------------------------------------
      real(8) :: pfl_ion(ndmfl) = 0.0_8, pfl_src(ndmfl) = 0.0_8
     >    , pfl_man(ndmfl) = 0.0_8, pfl_abs(ndmfl) = 0.0_8
     >    , pfl_pmp(ndmfl) = 0.0_8, pfl_err(ndmfl) = 0.0_8
     >    , psm_ion = 0.0_8, psm_ntl = 0.0_8, psm_src = 0.0_8
     >    , psm_man = 0.0_8, psm_abs = 0.0_8, psm_pmp = 0.0_8
     >    , psm_err = 0.0_8, psm_Fl = 0.0_8, psm_Si = 0.0_8
     >    , psm_Dt = 0.0_8
! added 1 line pfl_ntl(vkflx(i)),vkflx(i)=0 in ntbgprf.f  by kamata 2021/08/18
      real(8) :: pfl_ntl(0:ndmfl) = 0.0d0

      integer :: n6_src = 0, i6_src = 0

!-----------
!::cntsrc
!-----------
      real(8), allocatable :: sn0(:,:,:), ssn(:,:,:), ssp(:,:,:)
     >    , swi(:,:), swe(:,:)
      real(8) :: tnexh = 0.0_8, tnpmp = 0.0_8, tnpuf = 0.0_8
     >    , flexh = 0.0_8, flpmp = 0.0_8, totsn(10,ndsp) = 0.0_8
     >    , totsp(10,ndsp) = 0.0_8, totwe(10) = 0.0_8, totwi(10) = 0.0_8
     >    , sdotn2 = 0.0_8, sdotn = 0.0_8, sumsn(10,ndsp) = 0.0_8
     >    , sumsp(10,ndsp) = 0.0_8, sumwe(10) = 0.0_8, sumwi(10) = 0.0_8
     >    , dotn2 = 0.0_8, dotn = 0.0_8

!::total source in hot core region
      real(8) :: totmsn(ndsp) = 0.0_8, totmsp(ndsp) = 0.0_8
     >    , totmwe = 0.0_8, totmwi = 0.0_8

!-----------------------------------------------------------------------
!::random seed
!-----------------------------------------------------------------------
      integer :: iseed(ndptl) = 0, nsmpmx = 0

      end module cntmnt
