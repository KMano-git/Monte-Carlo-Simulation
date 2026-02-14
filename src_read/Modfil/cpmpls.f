! added replace all include files with module files by kamata 2021/08/18
!----------------------------------------------------------------------
!::main plasma parameter
!----------------------------------------------------------------------
      module cpmpls
      use csize, only : ndsp
      implicit none

!::[MPI_Send  in lnkntl_ini]  cpmpls/cpmain/ (flxni,emrk) 10/04/21
!::[MPI_Recv  in lnkntl_ini]  cpmpls/cpmain/ (flxni,emrk) 10/04/21
!::[MPI_Bcast in ntl_ini]     cpmpls/cpmain/ (flxni,emrk) 10/04/21

      real(8) :: flxni = 0.0_8, flxqi = 0.0_8, flxqe = 0.0_8
      real(8) :: fflna(ndsp) = 0.0_8, prfna(10,ndsp) = 0.0_8
     >    , prfni(10) = 0.0_8, prfti(10) = 0.0_8, prfte(10) = 0.0_8
     >    , fprna(ndsp) = 0.0_8, r0mp = 0.0_8
     >    , ramp = 0.0_8
     
      real(8), allocatable :: armp(:), arhmp(:), romp(:), rohmp(:)
     >    , vlmp(:), sfmp(:), dvmp(:), drmp(:), wdmp(:), wfmp(:)
     >    , wlmp(:)

      real(8) :: xmd1 = 0.0_8, ymd1 = 0.0_8, xmd2 = 0.0_8, ymd2 = 0.0_8
      integer :: jmd1 = 0, imd1 = 0, jmd2 = 0, imd2 = 0
      integer :: lnedg = 0, ltedg = 0

      real(8) :: bniedg = 0.0_8

      real(8) :: rmds(2) = 0.0_8, dn0s(2) = 0.0_8, pufs = 0.0_8
      real(8) :: rmdc(2) = 0.0_8, dn0c(2) = 0.0_8, taup = 0.0_8
     >    , wlosi = 0.0_8, wlose = 0.0_8
      real(8), allocatable :: an0mp(:), asnmp(:)

      real(8) :: clsav = 0.0_8, clsxp = 0.0_8, clsfc = 0.0_8
     >    , clsni = 0.0_8, clsti = 0.0_8, clste = 0.0_8
      real(8), allocatable :: alsfr(:)

      real(8), allocatable ::  mdp_ro(:), mdp_rh(:)
     >    , mdp_ni(:), mdp_ti(:), mdp_te(:)

!::[MPI_Bcast in cpmpls]  cpmpls/com_plmdprf/ (mdp_ro,emrk)  10/04/21

!::diffusion in main plasma from topics code
      real(8), allocatable :: adamp(:,:), avpmp(:,:), ahpmp(:,:)
     >    , aximp(:), axemp(:)

!::radiation power constraint 
!  nty_ctl: IMPMC number to be controlled for radiation power, like Ar.
!  rad_const_step: force limit on radiant power from the start point to rad_const_step.
!  rad_const_fac: constraint = (flxqi+flxqe)*rad_const_fac
      integer :: nty_ctl = 1, rad_const_step = -1
      real*8 :: rad_const_fac = 1.0d0

      end module cpmpls
