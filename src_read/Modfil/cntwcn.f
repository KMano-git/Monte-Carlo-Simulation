! added replace all include files with module files by kamata 2021/08/18
      module cntwcn
      implicit none

      real(8) :: wflx = 0.0_8, wend(10) = 0.0_8, wtot = 0.0_8
      real(8), allocatable :: wssn(:,:), wssp(:,:), wswe(:), wswi(:)
     >    , wsbr(:,:), wden(:,:), weng(:,:), wvlp(:,:), wtion(:,:)
     >    , wgden(:,:), wgeng(:,:), wemt(:), wwal(:)
      real(8) :: whta(30,3) = 0.0_8, whtm(30,3) = 0.0_8
      real(8), allocatable :: wfhta(:), wfhtm(:), wehta(:), wehtm(:)
      real(8) :: wreg(0:10) = 0.0_8, wsum = 0.0_8, wnrm = 0.0_8
     >    , wion = 0.0_8, wabs = 0.0_8, wpmp = 0.0_8, werr = 0.0_8

! 3d neutral flow (weight), toku 2016.06.23
      real(8), allocatable :: wnfl0x(:,:), wnfl0y(:,:), wnfl0z(:,:)
     >    , wnflgx(:,:), wnflgy(:,:), wnflgz(:,:)
      real(8), allocatable :: weema(:), weemm(:) ! 200513KH power(we) of reflected(em) atom/mole
      real(8), allocatable :: wees(:) ! 201021KH power(we) of emitted atom from plasma surface (sol/prv)

!::[MPI_Reduce in ntmont]  cntwcn/cntvar/ (wflx,werr)  10/04/21
!::[MPI_Bcast  in ntmont]  cntwcn/cntvar/ (wflx,werr)  10/04/21

!----------------------------------------------------------------------
!xxx  real(8)    wwork(nwksz); equivalence (wflx,wwork(1))  KSFUJI
      integer :: wisrc = 0
      character :: wcsrc*6 = ' ', wcest*6 = ' '

!----------------------------------------------------------------------
      real(8) :: swreg(0:10) = 0.0_8, swsum = 0.0_8, swnrm = 0.0_8
     >    , swion = 0.0_8, swabs = 0.0_8, swpmp = 0.0_8, swerr = 0.0_8

      end module cntwcn
