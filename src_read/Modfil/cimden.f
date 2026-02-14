! added replace all include files with module files by kamata 2021/08/18
      module cimden
      use csize, only : ndgt, ndmis
      implicit none

      integer :: nsput = 0, nzmx = 0, ncmx = 0, npsput(5) = 0
      real(8) :: fsput(5) = 0.0_8, wtsput(5) = 0.0_8
      real(8) :: bkspt(ndgt) = 0.0_8
      real(8) :: eipot(0:ndmis) = 0.0_8
      real(8), allocatable :: twrd(:), tdnz(:,:), twci(:), twne(:)  ! Ne at IMPMC cal 
      real(8), allocatable :: tengz(:,:), tprz(:,:) ! energy and pressure of impurity

! SYamoto: Friction, thermal, flowvel, ion. src., rec. src, respectively
      real(8), allocatable :: tfrz(:,:), tthz(:,:), tvlz(:,:)
     >    , tionZ(:,:), trecZ(:,:)

! SYamoto: Cooling rate
      real(8), allocatable :: tradiZ(:,:), tradliZ(:,:), tradrZ(:,:)
      character :: csput(5)*4 = ' ', nwspt*4 = ' '

!::[MPI_Recv in lnkimp_cal]   cimden/cimden/ (nsput,emrk)  10/04/21
!::[MPI_Send in lnkimp_cal]   cimden/cimden/ (nsput,emrk)  10/04/21
!::[MPI_Bcast in lnkimp_cal]  cimden/cimden/ (nsput,emrk)  10/04/21

! only IMPMC_TD :
      integer :: nsrc_spt = 0, nsrc_spt_dummy = 0

      end module cimden
! Note : difference between IMPMC and IMPMC_TD
! only IMPMC_TD :
! nsrc_spt
