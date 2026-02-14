! added replace all include files with module files by kamata 2021/08/18
      module cplimp
!     use csize, only : ndgt, ndmc, ndis2L => ndmis, ndwp, nzsmx
      implicit none

!     integer, parameter :: ndis2L = ndmis
!      integer, parameter :: nzsmx = 3    ==> csize

!:: cplimp_imcom
! cimcom_1      
! modified 2/2 lines treat 4 or more impurities with IMPMC by kamata 2022/05/08
!ik   real(8), dimension(nzsmx) :: aimasL, azmasL, amiL, amzL
!ik   integer    ismaxL(nzsmx)
      real(8), allocatable :: aimasL(:), azmasL(:), amiL(:), amzL(:)
      integer, allocatable :: ismaxL(:)

! cimcom_8      
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/05/08
!ik   real(8)    cdifL(nzsmx), cdifrgL(10,nzsmx)
      real(8), allocatable :: cdifL(:) !, cdifrgL(:,:)

! cimcom_12b
! modified 2/2 lines treat 4 or more impurities with IMPMC by kamata 2022/05/08
!ik   real(8)    denFlxZtoWallL(0:ndis2L,ndwp,nzsmx)
!ik   real(8)    eneFlxZtoWallL(0:ndis2L,ndwp,nzsmx)
!     real(8), allocatable :: denFlxZtoWallL(:,:,:)
!    >    , eneFlxZtoWallL(:,:,:)

!:: cplimp_imden from cimden
! modified 9/7 lines treat 4 or more impurities with IMPMC by kamata 2022/05/08
!ik   integer    nsputL(nzsmx), nzmxL(nzsmx)
!ik   integer    ncmxL(nzsmx), npsputL(5,nzsmx)
!ik   real(8)    fsputL(5,nzsmx), wtsputL(5,nzsmx) 
!ik   real(8)    bksptL(ndgt,nzsmx)
!ik   real(8)    eipotL(0:ndis2L,nzsmx)
!ik   real(8)    twrdL(ndmc,nzsmx), tdnzL(0:ndis2L,ndmc,nzsmx)
!ik   real(8)    twciL(ndmc,nzsmx)
!ik   real(8)    twneL(ndmc,nzsmx)  ! Ne at IMPMC cal 
!ik   character  csputL(5,nzsmx)*4, nwsptL(nzsmx)*4
!     integer, allocatable :: ncmxL(:), npsputL(:,:), nsputL(:)
      integer, allocatable :: nzmxL(:)
!     real(8), allocatable :: bksptL(:,:), fsputL(:,:), twneL(:,:)
!    >    , twrdL(:,:), wtsputL(:,:)
      real(8), allocatable :: eipotL(:,:), tdnzL(:,:,:), twciL(:,:)
! twneL :  Ne at IMPMC cal 
!     character, allocatable :: csputL(:,:)*4, nwsptL(:)*4
      real(8), allocatable :: tengzL(:,:,:), tprzL(:,:,:) ! energy and pressure of impurity

! modified 8/3 lines treat 4 or more impurities with IMPMC by kamata 2022/05/08
!ik   real(8)    tfrzL(0:ndis2L,ndmc,nzsmx)
!ik   real(8)    tthzL(0:ndis2L,ndmc,nzsmx)
!ik   real(8)    tvlzL(0:ndis2L,ndmc,nzsmx)
!ik   real(8)    tionZL(0:ndis2L,ndmc,nzsmx)
!ik   real(8)    trecZL(0:ndis2L,ndmc,nzsmx)
!ik   real(8)    tradiZL(0:ndis2L,ndmc,nzsmx)
!ik   real(8)    tradliZL(0:ndis2L,ndmc,nzsmx)
!ik   real(8)    tradrZL(0:ndis2L,ndmc,nzsmx)
!     real(8), allocatable :: tfrzL(:,:,:), tionZL(:,:,:)
!    >    , tradiZL(:,:,:), tradliZL(:,:,:), tradrZL(:,:,:)
!    >    , trecZL(:,:,:), tthzL(:,:,:), tvlzL(:,:,:)

      integer :: wmc_nty = 0

      end module cplimp
