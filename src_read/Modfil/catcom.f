! added replace all include files with module files by kamata 2021/08/18
      module catcom
      use csize, only : ndmis
      implicit none
!xx      tbcatmz = (/"C",    "Ne",   "Ar",   "Kr",  "Xe"/)
!xx      tbnatmz = (/ 6,      10,     18,     36,    54/)

!      ndxrw =    4    :    rank+1  <sig>,Ne,Te,Iz
!      ndxkn =    4    :    ion, rec, cxr, rad, rad, rad
!      ndxx1 =   24    :    N(Ne)
!      ndxx2 =   30    :    N(Te)
!      ndxdt = 4320x3  :    N(Ne)*N(Te)*N(Iz) = 24*30*6 = 4320
!      ndxzi =    6    :    N(Iz) < 100 

!     integer, parameter :: ndxzi = ndmis  ! Argon case
      integer, parameter :: ndxrw =  4
      integer, parameter :: ndxkn =  6
      integer, parameter :: ndxx1 = 41  ! <== <24> 2009/01/27 14/01/22
      integer, parameter :: ndxx2 = 48  ! <== <30> 2009/01/27
!xx   integer, parameter :: ndxzi = 54  ! <== <18> 2012/01/10

! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   integer, parameter :: ndxdt = ndxx1 * ndxx2 * ndxzi * ndxkn
      integer, parameter :: ndxdt = ndxx1 * ndxx2 * ndmis * ndxkn

! modified 11/11 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   real(8)    xdaty(ndxdt), xdtx1(ndxx1,ndxkn), xdtx2(ndxx2,ndxkn)
!ik   real(8)    xdatz(ndxzi,ndxkn)
!ik   real(8)    xwmin(ndxrw,ndxkn), xwmax(ndxrw,ndxkn)
!ik  >         , xwmlt(ndxrw,ndxkn)
!ik   integer    ip_ion, ip_rec, ip_cxr, ip_rad, ip_plt, ip_prb, ip_prc
!ik   integer    nxs, nyp
!ik   integer    xdrnk(ndxkn), xwnum(ndxrw,ndxkn)
!ik   integer    xjsta(ndxkn), xjend(ndxkn), xjnum(ndxkn)
!ik   character  xdsn(ndxkn)*80, xdnam(ndxkn)*20, xdevl(ndxkn)*20
!ik   character  xwvar(ndxrw,ndxkn)*12, xwunt(ndxrw,ndxkn)*12
!ik   character  xwspc(ndxrw,ndxkn)*12
      real(8), allocatable :: xdaty(:,:), xdatz(:,:,:), xdtx1(:,:,:)
     >    , xdtx2(:,:,:), xwmax(:,:,:), xwmin(:,:,:), xwmlt(:,:,:)
      integer, allocatable :: ip_cxr(:), ip_ion(:), ip_plt(:), ip_prb(:)
     >    , ip_prc(:), ip_rad(:), ip_rec(:)
     >    , nxs(:), nyp(:)
     >    , xdrnk(:,:), xwnum(:,:,:)
     >    , xjsta(:,:), xjend(:,:), xjnum(:,:)
      character, allocatable :: xdsn(:,:)*80, xdnam(:,:)*20
     >    , xdevl(:,:)*20
     >    , xwvar(:,:,:)*12, xwunt(:,:,:)*12
     >    , xwspc(:,:,:)*12

! modified 2/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   integer    prk(ndxkn), pin(ndxrw,ndxkn), plg(ndxrw,ndxkn)
!ik   real(8)    pxi(2,ndxkn), pxa(2,ndxkn), pxs(2,ndxkn), pxd(2,ndxkn)
      integer, allocatable :: prk(:,:), pin(:,:,:), plg(:,:,:)
      real(8), allocatable :: pxi(:,:,:), pxa(:,:,:), pxs(:,:,:)
     >    , pxd(:,:,:)

!-----
! modified 4/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   character(80)    dratm
!ik   character(2)     catmz
!ik   integer          natmz
!ik   real(8)          aionz
      character, allocatable :: catmz(:)*2, dratm(:)*80
      integer,   allocatable :: natmz(:)
      real(8),   allocatable :: aionz(:)

      integer :: nfatm = 0
      character(3),  dimension(7)  :: rtyp = ' ', ched = ' '
      character(20), dimension(10) :: tbdsn = ' '

!::fine grid points for irregular spacing of C96 data 
      integer, parameter :: ndfgd = 1001
      integer :: pfg(ndfgd,2,ndxkn) = 0, pfn(2,ndxkn) = 0
      real(8) :: pfs(2,ndxkn) = 0.0_8, pfe(2,ndxkn) = 0.0_8
     >    , pfd(2,ndxkn) = 0.0_8

      end module catcom
! Note :
!     catcom_2, catcom_3 : The following data is not available
!       ndfgd, pfd, pfe, pfg, pfn, pfs
