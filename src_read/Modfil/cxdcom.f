! added replace all include files with module files by kamata 2021/08/18
      module cxdcom
      implicit none

      integer, parameter :: ndxrw = 3, ndxkn = 20, ndxdt = 80000

      integer :: nxs = 0, nyp = 0
      real(8) :: xdaty(ndxdt) = 0.0_8
      character :: xdnam(ndxkn)*20 = ' ', xdevl(ndxkn)*20 = ' '
      integer :: xdrnk(ndxkn) = 0, xwnum(ndxrw,ndxkn) = 0
      integer :: xjsta(ndxkn) = 0, xjend(ndxkn) = 0, xjnum(ndxkn) = 0
      character :: xwvar(ndxrw,ndxkn)*40 = ' '
     >    , xwunt(ndxrw,ndxkn)*12 = ' ', xwspc(ndxrw,ndxkn)*20 = ' '
      real(8) :: xwmin(ndxrw,ndxkn) = 0.0_8, xwmax(ndxrw,ndxkn) = 0.0_8
     >    , xwmlt(ndxrw,ndxkn) = 0.0_8
      integer :: prk(ndxkn) = 0, pin(3,ndxkn) = 0, plg(3,ndxkn) = 0
      real(8) :: pxi(2,ndxkn) = 0.0_8, pxa(2,ndxkn) = 0.0_8
     >    , pxs(2,ndxkn) = 0.0_8, pxd(2,ndxkn) = 0.0_8

!::[MPI_Bcast in elinit]   cxdcom (nxs,emrk)  10/04/21

      end module cxdcom
