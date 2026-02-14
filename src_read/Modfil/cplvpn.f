! added replace all include files with module files by kamata 2021/08/18
      module cplvpn
      use csize, only : ndeq
      implicit none

!::input
      real(8), allocatable :: vpn_sol(:), vpn_man(:), vlvp(:), vdvp(:,:)
      real(8), allocatable :: qfy_vp(:,:,:)
      real(8) :: qsx_vp(ndeq) = 0.0_8

      end module cplvpn
