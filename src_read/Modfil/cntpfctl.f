! added by kamata 2021/11/18
      module cntpfctl
      use cntcom, only : ndpf
      implicit none

      ! input parameter
      integer :: lpfctl = 0     ! # of controled puff, 0 is nothing to be done
      integer :: lpfctl_rst = 0
      real(8) :: nsep_pfctl = 0.0_8 ! target nsep at outer midplane
      real(8) :: fac_pfctl = 0.0_8  ! strength of feedback
      real(8) :: rst_pfctl = 0.0_8
      real(8) :: fpfctl(ndpf) = 0.0_8

      end module cntpfctl
