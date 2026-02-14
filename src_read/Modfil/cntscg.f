! added replace all include files with module files by kamata 2021/08/18
!-----------------------------------------------------------------------
!::scoreing f:fluid/t:track-length/c:collision estimator
!-----------------------------------------------------------------------
      module cntscg
      use cntcom, only : ndrmx
      implicit none

      real(8), allocatable :: sn_mt(:,:), sp_mt(:,:), we_mt(:), wi_mt(:)
     >    , sn_mc(:,:), sp_mc(:,:), we_mc(:), wi_mc(:)

!::total source in region(10) for each reactions(ndrmx)
      real(8), dimension(ndrmx,10) :: tsn_mt = 0.0_8, tsp_mt = 0.0_8
     >    , twe_mt = 0.0_8, twi_mt = 0.0_8
      real(8), dimension(ndrmx,10) :: tsn_mc = 0.0_8, tsp_mc = 0.0_8
     >    , twe_mc = 0.0_8, twi_mc = 0.0_8

!::size
      integer, parameter :: nwkmp_sc = ndrmx * 10 * 8

      end module cntscg
