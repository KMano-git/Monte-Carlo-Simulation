! added replace all include files with module files by kamata 2021/08/18
      module czwflx
      implicit none

! modified 4/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   real(8)    zw_pflx(0:ndis,ndwp), zw_eflx(0:ndis,ndwp)
!ik   real(8)    zw_pflx2(0:ndis,ndwp), zw_eflx2(0:ndis,ndwp)
!ik   integer    zw_ismax, zw_ismax2
!ik   character(2)  zw_catmz, zw_catmz2
      real(8), allocatable :: zw_pflx(:,:,:), zw_eflx(:,:,:)
      integer, allocatable :: zw_ismax(:)
      character(2), allocatable :: zw_catmz(:)

      end module czwflx
