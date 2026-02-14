! added replace all include files with module files by kamata 2021/08/18
      module cmeffz
      implicit none

!::impurity effect
      integer :: mion = 0
      real(8), allocatable :: wai(:), wza(:), wma(:), wfeb(:), wfab(:,:)
     >    , wkab(:,:), wfrab(:,:)

!::electron flux due to imputity
      real(8), allocatable :: vdnz(:,:,:,:), vdnz0(:,:,:,:), dzan(:,:)
     >    , vsez0(:,:), vsez(:,:), vsezz(:,:), vsezg(:,:), xzwtm(:,:)
     >    , xzwtp(:,:), xzflz(:,:)

      end module cmeffz
