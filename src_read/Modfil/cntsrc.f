! added replace all include files with module files by kamata 2021/08/18
      module cntsrc
      use csize, only : ndsp
      implicit none

!::local common (plntsr_plt.f)
      real(8), allocatable :: tden0(:,:), teng0(:,:), tvlp0(:,:)
     >    , tdeng(:,:), tengg(:,:), tssn(:,:), tssp(:,:), tswi(:)
     >    , tswe(:)
      real(8) :: trgsn(10,ndsp) = 0.0_8, trgsp(10,ndsp) = 0.0_8
     >    , trgwi(10) = 0.0_8, trgwe(10) = 0.0_8

      end module cntsrc
