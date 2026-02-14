! added replace all include files with module files by kamata 2021/08/18
      module cntxwk
      implicit none

      integer, parameter :: ndhfl = 20, ndhun = 50  ! smothing,(off smothing:ndhfl=1, nxksz=1)
      integer :: ihfl = 0, jhfl = 0
      integer :: nhno(ndhfl) = 0, nhtm(ndhun,ndhfl) = 0
      real(8), allocatable :: xmwork(:,:,:)

      end module cntxwk
