! added replace all include files with module files by kamata 2021/08/18
      module cntsms
      implicit none

!::local common
      integer, parameter :: ndsm = 1000
      integer :: i6_ntsms = 0
      integer, dimension(ndsm,9) :: icsm = 0
      real(8), dimension(ndsm,9) :: wtsm = 0.0_8
      integer, dimension(ndsm)   :: ncsm = 0, itsm = 0
      integer :: nsmax = 0

      end module cntsms
