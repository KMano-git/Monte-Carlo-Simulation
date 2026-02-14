! added replace all include files with module files by kamata 2021/08/18
      module cplhst
      implicit none

!-----------------------------------------------------------------------
!::history data
!-----------------------------------------------------------------------
      integer, parameter :: ndht = 500, ndhp = 20, ndhv = ndhp*4

      character :: cvnm(ndhp)*20 = ' '
      integer :: iphs(ndhp) = 0, jphs(ndhp) = 0, nhpn = 0, nhtm = 0
     >    , nhvl = 0
      real(4) :: hstt(ndht) = 0.0_4, hstv(ndht,ndhv) = 0.0_4

!-----------------------------------------------------------------------
!::profile data
!-----------------------------------------------------------------------
      integer, parameter :: ndpr = 10
      integer :: icpr(ndpr) = 0, ncpr = 0

      end module cplhst
