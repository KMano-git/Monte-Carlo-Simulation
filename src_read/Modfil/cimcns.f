! added replace all include files with module files by kamata 2021/08/18
      module cimcns
      implicit none

!::normal random number
      integer, parameter :: ndghf = 1000, ndgau = 2 * ndghf + 1
      integer :: ngau = 0
      real(8) :: gaus0(ndgau) = 0.0_8, fgau = 0.0_8

!::normal random number in r-space
      integer, parameter :: ndgar = 1001
      integer :: ngar = 0
      real(8) :: gaus1(ndgar) = 0.0_8, fgar = 0.0_8

!::tompson model (sputtered velocity)
      integer, parameter :: ndtom = 2001
      integer :: ntom = 0
      real(8) :: tomran(ndtom) = 0.0_8, tomeng(ndtom) = 0.0_8
     >    , tomvel(ndtom) = 0.0_8, ftom = 0.0_8

!::seed of random number
      integer, parameter :: ndseed = 3000
      integer :: iseed(ndseed) = 0, nseed = 0

      end module cimcns
