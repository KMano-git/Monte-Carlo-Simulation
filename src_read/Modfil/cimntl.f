      module cimntl
      implicit none

!::local variables of neutral particle
!      integer ndmp2; parameter (ndmp2=20000)  ! OLD PSI
!      integer ndmp2; parameter (ndmp2=40000)  ! NEW IAEA
      integer, parameter :: ndmp2 = 500000  ! NEW Timedep or NEW TD
      real(8) :: sx0(ndmp2) = 0.0_8, sy0(ndmp2) = 0.0_8
     >    , sz0(ndmp2) = 0.0_8
      real(8) :: svx(ndmp2) = 0.0_8, svy(ndmp2) = 0.0_8
     >    , svz(ndmp2) = 0.0_8
      real(8) :: srn(ndmp2) = 0.0_8, sint(ndmp2) = 0.0_8
     >    , stm(ndmp2) = 0.0_8, stb(ndmp2) = 0.0_8

!::scoreing   No use
      integer, parameter :: ndsc = 2001
      integer :: nsc = 0, scic(ndsc) = 0
      real(8) :: scdt(ndsc) = 0.0_8

      end module cimntl
