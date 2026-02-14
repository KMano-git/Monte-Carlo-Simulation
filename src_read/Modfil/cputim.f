! added replace all include files with module files by kamata 2021/08/18
      module cputim
      use csize, only : ndmsr
      implicit none

!::SONIC
      real(8) :: cpustp = 0.0_8, cpubtm = 0.0_8, cpugrp(5) = 0.0_8
     >    , cputot = 0.0_8

!::NEUT2D (GROUP-2)
      real(8) :: cntstp = 0.0_8, cntsnd = 0.0_8, cntrev = 0.0_8
     >    , cntpls = 0.0_8
      real(8) :: cntcal(0:ndmsr) = 0.0_8, cntsum(0:ndmsr) = 0.0_8

!::init, exec  EL : elapsed time
      real(8) :: eltm00 = 0.0_8, eltm0 = 0.0_8, eltm1 = 0.0_8
     >    , eltm2 = 0.0_8, eltm3 = 0.0_8, eltm4 = 0.0_8, eltm5 = 0.0_8

      end module cputim
