! added replace all include files with module files by kamata 2021/08/18
      module crand
      implicit none

      integer :: iur0 = 0
      integer :: ia = 32771, ic = 1234567891, im = 2147483647
     >         , iur = 12345
      real(8) :: rm  = 4.65661287307739257812500000000d-10
      real(8) :: ur = 0.0_8
! ia : 2**15+3
! im : 2**31-1
! rm : 1/2**31

      end module crand
