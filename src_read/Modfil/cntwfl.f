! added replace all include files with module files by kamata 2021/08/18
      module cntwfl
      use csize, only : ndmsr
      implicit none

      real(8), allocatable :: xfhti(:,:), xfhta(:,:), xfhtm(:,:)
      real(8) :: xtfi(0:ndmsr) = 0.0_8, xtfa(0:ndmsr) = 0.0_8
     >    , xtfm(0:ndmsr) = 0.0_8
      real(8), allocatable :: xare(:), xdeg(:)
      integer :: xopi(10) = 0, xopa(10) = 0, xopm(10) = 0, xsno = 0
      character :: xcfl(0:ndmsr)*3 = ' '
      real(8), allocatable :: xehta(:,:), xehtm(:,:), xeema(:,:)
     >    , xeemm(:,:), xees(:,:)
      real(8) :: xtea(0:ndmsr,2) = 0.0_8, xtem(0:ndmsr,2) = 0.0_8
     >    , xtees(4) = 0.0_8

!::summation for neutral source
      real(8), allocatable :: dwntl(:,:)
      integer, allocatable :: jxwl(:), iywl(:)

      end module cntwfl
