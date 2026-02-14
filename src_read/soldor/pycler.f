!**********************************************************************
      subroutine pycler(nj)
!**********************************************************************
      use cplcom, only : cc, dd, ee, ff, gg, nion, ss
      use csize,  only : ndxy
      implicit none

!ik   integer nj
      integer, intent(in) :: nj ! dummy

      integer nequ, k1, k2, k3
!
!::/cslcde/
      nequ = 2*nion + 2
!
      do k3 = 1, ndxy
      do k2 = 1, nequ
      do k1 = 1, nequ
        cc(k1,k2,k3) = 0.0d0
        dd(k1,k2,k3) = 0.0d0
        ee(k1,k2,k3) = 0.0d0
        ss(k1,k2,k3) = 0.0d0
      enddo
      enddo
      enddo
!
      do k3 = 1, ndxy
      do k2 = 1, nequ
        ff(k2,k3) = 0.0d0
        gg(k2,k3) = 0.0d0
      enddo
      enddo
!
      return
      end
