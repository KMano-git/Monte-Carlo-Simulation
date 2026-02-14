!**********************************************************************
      subroutine pxcler(nt)
!**********************************************************************
      use cplcom, only : cc, dd, ee, ff, gg, fl1a, fl2a, fl3, fl4, nion
     >    , ss
      use csize, only : ndx, ndxy
      implicit none
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nt
      integer, intent(in) :: nt ! dummy

      integer  k1, k2, k3, nequ
!
!::/cslflx/
      do k2=1,nion
      do k1=1,ndx
        fl1a(k1,k2) = 0.0d0
        fl2a(k1,k2) = 0.0d0
      enddo
      enddo
      do k1=1,ndx
        fl3(k1) = 0.0d0
        fl4(k1) = 0.0d0
      enddo
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
      do k3 = 1, ndxy
      do k2 = 1, nequ
        ff(k2,k3) = 0.0d0
        gg(k2,k3) = 0.0d0
      enddo
      enddo
!
      return
      end
