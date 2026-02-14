!*********************************************************************
      subroutine plconsvq
!**********************************************************************
      use cphcns, only : cev
      use cplcom, only : ama, aza, nion, q1a, q2a, q3, q4, qw1a, qw2a
     >    , qw3, qw4, vna, vne, vte, vti, vva
      use cplmet, only : icel, itmax, itpve, itsls, jcel, jtmax, jtmin
      implicit none
!
!::local vairables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: it, jt, jtst, jten, jc, ic, ia, kbc
      integer :: it, jt, jc, ic, ia, kbc
      real(8) :: zne, zq3
!
!-----------------------------------------------------------------------
!::conservative variables & vne
!-----------------------------------------------------------------------
      do 510 it = 1, itmax
      do 520 jt = jtmin(it), jtmax(it)
      jc  = jcel(jt,it)
      ic  = icel(jt,it)
      zne = 0.0d0
      zq3 = 0.0
      do 530 ia = 1,nion
      q1a(jc,ic,ia) = ama(ia)*vna(jc,ic,ia)
      q2a(jc,ic,ia) = ama(ia)*vna(jc,ic,ia)*vva(jc,ic,ia)
      zne = zne + aza(ia)*vna(jc,ic,ia)
      zq3 = zq3 + 1.5d0*vna(jc,ic,ia)*vti(jc,ic)*cev
     >          + 0.5d0*ama(ia)*vna(jc,ic,ia)*vva(jc,ic,ia)**2
 530  continue
      vne(jc,ic) = zne
      q3(jc,ic)  = zq3
      q4(jc,ic)  = 1.5d0*vne(jc,ic)*vte(jc,ic)*cev
 520  continue
 510  continue
!
!-----------------------------------------------------------------------
!::boundary condition at sol-wall and prv-wall
!-----------------------------------------------------------------------
      kbc = 1
      it  = itsls
      jt  = (jtmin(it)+jtmax(it))/2
      jc = jcel(jt,it)
      ic = icel(jt,it)
      do ia = 1, nion
      qw1a(ia,kbc) = q1a(jc,ic,ia)
      qw2a(ia,kbc) = q2a(jc,ic,ia)
      enddo
      qw3(kbc) = q3(jc,ic)
      qw4(kbc) = q4(jc,ic)
!
      kbc = 2
      it = itpve
      jt = (jtmin(it)+jtmax(it))/2
      jc = jcel(jt,it)
      ic = icel(jt,it)
      do ia = 1, nion
      qw1a(ia,kbc) = q1a(jc,ic,ia)
      qw2a(ia,kbc) = q2a(jc,ic,ia)
      enddo
      qw3(kbc) = q3(jc,ic)
      qw4(kbc) = q4(jc,ic)
!
      return
      end
