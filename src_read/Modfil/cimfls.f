!::cimfls   impurity in/out-flux into the core   IMPMC/IMPACT
!
!    lfls       : option of calculating impurity flux
!    fls_roh(5) : roh values (0<roh<1.0)
!   
!    Note  ndis : iz charge state  ndmp : ip sample particle
!
!    fls_recI(ndmp,5)   : record when sample particles pass the surface
!    fls_recO(ndmp,5)   : record when sample particles pass the surface
!
!    fls_flxI(0:ndis,5) : inward  flux across surface with fls_roh(i)
!    fls_flxO(0:ndis,5) : outward flux across surface with fls_roh(i)
!
!    fls_dnz(0:ndis,5)  : averaged impurity density at roh = fls_roh(1:5)
!
!--------------------------------------------------------------------------
! added replace all include files with module files by kamata 2021/08/18
      module cimfls
      use csize, only : ndmis
      implicit none

      real(8) :: fls_roh(5) = 0.0_8
      integer :: lfls = 0

      real(8) :: fls_flxI(0:ndmis,5) = 0.0_8
      real(8) :: fls_flxO(0:ndmis,5) = 0.0_8
      real(8) :: fls_dnz(0:ndmis,5) = 0.0_8

      end module cimfls
