!***********************************************************************
      subroutine imsumry(kk)
!***********************************************************************
      use cimcom, only : hcnpt, idcal, idemt, iderr, idetc, idstk
      use cunit,  only : n6
      implicit none
!
!:argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  kk
      integer, intent(in) :: kk ! dummy
!
!::local variables
      real(8) :: zero = 1.0d-8
!
      write(n6,'(2x,"*** imsumry ****  emt =",i6,"  cal =",i6,
     >   "  stk =",i6,"  err =",i6,"  etc =",i6)')
     >  int(hcnpt(idemt)+zero), int(hcnpt(idcal)+zero),
     >  int(hcnpt(idstk)+zero), int(hcnpt(iderr)+zero),
     >  int(hcnpt(idetc)+zero)
!
      return
      end
