!**********************************************************************
      subroutine gdpcnv(kk,gxp,gyp,cxp,cyp)
!**********************************************************************
!
!    kk = 1  axis (gxp,gyp) ==> paper (cxp,cyp)
!    kk = 2  axis (gxp,gyp) <== paper (cxp,cyp)
!
!---------------------------------------------------------------------
!xx   implicit none
      include '#glinc'
!
!::local common
!xx   real(4) :: fsmi(0:2), fsma(0:2)
!xx   real(4) :: cpx1, cpx2, cpy1, cpy2
!
!::argument
      integer :: kk
      real(4) :: gxp, gyp
      real(4) :: cxp, cyp
!
!::local varibales
!
!::kk=1  (gxp,gyp) ==> (cxp,cyp)
      if( kk.eq.1 ) then
      cxp = (cpx2-cpx1)/(fsma(0)-fsmi(0))*(gxp-fsmi(0)) + cpx1
      cyp = (cpy2-cpy1)/(fsma(1)-fsmi(1))*(gyp-fsmi(1)) + cpy1 
!
!::kk = 2 (gxp,gyp) <== (cxp,cyp)
      elseif( kk.eq.2 ) then
      gxp = (fsma(0)-fsmi(0))/(cpx2-cpx1)*(cxp-cpx1) + fsmi(0)
      gyp = (fsma(1)-fsmi(1))/(cpy2-cpy1)*(cyp-cpy1) + fsmi(1)
!
!::error
      else
      write(6,'(/2x,"stop at sub. gdpcnv   kk =",i2)') kk
      call gdmod(0)
      stop
      endif
!
      write(6,'(2x,"gdpcnv  ",i3,2f12.4,2x,2f12.4)')
     >  kk,gxp,gyp,cxp,cyp
!
      end
