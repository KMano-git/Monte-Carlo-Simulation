!***********************************************************************
      subroutine set_corona(swrd)
!***********************************************************************
      use csize
      use cntcom
      use cplcom
      use cplwrd
      implicit none
!
!::argument
      real(8), dimension(ndmc) :: swrd
!
!::local variables
      integer :: i, ic, iz, ir
      character(80) :: clin
      character(4)  :: cnam
      real(8) :: cimz(10), wrd(ndmc), wrsm, wrgn(ndmc)

!::corona-model for background impurity
!xxx  call plwrdcr(wcr1,wcr2,ndmc)

      call plwrdcrG(swrd)

!::debug write
      do i = 1, wcr_nty
        iz = i
        cimz(1:10) = wcr_cnc(1:10,iz)
        cnam = wcr_typ(iz)
        wrd(1:ndmc) = wcr_wrd(1:ndmc,iz)

        call wrintg(wrd,wrsm,wrgn)
        write(clin,'(2x,"wrdcrG_",a2,2x,1p8e10.2,"  cimz=",1pe10.2)')
     >      cnam(1:2),wrsm,(wrgn(ir),ir=1,7),cimz(1)
        call gdput(trim(clin))
      enddo
!
      return
      end
      
