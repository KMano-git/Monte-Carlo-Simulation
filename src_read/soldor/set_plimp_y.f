!*******************************************************************************
      subroutine set_plimp_1y( kis )
!***********************************************************************
!     set cplimp_cimcom_1 in soldor from cimden_1 in IMPMC (IMP_2)
!         aimas, azmas, ami, amz, ismax  
!     set cplimp_cimcom_8 in soldor from cimden_8 in IMPMC (IMP_2)
!         cdif, cdifrg
!
!
      use csize
      use cimcom
      use cimden
      use cplimp
      implicit none

      integer  :: kis

      if(kis.eq.1) wmc_nty=0
      wmc_nty = wmc_nty + 1

      if(wmc_nty .gt. nzsmx)then
        call wexit("set_plimp","wmc_nty.gt.nzsmx")
      endif

      aimasL(kis) = aimas
      azmasL(kis) = azmas
      amiL(kis)   = ami
      amzL(kis)   = amz
      ismaxL(kis) = ismax

      cdifL(kis)  = cdif

      return
      end


!*******************************************************************************
      subroutine set_plimp_2y( kis )
!***********************************************************************
!     set cplimp_imden in soldor from cimden in IMPMC (IMP_2)
!
      use csize
      use cimcom
      use cimden
      use cplimp
      use cplimp_plimp
      implicit none

      integer  :: kis
      integer, save :: icnt = 0
     

!:: IMPMC/inc/cimden  --> soldor/inc/cplimp 
       nzmxL(kis)      = nzmx 
       tfrzL(:,:, kis) = tfrz
       tthzL(:,:, kis) = tthz 
       tvlzL(:,:, kis) = tvlz 
       tionZL(:,:,kis) = tionZ
       trecZL(:,:,kis) = trecZ
!
      return
      end
