!*******************************************************************************
      subroutine set_plimp_1( kis )
!***********************************************************************
!     set cplimp_cimcom_1 in soldor from cimden_1 in IMPMC (IMP_2)
!         aimas, azmas, ami, amz, ismax
!     set cplimp_cimcom_8 in soldor from cimden_8 in IMPMC (IMP_2)
!         cdif, cdifrg
!
!
      use cimcom, only : aimas, ami, amz, azmas, cdif, ismax
      use cplimp, only : aimasl, amil, amzl, azmasl, cdifl, ismaxl
     >    , wmc_nty
      use csize,  only : nzsmx
      implicit none

      integer, intent(in) :: kis

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
      subroutine set_plimp_2( kis )
!***********************************************************************
!     set cplimp_imden in soldor from cimden in IMPMC (IMP_2)
!
      use cimden, only : eipot, nzmx, tdnz, twci, tengz, tprz
      use cplimp, only : eipotl, nzmxl, tdnzl, twcil, tengzL, tprzL
      implicit none

      integer, intent(in) :: kis


!:: IMPMC/inc/cimden  --> soldor/inc/cplimp
       nzmxL(kis)      = nzmx
       eipotL(:, kis ) = eipot
       tdnzL(:,:, kis) = tdnz
       twciL(:, kis)   = twci
      tengzL(:,:,kis) = tengz
      tprzL(:,:,kis)  = tprz
!
      return
      end
