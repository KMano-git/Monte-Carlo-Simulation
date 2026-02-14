! added replace all include files with module files by kamata 2021/08/18
      module cimcxr
      use csize, only : ndcxrZ => ndmis
      implicit none
  
      integer, parameter :: ndcxrE = 201
      real(8) :: cxr_emin = 0.0_8, cxr_emax = 0.0_8
      real(8) :: cxr_edlt = 0.0_8, cxr_erel(ndcxrE) = 0.0_8
     >    , cxr_vrel(ndcxrE) = 0.0_8
      real(8) :: cxr_sigv(ndcxrE,ndcxrZ) = 0.0_8
      integer :: cxr_nemx = 0, cxr_nzmx = 0
  
      end module cimcxr
  