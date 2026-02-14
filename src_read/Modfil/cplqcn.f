! added replace all include files with module files by kamata 2021/08/18
!----------------------------------------------------------------------
!::conservation
!----------------------------------------------------------------------
      module cplqcn
      use csize, only : ndeq, ndmfl, ndsp

      implicit none

      real(8), allocatable :: qfx_cv(:,:,:), qfx_df(:,:,:)
     >    , qfy_df(:,:,:), qfx_cd(:,:,:), qfy_cd(:,:,:), qvl_pe(:,:,:)
     >    , qvl_cl(:,:,:), qvl_sc(:,:,:), qvl_al(:,:,:), qvl_dt(:,:,:)
     >    , qvl_pi(:,:,:), qfx_vh(:,:), qfy_vh(:,:)

      real(8) :: qsx_cv(ndeq) = 0.0_8, qsx_df(ndeq) = 0.0_8
     >    , qsy_df(ndeq) = 0.0_8, qsx_cd(ndeq) = 0.0_8
     >    , qsy_cd(ndeq) = 0.0_8, qsv_pe(ndeq) = 0.0_8
     >    , qsv_cl(ndeq) = 0.0_8, qsv_sc(ndeq) = 0.0_8
     >    , qsv_dt(ndeq) = 0.0_8, qsv_pi(ndsp) = 0.0_8
     >    , qsx_vh = 0.0_8, qsy_vh = 0.0_8

      real(8) :: pcn_tim = 0.0_8, pcn_dtm = 0.0_8
      real(8) :: pcn_ptl = 0.0_8, pcn_ptb = 0.0_8, pcn_pdt = 0.0_8
      real(8) :: pcn_ptl2 = 0.0_8, pcn_ptb2 = 0.0_8, pcn_pdt2 = 0.0_8
      real(8) :: pcn_psi = 0.0_8, pcn_ptf = 0.0_8, pcn_rhs = 0.0_8
     >    , pcn_pfl(ndmfl) = 0.0_8
      integer :: pcn_itm = 0, pcn_nlp = 0

      integer :: itqcn(10) = 0
      integer, allocatable :: mrgnp(:)

      end module cplqcn
