!**********************************************************************
      subroutine read_upldif(nft,ierr)
!**********************************************************************
      use cplcom, only : vlda, vlet, vlxe, vlxi
     >    , lerp_points_vldar, lerp_points_vletr
     >    , lerp_points_vlxir, lerp_points_vlxer
     >    , lerp_r_vldar, lerp_diff_vldar 
     >    , lerp_r_vletr, lerp_diff_vletr
     >    , lerp_r_vlxir, lerp_diff_vlxir
     >    , lerp_r_vlxer, lerp_diff_vlxer
      use cimcom, only : cdifrg_imp, cdif_imp, match_solimpdiff
     >    , lerp_points_cdif,lerp_r_cdif,lerp_diff_cdif
      implicit none
!::argument
! nft:: number for read namelist file
      integer, intent(in) :: nft
      integer, intent(out) :: ierr
!
!::local variable
!
! diffusion coefficient of soldor and IMPMC
      namelist /upldif/
! following is diffusion coefficient for soldor. IMPMC may use this.
     >  lerp_points_vldar,lerp_r_vldar,lerp_diff_vldar
     > ,lerp_points_vletr,lerp_r_vletr,lerp_diff_vletr
     > ,lerp_points_vlxir,lerp_r_vlxir,lerp_diff_vlxir
     > ,lerp_points_vlxer,lerp_r_vlxer,lerp_diff_vlxer
     > ,vlda,vlet,vlxi,vlxe
! following is diffusion coefficient for IMPMC, not use in soldor.
     > ,cdif_imp,cdifrg_imp
     > ,lerp_points_cdif,lerp_r_cdif,lerp_diff_cdif
     > ,match_solimpdiff

      rewind(nft)
      read(nft,upldif,iostat=ierr)
      end subroutine