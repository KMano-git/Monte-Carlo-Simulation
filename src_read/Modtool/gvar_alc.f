! modified dynamic allocation of arrays by kamata 2022/05/29
! added treat 4 or more impurities with IMPMC by kamata 2022/04/21
! allocate global variables
      subroutine gvar_alc( codenm, kk )
      use mod_externalgrid, only:use_exdata
      implicit none
! arguments
      character, intent(in) :: codenm*(*)
      integer,   intent(in) :: kk
! kk     : flag = 1: allocate,    = 2: deallocate
! codenm : code name ( 'master', 'soldor', 'neut2d', 'IMPMC', 'IMPMC2', 'IMPMC3' )

! allocate or deallocate module variables
      if( codenm(1:6) == 'master' ) then
      elseif( codenm(1:6) == 'soldor' ) then
        call cmeffz_alc( kk )
        call cntscg_alc( kk )
        call cntsrc_alc( kk )
        call cntwcn_alc( kk )
        call cplimp_alc( kk )
        call cplpst_alc( kk )
        call cplqcn_alc( kk )
        call cplwrd2_alc( kk )
      elseif( codenm(1:6) == 'neut2d' ) then
        call cntscg_alc( kk )
        call cntwcn_alc( kk )
        call cntxwk_alc( kk )
      elseif( codenm(1:5) == 'IMPMC' ) then
        call cimp_loc_alc( kk )
        call cimpuf_alc( kk )
        call cntsrc_alc( kk )
      elseif( codenm(1:5) == 'plimp' ) then
        call com_pth_line_alc()
        call com_wimp_alc()
        call cmeffz_alc( kk )
        call cplimp_alc( kk )
        call cntwcn_alc( kk )
        call cntxwk_alc( kk )
        call cplqcn_alc( kk )
        call cntsrc_alc( kk )
        call cimp_loc_alc( kk )
        call cimpuf_alc( kk )
        call cplwrd2_alc( kk )
        call cplpst_alc( kk )
        call cplimp_plimp_alc( kk )
      endif

! common setting
      call cplvpn_alc( kk )
      call catcom_alc( kk )
      call cgdcom_alc( kk )
      call cimcom_alc( kk )
      call cimden_alc( kk )
      call clocal_alc( kk )
      call cntcom_alc( kk )
      call cntmnt_alc( kk )
      call cntpls_alc( kk )
      call cntwfl_alc( kk )
      call cplcom_alc( kk )
      call cplmet_alc( kk )
      call cplwrd_alc( kk )
      call cpmpls_alc( kk )
      call czwflx_alc( kk )
      if(use_exdata) call exgrid_alc( kk )

      return
      end
