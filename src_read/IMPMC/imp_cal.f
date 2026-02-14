!***********************************************************************
      subroutine imp_cal
!***********************************************************************
      use cimctl, only : icalz, lstdy
      use cimpuf, only : bk_npt, bk_nty, lbkstw
      use csonic, only : itim
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer :: sv_bk_npt
!
!::cal. bk_mag
      if( bk_nty.gt.0 ) then
        sv_bk_npt = bk_npt
        if( mod(icalZ,8).ne.1 ) bk_npt = 0
!-------
        write(n6,'(2x,"SS-impcal imbkflx ",3i6)') itim, icalZ, bk_npt  
!-------
        call imbkflx
        bk_npt = sv_bk_npt
      endif
!
!::monte carlo calculation
      lbkstw = 0
      if( lstdy.eq.1 ) lbkstw = 1
!-------
        write(n6,'(2x,"SS-impcal immont  ",3i6)') itim, icalZ, lbkstw
!-------
      call immont
!
      return
      end
