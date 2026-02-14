!**********************************************************************
      subroutine dtypset(dtyp,mtyp)
!**********************************************************************
      use catcom,       only : catmz
      use cimcom,       only : ismax, lwrad, weflx, wpflx
      use cimden,       only : nsput
      use cplcom,       only : mdl_wrd
      use cplwrd,       only : wfac
      use cunit,        only : n6
      use czwflx,       only : zw_catmz, zw_eflx, zw_ismax, zw_pflx
      use mod_sizedef,  only : lnnam
      implicit none

      character(*), intent(in) :: dtyp
      integer, intent(in) :: mtyp ! 1:pre-Bcst,2:post-bcast,3:Send,4:Recv
      character(lnnam) :: ctyp 

! added 6 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
! local variables
      integer         mj, nc
      character(2)    np
! mj : number of characters
! nc : code group number
! np : processing number
      integer i

      ctyp = trim(dtyp)

!:: PLS_1 Recv
      if(ctyp == "PLS_1" .and. mtyp == 4) then
!::     lwrad <= mdl_wrd
        lwrad = mdl_wrd
      endif

! added 22 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
! IMP(nc)_2B send
      if( ctyp(1:3) == 'IMP' .and. mtyp == 3 ) then
! get processing number and code group number
        call getncnp( ctyp, nc, np )

! execution
        if( nc > 0 .and. np == '2B' ) then
          zw_catmz(nc) = catmz(1)
          zw_ismax(nc) = ismax
          zw_pflx(:,:,nc) = 0.0d0
          zw_eflx(:,:,nc) = 0.0d0
          do i = 1, nsput
            zw_pflx(1:ismax,:,nc) = zw_pflx(1:ismax,:,nc)
     >                            + wpflx(1:ismax,:,i)
            zw_eflx(1:ismax,:,nc) = zw_eflx(1:ismax,:,nc)
     >                            + weflx(1:ismax,:,i)
          enddo
          if(nc == 1)then !wfac is only applied to nc=1
          write(n6,*)"  dtypset: zw_flx is multiplied by wfac = ", wfac
          zw_pflx(:,:,nc) = zw_pflx(:,:,nc) * wfac
          zw_eflx(:,:,nc) = zw_eflx(:,:,nc) * wfac
          endif
        endif
      endif
      return
      end
