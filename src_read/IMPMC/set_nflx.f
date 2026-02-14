!***********************************************************************
      subroutine set_nflx(kk)
!***********************************************************************
!
!    soldor        neut2d      impmc
!    pfl_ion  ==>  vflux  ==>  pfl_ion
!
!-----------------------------------------------------------------------
      use cntmnt, only : mfmax, pfl_ion, pfl_ntl, vcsrc, vflux, visrc
     >    , vkflx, vnsrc
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  kk   ! no use
      integer, intent(in) :: kk   ! dummy
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  ia, i, isrc, iflx
      integer  i, isrc, iflx
!
      write(n6,'(/2x,"*** set_nflx ***  vflux => pfl_ion")')
!
!::soldor  flux from core edge  pfl_ion
!xx      call plwtfy
!xx      call plpflx
!xx      call plmflx
!
!::monte   vflux => pfl_ion
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   ia = 1
      do i = 1, mfmax
      if( i.eq.5 ) cycle
      pfl_ion(i) = 0.0d0
      pfl_ntl(i) = 0.0d0
      enddo
!
! modified 2/1 lines organize local variables and include files by kamata 2021/06/28
!ik   do i = 1, vnsrc
!ik   isrc = i
      do isrc = 1, vnsrc
      iflx = vkflx(isrc)
      visrc(iflx) = isrc
      if( vcsrc(isrc)(1:3).ne."vol" ) then
      pfl_ion(iflx) = vflux(isrc)
      pfl_ntl(iflx) = vflux(isrc)
      else
      pfl_ntl(7) = vflux(isrc)
      pfl_ion(8) = vflux(isrc)
      endif
      enddo
!
!::debug write
 100  continue
!xx      do i = 1, mfmax
!xx      iflx = i
!xx      isrc = visrc(iflx)
!xx      write(n6,'(i3,2x,i3,2x,a,2x,1p2e14.6)')
!xx     >  iflx, isrc, mcflx(iflx), pfl_ion(iflx), pfl_ntl(iflx)
!xx      enddo
!
      return
      end
