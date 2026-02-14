!***********************************************************************
      subroutine set_flx(kk)
!***********************************************************************
!
!        kk = 1     calculation flux in soldor  
!        kk = 2     read disk data of ntdisk
!
!-----------------------------------------------------------------------
      use csize
      use cntcom
      use cntmnt
      use csonic
      use cunit
      implicit none
!
!::argument
      integer  kk
!
!::local variables
      integer  ia, i, isrc, iflx
!
      write(n6,'(/2x,"*** set_flx ***  kk =",i2)') kk
!
!::flux from core edge
      call plwtfy
      call plpflx
      call plmflx
!
!::soldor
      if( kk.eq.1 ) then
        write(n6,'(2x,"   flux calculated in soldor-code")')
        goto 100
      endif
!
!::ntdisk
      write(n6,'(2x,"   flux from  ntdisk")')
      ia = 1
      do i = 1, mfmax
      if( i.eq.5 ) cycle
      pfl_ion(i) = 0.0d0
      pfl_ntl(i) = 0.0d0
      enddo
!
      do i = 1, vnsrc
      isrc = i
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
      do i = 1, mfmax
      iflx = i
      isrc = visrc(iflx)
      write(n6,'(i3,2x,i3,2x,a,2x,1p2e14.6)') 
     >  iflx, isrc, mcflx(iflx), pfl_ion(iflx), pfl_ntl(iflx)
      enddo
!
      return
      end
