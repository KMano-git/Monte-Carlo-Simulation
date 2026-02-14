!**********************************************************************
      subroutine ntdens(tflux)
!**********************************************************************
!
!       source terms of plasma
!
!----------------------------------------------------------------------
      use cntcom, only : den0, deng, eng0, engg, ncmax2
     >    , ngas, rmas, vlp0, volm, nfl0x, nfl0y, nfl0z
      use cntmnt, only : dotn
      use cntwcn, only : swnrm, wden, weng, wgden, wgeng, wnfl0x,
     >   wnfl0y, wnfl0z, wnflgx, wnflgy, wnflgz, wvlp
      use cphcns, only : cev
      use csize,  only : ndgs, ndmc
      implicit none
!
!::argument
      real(8), intent(in) :: tflux
!
!::local variables
      integer nsiz, ig, ic
      real*8  zfce, zfcn, zmas
!
!-----------------------------------------------------------------------
!::flux
!-----------------------------------------------------------------------
      dotn = tflux
!
!-----------------------------------------------------------------------
!::den0,eng0 for each source
!-----------------------------------------------------------------------
      nsiz = (ndmc+1)*ndgs
      call setd( den0, nsiz, 0.0d0 )
      call setd( eng0, nsiz, 0.0d0 )
      call setd( vlp0, nsiz, 0.0d0 )
      nsiz = (ndmc+1)*2
      call setd( deng, nsiz, 0.0d0 )
      call setd( engg, nsiz, 0.0d0 )
!
      do ig = 1, ngas
        zfce = 2.0d0/3.0d0*0.5d0*rmas(ig)/cev
        do ic = 1, ncmax2
          zfcn = dotn/swnrm/volm(ic)
          den0(ic,ig) = wden(ic,ig)*zfcn
          if( wden(ic,ig).gt.0.0d0 ) then
            eng0(ic,ig) = weng(ic,ig)/wden(ic,ig)*zfce
            vlp0(ic,ig) = wvlp(ic,ig)/wden(ic,ig)
! neutral flow 160623 toku
            nfl0x(ic,ig)=wnfl0x(ic,ig)/wden(ic,ig)
            nfl0y(ic,ig)=wnfl0y(ic,ig)/wden(ic,ig)
            nfl0z(ic,ig)=wnfl0z(ic,ig)/wden(ic,ig)
          endif
        enddo
      enddo
!
      do ig = 1, 2
        zmas = 2.0d0*rmas(1)
        zfce = 2.0d0/3.0d0*0.5d0*zmas/cev
        do ic = 1, ncmax2
          zfcn = dotn/swnrm/volm(ic)
          deng(ic,ig) = wgden(ic,ig)*zfcn
          if( wgden(ic,ig).gt.0.0d0 ) then
            engg(ic,ig) = wgeng(ic,ig)/wgden(ic,ig)*zfce
          endif
        enddo
      enddo
!
      return
      end
