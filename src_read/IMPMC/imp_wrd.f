!***********************************************************************
      subroutine imp_wrd
!***********************************************************************
!)
!)    imp_end.f  ==> imp_wrd.f  time averaged radiation loss
!)
!)    sub calwrd("ionrec", nv,wrd,wci,dnz)  wrd for a sputtering
!)    do k = 1, nspt; twrd = twrd + wrd
!)    twrd = (1-ftav0)*twrd + ftab0*wrd0
!)
!)----------------------------------------------------------------------
      use cimcom, only : ismax, ndis, phdnz, phfrz, phionz, phrecz
     >    , phthz, phvlz, phwci, phwrd, sflux, sptyc, swtot
      use cimden, only : bkspt, csput, fsput, nsput, tdnz, tfrz, tionz
     >    , trecz, tthz, tvlz, twci, twrd, wtsput
      use cimpuf, only : bk_mag
      use cntcom, only : ncmax2, nogt
      use csize,  only : ndmc
      use csonic, only : limp
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer   :: nv, k
      real(8)   :: zsum
      zsum = 0.0d0
      nv = ncmax2
      if( limp.ne.3 ) return

!::previous results
!  Note  do not calculate twrd in other routine
!::clear twrd, twci and tdnz in sub. imp_ini
      twrd(1:ndmc) = 0.0d0
      twci(1:ndmc) = 0.0d0
      tdnz(0:ndis,1:ndmc) = 0.0d0
      tfrz(0:ndis,1:ndmc) = 0.0d0
      tthz(0:ndis,1:ndmc) = 0.0d0
      tvlz(0:ndis,1:ndmc) = 0.0d0

!::sum up for generations
      do k = 1, nsput
      sptyc = csput(k)
      sflux = fsput(k)
      swtot = wtsput(k)
!
      write(n6,'(/2x,a,2x,i3,2x,a,2x,"sflux =",1pe12.4,"  swtot =",
     >   1pe12.4)') "imp_wrd", k, sptyc, sflux, swtot

      call imsave(sptyc,"r",k)
      call imsave(sptyc,"z",k)

      twrd(1:nv) = twrd(1:nv) + phwrd(1:nv)
      twci(1:nv) = twci(1:nv) + phwci(1:nv)
      tdnz(0:ismax,1:nv) = tdnz(0:ismax,1:nv) + phdnz(0:ismax,1:nv)
      tfrz(0:ismax,1:nv) = tfrz(0:ismax,1:nv) + phfrz(0:ismax,1:nv)
      tthz(0:ismax,1:nv) = tthz(0:ismax,1:nv) + phthz(0:ismax,1:nv)
      tionz(0:ismax,1:nv) = tionz(0:ismax,1:nv) + phionz(0:ismax,1:nv)
      trecz(0:ismax,1:nv) = trecz(0:ismax,1:nv) + phrecz(0:ismax,1:nv)
      tvlz(0:ismax,1:nv) = tvlz(0:ismax,1:nv) + phvlz(0:ismax,1:nv)
      zsum = sum(tvlz(:,:))
      enddo   !  loop(sputtering)
      write(n6,*) 'DBG: imp_wrd tvlz_total =', zsum
!
!::new variables to soldor
      bkspt(1:nogt) = bk_mag(1:nogt)
!
      return
      end
