!***********************************************************************
      subroutine calwrd(ctyp,nv,wrd,wci,dnz,engz)
!***********************************************************************
!
!         wrd(ic)     : radiation loss
!         wci(ic)     : cold ion due to scattering process
!         dnz(iz,ic)  : impurity dentity with charge state of iz
!         engz(iz,ic) : impurity energy with charge state of iz
!
!-----------------------------------------------------------------------
      use cimcom, only : amz, cdi, cdlz_cx, cdlz_i, cdlz_r, denz, ismax
     >    , ndis, sflux, sptyc, swtot, wsct, temZ
      use cimden, only : eipot
      use cntcom, only : ncmax2, volm
      use cntpls, only : dene
      use cphcns, only : cev
      use cplcom, only : tn0
      use csize,  only : ndmc, ndy
      use cunit,  only : n6
      implicit none

!::argument
      character, intent(in)  :: ctyp*(*)
      integer,   intent(out) :: nv
      real(8),   intent(out) :: wrd(ndmc), wci(ndmc), dnz(0:ndis,ndmc)
     > , engz(0:ndis,ndmc)
!
!::local variables
      integer ic, iz
      real*8  flnrm, twr_i, twr_l, twr_r, twr_cx
      real*8  twr, twr_rg(10), twr_it(ndy)
      real*8  twc, twc_rg(10), twc_it(ndy)
      real*8  zne, zsum, znz, zn0, zwc, ftz

!::clear
      nv = ncmax2
      wrd = 0.0d0
      wci = 0.0d0
      dnz = 0.0d0
      engz = 0.0d0
      twr_i  = 0.0d0
      twr_r  = 0.0d0
      twr_cx = 0.0d0
      twr_l  = 0.0d0
      twr    = 0.0d0

!::debug write
      write(n6,'(2x,"*** calwrd ***",2x,a)') ctyp
      if(swtot.le.0.0d0)then
        write(n6,'(4x,"skip calwrd due to swtot <= 0.0d0")')
        return
      endif
!
      flnrm = sflux/swtot  ! see sub. imemit & imatom
!
!::absolue density
      ftz = 0.5d0*amz/cev*2.0d0/3.0d0
      do ic = 1, nv
        do iz = 0, ismax
          dnz(iz,ic) = flnrm*denZ(iz,ic)/volm(ic)
          if(denZ(iz,ic) .ne. 0) then
            engz(iz,ic) = ftz*temZ(iz,ic)/denZ(iz,ic)
          else
            engz(iz,ic) = 0.0d0
          endif
        enddo
      enddo
!
!::ionization
      if( index(ctyp,"ion").gt.0 ) then
      do ic = 1, nv
      zne  = dene(ic)
      zsum = 0.0d0
      do iz = 0, ismax
      znz  = flnrm*denZ(iz,ic)/volm(ic)
      zsum = zsum + zne*znz*cdLz_i(iz,ic)
      enddo
      wrd(ic) = wrd(ic) + zsum
      twr_i   = twr_i + zsum*volm(ic)
      enddo
!
!::Note    cdi = Ne*<sigv>ion
      do ic = 1, nv
      zne  = dene(ic)
      zsum = 0.0d0
      do iz = 0, ismax
      znz  = flnrm*denZ(iz,ic)/volm(ic)
      zsum = zsum + znz*cdi(iz,ic)*eipot(iz)*cev
      enddo
      wrd(ic) = wrd(ic) + zsum
      twr_l   = twr_l + zsum*volm(ic)
      enddo
      endif
!
!::recombination
      if( index(ctyp,"rec").gt.0 ) then
      do ic = 1, nv
      zne  = dene(ic)
      zsum = 0.0d0
      do iz = 0, ismax
      znz  = flnrm*denZ(iz,ic)/volm(ic)
      zsum = zsum + zne*znz*cdLz_r(iz,ic)
      enddo
      wrd(ic) = wrd(ic) + zsum
      twr_r   = twr_r + zsum*volm(ic)
      enddo
      endif
!
!::cxr  energy loss for H0 (Not H+)
      if( index(ctyp,"cxr").gt.0 ) then
      call wexit("calwrd","under development")
      do ic = 1, nv
      zn0  = tN0(ic,1)
      zsum = 0.0d0
      do iz = 0, ismax
      znz  = flnrm*denZ(iz,ic)/volm(ic)
      zsum = zsum + zn0*znz*cdLz_cx(iz,ic)
      enddo
      wrd(ic) = wrd(ic) + zsum
      twr_cx  = twr_cx + zsum*volm(ic)
      enddo
      endif
!
!::energy exchange with plasma ions
      do ic = 1, nv
      zwc = flnrm*0.5d0*amz*wsct(ic)/volm(ic)
      wci(ic) = zwc
      enddo
!
!::integration
! twr_i, twr_l, and twr_r include core radiation SYamoto
      call totrgn(wrd, twr, twr_rg, twr_it)
      call totrgn(wci, twc, twc_rg, twc_it)
!
      write(n6,'(2x,a,2x,"sflux =",1pe12.4,"  swtot =",1pe12.4,
     >   "  flnrm =",1pe12.3)') sptyc, sflux, swtot, flnrm
      write(n6,'(2x,"totwr =",1pe12.3,"  wr_i =",1pe12.3,"  wr_l =",
     >  1pe12.3,"  wr_r =", 1pe12.3,"  wr_cx =",1pe12.3,"  wci =",
     >  1pe12.3)')  twr, twr_i, twr_l, twr_r, twr_cx, twc
      write(n6,'(2x,"wr_rg =",1p10e12.3)') twr_rg
      return
      end
