!***********************************************************************
      subroutine calflw(ctyp,nv,frz,thz,vlvlz,ldenZ,lionZ,lrecZ,dnz,
     > radi,radli,radr)
!***********************************************************************
!
!         wrd(ic)    : radiation loss
!         wci(ic)    : cold ion due to scattering process
!         dnz(iz,ic) : impurity dentity with charge state of iz
!
!-----------------------------------------------------------------------
      use cimcom, only : cdi, cdlz_i, cdlz_r, denz, friz, ionz, ismax
     >    , ndis, recz, sflux, swtot, thfz, vzpz
      use cimden, only : eipot
      use cntcom, only : ncmax2, volm
      use cntpls, only : dene
      use cphcns, only : cev
      use csize,  only : ndmc
      use cunit,  only : n6
      implicit none

!::argument
      character, intent(in)  :: ctyp*(*)
      real(8),   intent(in)  :: dnz(0:ndis,ndmc)
      integer,   intent(out) :: nv
      real(8),   intent(out) :: frz(0:ndis,ndmc), thz(0:ndis,ndmc)
      real(8),   intent(out) :: vlvlz(0:ndis,ndmc)
      real(8),   intent(out) :: ldenZ(0:ndis,ndmc)
      real(8),   intent(out) :: lionZ(0:ndis,ndmc), lrecZ(0:ndis,ndmc)
      real(8),   intent(out) :: radi(0:ndis,ndmc), radli(0:ndis,ndmc)
      real(8),   intent(out) :: radr(0:ndis,ndmc)
!
!::local variables
      integer ic, iz
      real(8)  flnrm, zne, znz
!
!::clear
      nv = ncmax2
      frz = 0.0d0
      thz = 0.0d0
      vlvlz = 0.0d0
      ldenZ = 0.0d0
      lionZ = 0.0d0
      lrecZ = 0.0d0
      radi = 0.0d0
      radli = 0.0d0
      radr = 0.0d0

!::debug write
      write(n6,'(2x,"*** calflw ***",2x,a)') ctyp
      if(swtot.le.0.0d0)then
        write(n6,'(4x,"skip calwrd due to swtot <= 0.0d0")')
        return
      endif
!
      flnrm = sflux/swtot  ! see sub. imemit & imatom
!
!
!::absolue density
      do ic = 1, nv
      do iz = 0, ismax
      frz(iz,ic) = friZ(iz,ic)
      thz(iz,ic) = thfZ(iz,ic)
      vlvlz(iz,ic) = vzpZ(iz,ic)
      ldenZ(iz,ic) = denZ(iz,ic)
      lionZ(iz,ic) = flnrm*ionZ(iz,ic)/volm(ic)
      lrecZ(iz,ic) = flnrm*recZ(iz,ic)/volm(ic)
      enddo
      enddo

      if( index(ctyp,"ion").gt.0 ) then
      do ic = 1, nv
      zne  = dene(ic)
      do iz = 0, ismax
      znz  = dnz(iz,ic)
      radi(iz,ic)=zne*znz*cdLz_i(iz,ic)
      enddo
      enddo
!
!::Note    cdi = Ne*<sigv>ion
      do ic = 1, nv
      zne  = dene(ic)
      do iz = 0, ismax
      znz  = dnz(iz,ic)
      radli(iz,ic) = znz*cdi(iz,ic)*eipot(iz)*cev
      enddo
      enddo
      endif
!
!::recombination
      if( index(ctyp,"rec").gt.0 ) then
      do ic = 1, nv
      zne  = dene(ic)
      do iz = 0, ismax
      znz  = dnz(iz,ic)
      radr(iz,ic) = zne*znz*cdLz_r(iz,ic)
      enddo
      enddo
      endif
      end
