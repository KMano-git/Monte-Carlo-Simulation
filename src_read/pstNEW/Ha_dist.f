!> \brief calculate Da distribution from densiy of H and electorons
!**********************************************************************
      subroutine Ha_dist(n_pvar, pvar, nv)
!      subroutine pth_load(ctyp,pvar)
!**********************************************************************
!
!       HAtot : Ha_ion +  Ha_rec
!       HAion : Ha_ion
!       HArec : Ha_rec
!
!         svhion :   Si = <sig*v>(Te,ne) [m3/s]
!         emhion :   Eps(Te,ne) = Si/ndot(H_alpha)  [event/photon]
!
!----------------------------------------------------------------------
      use csize
      use cntcom
      use cntpls
      use cntsrc
      use csonic
      use cunit
      implicit none

!         integer, parameter :: ndmc (sonic/inc/csize)
!         /cntpla/ dene, deni, teme  (monte/inc/cntsrc)
!         /cntfig/ tden0             (monte/inc/cntsrc)
!         /cntgrd/ ncmax             (monte/inc/cntcom)
!
!::argument
      !!!!character*(*), intent(in)  :: ctyp
      integer,       intent(in)  :: n_pvar
      real*8,        intent(out) :: pvar(n_pvar)
      integer,       intent(out) :: nv
!
!::local variables
      integer :: ic, ig, ia
      real*8 ::  ane, ani, an0, ate
      !!!!real*8 ::  aha
      !!!!real*8 ::  aha_ion, aha_rec
      real*8 ::  aha_tot
      real*8 :: emhion, emhrcm
!          real*8 functions emhion, enhrcm  \see{dtdedge/atomhy.f}
      real*8 ::  pvmn, pvmx
!
      !!!!write(n6,'(/2x,"*** Da_dist ***  [",a,"]" )') trim(ctyp)
      write(n6,'(/2x,"*** Da_dist ***  [",a,"]" )') 'Ha_tot'
!
      nv=ncmax
      pvmn =  1.0d20
      pvmx = -1.0d20
      do ic = 1, ncmax
          pvar(ic) = 0.0d0
      enddo    
!
!::HA
      !!!!if( ctyp(1:2).eq."HA" ) then
      ig = 1
      ia = 1
      do ic = 1, ncmax
          ane = dene(ic)
          ani = deni(ic,ia) 
          an0 = tden0(ic,ig)
          ate = teme(ic)
          if(ane .eq. 0.0d0 .or. ate .eq. 0.0d0 ) cycle

!           write(*,'(1x, "Ha: ",i, 4e12.5)') ic, ane,ani,an0,ate
!
          aha_tot = ane*(an0*emhion(ate,ane) +ani*emhrcm(ate,ane))
          !!!!aha_ion = ane*an0*emhion(ate,ane)
          !!!!aha_rec = ane*ani*emhrcm(ate,ane)
          !!!!aha_tot = aha_ion + aha_rec
!
          !!!! aha = 0.0d0
          !!!!if( ctyp(1:6).eq."HA_tot" ) aha = aha_tot
          !!!!if( ctyp(1:6).eq."HA_ion" ) aha = aha_ion
          !!!!if( ctyp(1:6).eq."HA_rec" ) aha = aha_rec
!
          pvar(ic) = aha_tot
          pvmn = dmin1( pvmn, aha_tot )
          pvmx = dmax1( pvmx, aha_tot )
          !!!!pvar(ic) = aha
          !!!!pvmn = dmin1( pvmn, pvar(ic) )
          !!!!pvmx = dmax1( pvmx, pvar(ic) )
      enddo
      !!!!endif
!
!::debug write
      write(n6,'(2x,"pvmn =",1pe12.3,"  pvmx =",1pe12.3)') pvmn, pvmx
!
      return
      end
