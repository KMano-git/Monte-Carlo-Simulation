!***********************************************************************
      subroutine setv2d_mc(cnam,scl,nmc,vmc,var,kcn)
!***********************************************************************
!
!       monte variables (deni,temi,twrd) ==> var
!
!      real*8    dene(0:ndmc),deni(0:ndmc,ndgs)
!     >         ,teme(0:ndmc),temi(0:ndmc),vflw(0:ndmc,ndgs)
!
!     nmc  = ncmax,  ncmax2 (when include Vac-region)
!     scl  = 1.0d0, Ne 1.0d-19, twrd 1.0d-6  vwrd -1.0d-6
!         
!-----------------------------------------------------------------------
      use csize
      use csonic
      use cunit
      implicit none
!
!::argument
      character(*), intent(IN) :: cnam
      real(8),      intent(IN) :: scl
      integer,      intent(IN) :: nmc
      real(8), dimension(ndmc), intent(IN) :: vmc
!
      real(4), dimension(ndmc), intent(OUT) :: var
      integer,      intent(OUT) :: kcn
!
!::local variables
      real(8) :: vmin, vmax, tmp
      integer :: ii, ic

!
      vmax = -1.0d20
      vmin = +1.0d20
!
      var(1:ndmc) = 0.0d0
!
      ii = 0
      do ic = 1, nmc
      tmp = vmc(ic)*scl
      if( tmp.ne.0.0d0 ) then
      vmax = dmax1( vmax, tmp )
      vmin = dmin1( vmin, tmp )
      endif
      !if( tmp.eq.0.0d0 ) tmp = -1.0d30
      var(ic) = tmp
      if( tmp.gt.0.0d0 ) ii = ii + 1
      enddo
!
      kcn = 0
      if( ii.le.0 ) kcn = 1
!
!::min,max
      write(n6,'(/2x,"*** setv2d_mc ***",2x,a4,2x,2i7,2x,1p2e12.4)')
     >  trim(cnam), nmc, ii, vmin, vmax
!
      return
      end
