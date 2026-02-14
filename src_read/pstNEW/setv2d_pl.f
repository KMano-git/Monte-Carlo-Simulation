!***********************************************************************
      subroutine setv2d_pl(cnam,scl,nmc,vpl,var,kcn)
!***********************************************************************
!
!       soldor variables (vni,vte,vwrd) ==> var
!
!     nmc  = ncmax,  ncmax2 (when include Vac-region)
!     scl  = 1.0d0, Ne 1.0d-19, vwrd -1.0d-6
!         
!      real*8    vna(ndx,ndy,ndsp),vne(ndx,ndy),vni(ndx,ndy)
!     >         ,vnezef(ndx,ndy)     ! <== impurity effect 
!     >         ,vzf(ndx,ndy),vva(ndx,ndy,ndsp),vve(ndx,ndy)
!     >         ,vti(ndx,ndy),vte(ndx,ndy),vcs(ndx,ndy)
!
!-----------------------------------------------------------------------
      use csize
      use cntcom
      use csonic
      use cunit
      implicit none
!
!::argument
      character(*), intent(IN) :: cnam
      real(8),      intent(IN) :: scl
      integer,      intent(IN) :: nmc
      real(8), dimension(ndx,ndy), intent(IN) :: vpl
!
      real(4), dimension(ndmc), intent(OUT) :: var
      integer,      intent(OUT) :: kcn
!
!::local variables
      real(8) :: vmin, vmax, tmp
      integer :: ii, ic, j, i
!
      vmax = -1.0d20
      vmin = +1.0d20
!
      var(1:ndmc) = 0.0d0
!
      ii = 0
      do ic = 1, nmc
      j = iplx(ic)
      i = iply(ic)
      if( j.le.0 .or. i.le.0 ) cycle
      tmp = vpl(j,i)*scl
      if( tmp.ne.0.0d0 ) then
      vmax = dmax1( vmax, tmp )
      vmin = dmin1( vmin, tmp )
      endif
      var(ic) = tmp
      if( tmp.gt.0.0d0 ) ii = ii + 1
      enddo
!
      kcn = 0
      if( ii.le.0 ) kcn = 0
!
!::min,max
      write(n6,'(/2x,"*** setv2d_pl ***",2x,a4,2x,2i7,2x,1p2e12.4)')
     >  trim(cnam), nmc, ii, vmin, vmax
!
      return
      end
