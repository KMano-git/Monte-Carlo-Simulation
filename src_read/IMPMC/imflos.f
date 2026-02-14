!***********************************************************************
      subroutine imflos(ic, iz, weight, roh, rohb, kcnd)
!***********************************************************************
!
!       kcnd = ien(ip) : = kind see sub. imtrci.f
!                        = 6 (stick) = 7 (absorption) = 9 (error)
!
!       mrgn(ic)    : = 6 (core edge) = 7 (core)
!       iz = is(ip) :  charge state
!
!     inward flux
!       r    ri          rb
!     *<===|===========*   rec = 0  r-ri < 0  rb-ri > 0
!
!     outward flux
!       rb   ri          r
!       *====|==>=======>*   rec = 0  r-ri > 0  rb-ri < 0
!    0.0                  1.0
!
!-----------------------------------------------------------------------
      use cimfls, only : fls_flxi, fls_flxo, fls_roh
      use cntcom, only : mrgn
      implicit none
!
!::argument
      integer, intent(in) :: ic, iz, kcnd
      real(8), intent(in) :: weight, roh, rohb

!::local variables
      integer :: n
      real(8) :: zroi


!::not subject to flux loss at the core edge
      if( kcnd >= 6 ) return
      if( mrgn(ic) /= 6 .and. mrgn(ic) /= 7 ) return

!::roh
      if( roh > 1.000d0 ) return

!::new type
      do n = 1, 5
        zroi = fls_roh(n)
!---------
        if( roh - zroi < 0.0d0 .and. rohb - zroi > 0.0d0 ) then
          fls_flxI(iz,n) = fls_flxI(iz,n) + weight
        endif
        if(  roh - zroi > 0.0d0 .and. rohb - zroi < 0.0d0 ) then
          fls_flxO(iz,n) = fls_flxO(iz,n) + weight
        endif
      enddo

      return
      end

!***********************************************************************
      subroutine lstflos(ctyp)
!***********************************************************************
!
!        1PE  :  particle number in 1PE
!        TPE  :  MPI_reduce
!        flux :  flux at core edge  (5 point)
!
!-----------------------------------------------------------------------
      use cimcom, only : ismax, npmax, sflux, sptyc, swtot
      use cimctl, only : icalz, nprg
      use cimfls, only : fls_flxi, fls_flxo, fls_roh
      use cplwrd, only : wfac
      use csonic, only : itim
      use cunit,  only : lmype, ufl6
      use mod_shexe, only : impmc_model
      use csize,  only : ndmis
      implicit none

      character(*), intent(in) :: ctyp

      integer :: iz
      real(8) :: fac
      real(8) :: ednz(0:ndmis,5)
      real(8) :: edne(5), edni(5), etne(5), etni(5)
      character(20) :: dsn
      character(20)    cbf
      integer :: i6 = 0
      save :: i6

      if( trim(ufl6) == "/dev/null" ) return
      if( i6 == 0 ) then
        i6 = 250000
        write(dsn,'(i5)') nprg
        write(cbf,'(i10)') lmype
        dsn = adjustl( dsn )
        cbf = adjustl( cbf )
        dsn = 'flos_' // trim( dsn ) // '_' // trim( cbf )
        open(unit=i6,file=trim(dsn))
      endif

      call caldnz(fls_roh,ednz)
      if( impmc_model == 0 ) then
        call caldne(fls_roh,edni,etni,1) !ion
        call caldne(fls_roh,edne,etne,2) !electron
      endif

      if( trim(ctyp) == "1PE" .or. trim(ctyp) == "TPE" ) then
        write(i6,'(/2x,"*** lstflos ***",2x,a)') trim(ctyp)
        write(i6,'(2x,"lmype =",i4,"  itim =",i6,"  icalZ =",i3,2x,a,
     > "  npmax =",i5,2x,a)') lmype,itim,icalZ,sptyc,npmax

        write(i6,'(2x," iz",1p5e10.3,2x,1p5e10.3)')
     >    fls_roh(1:5),fls_roh(1:5)
        do iz = 0, ismax
          if( mod(iz,3) /= 1 ) cycle
          write(i6,'(2x,i3,1p5e10.3,2x,1p5e10.3,2x,1p5e10.3)') iz,
     >     fls_flxI(iz,1:5), fls_flxO(iz,1:5)
        enddo
      endif

!::flux
      if( trim(ctyp) == "flux" ) then
        write(i6,'(/2x,"*** lstflos ***",2x,a)') trim(ctyp)
        write(i6,'(2x,"wfac =",f10.3,"  itim =",i6,"  icalZ =",i3, 2x,
     >   a,"  npmax =",i5,"  sflux =",1pe12.3,
     >   "  swtot =",1pe12.3,2x,a)')
     >   wfac, itim, icalZ, sptyc, npmax, sflux, swtot, trim(sptyc)

      if( impmc_model == 0 ) then
        write(i6,'(2x,4x,"dene",48x,"etne")')
        write(i6,'(2x,(2x,1p5e10.3,3x,1p5e10.3))')
     >      edne(1:5), etne(1:5)
        write(i6,'(2x,4x,"deni",48x,"etni")')
        write(i6,'(2x,(2x,1p5e10.3,3x,1p5e10.3))')
     >      edni(1:5), etni(1:5)
      endif

      write(i6,'(2x,4x,"flxI",48x,"flxO",48x,"denZ")')

      fac = wfac
      if( impmc_model == 0 ) fac = sflux / swtot * wfac

        write(i6,'(2x," iz",1p5e10.3,3(2x,1p5e10.3))')
     >    fls_roh(1:5),fls_roh(1:5),fls_roh(1:5)
        do iz = 0, ismax
          write(i6,'(2x,i3,1p5e10.3,3(2x,1p5e10.3))') iz,
     >      fac*fls_flxI(iz,1:5), fac*fls_flxO(iz,1:5),
     >      wfac*ednz(iz,1:5)
        enddo
      endif

!::<Nz(iz)> at core edge
      if( trim(ctyp) == "denZ" ) then
        write(i6,'(/2x,"*** lstflos ***",2x,a)') trim(ctyp)
        write(i6,'(2x,"wfac =",f10.3,"  itim =",i6,"  icalZ =",i3,2x,a,
     > "  npmax =",i5,"  sflux =",1pe12.3,"  swtot =",1pe12.3,2x,a)')
     >   wfac, itim, icalZ, sptyc, npmax, sflux, swtot, trim(sptyc)

        write(i6,'(2x,3x,1p5e10.3)') fls_roh(1:5)
        do iz = 0, ismax
          write(i6,'(2x,i3,1p5e10.3)') iz, wfac*ednz(iz,1:5)
        enddo
      endif

      return
      end
