!***********************************************************************
      subroutine test_eqdsk_efit
!***********************************************************************
      implicit none
      integer  nft
      character  dsn*80
!
      nft = 10
      dsn = "/home1/jt60u/MESH/jt60/G45330a/EFIT/EFIT_E45330.dat"
!
      open(unit=nft,file=dsn)
      call eqdsk_efit(nft,"read")
!
      stop
      end
!
!***********************************************************************
      subroutine eqdsk_efit( nft, cwr )
!***********************************************************************
!
!        original /analysis/src/selene/efitout.f    2009/2/18
!
!-----------------------------------------------------------------------
!  Variables
!-----------------------------------------------------------------------
!  CASE   : Identification character string.
!  NW     : Number of horizontal R grid points.
!  NH     : Number of vertical Z grid points.
!  RDIM   : Horizontal dimension in meter of computational box.
!  ZDIM   : Vertical dimension in meter of computational box.
!  RLEFT  : Minimun R in meter of rectanglar computational box.
!  ZMID   : Z of center of computational box in meter.
!  RMAXIS : R of magnetic axis in meter.
!  ZMAXIS : Z of magnetic axis in meter.
!  SIMAG  : Poloidal flux at magnetic axis in Weber/rad.
!  SIBRY  : Poloidal flux at the plasnma boundary in Weber/rad.
!  RCENTR : R in meter of vacuume troidal magnetic field BCENTR.
!  BCENTR : Vacuume troidal magnetic field in Tesla at RCENTRA.
!  CURRENT: Plasma current in Ampere.
!  FPOL   : Poloidal current function in m-T, F=R*Bt on flux grid.
!  PRES   : Plasma pressure in nt/m2 on uniform flux grid.
!  FFPRIM : FF'(psi) in {(mT)^2/(Weber/rad)} on uniform flux grid.
!  PPRIME : P'(psi) in {(mT)^2/(Weber/rad)} on uniform flux grid.
!  PSIZR  : Poloidal flux in Waber/rad on the rectangular grid points
!  QPSI   : Q values on uniform flux grid from axis to boundary
!  NBBBS  : Number of boundary points
!  LIMITR : Number of limiter points
!  RBBBS  : R of boundary points in meter
!  ZBBBS  : Z of boundary points in meter
!  RLIM   : R of surrounding limiter contour in meter
!  ZLIM   : Z of surrounding limiter contour in meter
!-----------------------------------------------------------------------
!  KVTOR  : Toroidal rotation switch
!  RVTOR  : Toroidal rotation characteristic major radius in m
!  NMASS  : Mass density switch
!  PRESSW : Rotational pressure Pomega(psi) in N/m^2
!  PWPRIM : P'omega(psi) in (N/m^2)/(Web/rad)
!  DMION  : Mass density on uniform poloidal flux grid
!  RHOVN  : Normalized toroidal flux on uniform poloidal flux grid
!-----------------------------------------------------------------------
!  Toroidal Current Density Jt
!     Jt(A/m^2) = R P'(psi) + FF'(psi)/R/u0
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  Edit points
!   (5)  Press(j) = PRV(j)/cmu0  cmu0 = 1.0d-6 (Bug) ==> 4*pi*1.0d-6
!   (9)  epsiw=(siw(nv)-siw(1))/nf (Bug) ==> deps = xxx/dfloat(nf-1)
!
!   cnmu = 4.0*pi*1.0d-7 = 4.0*pi/10.0d0*1.0d-6 = 1.2566d-6
!
!                                          2009/02/17  K. Shimizu
!-----------------------------------------------------------------------
      use com_eqdat,   only : dr, dz, ndr, ndz, nr, nz, paxs, psep, psi
     >    , raxs, rbt0, rg, rmax, rmin, rsep, zaxs, zg, zmax, zmin, zsep
      use cunit,       only : lmspe, lmype, n6
      implicit none
!
!::argument
      integer,   intent(in) :: nft
      character, intent(in) :: cwr*(*) ! dummy
!
!
!::local variables
      integer  :: i, j, k
      integer  :: ishot, it
      real*8   :: eqtime
      character :: xdbnm*8, xdate*8, xdver*8
!
!local EFIT side
      integer, parameter :: irdm = ndr, izdm = ndz
      integer, parameter :: ifdm = irdm
      integer  :: nw, nh, nbbbs, limitr
      real*8   :: rdim, zdim, rcentr, rleft, zmid
      real*8   :: rmaxis, zmaxis, simag, sibry, bcentr, current
!
      character  case(6)*8
      real*8, dimension(ifdm) :: fpol, pres, ffprim, pprime, qpsi
      real*8, dimension(irdm,izdm) :: psirz
      real*8   :: rbbbs(1000), zbbbs(1000)
      real*8   :: rlim(100), zlim(100)
      integer  :: idum, nf
      real*8   :: xdum
!
      if( lmype.eq.lmspe ) then
      write(n6,'(/2x,"*** eqdisk EFIT (G EQDISK)*** (read)")')
!-----------------------------------------------------------------------
!::EFIT-data
!-----------------------------------------------------------------------
      read(nft,2000) (case(i),i=1,6), idum, nw, nh
      nf = nw
!
      xdbnm = case(1)
      xdate = case(2)
      read(case(3),'(i8)') ishot
      read(case(4),'(i6,"ms")') it
      eqtime  = dfloat(it)*1.0d-3
      xdver = case(5)
!
!::dimension check
      write(n6,'(2x,"*** eqdsk_efit ***")')
      write(n6,'(2x,"dbnm = ",a,"  date = ",a,"  ver =",a)')
     >    xdbnm, xdate, xdver
      write(n6,'(2x,"ishot =",i6,"  time =",f9.3)') ishot, eqtime
!
      write(n6,'(2x,"nw =",i4," <= ",i4,"  nh =",i4," <= ",i4,
     >  "  nf =",i4," <= ",i4)')  nw, irdm, nh, izdm, nf, ifdm
      if( nw.gt.irdm .or. nh.gt.izdm .or. nf.gt.ifdm ) goto 910
!
      read(nft,2020)  rdim, zdim, rcentr, rleft, zmid
      read(nft,2020)  rmaxis, zmaxis, simag, sibry, bcentr
      read(nft,2020)  current, simag, xdum, rmaxis, xdum
      read(nft,2020)  zmaxis, xdum, sibry, xdum, xdum
!
      read(nft,2020) (fpol(i),i=1,nf)
      read(nft,2020) (pres(i),i=1,nf)
      read(nft,2020) (ffprim(i),i=1,nf)
      read(nft,2020) (pprime(i),i=1,nf)
      read(nft,2020) ((psirz(i,j),i=1,nw),j=1,nh)
      read(nft,2020) (qpsi(i),i=1,nf)
!
      read(nft,2022)  nbbbs,limitr
      write(n6,'(2x,"nbbbs  =",i4," <= ",i4,"  limitr =",i4," <= ",
     >  i4)')  nbbbs, 1000,limitr, 100
      if( nbbbs.gt.1000 .or. limitr.gt.100 ) goto 920
!
      read(nft,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
      read(nft,2020) (rlim(i),zlim(i),i=1,limitr)
!
!-----------------------------------------------------------------------
!::grid
!-----------------------------------------------------------------------
      nr = nw
      nz = nh
      write(n6,'(2x,"nr = ",i4," <= ",i4,"  nz = ",i4," <= ",i4)')
     >   nr, ndr, nz, ndz
      if( nr.gt.ndr .or. nz.gt.ndz ) goto 930
!
      dr = rdim/dfloat(nr-1)
      dz = zdim/dfloat(nz-1)
!
      rmin = rleft
      do i = 1, nr
      rg(i) = rmin + dr*dfloat(i-1)
      enddo
!
      zmin = zmid - 0.5d0*zdim
      do j = 1, nz
      zg(j) = zmin + dz*dfloat(j-1)
      enddo
!
!-----------------------------------------------------------------------
!::psi
!-----------------------------------------------------------------------
      do j = 1, nz
      do i = 1, nr
      k = (j-1)*nr+i
      psi(k) = psirz(i,j)
      enddo
      enddo
!
!-----------------------------------------------------------------------
!::others
!-----------------------------------------------------------------------
      rbt0 = rcentr*bcentr
      raxs = rmaxis
      zaxs = zmaxis
      paxs = simag
!
      rsep = 0.0d0
      zsep = 0.0d0
      psep = sibry
!
      endif

      rmin = rg(1)
      rmax = rg(nr)
      zmin = zg(1)
      zmax = zg(nz)

!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
!      write(n6,'(2x)')
!      write(n6,'(2x,"dr  =",1pe12.4,"  dz  =",1pe12.4)') dr, dz
!      write(n6,'(2x,"rg  =",1p3e12.4," ...",1p3e12.4)')
!     >  rg(1), rg(2), rg(3), rg(nr-2), rg(nr-1), rg(nr)
!      write(n6,'(2x,"zg  =",1p3e12.4," ...",1p3e12.4)')
!     >  zg(1), zg(2), zg(3), zg(nz-2), zg(nz-1), zg(nz)
!      write(n6,'(2x,"psi =",1p3e12.4," ...",1p3e12.4)')
!     >  psi(1), psi(2), psi(3), psi(nr*nz-2), psi(nr*nz-1), psi(Nr*nz)
!      write(n6,'(2x,"rbt0=",1pe12.4)') rbt0
!      write(n6,'(2x,"raxs=",1pe12.4,"  zaxs =",1pe12.4,"  paxs =",
!     >  1pe12.4)') raxs, zaxs, paxs
!      write(n6,'(2x,"rsep=",1pe12.4,"  zsep =",1pe12.4,"  psep =",
!     >  1pe12.4)') rsep, zsep, psep
!
      return
!
!::format
 2000 format(6a8,3i4)
 2020 format(5e16.9)
 2022 format(2i5)
 2024 format(i5,e16.9,i5)
!
 910  continue
 920  continue
 930  continue
      write(n6,'(2x,"dimension error")')
      end
