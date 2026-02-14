!***********************************************************************
      subroutine test_eqdisk
!xx      program test_eqdisk
!***********************************************************************
      implicit none
!
!::local variables
      character  dsn*20
      integer    nft
!
      nft = 25
      dsn = "dat.eqil"
!
      open( unit=nft, file=dsn, form="formatted", action="read" )
      call eqdisk( nft, "read" )
!
!x      nft = 26
!x      dsn = "@newEQ.dat"
!x      open( unit=nft, file=dsn, form="formatted", action="write" )
!x      call eqdisk( nft, "write" )
!
      stop
      end
!
!***********************************************************************
      subroutine eqdisk( nft, cwr )
!***********************************************************************
!
!       nft           : unit number of equil data
!       cwr           : "write" or "read"
!
!   include 'com_eqdat'
!
!       nr            : r-mesh number
!       nz            : z-mesh number
!       dr            : width of r-mesh [m]
!       dz            : width of z-mesh [m]
!       rg(ndr)       : r-grid  [m]
!       zg(ndz)       : z-grid  [m]
!       psi(ndr*ndz)  : poloidal flux function [Wb]
!       rbt0          : Rp*Bt   [m*T]
!       raxs          : r-poistion of axis [m]
!       zaxs          : z-position of axis [m]
!       paxs          : poloidal flux function at axis [Wb]
!       rsep          : r-position of X-point [m]
!       zsep          : z-position of X-point [m]
!       psep          : poloidal flux function at X-point [Wb]
!
!-----------------------------------------------------------------------
      use com_eqdat,   only : dr, dr0, dz, dz0, ndr, ndz, nr, nz, paxs
     >    , psep, psi, rg, raxs, rbt0, rmax, rmin, rsep, zaxs, zg, zmax
     >    , zmin, zsep
      use cunit,       only : lmspe, lmype, n6
      implicit none
!
      integer,   intent(in) :: nft
      character, intent(in) :: cwr*(*)
!
!::local variables
      character dsn*80
      integer i, imax, j, mj
!
      character cver*30, cfmtI*16, cfmtR*16, clin*120
      character cn(10)*1
      real*8    frax, fzax, frsp, fzsp
      data   cn/"1","2","3","4","5","6","7","8","9","0"/
      data cver  /"1.01 (03/12/10)"/; save cver
! function
      integer    lenx
!
      inquire( unit=nft, name=dsn )
!
!----------------------------------------------------------------------
!::write
!----------------------------------------------------------------------
      if( cwr(1:1).eq."w" ) then
      if( lmype.eq.lmspe ) return
      dr0 = dr
      dz0 = dz
      rmin = rg(1)
      rmax = rg(nr)
      zmin = zg(1)
      zmax = zg(nz)
!
      write(n6,'(/2x,"*** eqdisk *** (write)")')
!-----
      cfmtI = "(2i5)"
      cfmtR = "(1p5e15.7)"
      write(n6,'(2x,"version = ",a)') cver
      write(n6,'(2x,"formatI = ",a)') cfmtI
      write(n6,'(2x,"formatR = ",a)') cfmtR
!-----
      write(nft,'("* ",a)') cver
      write(nft,'("*",80a1)') (cn(mod(i-1,10)+1),i=2,75)
      write(nft,'("* formatI  ",a)') cfmtI
      write(nft,'("* formatR  ",a)') cfmtR
      write(nft,'("* data")')
      write(nft,cfmtI) nr, nz
      write(nft,cfmtR) rmin, rmax, zmin, zmax, dr0, dz0
      write(nft,cfmtR) rbt0
      write(nft,cfmtR) raxs, zaxs, paxs
      write(nft,cfmtR) rsep, zsep, psep
      write(nft,cfmtR) (psi(i),i=1,nr*nz)
      close(nft)
!
!::debug write
      write(n6,'(2x,"nft  =",i3,"  dsn =",a)') nft, dsn(1:lenx(dsn))
      write(n6,'(2x,"nr   =",2i5,"  nz =",2i5)') nr, ndr, nz, ndz
      return
      endif
!
!----------------------------------------------------------------------
!::read
!----------------------------------------------------------------------
      if( cwr(1:1).eq."r" ) then
      if( lmype.eq.lmspe ) then
      write(n6,'(/2x,"*** eqdisk *** (read)")')
!-----
      do i = 1, 20
      read(nft,'(a)') clin
      write(n6,'(2x,"[",a,"]")') trim(clin)

      call lincut(clin,"version",mj)
      if( mj > 0 ) cver = clin(mj:)
      call lincut(clin,"formatI",mj)
      if( mj > 0 ) cfmtI = clin(mj:)
      call lincut(clin,"formatR",mj)
      if( mj > 0 ) cfmtR = clin(mj:)
      if( index(clin,"data") > 0 ) exit
      enddo

      cfmtI = "("//cfmtI(1:lenx(cfmtI))//")"
      cfmtR = "("//cfmtR(1:lenx(cfmtR))//")"
      write(n6,'(2x,"version = ",a)') cver
      write(n6,'(2x,"formatI = ",a)') cfmtI
      write(n6,'(2x,"formatR = ",a)') cfmtR
!-----
!
      read(nft,'(2i5)') nr, nz
      read(nft,cfmtR) rmin, rmax, zmin, zmax, dr0, dz0
      read(nft,cfmtR) rbt0
      read(nft,cfmtR) raxs, zaxs, paxs
      read(nft,cfmtR) rsep, zsep, psep
      read(nft,cfmtR) (psi(i),i=1,nr*nz)
      close(nft)
!
      dr = (rmax-rmin)/dfloat(nr-1)
      dz = (zmax-zmin)/dfloat(nz-1)
      do i = 1, nr
      rg(i) = rmin + dr*dfloat(i-1)
      enddo
      do j = 1, nz
      zg(j) = zmin + dz*dfloat(j-1)
      enddo
      endif

      rmin = rg(1)
      rmax = rg(nr)
      zmin = zg(1)
      zmax = zg(nz)
!
      endif
!
!::debug write
      write(n6,'(2x,"nft  =",i3,"  dsn =",a)') nft, dsn(1:lenx(dsn))
      write(n6,'(2x,"nr   =",2i5,"  nz =",2i5)') nr, ndr, nz, ndz
      if( nr.gt.ndr .or. nz.gt.ndz ) goto 910
      write(n6,'(2x,"rmin =",1pe12.4,"  rmax =",1pe12.4,"  zmin =",
     >  1pe12.4,"  zmax =",1pe12.4)') rmin, rmax, zmin, zmax
!
      imax = nr*nz
      frax = (raxs-rmin)/dr+1.0d0
      fzax = (zaxs-zmin)/dz+1.0d0
      frsp = (rsep-rmin)/dr+1.0d0
      fzsp = (zsep-zmin)/dz+1.0d0
      if( raxs.eq.99.0 .or. zaxs.eq.99.0 ) then
        frax = 0; fzax = 0; endif
      if( rsep.eq.99.0 .or. zsep.eq.99.0 ) then
        frsp = 0; fzsp = 0; endif
!
      write(n6,'(2x,"dr   =",1p2e12.4,"  eps-r =",1pe12.4)')
     >      dr,dr0,dr-dr0
      write(n6,'(2x,"dz   =",1p2e12.4,"  eps-z =",1pe12.4)')
     >      dz,dz0,dz-dz0
      write(n6,'(2x,"rbt0 =",1pe12.4)') rbt0
      write(n6,'(2x,"raxs =",1pe12.4,"  zaxs =",1pe12.4,"  paxs =",
     >   1pe12.4,"  ir,iz =",0p2f9.4)') raxs, zaxs, paxs, frax, fzax
      write(n6,'(2x,"rsep =",1pe12.4,"  zsep =",1pe12.4,"  psep =",
     >   1pe12.4,"  ir,iz =",0p2f9.4)') rsep, zsep, psep, frsp, fzsp
      write(n6,'(2x,"psi(i) i = 1,2,  =",5(1x,i6,1pe12.4))')
     >   (i,psi(i),i=1,5)
      write(n6,'(2x,"psi(i) i = (imax)=",5(1x,i6,1pe12.4))')
     >   (i,psi(i),i = imax-4, imax)
!
      return

 910  continue
      call wexit("eqdisk","dimension error")
      end
