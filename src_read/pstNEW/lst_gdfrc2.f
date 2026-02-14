!**********************************************************************
      subroutine lst_gdfrc2(nty, it)
!**********************************************************************
      use csize
      use cimcom
      use cimden
      use cntcom
      use cntpls
      use cplcom
      use cplmet
      use cphcns
      use cunit
      use csonic
!      use com_gdfrc, only : dlnx
      implicit none
!
!::argument
      integer  it, nty
!
!::local varibales
      character  dsn*80
      real*8   frez(ndx),frte(ndx),frti(ndx),flge(ndx),flgi(ndx)
      integer  jt, j, i, jte, ic, imox, imoy, i6
      real*8   fjt
      integer  ia, iz, k, kmax
      real*8   zi, cftz, tauz, ftot, ffmx, zvz, rmfp
      real*8   hlge, hlgi, hgte, hgti, rmfpe, rmfpi, rlnti, rlnte
      real(8) :: grdti, grdte, zcti
      real(8) :: alnx(ndx)
      real(8) :: zlmx, sgnti, sgnte
      real(8), dimension(ndx, 0:ndmis) :: hvnz
      real(8), dimension(ndx) :: hvtz
      integer  jtp, jtm, jp, jm, icp, icm
      real(8) :: gdpzp, gdpzm
!
      character(80) :: drimp, drFRC
      character(10) :: ctub, cwxdr
      integer  ipwxdr, ipsla
      real(8), dimension(150) :: Wtauz, Wftot, Wffmx
      integer, dimension(150) :: Wiz
      integer :: lp = 0
      save :: lp, kmax, Wiz, drimp
!
!::input data
!      n6 = 6
      lp = lp + 1
      if( lp.eq.1 ) then
      write(n6,'(/2x,"*** gdfrc ***  set Wiz")')
      kmax = 25
      if( ismax.le.kmax ) kmax = ismax
      do i = 1, kmax
      Wiz(i) = i
      enddo
      write(n6,'(2x,"Wiz =",10i4)') (Wiz(i),i=1,kmax)
!
      call getenv("IMP1D",drimp)
!x      write(n6,'(4x,"drimp  = ",a)') trim(drimp)
      drFRC=trim(drimp) // "/FRC"
      call system("mkdir "//trim(drFRC))
      endif
!::constant
      jte = jtmax(it)
      ia  = 1
!xx   iz  = 1
!xx   zi  = iz
      cftz = 1.0d0/dsqrt(aimas)*azmas/(azmas+aimas)
      call plenc(it,alnx)
      zlmx = alnx(jte)
!::define hvnz(j,iz) it fixed
      do jt=1,jte
        j = jcel(jt,it)
        i = icel(jt,it)
        ic = mcel(j,i)
        if( ic .eq. 0) cycle
        hvnz(j,0:ismax) = tdnz(0:ismax,ic)
        hvtz(j) = vti(j,i)
      enddo
!
!::charcteristic length
!      call mfpath(j,i,rmfpe,rmfpi)
!      rlnti = vti(j,i)/dabs(frti(jt))
!      rlnte = vte(j,i)/dabs(frte(jt))
!
!::charge state
      do k = 1, kmax
      iz = Wiz(k) 
      zi = iz
!
!::file name
      call getenv("IMP1D",drimp)
!x      write(n6,'(4x,"drimp  = ",a)') trim(drimp)
      write(ctub,'(i2.2)') it
      ipwxdr = index(drimp,"wxdr_")
      ipsla  = index(drimp(ipwxdr:80),"/")
      cwxdr=drimp(ipwxdr+5:ipsla)
!      dsn = trim(drimp) // "/it" // trim(ctub) // "FRC,txt"
!      write(dsn,'(a,a)')
!     >  "./FRC", cwxdr 
!      call system("mkdir " // trim(dsn) ) 
!      write(dsn,'(a,a,a,a)')
!     >  "./FRC", cwxdr , "/it", trim(ctub) 
!      call system("mkdir " // trim(dsn) ) 
!      write(dsn,'(a,a,a,a,a,a,a,a,i2.2,a,i2.2,a)')
!     >  trim(drimp), "/FRC", cwxdr , "/it", trim(ctub), 
!     >  "wx", cwxdr, "it", it, "zz", k, ".txt"
      write(dsn,'("/frc_it",i2.2,"_zz",i2.2,".txt")')
     >   it,k
      call delspc(dsn)
!
!      write(n6,'(4x,"dsn = ",a)') dsn
!      call flush(n6)
!      return
!
      i6 = 91
      open(unit=i6,file=trim( trim(drFRC) // trim(dsn)))
!
!::header
      write(i6,'(2x,"*** gdfrc2 ***   it =",i3,"  time =",
     >  1pe16.8)') it, time
!      write(i6,'(2x,"cfez  =",10(2x,i2.2,1pe10.3))')
!     >   (Wiz(k),cfez(Wiz(k)),k=1,kmax)
!      write(i6,'(2x,"cfgTe =",10(2x,i2.2,1pe10.3))')
!     >   (Wiz(k),cfgte(Wiz(k)),k=1,kmax)
!      write(i6,'(2x,"cfgTi =",25(2x,i2.2,1pe10.3))')
!     >   (Wiz(k),cfgti(Wiz(k)),k=1,kmax)
!
      write(i6,'(4x,"jt",3x,"j",3x,"i",2x,"Lp",9x,"Lp2",8x,
     >  "Ni",9x,"Ti",9x,"Ne",9x,"Te",9x,"Nz",9x,
     >  "mfpi",7x,"LTi",8x,"mfpe",7x,"LTe",8x,
     >  "sgTi",1x,"grdTi",6x,"gdTi",7x,
     >  "sgTe",1x,"grdTe",6x,"gdTe",7x,
     >  "tauz",i2.2,5x,"T",i2.2,2x,"ftot",i2.2,5x,
     >                 "F",i2.2,2x,"ffmx",i2.2,4x,
     >                 "Pp",i2.2,2x,"gdpzp",i2.2,3x,
     >                 "Pm",i2.2,2x,"gdpzm",i2.2,4x)')
     >  k,k,k,k,k,k,k,k,k
!
!::output
      do jt=1,jte
      j  = jcel(jt,it)
      i  = icel(jt,it)
      ic  = mcel(j,i)
      jtp = min0(jt+1,jte)
      jtm = max0(jt-1,1)
      jp  = jcel(jtp,it)
      jm  = jcel(jtm,it)
      icp = mcel(jp,i)
      icm = mcel(jm,i)
      if(ic.eq.0) cycle
      fjt = jt
!::thermal force & its limit
      call gdfrc(it,frez,frte,frti,flge,flgi)
!
!::new limit of temperature gradient
      zcti  = cfgti(iz)
      grdti = dmin1( glti(ic), dabs(gdti(ic)) )
      grdti = dsign( grdti, gdti(ic) )
      grdte = dmin1( glte(ic), dabs(gdte(ic)) )
      grdte = dsign( grdte, gdte(ic) )
!
!::thermal force
      ftot  = cfez(iz)*gdez(ic)+cfgte(iz)*grdte+zcti*grdti
!xx   ftot  = cfez(iz)*gdez(ic)+cfgte(iz)*gdte(ic)+cfgti(iz)*gdti(ic)
!
!::frctition force
      tauz = 1.0d0/(zi**2*slw0(ic,2)*slnv(ic,2))*cftz
      ffmx = amz*vva(j,i,ia)/tauz/cev
!
!::save
!      Wtauz(jt) = tauz
!      Wftot(jt) = ftot
!      Wffmx(jt) = ffmx
!
!::new limit of temperature gradient
      zcti  = cfgti(iz)
      grdti = dmin1( glti(ic), dabs(gdti(ic)) )
      grdti = dsign( grdti, gdti(ic) )
      grdte = dmin1( glte(ic), dabs(gdte(ic)) )
      grdte = dsign( grdte, gdte(ic) )
      call mfpath(j,i,rmfpe,rmfpi)
      rlnti = vti(j,i)/dabs(frti(jt))
      rlnte = vte(j,i)/dabs(frte(jt))
      sgnti = dsign(1.0d0,grdti)
      sgnte = dsign(1.0d0,grdte)
!::grad(pz) = grad(Nz*Ti)
! delete 1 line, if state because dlnx = 0 always
!      if( dlnx(jt).eq.0.0 .or. hvnz(j,iz).eq.0.0 ) then
        gdpzp = 0.0d0
        gdpzm = 0.0d0
! delete 1 line, if state because dlnx = 0 always
!      else
!        gdpzp = -(hvnz(jp,iz)*hvtz(jp)-hvnz(j,iz)*hvtz(j))/dlnx(jt)
!     >          /hvnz(j,iz)
!        gdpzm = -(hvnz(j,iz)*hvtz(j)-hvnz(jm,iz)*hvtz(jm))/dlnx(jt)
!     >          /hvnz(j,iz)
!      endif
      write(i6,'(2x,3i4,1p11e11.3,0pf5.1,1p2e11.3,0pf5.1,1p2e11.3,
     >  1pe11.3,2(0pf5.1,1pe11.3,0pf5.1,1pe11.3))')
     >  jt, j, i, alnx(jt), zlmx-alnx(jt),
     >  vni(j,i), vti(j,i), vne(j,i), vte(j,i), hvnz(j,iz),
     >  rmfpi, rlnti, rmfpe, rlnte, 
     >  sgnti,  dabs(grdti), dabs(gdti(ic)),
     >  sgnte,  dabs(grdte), dabs(gdte(ic)),
     >  tauz, dsign(1.0d0, ftot),      dabs(ftot),
     >        dsign(1.0d0, ffmx),      dabs(ffmx),
     >        dsign(1.0d0, gdpzp)    , dabs(gdpzp),
     >        dsign(1.0d0, gdpzm)    , dabs(gdpzm)
!
      enddo  ! jt
      close(i6)
      enddo  ! k
!
      return
      end
