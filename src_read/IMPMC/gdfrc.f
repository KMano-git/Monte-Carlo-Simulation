!**********************************************************************
!      subroutine test_gdfrc
!**********************************************************************
!
!      call set_forc
!      call set_gcoef
!
!      do it = itsls+2, itmps+5
!      if( it.eq.itpvs ) cycle
!      call lst_frc(it)
!      enddo
!
!      return
!      end
!
!**********************************************************************
      subroutine gdfrc(it,frez,frte,frti,flge,flgi)
!**********************************************************************
!
!        it   :  tube number
!        frez :  electric field
!        frte :  Te gradient         [eV/m]
!        frti :  Ti gradient         [eV/m]
!        flge :  limit of Te gradient      0.3*Te/(Vthe*Taue)
!        flgi :  limit of Ti gradient      0.3*Ti/(Vthi*Taui)
!
!     Limit of Thermal Force for long mean free path
!
!       rmfp < 0.3*LTi = 0.3*Ti/(dTi/ds)
!       dTi/ds < 0.3*Ti/rmfp  [eV/m]      rmfp = vthi*Taui
!
!       gradTi = dmin1( dTi/ds,   0.3*Ti/(Vthi*Taui) )
!       gradTi = dmin1( gdti(ic), glti(ic) )
!
!----------------------------------------------------------------------
      use cplcom, only : vne, vte, vti
      use cplmet, only : hwtm, hwtp, icel, itmax, itmps, jcel, jtmax
      use csize,  only : ndx, ndy
      implicit none
!
!::argument
      integer, intent(in)  :: it
      real(8), intent(out) :: frez(ndx), frte(ndx), frti(ndx)
      real(8), intent(out) :: flge(ndx), flgi(ndx)
!
!::local common
      real*8  dlnx(ndx), zpes(ndx), ztes(ndx), ztis(ndx)
!
!::local
      real*8    gpe, gte, gti, zne, rmfpe, rmfpi
      integer  jte, jt, j, i, jp
!
!::length in cell along the magnetic field line
      call gdlen(it,dlnx)
!
!::efield,thermal force
      frez(1:ndx) = 0.0d0
      frte(1:ndx) = 0.0d0
      frti(1:ndx) = 0.0d0
      flge(1:ndx) = 0.0d0
      flgi(1:ndy) = 0.0d0
!
      if( it.le.0 .or. it.gt.itmax ) return
!
!::Ne,Te,Ti at cell boundary
      zpes(1:ndx) = 0.0d0
      ztes(1:ndx) = 0.0d0
      ztis(1:ndx) = 0.0d0
!
      jte = jtmax(it)
      do jt = 1, jte-1
      j  = jcel(jt,it)
      i  = icel(jt,it)
      jp = jcel(jt+1,it)
      zpes(jt) = hwtp(j,i)*vne(j,i)*vte(j,i)
     >          +hwtm(j,i)*vne(jp,i)*vte(jp,i)
      ztes(jt) = hwtp(j,i)*vte(j,i)+hwtm(j,i)*vte(jp,i)
      ztis(jt) = hwtp(j,i)*vti(j,i)+hwtm(j,i)*vti(jp,i)
!           vte[joule/cev = eV], vti[eV], zpes(jt) [pressure/cev]
!          [pressure/cev] = [m^(-3) eV]
      enddo
!
!::gradient at cell center
      do jt = 1, jte
      j  = jcel(jt,it)
      i  = icel(jt,it)
      zne = vne(j,i)
!
      if( jt.eq.1 .or. jt.eq.jte ) then
      else
      gpe = (zpes(jt)-zpes(jt-1))/dlnx(jt)
!           gpe [ (pressure/cev) /m]
      gte = (ztes(jt)-ztes(jt-1))/dlnx(jt)
      gti = (ztis(jt)-ztis(jt-1))/dlnx(jt)
!           gte [eV/m], gti [eV/m]
      frez(jt) = - gpe/zne - 0.71*gte
!     cf. Eq. (2.2-11) in the file
!           frez is steady paralle electric field
!     derived from the electron momentum equation
!     under assumption that
!     (1) steady state d/dt = 0,
!     (2) mass m_e -> 0 limit.
      frte(jt) = gte
      frti(jt) = gti
!           frte [eV/m], frti [eV/m]
!
!::Maximum gradients based on e-e and i-i MFP
      call mfpath(j,i,rmfpe,rmfpi)
      flge(jt) = 0.3d0*vte(j,i)/rmfpe
      flgi(jt) = 0.3d0*vti(j,i)/rmfpi
      endif
      enddo
!
!::cyclic condition
      if( it.ge.itmps ) then
      frez(1) = frez(jte-1)
      frte(1) = frte(jte-1)
      frti(1) = frti(jte-1)
      flge(1) = flge(jte-1)
      flgi(1) = flgi(jte-1)
!
      frez(jte) = frez(2)
      frte(jte) = frte(2)
      frti(jte) = frti(2)
      flge(jte) = flge(2)
      flgi(jte) = flgi(2)
      endif
!
      return
      end
!
!**********************************************************************
      subroutine lst_gdfrc(it)
!**********************************************************************
      use cimcom, only : aimas, amz, azmas, cfez, cfgte, cfgti, gdez
     >    , gdte, gdti, glte, glti, ismax, lfgtzf, slnv, slw0
      use cntcom, only : mcel
      use cntpls, only : zefm
      use cphcns, only : cev
      use cplcom, only : vni, vte, vti, vva
      use cplmet, only : icel, jcel, jtmax
      use csize,  only : ndx
      use csonic, only : time
      use cunit,  only : n6
      implicit none
!
!::argument
      integer, intent(in) :: it
!
!::local varibales
      character  dsn*80
      real*8   frez(ndx),frte(ndx),frti(ndx),flge(ndx),flgi(ndx)
      integer  jt, j, i, jte, ic, i6
      integer  ia, iz, k, kmax
      real*8   cftz, tauz, ftot, ffmx
      real*8   rmfpe, rmfpi, rlnti, rlnte
      real(8) :: grdti, grdte
      real(8) :: zfca, zcti, zcte
      real(8) :: alnx(ndx)
      real(8) :: zlmx, sgnti, sgnte
!
      character(80) :: drimp
      character(10) :: ctub
      real(8), dimension(50) :: Wtauz, Wftot, Wffmx
      integer, dimension(50) :: Wiz
      integer :: lp = 0
      save :: lp, kmax, Wiz, drimp
!
!::input data
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
      write(n6,'(4x,"drimp  = ",a)') trim(drimp)
      endif
!
!::file name
      write(ctub,'(i2)') it
      dsn = trim(drimp) // "/it" // trim(ctub) // "FRC.txt"
      call delspc(dsn)
!
      i6 = 21
      open(unit=i6,file=dsn)
!
!::thermal force & its limit
      call gdfrc(it,frez,frte,frti,flge,flgi)
!
!::header
      write(i6,'(2x,"*** gdfrc ***   it =",i3,"  time =",
     >  1pe16.8)') it, time
      write(i6,'(2x,"cfez  =",10(2x,i2.2,1pe10.3))')
     >   (Wiz(k),cfez(Wiz(k)),k=1,kmax)
      write(i6,'(2x,"cfgTe =",10(2x,i2.2,1pe10.3))')
     >   (Wiz(k),cfgte(Wiz(k)),k=1,kmax)
      write(i6,'(2x,"cfgTi =",25(2x,i2.2,1pe10.3))')
     >   (Wiz(k),cfgti(Wiz(k)),k=1,kmax)
!
      write(i6,'(4x,"jt",3x,"j",3x,"i",3x,"Lp",9x,"Lp2",8x,
     >  "Ni",9x,"Ti",9x,"Te",9x,
     >  "mfpi",7x,"LTi",8x,"mfpe",7x,"LTe",8x,
     >  "sgTi",1x,"grdTi",6x,"gdTi",7x,
     >  "sgTe",1x,"grdTe",6x,"gdTe",7x,
     >  25("tauz",i2.2,5x,"T",i2.2,2x,"ftot",i2.2,5x,
     >                    "F",i2.2,2x,"ffmx",i2.2,5x))')
     >  (Wiz(k),Wiz(k),Wiz(k),Wiz(k),Wiz(k),k=1,kmax)
!
!::constant
      jte = jtmax(it)
      ia  = 1
      cftz = 1.0d0/dsqrt(aimas)*azmas/(azmas+aimas)
      call plenc(it,alnx)
      zlmx = alnx(jte)
!
!::loop (jt)
      do jt = 1, jte
      j  = jcel(jt,it)
      i  = icel(jt,it)
      ic  = mcel(j,i)
      if(ic.eq.0) cycle
!
!::charcteristic length
      call mfpath(j,i,rmfpe,rmfpi)
      rlnti = vti(j,i)/dabs(frti(jt))
      rlnte = vte(j,i)/dabs(frte(jt))
!
!::charge state
      do k = 1, kmax
      iz = Wiz(k)
      zfca = 1.0d0
      if(lfgtzf.eq.1)  zfca = 1.0d0/iz**2*(iz/zefm(ic)-1.0d0)*iz
      zcti = cfgti(iz)*zfca
      zcte = cfgte(iz)*zfca
!
!::new limit of temperature gradient
      grdti = dmin1( glti(ic), dabs(gdti(ic)) )
      grdti = dsign( grdti, gdti(ic) )
      grdte = dmin1( glte(ic), dabs(gdte(ic)) )
      grdte = dsign( grdte, gdte(ic) )
!
!::thermal force
      ftot  = cfez(iz)*gdez(ic)+zcte*grdte+zcti*grdti
!
!::frictition force
      tauz = 1.0d0/(iz**2*slw0(ic,2)*slnv(ic,2))*cftz
      ffmx = amz*vva(j,i,ia)/tauz/cev
!
!::save
      Wtauz(k) = tauz
      Wftot(k) = ftot
      Wffmx(k) = ffmx
      enddo
!
!::output
      sgnti = dsign(1.0d0,grdti)
      sgnte = dsign(1.0d0,grdte)
      write(i6,'(2x,3i4,1p9e11.3,0pf5.1,1p2e11.3,0pf5.1,1p2e11.3,
     >  25(1pe11.3,0pf5.1,1pe11.3,0pf5.1,1pe11.3))')
     >  jt, j, i, alnx(jt), zlmx-alnx(jt),
     >  vni(j,i), vti(j,i), vte(j,i),
     >  rmfpi, rlnti, rmfpe, rlnte,
     >  sgnti,  dabs(grdti), dabs(gdti(ic)),
     >  sgnte,  dabs(grdte), dabs(gdte(ic)),
     >  (Wtauz(k), dsign(1.0d0, Wftot(k)), dabs(Wftot(k)),
     >             dsign(1.0d0, Wffmx(k)), dabs(Wffmx(k)), k=1,kmax)
      enddo
!
      close(i6)
      return
      end
