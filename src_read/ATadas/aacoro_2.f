!***********************************************************************
      program aacoro
!***********************************************************************
!
!      corona-Lz(Te) using ADAS-data
!          include CXR
!
!-----------------------------------------------------------------------
      implicit none
      include 'catcom'
      include 'csonic'
!
!::local common
      real*8  cpi, cev, cmp, cme, cmz
      common /com_phconst/ cpi, cev, cmp, cme, cmz
!
!::local variables
      integer nds, ns
      parameter (nds=ndatm)
      real(8) ane0, ate0, ane, ate
      real(8), dimension(0:nds) ::  alf, bet, gam, plt, anz, alz
      real(8) ::  savz, szef, stlz, rlos
      character(80) :: cfmt1, cfmt2, cfmt3
!
      character(2), dimension(9) :: tbcatmz
      integer, dimension(9)      :: tbnatmz
      real(8), dimension(9)      :: tbaionz
!
      character(20) :: cinp
      character(80) :: cdir
      integer :: nf_ion, nf_rec, nf_plt, nf_wlz, nf_cxr
      integer :: nmax, i, j, is, is1, lold, ityp
      real(8) :: atemin, atemax, dte, zero
      real(8) :: tbne(201), tbte(201)
      real(8) :: acof(nds)
      real(8) :: funLz_Ar
      real(8) :: fn0, ati, vrl2, vrl, erl
!
!::MPI
      n6 = 10
      lnope = 1
      lmspe = 0
      lmype = 0
!
!::dir => catmz
      write(n6,'(2x,"*** aacoro_2 ***  2012/01/18")')
      call getenv("drads", cdir)
      write(n6,'(2x,"setenv drads  cdir = ",a)') trim(cdir)
      call atm_set(cdir,catmz)
!
!::choice atom
      tbcatmz = (/"C ",  "Ne",    "Ar",    "Kr",   "Xe",   
     &            "He",  "Be",     "W",     "N"/)
      tbnatmz = (/   6,    10,      18,      36,     54, 
     &               2,     4,      74,       7/)
      tbaionz = (/12.011d0,  20.179d0,  39.95d0,  83.80d0,  131.3d0,
     >             4.003d0,   9.012d0,  183.8d0, 14.01d0/)
!
      write(n6,'(2x,"catmz = ",a,"  tbcatmz = ",5(a,2x))')
     >   trim(catmz), (tbcatmz(i),i=1,9)
      ityp = 0
      do i = 1, 9
      if( trim(catmz).eq.tbcatmz(i)) then
        ityp = i
        exit
      endif
      enddo
      if( ityp.eq.0 ) then
      write(n6,'(2x,"No found catmz in tbcatmz  ",a)') trim(catmz)
      stop
      endif
      write(*,'(2x,"ityp = ",i3)') ityp
!
      catmz = tbcatmz(ityp)
      natmz = tbnatmz(ityp)
      aionz = tbaionz(ityp)
!
!::input data
      rlos = 0.0d0
      write(n6,'(2x)')
      write(*,'(2x,"rlos = 1/(Ne*TauI)")')
      write(*,'(2x,">> enter rlos (0:coro/1.0d-17/1.0d-15) ==> ",$)')
      read(*,'(a)') cinp
      read(cinp,*) rlos
      write(*,'(2x,">> enter Ne (1.0d20) ==> ",$)')
      read(*,'(a)') cinp
      read(cinp,*) ane
      write(*,'(2x,">> enter N0/Ne (1.0d-2) ==> ",$)')
      read(*,'(a)') cinp
      read(cinp,*) fn0
!
      write(*,'(2x,"catmz = ",a,"  natmz =",i3,"  aionz =",f8.3)')
     > trim(catmz), natmz, aionz
      write(*,'(2x,"rlos =",1pe12.3,"  Ne =",1pe12.3,"  fN0 =",
     >  1pe12.3)') rlos, ane, fn0
!
      lold = 0
      write(n6,'(2x,"catmz = ",a,"  natmz =",i3,"  aionz =",f7.3,
     >  "  lold =",i2 )') catmz, natmz, aionz, lold
      write(n6,'(2x,"rlos =",1pe12.3,"  Ne =",1pe12.3,"  fN0 =",
     >  1pe12.3)') rlos, ane, fn0
!
!::const
      cpi = 4.0d0*datan(1.0d0)
      cev = 1.60210e-19
      cmp = 1.67252e-27
      cme = 9.10908e-31
      cmz = cmp*aionz   ! old-version  12.0d0 => aionz
!
!::read data
      if( lold.eq.0 ) call set_adas      
      if( lold.eq.1 ) call set_adas_old
!
      ane0 = 1.0d19
      ate0 = 1.0d0
      if( lold.eq.0 ) call atm_eval2(ip_ion,ane0,ate0,acof,ns)
      if( lold.eq.1 ) call atm_eval2_old(ip_ion,ane0,ate0,acof,ns)
      write(n6,'(/2x,"ns =",i3,"  natmz =",i3," < ndatm =",i3)')
     >   ns, natmz, ndatm
      write(n6,'(2x,"aionz =",f7.3,"  cmz =",1pe12.3,2x,
     >   "  cmz/cmp =",1pe12.3)')  aionz, cmz, cmz/cmp
      if( ns.ne.natmz .or. natmz.gt.ndatm .or.
     >   dabs(aionz-cmz/cmp).gt.0.001d0 ) then
        write(n6,'(2x,"invalid  ns.ne.natmz .or. aionz.ne.cmz/cmp")')
        stop
      endif
!
!::check
!xx   call chk_adas
      write(n6,'(/2x,"ip_ion =",i2,"  ip_rec =",i2,
     >  "  ip_cxr(option) =",i2,"  ip_plt =",i2)')
     >  ip_ion, ip_rec, ip_cxr, ip_plt
      if( ip_ion.eq.0 .or. ip_rec.eq.0 .or. ip_plt.eq.0 ) then
      call wexit("aacoro_2","file allocation errorr")
      endif
!
      zero = 1.0d-70
      atemin = 1.0d0
      atemax = 1.0d4
      nmax = 101  ! 41  < 201
      dte = dlog(atemax/atemin)/dfloat(nmax-1)
      do i = 1, nmax
      tbne(i) = ane
      tbte(i) = atemin*dexp(dte*dfloat(i-1))
      enddo
!
!::file
      nf_ion = 61
      nf_rec = 62
      nf_cxr = 63
      nf_plt = 64
      nf_wlz = 65
      if( ip_cxr.eq.0 ) nf_cxr = 0
!
      cfmt1 = '(4x,"i",3x,"Ne",10x,"Te",10x,30(:"sv",i2.2,8x))'
      write(cfmt1(30:31),'(i2)') ndatm+1
      cfmt2 = '(4x,"i",3x,"Ne",10x,"Te",10x,30(:"Lz",i2.2,8x))'
      write(cfmt2(30:31),'(i2)') ndatm+1
!
      open(unit=nf_ion,file="dtion.txt")
      write(nf_ion,'(2x,"=== ionization <sv> [m3/s] ===")')
      write(nf_ion,trim(cfmt1)) (is,is=0,ns)
!
      open(unit=nf_rec,file="dtrec.txt")
      write(nf_rec,'(2x,"=== recombination <sv> [m3/s] ===")')
      write(nf_rec,trim(cfmt1)) (is,is=0,ns)
!
      if( nf_cxr.gt.0 ) then
      open(unit=nf_cxr,file="dtcxr.txt")
      write(nf_cxr,'(2x,"=== CXR <sv> [m3/s] ===")')
      write(nf_cxr,trim(cfmt1)) (is,is=0,ns)
      endif
!
      open(unit=nf_plt,file="dtplt.txt")
      write(nf_plt,'(2x,"=== line rad.  Lz(Te) [W/m3] ===")')
      write(nf_plt,trim(cfmt2)) (is,is=0,ns)
!
      open(unit=nf_wlz,file="dtwlz.txt")
      write(nf_wlz,'(2x,"=== corona_3  WLz(Te) [W/m3] ===  rlos =",
     >  1pe12.3,"  fn0 =",1pe12.3,"  Ne =",1pe12.3)') rlos, fn0, ane
      write(nf_wlz,'(4x,"i",3x,"Ne",10x,"Te",10x,"<z>",9x,"<z2>",8x,
     > "tLz",9x,"fit_Lz")')
!
!::loop
      write(6,'(2x)')
      do i = 1, nmax
      ane = tbne(i)
      ate = tbte(i)
      ati = ate
      vrl2= 2.0d0*ati*cev/cmz
      vrl = sqrt(vrl2)
      erl = 0.5d0*cmp*vrl2/cev
!
!::ionization rate <sv> [m3/s]
      if( lold.eq.0 ) call atm_eval2(ip_ion,ane,ate,acof,ns)
      if( lold.eq.1 ) call atm_eval2_old(ip_ion,ane,ate,acof,ns)
      do is = 1, ns
      is1 = is-1
      alf(is1) = acof(is)
      enddo
      alf(ns) = 0.0d0

      cfmt3 = '(2x,i3,1p2e12.3,1p30e12.3)'
      write(cfmt3(19:20),'(i2)') ndatm+1

      write(nf_ion,trim(cfmt3))
     >  i, ane, ate, (dmax1(alf(is),zero),is=0,ns)
!
!::recombination rate <sv> [m3/s]
      if( lold.eq.0 ) call atm_eval2(ip_rec,ane,ate,acof,ns)
      if( lold.eq.1 ) call atm_eval2_old(ip_rec,ane,ate,acof,ns)
      do is = 1, ns
      bet(is) = acof(is)
      enddo
      bet(0) = 0.0d0
      write(nf_rec,trim(cfmt3)) 
     >  i, ane, ate, (dmax1(bet(is),zero),is=0,ns)
!
!::CXR rate <sv> [m3/s]   (ane,erl) ==> (ane,ati)
!       1/2*mH*Vrel^2 = 3/2*Ti + 1/2*8/pi*mH/mZ*ATz ~ 3/2*Ti
!
      if( ip_cxr.gt.0 ) then
      if( lold.eq.0 ) call atm_eval2(ip_cxr,ane,ati,acof,ns)
      if( lold.eq.1 ) call atm_eval2_old(ip_cxr,ane,ati,acof,ns)
      else
      do is = 1, ns
      acof(is) = 0.0d0
      enddo
      endif
!
      do is = 1, ns
      gam(is) = acof(is)
      enddo
      gam(0) = 0.0d0
      if( nf_cxr.gt.0 ) then
      write(nf_cxr,trim(cfmt3))
     >  i, ane, ate, (dmax1(gam(is),zero),is=0,ns)
      endif
!
!::line radiation Lz      
      if( lold.eq.0 ) call atm_eval2(ip_plt,ane,ate,acof,ns)
      if( lold.eq.1 ) call atm_eval2_old(ip_plt,ane,ate,acof,ns)
      do is = 1, ns
      is1 = is - 1
      plt(is1) = acof(is)
      enddo
      plt(ns) = 0.0d0
      write(nf_plt,trim(cfmt3))
     >  i, ane, ate, (dmax1(plt(is),zero),is=0,ns)
!
!::corona radiation
      call corona_3(rlos,fn0,ns,alf,bet,gam,plt,anz,alz,savz,szef,stlz)
!
!::debug write
      write(nf_wlz,'(2x,i3,1p5e12.3)') 
     >   i, ane, ate, savz, szef, dmax1(stlz,zero)
      enddo
!
      stop
      end
!
!***********************************************************************
      function funLz_Ar(ate)
!***********************************************************************
!
!      clooling rate of Ar
!
!           Ne = 1.0d20 1/m3
!           1 eV  <  Te < 1.0 keV
!           rlos = 1/(Ne*TauI) = 1.0d-16
!                                                2007/06/13
!
!-----------------------------------------------------------------------
      implicit none
!
      real*8  ate, funLz_Ar
      real*8  zte, zx, zy
!
      zte = ate
      if( zte.gt.1.0d3 ) zte = 1.0d3
!
      zx = dlog10(ate)
      zy = -36.664d0 + 11.499d0*zx - 9.9541d0*zx**2
     >               + 7.2144d0*zx**3 - 4.2755d0*zx**4
     >               + 1.3616d0*zx**5 - 0.16327d0*zx**6
!
      funLz_Ar = 10.0d0**zy
      return
      end
