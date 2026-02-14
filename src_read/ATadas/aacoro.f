!***********************************************************************
      program aacoro  
!***********************************************************************
!
!      corona-Lz(Te) using ADAS-data
!          include CXR
!
!      aacoro_1.f  aacoro_2.f  aacoro.f     2014/01/15
!
!      lold = 0
!      nds = ndatm => nds = ndmis  (csize)
!      atm_set(cdir,catmz) => set_atomz(drnm)
!      cmz => amz
!
!    setenv DRADAS "/home/shimizu/PRJ/sonicV4/adas/Ar_1220"
!
!    ip_ion, ip_rec, ip_cxr, ip_rad, ip_plt, ip_prb, ip_prc  &
!
!-----------------------------------------------------------------------
      use csize, only:ndmis,ndmc,ndmg,ndx,ndy,nvyd,ndgtp
      use catcom, only:catmz,natmz,aionz
     > ,ip_ion,ip_rec,ip_cxr,ip_plt,ip_prb
      use cplcom, only:aion
      use cimden, only:eipot
      use cphcns, only:cmp,cev
      use cunit, only:n6,lnope,lmspe,lmype
      use cimcom, only:amz,ane,ate
      use cimctl, only:nprg
      implicit none
!
!::local variables
      integer nds, ns
      parameter (nds=ndmis)
      real(8) ane0, ate0
      real(8), dimension(0:nds) :: alf, bet, gam, anz
      real(8), dimension(0:nds) :: pei, plt, prb
      real(8) :: lzei, lzln, lzrc, tlz1, tlz2, tlz3
      real(8) :: zero = 1.0d-60

      real(8) :: savz, szef, rlos
      character(80) :: cfmt1, cfmt2, cfmt3
!
      character(20) :: cinp
      character(80) :: cdir
      integer :: nf_ion, nf_rec, nf_plt, nf_wlz, nf_cxr
      integer :: nmax, i, j, iz
      real(8) :: atemin, atemax, dte
      real(8) :: tbne(201), tbte(201)
      real(8) :: acof(nds)
      real(8) :: fn0

      character(80) :: DRADAS
      character(*), parameter :: namelist_name = "corona_input"
      integer IMPMC_num

      IMPMC_num = 1
      DRADAS = "../../../adas/Ar_1220"

!::set grobal value
      namelist /corona/
     > DRADAS, IMPMC_num
     > ,ndmc,ndmg,ndx,ndy,nvyd,ndgtp
      open(10,file=namelist_name)
      read(10,corona)
      close(10)

      n6 = 10
      call git_mesg(n6) 
      write(n6,'(2x,"*** aacoro ***  2024/06/28")')

      call catcom_alc(1)
      lnope = 1
      lmspe = 0
      lmype = 0
      nprg = IMPMC_num
      call phcnst
      aion(1) = 2.0d0  !  D+ plasma
!
!::dir => catmz
      cdir = DRADAS
      if( len_trim(cdir) <= 0 ) then
        call wexit("aacoro","No DRADAS in namelist")
      endif
      write (n6,'(/2x,"DRADAS in namelist",a)') trim(cdir)
      call set_atomz_x(cdir,1)
!
!::choice atom  see sub. IMPMC/set_atomz.f
      write(n6,'(2x,"catmz = [",a,"]","  natmz =",i3,"  aionz =",f8.3
     >  ,"  amz =",1pe11.3)') catmz(1), natmz(1), aionz(1), amz

!::read data
      call set_adas

!::input data
      rlos = 0.0D0
      write (n6,'(2x)')
      write (*,'(2x,"rlos = 1/(Ne*TauI)")')
      write (*,'(2x,">> enter rlos (0:coro/1.0d-17/1.0d-15) ==> ",$)')
      read (*,'(a)') cinp
      read (cinp,*) rlos
      write (*,'(2x,">> enter Ne (1.0d20) ==> ",$)')
      read (*,'(a)') cinp
      read (cinp,*) ane
      write (*,'(2x,">> enter N0/Ne (1.0d-2) ==> ",$)')
      read (*,'(a)') cinp
      read (cinp,*) fn0
!
      write (*,'(2x,"catmz = ",a,"  natmz =",i3,"  aionz =",f8.3)')
     >   trim(catmz(1)), natmz(1), aionz(1)
      write (*,'(2x,"rlos =",1pe12.3
     > ,"  Ne =",1pe12.3,"  fN0 =", 1pe12.3)')
     >  rlos, ane, fn0
!
      write (n6,'(2x,"catmz = ",a,"  natmz =",i3,"  aionz =",f7.3)')
     > catmz(1), natmz(1), aionz(1)
      write (n6,'(2x,"rlos =",1pe12.3
     > ,"  Ne =",1pe12.3,"  fN0 =", 1pe12.3)')
     >  rlos, ane, fn0

      ane0 = 1.0D19
      ate0 = 1.0D0
      call atm_eval2(ip_ion,ane0,ate0,acof,ns,nds,IMPMC_num)
      write (n6,'(/2x,"ns =",i3,"  natmz =",i3," < nds =",i3)')
     >  ns, natmz(1), nds
      write (n6,'(2x,"aionz =",f7.3,"  amz =",1pe12.3,2x, "  amz/cmp =",
     >  1pe12.3)') aionz(1), amz, amz/cmp
      if ( ns/=natmz(1) .or. natmz(1) > nds ) then
        write (n6,'(2x,"invalid  ns.ne.natmz .or. aionz.ne.amz/cmp")')
        stop
      endif

!::check
      write (n6,'(/2x,"ip_ion =",i2,"  ip_rec =",i2
     > , "  ip_cxr(option) =",i2,
     >  "  ip_plt =",i2)') ip_ion(1), ip_rec(1), ip_cxr(1), ip_plt(1)
      if ( ip_ion(1)==0 .or. ip_rec(1)==0 .or. ip_plt(1)==0 ) then
        call wexit("aacoro_2","file allocation errorr")
      endif

      zero = 1.0D-70
      atemin = 1.0D-1
      atemax = 1.0D2
      nmax = 200      ! 41  < 201
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
      if ( ip_cxr(1)==0 ) nf_cxr = 0

      cfmt1 = '(4x,"i",3x,"Ne",10x,"Te",10x,30(:"sv",i2.2,8x))'
      write (cfmt1(30:31),'(i2)') nds + 1
      cfmt2 = '(4x,"i",3x,"Ne",10x,"Te",10x,30(:"Lz",i2.2,8x))'
      write (cfmt2(30:31),'(i2)') nds + 1

      open (unit=nf_ion,file="dtion.txt")
      write (nf_ion,'(2x,"=== ionization <sv> [m3/s] ===")')
      write (nf_ion,trim(cfmt1)) (j,j=0,ns)

      open (unit=nf_rec,file="dtrec.txt")
      write (nf_rec,'(2x,"=== recombination <sv> [m3/s] ===")')
      write (nf_rec,trim(cfmt1)) (j,j=0,ns)

      if ( nf_cxr>0 ) then
        open (unit=nf_cxr,file="dtcxr.txt")
        write (nf_cxr,'(2x,"=== CXR <sv> [m3/s] ===")')
        write (nf_cxr,trim(cfmt1)) (j,j=0,ns)
      endif

      open (unit=nf_plt,file="dtplt.txt")
      write (nf_plt,'(2x,"=== line rad.  Lz(Te) [W/m3] ===")')
      write (nf_plt,trim(cfmt2)) (j,j=0,ns)

      open (unit=nf_wlz,file="dtwlz.txt")
      write (nf_wlz,'(2x,"=== corona_3  WLz(Te) [W/m3] ===  rlos ="
     > , 1pe12.3,
     >  "  fn0 =",1pe12.3,"  Ne =",1pe12.3)') rlos, fn0, ane
      write (nf_wlz,'(4x,"i",3x,"Ne",10x,"Te"
     > ,10x,"<z>",9x,"<z2>",8x, "Lz_ei",7x,
     >  "Lz_ln",7x,"Lz_rc")')

!::loop
      write (6,'(2x)')
      do i = 1, nmax
        ane = tbne(i)
        ate = tbte(i)

!::ionization rate <sv> [m3/s]
        call atm_eval2(ip_ion,ane,ate,acof,ns,nds,IMPMC_num)
        do j = 1, ns
          iz = j - 1
          alf(iz) = acof(j)
        enddo
        alf(ns) = 0.0D0

        cfmt3 = '(2x,i3,1p2e12.3,1p30e12.3)'
        write (cfmt3(19:20),'(i2)') nds + 1

        write (nf_ion,trim(cfmt3)) i, ane, ate
     >   , (dmax1(alf(iz),zero),iz=0,ns)

!::recombination rate <sv> [m3/s]
        call atm_eval2(ip_rec,ane,ate,acof,ns,nds,IMPMC_num)
        do j = 1, ns
          iz = j
          bet(iz) = acof(j)
        enddo
        bet(0) = 0.0D0
        write (nf_rec,trim(cfmt3)) i, ane, ate
     >   , (dmax1(bet(iz),zero),iz=0,ns)

!::CXR rate <sv> [m3/s]   (ane,erl) ==> (ane,ate)
!       1/2*mH*Vrel^2 = 3/2*Ti + 1/2*8/pi*mH/mZ*ATz ~ 3/2*Ti

        if ( ip_cxr(1) > 0 ) then
          call atm_eval2(ip_cxr(1),ane,ate,acof,ns,nds,IMPMC_num)
        else
          do j = 1, ns
            acof(j) = 0.0D0
          enddo
        endif

        do j = 1, ns
          iz = j
          gam(iz) = acof(j)
        enddo
        gam(0) = 0.0D0
        if ( nf_cxr > 0 ) then
          write (nf_cxr,trim(cfmt3)) i, ane, ate
     >     , (dmax1(gam(iz),zero),iz=0,ns)
        endif

!::corona ratio of Nz
        call corona_3(rlos,fn0,ns,alf,bet,gam,anz)

!::ele. ionization
        call atm_eval2(ip_ion,ane,ate,acof,ns,nds,IMPMC_num)
        do j = 1, ns
          iz = j  - 1
          pei(iz) = acof(j)*eipot(iz)*cev
        enddo
        pei(ns) = 0.0d0

!::line radiation Lz
        call atm_eval2(ip_plt,ane,ate,acof,ns,nds,IMPMC_num)
        do j = 1, ns
          iz = j  - 1
          plt(iz) = acof(j)
        enddo
        plt(ns) = 0.0D0
        write (nf_plt,trim(cfmt3)) i, ane, ate
     >   , (dmax1(plt(iz),zero),iz=0,ns)

!::recombination & Bremsstrahlung Lz
      call atm_eval2(ip_prb,ane,ate,acof,ns,nds,IMPMC_num)
      do j = 1, ns  
        iz = j
        prb(iz) = acof(j)  
      enddo
      prb(0) = 0.0d0

!::<Z>, Zeff, Lz
      savz = 0.0D0
      szef = 0.0D0
      lzei = 0.0D0    ! ele. ionization
      lzln = 0.0d0    ! line radiation
      lzrc = 0.0d0    ! recombinaion 
      do iz = 0, ns
        savz = savz + dfloat(iz)*anz(iz)
        szef = szef + dfloat(iz)**2*anz(iz)
        lzei = lzei + pei(iz)*anz(iz)
        lzln = lzln + plt(iz)*anz(iz)
        lzrc = lzrc + prb(iz)*anz(iz)
      enddo

      if(ate.ge.55.0d0.and.ate.lt.120.0d0)then
      do iz=0,ns
      write(334,*) ate, iz, anz(iz)
      end do
      endif
      write(434,*) ate, pei(0)

      tlz1 = lzei
      tlz2 = tlz1 + lzln
      tlz3 = tlz2 + lzrc
      write (nf_wlz,'(2x,i3,1p8e12.3)') i, ane, ate, savz, szef,
     > dmax1(tlz1,zero), dmax1(tlz2,zero), dmax1(tlz3,zero) 

      enddo  ! loop  Te
!
      close(nf_ion)
      close(nf_rec)
      close(nf_cxr)
      close(nf_plt)
      close(nf_wlz)
      end program aacoro
