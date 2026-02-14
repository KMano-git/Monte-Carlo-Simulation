!***********************************************************************
      subroutine atm_copy(cxsp,ixsp)
!***********************************************************************
!
!         cxsp  (i)   "ATOM_ION"   ionization of carbon
!                     "ATOM_REC"   recombination of carbon
!                     "ATOM_RAD"   radizaion of carbon   Lz
!                     "ATOM_CXR"   charge exchange recombination
!         dsn   (o)   name of data file
!
!     iclass  dsn           meaning
!  ----------------------------------------------------
!      (1)  acd96_c.data : recombination rate
!      (2)  scd96_c.data : ionization rate
!      (3)  ccd96_c.data : CX-transfer
!      (4)  prb96_c.data : Power recombination
!      (5)  prc96_c.data : Power continuum
!      (6)  qcd96_c.data : ??
!      (7)  xcd96_c.data : ??
!      (8)  plt96_c.data : Power ionization
!      (9)  pls96_c.data : Power charge transfer
!  ----------------------------------------------------
!
! Header of data
!    6   24   30    1    6     /CARBON             /GCR PROJECT
!
! Header of data  acd96r_he.dat
!    2   24   30    1    2     /HELIUM             /GCR PROJECT
!    2    1    1          meta-stable
!-----------------------------------------------------------------------
      use catcom, only : ndxdt, ndxkn, ndxx1, ndxx2, nxs, nyp
     >    , xdaty, xdatz, xdevl, xdnam, xdrnk, xdsn, xdtx1, xdtx2, xjend
     >    , xjnum, xjsta, xwmax, xwmin, xwmlt, xwnum, xwspc, xwunt
     >    , xwvar
      use cunit,  only : n6
      use csize, only : ndmis
      implicit none
!
!::argument
      character, intent(in)  :: cxsp*(*)
      integer,   intent(out) :: ixsp
!
!::local variables for atm_read
      integer isdimd, izdimd, itdimd, iddimd
      parameter ( isdimd = ndmis, izdimd = isdimd )
      parameter ( itdimd = ndxx2, iddimd = ndxx1 )
!
      character dsninc*120
      integer   ifail, iz0
      integer   ismaxd,izmaxd,itmaxd,idmaxd
      logical   lpart
!
!::local variables
      integer n, i, it, nw, ynum, ynum2, ismx
      integer ndty, m1, m2
      character clin*120, cttl(10)*80, cmsg*80, dsn*80, cout*12
      character cknd*3, catm*2
      integer   nft, nknd, nymx
      integer   i6, jd, mj, id, isel
      real*8    zero, zmin, zmax
      integer nfl; data nfl/0/; save nfl
!
      integer mjs(20),mje(20),nk
      integer idszx, itszx, izmnx, izmxx, nszx
      character catmx*10, cspc*120
! function
      integer    indxr, lenx
!
!::see sub. atm_read
!
      data cttl(1)/"recombination rate"/
      data cttl(2)/"ionization rate"/
      data cttl(3)/"CX-transfer"/
      data cttl(4)/"Power-recombination"/     ! prb  2004/11/12
      data cttl(5)/"Power-continuum"/         ! prc  2004/11/12
      data cttl(6)/"qcd"/
      data cttl(7)/"xcd"/
      data cttl(8)/"Power-ionizaion"/         ! plt
      data cttl(9)/"Power-charage transfer"/  ! pls  2004/11/12
!
!-----------------------------------------------------------------------
!::initial
!-----------------------------------------------------------------------
!
!::clear
      if( index(cxsp,"clear").gt.0 ) then
      write(n6,'(/2x,"*** atm_copy ***   Adas-data are cleared. ",
     >  "  nfl =",i3)') nfl
      nfl = 0
      return
      endif
!
      call atm_fnam(cxsp,dsn)
!
!::no found file
      if( len_trim(dsn).eq.0 ) then
        ixsp = 0
        return
      endif
!
!::clear counter variables & data
      if( nfl.ne.0 .and. nxs(1).eq.0 ) then
      write(n6,'(/2x,"*** atm_copy ***  Warning !")')
      write(n6,'(4x,"Adas-data are cleared. nfl =",i3)') nfl
      nfl = 0
      endif
!
      nfl = nfl + 1
      if( nfl.eq.1 ) then
        nxs(1) = 0;  nyp(1) = 0
        xdaty(1:ndxdt,1) = 0.0_8
        xdtx1(1:ndxx1,1:ndxkn,1) = 0.0_8
        xdtx2(1:ndxx2,1:ndxkn,1) = 0.0_8
        xdatz(1:ndmis,1:ndxkn,1) = 0.0_8
      endif
!
!-----------------------------------------------------------------------
!::data type and atom
!-----------------------------------------------------------------------
      nknd = 0
      cknd = cxsp(6:8)
      m1 = index( dsn, "_", .true. )
      m2 = index( dsn, ".", .true. )
      catm = dsn(m1+1:m2-1)
      if( cknd.eq."REC" ) nknd = 1
      if( cknd.eq."ION" ) nknd = 2
      if( cknd.eq."CXR" ) nknd = 3
      if( cknd.eq."RAD" ) nknd = 8
!-----
      if( cknd.eq."PRB" ) nknd = 4
      if( cknd.eq."PRC" ) nknd = 5
      if( cknd.eq."PLT" ) nknd = 8
      if( cknd.eq."PLS" ) nknd = 9
!-----
      write(n6,'(/2x,"*** atm_copy ***",
     >  3x,a,"  knd =",i3,2x,a,"   atm =",a)')
     >  cxsp(1:lenx(cxsp)),nknd,cknd,trim(catm)
      if( nknd.le.0 ) then
        call wexit("atm_copy","incorrect nknd")
      endif
      write(n6,'(2x,"dsn =",a)') dsn(1:lenx(dsn))
!
!::file open
      nft = 21
      open( unit=nft, file=dsn, form="formatted", action="read" )
!
!::header & ismaxd
      ismx = 0
      do i = 1, 200000
      read(nft,'(a)',end=120) clin
      if( i.eq.1 ) cspc = clin
      if( i.le.17 ) write(n6,'(2x,a)') clin
      if( index(clin,"IGRD").ne.0 .or. index(clin,"IPRT").ne.0 )
     >  ismx = ismx + 1
      enddo
 120  continue
      rewind nft
!
!-----------------------------------------------------------------------
!::data size
!-----------------------------------------------------------------------
      call linsep(cspc," /",nk,mjs,mje,20)
      read( cspc(mjs(1):mje(1)), * ) iz0
      read( cspc(mjs(2):mje(2)), * ) idszx
      read( cspc(mjs(3):mje(3)), * ) itszx
      read( cspc(mjs(4):mje(4)), * ) izmnx
      read( cspc(mjs(5):mje(5)), * ) izmxx
      catmx = cspc(mjs(6):mje(6))
      nszx  = idszx*itszx*ismx
!
      write(n6,'()')
      write(n6,'(2x,"atom =",a,"  ismx =",i3)')
     >  catmx(1:lenx(catmx)),ismx
      write(n6,'(2x,"iz0 =",i3,"   idsz =",i3,"   itsz =",i3,
     >  "   izmn =",i3,"   izmx =",i3,"  nszx =",i7)')
     >  iz0, idszx, itszx, izmnx, izmxx, nszx
!
!::pointer
      nymx   = nyp(1) + nszx
      nxs(1) = nxs(1) + 1
      nyp(1) = nyp(1) + 1
!
!::dimension error
      if( nxs(1).gt.ndxkn ) then
      write(cmsg,'("too many variables  nxs > ndxkn ",2i5)')
     >    nxs(1),ndxkn
      call wexit("atm_copy",cmsg)
      endif
      if( nymx.gt.ndxdt ) then
      write(cmsg,'("too many data  nyp+nszx > ndxdt ",2i8)') nymx,ndxdt
      call wexit("atm_copy",cmsg)
      endif
!
      if( iz0 .gt.isdimd ) then
        write(n6,'(2x,"iz0.gt.isdimd  ",2i5)') iz0, isdimd
        call wexit("atm_copy","iz0.gt.isdimd")
      endif
      if( itszx.gt.itdimd ) then
        write(n6,'(2x,"itszx.gt.itdimd ",2i5)') itszx, itdimd
        call wexit("atm_copy","itsz.gt.itdimd")
      endif
      if( idszx.gt.itdimd ) call wexit("atm_copy","idsz.gt.itdimd")
      if( ismx.gt.izdimd ) call wexit("atm_copy","izsz.gt.izdimd")
      if( idszx.gt.ndxx1  ) then
        write(n6,'(2x,"idszx.gt.ndxx1 ",2i5)') idszx, ndxx1
        call wexit("atm_copy","idsz.gt.ndxx1")
      endif
      if( itszx.gt.ndxx2  ) call wexit("atm_copy","itsz.gt.ndxx2")
!
!-----------------------------------------------------------------------
!::adas-input data
!-----------------------------------------------------------------------
!-INPUT
      dsninc = dsn     ! master file name
      lpart  = .false. ! .T. partial(resolved) .F. unresolved master dt.
      mj = indxr(dsn,"_")
      if( dsn(mj-1:mj-1).eq."r" ) lpart = .true.
!--
!     isdimd = maximum number of (charge, parent, gound)
!     izdimd = maximum number of charge states
!     itdimd = maximum number of temp or dens values
!
!--OUTPUT
!     ifail  = 0 (sucess); 1 (open error)
!     ismaxd = number of (charge, parent, gound)
!     izmaxd = number of charge state (zdata)
!     itmaxd = number of Te (dtevd)
!     idmaxd = number of Ne (ddensd)
!     dtevd  = dlog10(Te)    [eV]
!     ddensd = dlog10(Ne)    [1/cm^3]
!     drcofd = dlog10(<sig*v>(Is,Te,Ne))   [cm^3/s]
!     zdata  = charge + 1 for ions
!
!::read master-file
      ndty = ndxdt - (nyp(1) - 1)
      if( nszx.gt.ndty ) then
      write(cmsg,'("dimension error  nszx,ndty,nyp,ndxdt ",4i7)')
     >   nszx,ndty,nyp(1),ndxdt
      call wexit("atm_copy",cmsg)
      endif
!
      call atm_read( dsninc, nft, iz0, nknd, lpart,  ifail,
     >  isdimd, izdimd, itdimd,  iddimd,
     >  ismaxd, izmaxd, itmaxd,  idmaxd,
     >  xdatz(1,nxs(1),1), xdtx2(1,nxs(1),1), xdtx1(1,nxs(1),1)
     >    , xdaty(nyp(1),1), ndty )
      close(nft)
!
      write(n6,'(2x,"isdimd, izdimd, itdimd, iddimd ",4i8)')
     >  isdimd, izdimd, itdimd, iddimd
      write(n6,'(2x,"ismaxd, izmaxd, itmaxd, idmaxd ",4i8)')
     >  ismaxd, izmaxd, itmaxd, idmaxd
!
      if( ifail.ne.0 ) then
      call wexit("atm_copy","no found file")
      endif
!
!::assumption
      if( iz0.ne.ismaxd ) then
      write(n6,'(2x,"Warning !  iz0 =",i3," .ne. ismaxd =",i3)')
     >   iz0, ismaxd
      endif
!
!-----------------------------------------------------------------------
!:save
!-----------------------------------------------------------------------
      ixsp = nxs(1)
      n    = nxs(1)
      xdsn(n,1)  = dsn
      xdnam(n,1) = cxsp
      xdevl(n,1) = "table_calc"
      xdrnk(n,1) = 3
!-----
      xwvar(1,n,1) = "<sig*v>"
      xwvar(2,n,1) = "Ne"; xwvar(3,n,1) = "Te"; xwvar(4,n,1) = "Iz"
      xwnum(1,n,1) = nszx
      xwnum(2,n,1) = idmaxd; xwnum(3,n,1) = itmaxd
      xwnum(4,n,1) = ismaxd
      xwmin(1,n,1) = 0.0d0;
      xwmin(2,n,1) = xdtx1(1,n,1); xwmin(3,n,1) = xdtx2(1,n,1)
      xwmin(4,n,1) = xdatz(1,n,1)
      xwmax(1,n,1) = 0.0d0;
      xwmax(2,n,1) = xdtx1(idmaxd,n,1)
      xwmax(3,n,1) = xdtx2(itmaxd,n,1)
      xwmax(4,n,1) = xdatz(ismaxd,n,1)
      xwmlt(1,n,1) = 1.0d-6
      xwmlt(2,n,1) = 1.0d0; xwmlt(3,n,1) = 1.0d0;  xwmlt(4,n,1) = 1.0d0
      xwunt(1,n,1) = "cm^3/s"
      xwunt(2,n,1) = "1/cm^3"; xwunt(3,n,1) = "eV"; xwunt(4,n,1) = "[]"
      xwspc(1,n,1) = "log10"
      xwspc(2,n,1) = "log10"; xwspc(3,n,1) = "log10"
      xwspc(4,n,1) = "--"
      nw = xdrnk(n,1) + 1
      ynum = xwnum(1,n,1)
      ynum2 = 1
      do i = 2, nw
      ynum2 = ynum2*xwnum(i,n,1)
      enddo
      xjsta(n,1) = nyp(1)
      xjend(n,1) = nyp(1) + ynum - 1
      xjnum(n,1) = ynum
!
      nyp(1) = xjend(n,1)
!
      return
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      mj = indxr(dsn,"/")
      write(n6,'(/2x,"*** xdt_copy ***  [",a,"]  nxs =",i3,"  cvar =",a
     > ,2x,"devl =",a12,"  drnk =",i2)')  trim( xdnam(n,1) )
     > ,nxs(1),dsn(mj+1:lenx(dsn)),xdevl(n,1),xdrnk(n,1)
      do i = 1, nw
      zmin = xwmin(i,n,1); zmax = xwmax(i,n,1)
      write(n6,'(4x,"wvar =",a,"  wnum =",i6,"  wmin =",1pe10.3,
     >  "  wmax =",1pe10.3,"  wmlt =",1pe10.3,"  wunt =",a8,
     >  "  wspc =",a8)') xwvar(i,n,1),xwnum(i,n,1),zmin,zmax
     > ,xwmlt(i,n,1),xwunt(i,n,1),xwspc(i,n,1)
      enddo
      write(n6,'(4x,"ynum =",2i8,
     >     "  daty =",i8,1x,1pe10.3,"  daty =",i8,1x,1pe10.3
     >  ,"  dmsz =",2i8)')
     >   ynum,ynum2,xjsta(n,1),xdaty(xjsta(n,1),1)
     >   ,xjend(n,1),xdaty(xjend(n,1),1), nyp(1), ndxdt
!
      write(n6,'(2x," Ne =",1p10e12.3)')
     >    (10.0d0**xdtx1(i,n,1),i=1,idmaxd)
      write(n6,'(2x," Te =",1p10e12.3)')
     >    (10.0d0**xdtx2(i,n,1),i=1,itmaxd)
      if( nyp(1).gt.ndxdt ) then
      call wexit("atm_copy","too many data  nyp > ndxdt")
      endif
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      i6 = 81
      jd = 17     ! <== ne = 1.0e13
!
      cout = "data_"//cknd//"_"//catm
      open( unit=i6, file=cout(1:lenx(cout)) )
!
      zero = 1.0d-70
      write(i6,'(/2x,"*** <atom> "a," ***")')
     >   cttl(nknd)(1:lenx(cttl(nknd)))
      write(i6,'(2x,"natm =",i3,"  rtyp =",a)')   iz0,cknd
      write(i6,'(2x,"unit =",i3,"  dsn =",a)') nft,dsn(1:80)
      write(i6,'(2x,"iz0  =",i3,"  iclass =",i3,"  lpart =",l4,
     >  "  ifail =",i3)') iz0,nknd,lpart,ifail
      write(i6,'(2x,"ismaxd =",i3,"  izmaxd =",i3,"  itmaxd =",i3,
     >  "  idmaxd =",i3)') ismaxd,izmaxd,itmaxd,idmaxd
      write(i6,'(2x,"dne  =",1p10e11.2)')
     >  (10.0**xdtx1(id,n,1),id=1,idmaxd)
      write(i6,'(2x,"dte  =",1p10e11.2)')
     >  (10.0**xdtx2(it,n,1),it=1,itmaxd)
      write(i6,'(2x,"dne =",i5,2x,1pe11.2)') jd,10.0**xdtx1(jd,n,1)
      write(i6,'(7x,"Te",5x,10(3x,i2,6x))')
     >    (int(xdatz(isel,n,1)+0.001),isel=1,ismaxd)
      do it = 1, itmaxd
      write(i6,'(2x,1p12e11.2)') 10.0**xdtx2(it,n,1),(dmax1(zero,
     >  10.0_8**xdaty(xjsta(n,1)-1
     >                +jd+(it-1)*idmaxd+(isel-1)*idmaxd*itmaxd,1))
     > ,isel=1,ismaxd)
      enddo
!
      return
      end
