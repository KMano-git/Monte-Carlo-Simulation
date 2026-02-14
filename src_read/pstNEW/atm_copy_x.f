!***********************************************************************
      subroutine atm_copy_x(cxsp,ixsp,kis)
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
! modified 1/2 lines organize local variables and include files by kamata 2021/06/28
!ik   character cxsp*(*) ; integer ixsp
      character, intent(in)  :: cxsp*(*)
      integer,   intent(out) :: ixsp
      integer,   intent(in) :: kis
!
!::local variables for atm_read
      integer isdimd, izdimd, itdimd, iddimd
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   parameter ( isdimd = ndxzi, izdimd = isdimd )
      parameter ( isdimd = ndmis, izdimd = isdimd )
      parameter ( itdimd = ndxx2, iddimd = ndxx1 )
!
      character dsninc*120
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer   ifail, iz0, iclass
      integer   ifail, iz0
      integer   ismaxd,izmaxd,itmaxd,idmaxd
      logical   lpart
!
!xx   real*8    dtevd(itdimd), ddensd(itdimd), zdata(isdimd)
!xx   real*8    ddensd(itdimd) ==> xdtx1(ndxx1,ndxknd)
!xx   real*8    dtevd(itdimd)  ==> xdtx2(ndxx2,ndxknd)
!xx   real*8    zdata(isdimd)  ==> xdatz(ndxzi,ndxknd)
!xx   real*8    drcofd(isdimd,itdimd,itdimd) ==> xdaty(ndxdt)
!
!::local variables
! modified 2/2 lines organize local variables and include files by kamata 2021/06/28
!ik   integer n, jst, i, j, is, in, it, nw, ynum, ynum2, ismx, j0
!ik   integer ndrw, ndty, nsiz, m1, m2
      integer n, i, it, nw, ynum, ynum2, ismx
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   integer ndty, nsiz, m1, m2
      integer ndty, m1, m2
      character clin*120, cttl(10)*80, cmsg*80, dsn*80, cout*12
      character cknd*3, catm*2
! modified 3/3 lines organize local variables and include files by kamata 2021/06/28
!ik   integer   lenx, nft, natm, nknd, nymx, indxr
!ik   integer   iunt12, i6, jd, mj, id, isel
!ik   real*8    eps, zero, zmin, zmax, zsta, zend
      integer   nft, nknd, nymx
      integer   i6, jd, mj, id, isel
      real*8    zero, zmin, zmax
      integer nfl; data nfl/0/; save nfl
!
      integer mjs(20),mje(20),nk
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer iz0x, idszx, itszx, izmnx, izmxx, izszx, nszx
      integer idszx, itszx, izmnx, izmxx, nszx
      character catmx*10, cspc*120
! added 2 lines organize local variables and include files by kamata 2021/06/28
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
      call atm_fnam_x(cxsp,dsn,kis)
!
!::no found file
      if( len_trim(dsn).eq.0 ) then
        ixsp = 0
        return
      endif
!
!::clear counter variables & data
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   if( nfl.ne.0 .and. nxs.eq.0 ) then
      if( nfl.ne.0 .and. nxs(kis).eq.0 ) then
      write(n6,'(/2x,"*** atm_copy ***  Warning !")')
      write(n6,'(4x,"Adas-data are cleared. nfl =",i3)') nfl
      nfl = 0
      endif
!
      nfl = nfl + 1
      if( nfl.eq.1 ) then
! modified 5/5 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik     nxs = 0;  nyp = 0
!ik     nsiz = ndxdt;        call setd( xdaty, nsiz, 0.0d0 )
!ik     nsiz = ndxx1*ndxkn;  call setd( xdtx1, nsiz, 0.0d0 )
!ik     nsiz = ndxx2*ndxkn;  call setd( xdtx2, nsiz, 0.0d0 )
!ik     nsiz = ndxzi*ndxkn;  call setd( xdatz, nsiz, 0.0d0 )
        nxs(kis) = 0
        nyp(kis) = 0
        xdaty(1:ndxdt,kis) = 0.0_8
        xdtx1(1:ndxx1,1:ndxkn,kis) = 0.0_8
        xdtx2(1:ndxx2,1:ndxkn,kis) = 0.0_8
        xdatz(1:ndmis,1:ndxkn,kis) = 0.0_8
      endif
!
!-----------------------------------------------------------------------
!::data type and atom
!-----------------------------------------------------------------------
      nknd = 0
      cknd = cxsp(6:8)
!xx   catm = cxsp(10:)
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
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   read( cspc(mjs(1):mje(1)), * ) iz0x
      read( cspc(mjs(1):mje(1)), * ) iz0
      read( cspc(mjs(2):mje(2)), * ) idszx
      read( cspc(mjs(3):mje(3)), * ) itszx
      read( cspc(mjs(4):mje(4)), * ) izmnx
      read( cspc(mjs(5):mje(5)), * ) izmxx
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   natm  = iz0x
      catmx = cspc(mjs(6):mje(6))
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   izszx = ismx
      nszx  = idszx*itszx*ismx
!
      write(n6,'()')
      write(n6,'(2x,"atom =",a,"  ismx =",i3)')
     >  catmx(1:lenx(catmx)),ismx
      write(n6,'(2x,"iz0 =",i3,"   idsz =",i3,"   itsz =",i3,
     >  "   izmn =",i3,"   izmx =",i3,"  nszx =",i7)')
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik  >  iz0x, idszx, itszx, izmnx, izmxx, nszx
     >  iz0, idszx, itszx, izmnx, izmxx, nszx
!
!::pointer
! modified 3/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   nymx = nyp + nszx
!ik   nxs  = nxs + 1
!ik   nyp  = nyp + 1
      nymx   = nyp(kis) + nszx
      nxs(kis) = nxs(kis) + 1
      nyp(kis) = nyp(kis) + 1
!
!::dimension error
! modified 2/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   if( nxs.gt.ndxkn ) then
!ik   write(cmsg,'("too many variables  nxs > ndxkn ",2i5)') nxs,ndxkn
      if( nxs(kis).gt.ndxkn ) then
      write(cmsg,'("too many variables  nxs > ndxkn ",2i5)')
     >    nxs(kis),ndxkn
      call wexit("atm_copy",cmsg)
      endif
      if( nymx.gt.ndxdt ) then
      write(cmsg,'("too many data  nyp+nszx > ndxdt ",2i8)') nymx,ndxdt
      call wexit("atm_copy",cmsg)
      endif
!
! modified 2/2 lines organize local variables and include files by kamata 2021/06/28
!ik   if( iz0x .gt.isdimd ) then
!ik     write(n6,'(2x,"iz0x.gt.isdimd  ",2i5)') iz0x, isdimd
      if( iz0 .gt.isdimd ) then
        write(n6,'(2x,"iz0.gt.isdimd  ",2i5)') iz0, isdimd
        call wexit("atm_copy","iz0.gt.isdimd")
      endif
      if( itszx.gt.itdimd ) then
        write(n6,'(2x,"itszx.gt.itdimd ",2i5)') itszx, itdimd
        call wexit("atm_copy","itsz.gt.itdimd")
      endif
      if( idszx.gt.itdimd ) call wexit("atm_copy","idsz.gt.itdimd")
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   if( izszx.gt.izdimd ) call wexit("atm_copy","izsz.gt.izdimd")
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
! deleted 2 lines organize local variables and include files by kamata 2021/06/28
!ik   iunt12 = nft     ! unit number
!ik   iz0    = natm    ! nuclear charge
      lpart  = .false. ! .T. partial(resolved) .F. unresolved master dt.
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   iclass = nknd    ! class of data (1-9)
!
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
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   jst  = nyp
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   ndty = ndxdt - (nyp - 1)
      ndty = ndxdt - (nyp(kis) - 1)
      if( nszx.gt.ndty ) then
      write(cmsg,'("dimension error  nszx,ndty,nyp,ndxdt ",4i7)')
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik  >   nszx,ndty,nyp,ndxdt
     >   nszx,ndty,nyp(kis),ndxdt
      call wexit("atm_copy",cmsg)
      endif
!
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   call atm_read( dsninc, iunt12, iz0, iclass, lpart,  ifail,
      call atm_read( dsninc, nft, iz0, nknd, lpart,  ifail,
     >  isdimd, izdimd, itdimd,  iddimd,
     >  ismaxd, izmaxd, itmaxd,  idmaxd,
! modified 2/2 lines organize local variables and include files by kamata 2021/06/28
!ik  >  xdatz(1,nxs), xdtx2(1,nxs), xdtx1(1,nxs), xdaty(jst), ndty )
!ik   close(iunt12)
! modified 1/2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik  >  xdatz(1,nxs), xdtx2(1,nxs), xdtx1(1,nxs), xdaty(nyp), ndty )
     >  xdatz(1,nxs(kis),kis), xdtx2(1,nxs(kis),kis), 
     >  xdtx1(1,nxs(kis),kis), xdaty(nyp(kis),kis), ndty )
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
!xx   write(clin,'("iz0 =",i3," .ne. ismaxd =",i3)') iz0,ismaxd
!xx   call wexit("atm_copy",clin)
      endif
!
!-----------------------------------------------------------------------
!:save
!-----------------------------------------------------------------------
! modified 2/2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   ixsp = nxs
!ik   n    = nxs
      ixsp = nxs(kis)
      n    = nxs(kis)
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   ndrw = ndxrw
! modified 4/4 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   xdsn(n)  = dsn
!ik   xdnam(n) = cxsp
!ik   xdevl(n) = "table_calc"
!ik   xdrnk(n) = 3
      xdsn(n,kis)  = dsn
      xdnam(n,kis) = cxsp
      xdevl(n,kis) = "table_calc"
      xdrnk(n,kis) = 3
!-----
!xx      write(n6,'(2x,8f10.5)') (xdatz(i,n),i=1,ismaxd)
!xx      write(n6,'(2x,8f10.5)') (xdaty(i),i=jst,jst+55)
!-----
! modified 18/21 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   xwvar(1,n) = "<sig*v>"
!ik   xwvar(2,n) = "Ne";    xwvar(3,n) = "Te";   xwvar(4,n) = "Iz"
!ik   xwnum(1,n) = nszx
!ik   xwnum(2,n) = idmaxd;  xwnum(3,n) = itmaxd; xwnum(4,n) = ismaxd
!ik   xwmin(1,n) = 0.0d0;
!ik   xwmin(2,n) = xdtx1(1,n); xwmin(3,n) = xdtx2(1,n);
!ik  >                         xwmin(4,n) = xdatz(1,n)
!ik   xwmax(1,n) = 0.0d0;
!ik   xwmax(2,n) = xdtx1(idmaxd,n); xwmax(3,n) = xdtx2(itmaxd,n);
!ik  >                              xwmax(4,n) = xdatz(ismaxd,n)
!ik   xwmlt(1,n) = 1.0d-6
!ik   xwmlt(2,n) = 1.0d0;   xwmlt(3,n) = 1.0d0;  xwmlt(4,n) = 1.0d0
!ik   xwunt(1,n) = "cm^3/s"
!ik   xwunt(2,n) = "1/cm^3"; xwunt(3,n) = "eV"; xwunt(4,n) = "[]"
!ik   xwspc(1,n) = "log10"
!ik   xwspc(2,n) = "log10";  xwspc(3,n) = "log10"; xwspc(4,n) = "--"
!ik   nw = xdrnk(n) + 1
!ik   ynum = xwnum(1,n)
      xwvar(1,n,kis) = "<sig*v>"
      xwvar(2,n,kis) = "Ne"
      xwvar(3,n,kis) = "Te"
      xwvar(4,n,kis) = "Iz"
      xwnum(1,n,kis) = nszx
      xwnum(2,n,kis) = idmaxd
      xwnum(3,n,kis) = itmaxd
      xwnum(4,n,kis) = ismaxd
      xwmin(1,n,kis) = 0.0d0;
      xwmin(2,n,kis) = xdtx1(1,n,kis)
      xwmin(3,n,kis) = xdtx2(1,n,kis)
      xwmin(4,n,kis) = xdatz(1,n,kis)
      xwmax(1,n,kis) = 0.0d0;
      xwmax(2,n,kis) = xdtx1(idmaxd,n,kis)
      xwmax(3,n,kis) = xdtx2(itmaxd,n,kis)
      xwmax(4,n,kis) = xdatz(ismaxd,n,kis)
      xwmlt(1,n,kis) = 1.0d-6
      xwmlt(2,n,kis) = 1.0d0
      xwmlt(3,n,kis) = 1.0d0
      xwmlt(4,n,kis) = 1.0d0
      xwunt(1,n,kis) = "cm^3/s"
      xwunt(2,n,kis) = "1/cm^3"
      xwunt(3,n,kis) = "eV"
      xwunt(4,n,kis) = "[]"
      xwspc(1,n,kis) = "log10"
      xwspc(2,n,kis) = "log10"
      xwspc(3,n,kis) = "log10"
      xwspc(4,n,kis) = "--"
      nw = xdrnk(n,kis) + 1
      ynum = xwnum(1,n,kis)
      ynum2 = 1
      do i = 2, nw
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   ynum2 = ynum2*xwnum(i,n)
      ynum2 = ynum2*xwnum(i,n,kis)
      enddo
! modified 2/2 lines organize local variables and include files by kamata 2021/06/28
!ik   xjsta(n) = jst
!ik   xjend(n) = jst + ynum - 1
! modified 3/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   xjsta(n) = nyp
!ik   xjend(n) = nyp + ynum - 1
!ik   xjnum(n) = ynum
      xjsta(n,kis) = nyp(kis)
      xjend(n,kis) = nyp(kis) + ynum - 1
      xjnum(n,kis) = ynum
!
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   nyp = xjend(n)
      nyp(kis) = xjend(n,kis)
!
      return
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      mj = indxr(dsn,"/")
! modified 3/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   write(n6,'(/2x,"*** xdt_copy ***  [",a,"]  nxs =",i3,"  cvar =",a
!ik  > ,2x,"devl =",a12,"  drnk =",i2)')  xdnam(n)(1:lenx(xdnam(n)))
!ik  > ,nxs,dsn(mj+1:lenx(dsn)),xdevl(n),xdrnk(n)
      write(n6,'(/2x,"*** xdt_copy ***  [",a,"]  nxs =",i3,"  cvar =",a
     > ,2x,"devl =",a12,"  drnk =",i2)')  trim( xdnam(n,kis) )
     > ,nxs(kis),dsn(mj+1:lenx(dsn)),xdevl(n,kis),xdrnk(n,kis)
      do i = 1, nw
! modified 5/5 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   zmin = xwmin(i,n); zmax = xwmax(i,n)
!ik   write(n6,'(4x,"wvar =",a,"  wnum =",i6,"  wmin =",1pe10.3,
!ik  >  "  wmax =",1pe10.3,"  wmlt =",1pe10.3,"  wunt =",a8,
!ik  >  "  wspc =",a8)') xwvar(i,n),xwnum(i,n),zmin,zmax
!ik  > ,xwmlt(i,n),xwunt(i,n),xwspc(i,n)
      zmin = xwmin(i,n,kis); zmax = xwmax(i,n,kis)
      write(n6,'(4x,"wvar =",a,"  wnum =",i6,"  wmin =",1pe10.3,
     >  "  wmax =",1pe10.3,"  wmlt =",1pe10.3,"  wunt =",a8,
     >  "  wspc =",a8)') xwvar(i,n,kis),xwnum(i,n,kis),zmin,zmax
     > ,xwmlt(i,n,kis),xwunt(i,n,kis),xwspc(i,n,kis)
      enddo
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   zsta = xdaty(xjsta(n)); zend = xdaty(xjend(n))
      write(n6,'(4x,"ynum =",2i8,
     >     "  daty =",i8,1x,1pe10.3,"  daty =",i8,1x,1pe10.3
     >  ,"  dmsz =",2i8)')
! modified 2/2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik  >   ynum,ynum2,xjsta(n),xdaty(xjsta(n))
!ik  >   ,xjend(n),xdaty(xjend(n)), nyp, ndxdt
     >   ynum,ynum2,xjsta(n,kis),xdaty(xjsta(n,kis),kis)
     >   ,xjend(n,kis),xdaty(xjend(n,kis),kis), nyp(kis), ndxdt
!x      write(n6,'(2x,"LNe =",1p10e12.3)') (xdtx1(i,n),i=1,idmaxd)
!x      write(n6,'(2x,"LTe =",1p10e12.3)') (xdtx2(i,n),i=1,itmaxd)
!
!x      do is = 1, ismx
!x      j0 = xjsta(n) + xwnum(2,n)*xwnum(3,n)*(is-1)-1
!x      write(n6,'(2x,"LSV =",i3)') int(xdatz(is,n)+0.001d0)
!x      write(n6,'(2x,5x,3x,2x,8f10.5)') (xdaty(j0+j),j=1,xwnum(2,n))
!x      enddo
!
! modified 3/5 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   write(n6,'(2x," Ne =",1p10e12.3)') (10.0d0**xdtx1(i,n),i=1,idmaxd)
!ik   write(n6,'(2x," Te =",1p10e12.3)') (10.0d0**xdtx2(i,n),i=1,itmaxd)
!ik   if( nyp.gt.ndxdt ) then
      write(n6,'(2x," Ne =",1p10e12.3)')
     >    (10.0d0**xdtx1(i,n,kis),i=1,idmaxd)
      write(n6,'(2x," Te =",1p10e12.3)')
     >    (10.0d0**xdtx2(i,n,kis),i=1,itmaxd)
      if( nyp(kis).gt.ndxdt ) then
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
! modified 4/4 lines organize local variables and include files by kamata 2021/06/28
!ik   write(i6,'(2x,"natm =",i3,"  rtyp =",a)')   natm,cknd
!ik   write(i6,'(2x,"unit =",i3,"  dsninc =",a)') iunt12,dsninc(1:80)
!ik   write(i6,'(2x,"iz0  =",i3,"  iclass =",i3,"  lpart =",l4,
!ik  >  "  ifail =",i3)') iz0,iclass,lpart,ifail
      write(i6,'(2x,"natm =",i3,"  rtyp =",a)')   iz0,cknd
      write(i6,'(2x,"unit =",i3,"  dsn =",a)') nft,dsn(1:80)
      write(i6,'(2x,"iz0  =",i3,"  iclass =",i3,"  lpart =",l4,
     >  "  ifail =",i3)') iz0,nknd,lpart,ifail
      write(i6,'(2x,"ismaxd =",i3,"  izmaxd =",i3,"  itmaxd =",i3,
     >  "  idmaxd =",i3)') ismaxd,izmaxd,itmaxd,idmaxd
      write(i6,'(2x,"dne  =",1p10e11.2)')
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik  >  (10.0**xdtx1(id,n),id=1,idmaxd)
     >  (10.0**xdtx1(id,n,kis),id=1,idmaxd)
      write(i6,'(2x,"dte  =",1p10e11.2)')
! modified 2/2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik  >  (10.0**xdtx2(it,n),it=1,itmaxd)
!ik   write(i6,'(2x,"dne =",i5,2x,1pe11.2)') jd,10.0**xdtx1(jd,n)
     >  (10.0**xdtx2(it,n,kis),it=1,itmaxd)
      write(i6,'(2x,"dne =",i5,2x,1pe11.2)') jd,10.0**xdtx1(jd,n,1)
      write(i6,'(7x,"Te",5x,10(3x,i2,6x))')
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik  >    (int(xdatz(isel,n)+0.001),isel=1,ismaxd)
     >    (int(xdatz(isel,n,kis)+0.001),isel=1,ismaxd)
!x      do it = 1, itmaxd
!x      write(i6,'(2x,1p12e11.2)') 10.0**dtevd(it),
!x     >  (dmax1(10.0**drcofd(isel,it,jd),zero),isel=1,ismaxd)
!x      enddo
      do it = 1, itmaxd
! modified 2/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   write(i6,'(2x,1p12e11.2)') 10.0**xdtx2(it,n),(dmax1(zero,
!ik  >  10.0**xdaty(xjsta(n)-1+jd+(it-1)*idmaxd+(isel-1)*idmaxd*itmaxd))
      write(i6,'(2x,1p12e11.2)') 10.0**xdtx2(it,n,kis),(dmax1(zero,
     >  10.0_8**xdaty(xjsta(n,kis)-1
     >                +jd+(it-1)*idmaxd+(isel-1)*idmaxd*itmaxd,kis))
     > ,isel=1,ismaxd)
      enddo
!
      return
      end
