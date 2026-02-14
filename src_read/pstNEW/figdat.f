!**********************************************************************
      subroutine plfgdt(cdat,cact)
!**********************************************************************
!
!       cdat = "DTPLV, DTNTL, DTIMP, DTMCP"   pls_end
!       cact = w : write  r : read
!
!---------------------------------------------------------------------
      use csize
      use cimden
      use cntcom
      use cntpls
      use cntsrc
      use cplcom
      use cplmet
      use cimi60
      use cunit
      use csonic
      implicit none

!
!::argument
      character(*) :: cdat,cact
!
!::local variables
!
      integer :: i60
      integer :: nskp = 100
      integer :: kio, nft, ia, ic, jc, n6sv, i, iz, ir, k
      integer :: imox, imoy, icx, icy
      integer :: lp = 0
      integer :: itimb = -999
      character(40) :: cdsn, ched, ctyp
      real*8  ::  hwrsm
      real*8, dimension(10) :: hwrgn
      real(8) ::  twr
      real(8), dimension(10) ::  twr_rg
      real(8), dimension(ndy) :: twr_it
      integer :: nk
      integer, dimension(20) :: mjs, mje
!
      i60 = 180000 + mype
!
      lp = lp + 1
!
      if( cact(1:1).eq."w" ) kio = 1
      if( cact(1:1).eq."r" ) kio = 2
      if( cact(1:1).eq."z" ) kio = 2
!
      if( kio.eq.1 .and. itim.eq.itimb ) then
      write(i60,'(/2x,"*** plfgdtxx ***  lp =",i3,"  itim =",i6,
     >  "  time =",1pe14.6,"  cact = ",a,2x,"  kio =",i2)')
     >  lp, itim, time, cact(1:1), kio
      write(i60,'(2x,"Not executed because of itim = itimB ",2i7)')
     >    itim, itimB
      close(i60)
      return
      endif
!
      write(i60,'(/2x,"*** plfgdt ***  lp =",i3,"  itim =",i6,
     >  "  time =",1pe14.6,"  cact = ",a,2x,"  kio =",i2)')
     >  lp, itim, time, cact(1:1), kio
      write(i60,'(2x,"cdat = ",a,"  limp =",i2)') trim(cdat), limp
!      itimb = itim
!
!::loop
      call linsep(cdat," ,",nk,mjs,mje,20)
      do k = 1, nk
      ctyp = cdat(mjs(k):mje(k))
!
!---------------------------------------------------------------------
!::DTPLS : conservative variables in soldor (q1a,q2a,q3,q4)
!::DTPLV : non-conservative variables       (vna,vva,vti,vte)
!---------------------------------------------------------------------
      if( index(ctyp,"DTPLV").gt.0 ) then
      nft = 21
      cdsn = "./" // trim(ctyp)
      if( kio.eq.2 ) call inqdsn(cdsn)
      call plvdsk(nft,kio,cdsn)
!
      write(i60,'(/2x,"*** plvdsk *** ",a,"  vne,vna,vti,vte")') 
     >   trim(cdsn)
      write(i60,'(2x,"time,qtim =",1p2e14.6)') time, qtim
!
      write(i60,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",2x,"vne",9x,"vne2",
     >  8x,"vni",9x,"vte",9x,"vti")')
      do ic = nskp, ncmax, nskp
        icx = iplx(ic)
        icy = iply(ic)
        ir  = mrgn(ic)
        if( icx.le.0 .or. icy.le.0 ) cycle
        write(i60,'(2x,i6,i5,i4,i3,1p5e12.4)')
     >   ic, icx, icy, ir, vne(icx,icy), vnezef(icx,icy), vni(icx,icy),
     >   vte(icx,icy), vti(icx,icy)
      enddo
      endif
!
!---------------------------------------------------------------------
!::DTNTL : neutral data
!---------------------------------------------------------------------
!   vtime, vflux, vmwork, vitim, vitnt, vnsrc, vsmty, vsmno
!   visrc, vkflx, vksty,  vnsmp, vcsrc
!  --> NTL_2
!
      if( index(ctyp,"DTNTL").gt.0 ) then
      n6sv = n6
      n6   = i60
      nft = 21 ! KH150906
!
      cdsn = "./" // trim(ctyp)
      if( kio.eq.2 ) call inqdsn(cdsn)
      write(i60,'(/2x,"*** ntdisk ***  ",a,"  vflux,vmwork")')
     >    trim(cdsn)
!
      call ntdisk(nft,kio,cdsn)
!
      n6 = n6sv
      endif
!
!---------------------------------------------------------------------
!::DTIMP/DTIPF : impurity data
!---------------------------------------------------------------------
!  nsizp, nsizs, nsizc, nsput, nzmx, ncmx, fsput, twrd, tdnz, csput
!
!  IMP_2: nsput, nzmx, ncmx, fsput, twrd, tdnz, csput
!  nsizp, nsizs, nsizc
!
      if( index( ctyp,"DTIMP" ).gt.0 ) then
      cdsn = "./" // trim(ctyp)
      if( kio.eq.2 ) call inqdsn(cdsn)
      ched = "*** imdisk ***"
      call imdisk(cdsn,cact)
!
      write(i60,'(/2x,a,2x,a,"  twrd,tdnz")') trim(ched), trim(cdsn)
      call totrgn(twrd, twr, twr_rg, twr_it)
      write(i60,'(2x,"nsput =",i2,"  nzmx =",i2,"  ncmx =",i5)')
     >   nsput, nzmx, ncmx
      write(i60,'(2x,5(a," =",1pe12.4,2x))') 
     >  (csput(i), fsput(i), i = 1, nsput)
      write(i60,'(2x,"totwr=",1pe14.6,2x,a)')   twr,
     >   " = Intg_twrd(1:6) + wrd_core"
      call wrintg(twrd,hwrsm,hwrgn)
      write(i60,'(2x,"Intg-twrd =",1pe14.6,2x,1p10e12.4)')
     >   hwrsm,(hwrgn(ir),ir=1,7)
      write(i60,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",2x,"twrd",8x,
     >    6("tdnz_",i1,6x))') (iz,iz=0,5)
      do ic = nskp, ncmax, nskp
      ir = mrgn(ic)
      write(i60,'(2x,i6,i5,i4,i3,1p10e12.4)')
     >  ic, iplx(ic), iply(ic), ir, twrd(ic), (tdnz(iz,ic),iz=0,5)
      enddo
      endif

      if( index( ctyp,"IMPSY" ).gt.0 ) then
      cdsn = "./" // trim(ctyp)
      if( kio.eq.2 ) call inqdsn(cdsn)
      ched = "*** imdisk ***"
      call imdisk(cdsn,cact)
      write(i60,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",8x,
     >    6("tthz_",i1,6x))') (iz,iz=1,6)
      do ic = nskp, ncmax, nskp
      ir = mrgn(ic)
      write(i60,'(2x,i6,i5,i4,i3,1p10e12.4)')
     >  ic, iplx(ic), iply(ic), ir, (tthz(iz,ic),iz=1,6)
      enddo

      write(i60,'(6x,"ic",2x,"icx",1x,"icy",1x,"ir",8x,
     >    6("tvlz_",i1,6x))') (iz,iz=1,6)
      do ic = nskp, ncmax, nskp
      ir = mrgn(ic)
      write(i60,'(2x,i6,i5,i4,i3,1p10e12.4)')
     >  ic, iplx(ic), iply(ic), ir, (tvlz(iz,ic),iz=1,6)
      enddo
      endif
!
!
!---------------------------------------------------------------------
!::DTWRD : radiation in soldor
!---------------------------------------------------------------------
!     wtime, wfac, wcre, wmce, wime  (wmci, wimi)
!
      if( index(ctyp,"DTWRD").gt.0 ) then
      nft = 21
      cdsn = "./" // trim(ctyp)
      if( kio.eq.2 ) call inqdsn(cdsn)
      open(unit=nft,file=cdsn,form="unformatted")
      call plwrdt(nft,kio,cdsn)   ! close
      endif
!
!::loop
      enddo
!
      return
      end
!
!**********************************************************************
      subroutine inqdsn(cdsn)
!**********************************************************************
      use csize
      use cimden
      use cntcom
      use cntpls
      use cntsrc
      use cplcom
      use cplmet
      use cimi60
      use cunit
      implicit none

      character(*) :: cdsn
      logical :: lex
!
      inquire(file=trim(cdsn),exist=lex)
      write(n6,'(2x,"inqdsn  dsn =",a,2x,l2)') trim(cdsn), lex
!
      if(.not.lex) then
      write(n6,'(/2x,"No found file when read ",a)') trim(cdsn)
      call wexit("inqdsn","No found file")
      endif
!
      return
      end

!***********************************************************************
      subroutine chk_tden0(cmsg)
!***********************************************************************
      use csize
      use cimden
      use cntcom
      use cntpls
      use cntsrc
      use cplcom
      use cplmet
      use cimi60
      implicit none

      character(*) :: cmsg
      integer :: ic

!::debug write
      write(6,'(a,2x,"tden0 =",5(2x,1p10e12.3))') trim(cmsg),
     >  (tden0(ic,1),ic=1,ncmax,200)

      return
      end
