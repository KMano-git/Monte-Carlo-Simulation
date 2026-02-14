!***********************************************************************
      subroutine imstpbk(tst,dtm)
!***********************************************************************
!
!       kemt /0:start/1:gen/2:rec/3:ion/4:wln/5:wli/6:wab/9:err
!            /7:genI/
!
!       kcnd /0:finish trace/kemt
!
!       set is(ip) = -1  generation of neutral impurity
!                                                   at imp_cal
!
!      istep : step size not dtimZ = 1.0d-6
!
!       emit => ntl => ion => ntl ... => wal => emit ..  => exit
!        1       2      3      4          5      6          128
!                     20 msec                    30 msec => 0.1 sec
!                                                           dtout
!-----------------------------------------------------------------------

!)
!)      dtimZ   dtemt    dtscr    dtprf  dtout
!)      10^-6   10^-4   5*10^-3   10^-1  ???
!)
!)     tst = 0.0d0
!)     dtm = dtscr
!)     call imp_step(tst,dtm)
!)        dtsiz = dtm
!)        denz(iz,ic) = denz(iz,ic)/dtsiz
!)
!)              ___
!)           __|
!)        __|--  emitted   -------------->
!) npext |_________    particle_
!) npmax |
!)   2   |  already emitted  ------------>
!)   1   |__|__|__|__|__|____|__|__|__|__|__|__|
!)       *           A                   *
!)      tst          |                 tmend
!)                 tnow                 tst = tmend
!)          <-> dtemt
!)       <---------------  dtm  --------->
!)                   <-- dtscr -->
!)
!)  prg_prof : dtscr => sub. imp_step(dtm) => dtsiz inc/cimcom
!)  prg_exec : dtexe => sub. imp_step(dtm) => dtsiz inc/cimcom
!)
!)  loop for already emitted partcile
!)   if( npmax > 0 )
!)     do ip = 1, npmax  check  tt(ip) = tst
!)       200 continue
!)       if( is(ip) = 0 ) call imtrcn(ip,tmend,kcnd)
!)       if( is(ip) > 0 ) call imtrci(ip,tmend,kcnd)
!)       if( kcnd == 0 .or. error ) goto 300    when tt(ip) == tmend
!)       goto 200
!)       300 continue
!)     enddo
!)   endif
!)
!)  loop for emitted particle
!)     tnow  = tst
!)     300 continue
!)     npext = npmax + 1
!)     tnow = tnow + dtemt
!)     temt = tnow
!)     call imp_lnch(temt,tmend)  npmax = npmax + npemt  tt(ip) = temt
!)     do ip = npext, npmax
!)       400 continue
!)       if( is(ip) = 0 ) call imtrcn(ip,tmend,kcnd)
!)       if( is(ip) > 0 ) call imtrci(ip,tmend,kcnd)
!)       if( kcnd == 0 .or. error ) goto 500    when tt(ip) == tmend
!)       goto 400
!)     enddo
!)     500 continue
!)     if( tnow < tmend ) goto 300
!)
!)     return
!)
!)::flag to list trace
!x      iswtr = 0
!x      if( i6_trac.gt.0 .and. mype.eq.mpetr ) then
!x        if( ipcm.ge.ip1tr .and. ipcm.le.ip2tr .and.
!x     >      lpcm.ge.lp1tr .and. lpcm.le.lp2tr ) iswtr = 1
!x      endif
!)
!)::conservation check  mypcn
!)     11 mypcn(igtst) : pemt(ip)*weit(ip) at tst
!)     12 mypcn(igten) : pemt(ip)*weit(ip) at ten
!)     13 mypcn(igtav) : Integ[jc] Nzav*dvol
!)     14 mypcn(igtbl) : Wtbl = Wtst + Wemt - Wexh
!)
!)     mypcn(igtbl) =  mypcn(igtst) + mypcn(igemt) - mypcn(igexh)
!)check
!)     mypcn(igten) = mypcn(igtbl)
!)     mypcn(igtav) ~ (mypcn(igtst)+mypcn(igten))/2
!)
!)-------------------------------------------------------------------
      use cimcom, only : bkfrac, dtemt, dtimz, dtout, dtscr, dtsiz
     >    , flximpamp, fwcal, fwstp, hcpum, hmype, hpmax, hstep, htout
     >    , ichcm, ichmx, ien, igtav, igten, ipcm, is, istcm, itimz
     >    , lgstk, lpcm, lprf, lpsct, ltmax, lwstk, mypcn, myptl, mywgt
     >    , npemt, npmax, npmaxb, npmaxl, pemt, sitmz, sptyc, stimz
     >    , timez, tout, wght, wstp, nbz
      use cimden, only : csput, npsput, nsput, nsrc_spt
      use cimpuf, only : bk_mag, bk_npt, bk_nty, lbkstw, pf_mag, pf_nty
      use cntcom, only : chwl, nhwl, nogt
      use csize,  only : ndwh
      use csonic, only : itim
      use cunit,  only : lmype, n6
      use mpi!,    only : mpi_wtime
      implicit none

!::argument
      real(8), intent(in) :: tst, dtm
!
!::local variables
      real(8) :: dtunt, tmend, zwgt
      real(8) :: ptz
      real(8) :: tstbk
      integer :: i, ih, itm, ihst, itmmx, ip, npext
      integer :: kend
      integer :: kcnd
      integer :: ltmx
      real(8) :: cptm0, cptm1
      integer :: sv_nsput
      integer, dimension(5) :: sv_npsput
      integer :: bkitimZ
      real(8) :: bktimeZ
      integer :: sv_lgstk,  sv_ltmax
      real(8) :: sv_fwstp
      integer, dimension(ndwh) :: sv_lwstk
      integer :: sv_itm, sv_lpsct, sv_npemt, sv_nsrc_spt
      real(8) :: sv_dtunt, sv_dtout, sv_tout
      integer :: sv_itimZ
      real(8) :: sv_timeZ
      real(8) :: sv_hmype, sv_htout, sv_hstep, sv_hcpum, sv_hpmax
      integer :: nbks, nbke
      real(8) :: dttrak
!
      real(8) :: totpf
      character(20) :: cmsg
!
      dtunt = 0.0d0
      cptm0 = MPI_WTIME()
      ichmx = 4000
      kend = 0 ! only for BK model
!:: time (from IMPV5/imp_step.f)
      tstbk = tst
      dttrak = 10.d0
      tmend = tstbk + dtm
      itm = 0
      itmmx = int(dtm/dtemt + 1.0d-9)
      dtsiz = dtm
      ihst = 0
!::debug write (from IMPV5)
      write(n6,'(2x,a,2x,a,"  tst/dtm/dtemt =",1p3e11.4,
     > "  itmmx =",i5,"  npmax/npemt =",2i5,"  mdlgn/mdlem =",2i2,
     > "  e0emt =",1pe11.4)') "mdl", sptyc,
     >   tstbk, dtm, dtemt, itmmx, npmax, npemt!, mdlgn, mdlem, e0emt

!
!:: imclear routines are moved from immont.f because we do not use immont
!:: in time-dependent calculation 09.10.2019 SY
!:: imclear("particle") is no longer needed
!:: clear
      call imclear("score")
      call imclear("floss")

!::variables for output
      bk_nty = nogt
      totpf = 0.0d0
      if(csput(1).eq."Arpf") then
        do i = 1, pf_nty
        totpf = totpf + pf_mag(i)
        enddo
      elseif(csput(1).eq."Cedg")then
        totpf = flxImpAmp
      endif
!
!::input data
      if( bk_npt.le.0 ) then
      cmsg = "inpt"
      write(n6,'(2x,"imbkflx  itim =",i7,2x,a,"  totpf =",1pe12.4,
     >  "  bk_mag =",1p10e12.4)') itim, cmsg, totpf, bk_mag(1:bk_nty)
      return
      endif
!
!::input: save parameters
      bkitimZ = itimZ
      bktimeZ = timeZ

!SY  save everything just in case
      sv_nsput = nsput
      sv_npsput(1:5) = npsput(1:5)
      sv_lgstk  = lgstk
      sv_fwstp  = fwstp
      sv_ltmax  = ltmax
      sv_lwstk(1:ndwh) = lwstk(1:ndwh)
      sv_itm    = itm
      sv_dtunt  = dtunt
      sv_dtout  = dtout
      sv_tout   = tout
      sv_lpsct  = lpsct
      sv_hmype  = hmype
      sv_htout  = htout
      sv_hstep  = hstep
      sv_hcpum  = hcpum
      sv_hpmax  = hpmax
      sv_itimZ  = itimZ
      sv_timeZ  = timeZ
      sv_npemt  = npemt
      sv_nsrc_spt = nsrc_spt
!
      nsput = 1
      npsput(1) = bk_npt
      nsrc_spt = 5 !5 reserved for BK: should be OK
      lgstk  = 0
      lbkstw = 0
!      lbkstw = 1  ! debug write in tmwgtn
      lwstk(1:ndwh) = 0
      if( lgstk.eq.1 ) then
        do ih = 1, nhwl
          if( chwl(ih)(3:3).eq."g" ) lwstk(ih) = 1
          if( chwl(ih)(1:1).eq."g" ) lwstk(ih) = 1
        enddo
      endif
!
      fwstp  = 0.20   !  stop condition
      npemt = bk_npt
      itm = 0

!======== FORMER LT TIMELOOP ========
!Particle lauch
      npext = npmax + 1 !npmax before update
      if(nbz > 0) then
        call imp_lnch(0.d0,tmend)
      endif
      npmaxL(nsrc_spt) = npmax
      write(n6,*) '=== IMSTPBK CALLED ==='
      write(n6,*) 'number of test particles:', npmax-npext+1
      cptm1 = MPI_WTIME()
      write(n6,'(2x,a,f12.3)') 'cpu =', cptm1-cptm0

      do !time loop
        ihst = ihst + 1
        if (tstbk+1.0d-12 > dttrak) exit
        tmend = tstbk + dtm
        itm = itm + 1
        dtsiz = dtm
!::advance time till time = tout
        dtunt = dtimZ
        dtout = dtscr
        tout  = timeZ + dtout
        lpsct = 0
        ltmx = 0      ! maximum of step number (istep)
        npmaxB = npmax
!-----------------------------------------------------------------
!::emitted partcile  ext : extra particles
!-----------------------------------------------------------------

!::loop ip
        do ip = npext, npmax
          ipcm  = ip
          lpcm  = itm
          ichcm = 0   ! Ar0 <==> Ar+++
          istcm = 0   ! step number in imtrci/imtrcn

!::terminate trace
          if( ien(ip).ge.6 ) goto 300
          kcnd = 0
 200      continue
          if( ien(ip).ge.6 ) goto 300
          ichcm = ichcm + 1  ! change Ar0 <==> Ar+ at imsep
          if( ichcm.ge.ichmx ) then
            write(n6,'(2x,"Stop too many loop (ichncm at imp_step  ",
     >        2i6)') ichcm, ichmx
            call wexit("imp_step","ichcm > ichmx")
          endif

!::neutral impurity  (loop: istcm)  tnow => tmend
          if( is(ip).eq.0 ) then
            call imtrcn(ip,tmend,kcnd)
            if( kcnd.eq.0 .or. kcnd.ge.6 ) goto 300
          else
!::ion impurity  (loop: istcm)
            call imtrci(ip,tmend,kcnd)
            if( kcnd.eq.0 .or. kcnd.ge.6 ) goto 300
          endif

!::Ar0 <=> Ar+  Ar0/+ => wall  (loop: icnng)
          if( ien(ip).ge.6 ) goto 300
          call imshft(ip,tmend,kcnd)
          if( ien(ip).ge.6 ) goto 300
          goto 200

!::terminate
 300      continue
          ltmx = max( ltmx, ichcm )  ! istcm => ichcm
        enddo                     ! loop ip

!::loop (lt)
!::particle number mypcn at tmend
        do ip = npext, npmax
          zwgt = wght(ip)
          if( wstp(ip) > 0.0d0  ) then
            zwgt = wstp(ip)
            if( wght(ip) /= 0.0d0 ) then
              write(n6,'(2x,"imp_step  ihtst ip,wght,wstp =",
     >        i8,1pe12.4)') ip, wght(ip), wstp(ip)
              call wexit("imp_step","wght > 0 when wstp > 0")
            endif
          endif
!::Wten = Wsum
          myptl(igten) = myptl(igten) + 1
          mywgt(igten) = mywgt(igten) + zwgt
          mypcn(igten) = mypcn(igten) + pemt(ip)*zwgt
        enddo

!::total particle number calculated from avaraged density
        ptz = 0.0d0
        mypcn(igtav) = ptz
!::next step
        bkitimZ = bkitimZ + 1
        bktimeZ = bktimeZ + dtout
!
!::KSFUJI  hist variables  see New imhist
!
!::variables for imstep
        hmype = lmype
        htout = tout
        hstep = ltmx
        hcpum = cptm1-cptm0
        hpmax = npmax
!
!::scatter test partciles
        tstbk = tstbk + dtm
        stimz = tst
        sitmz = itm
        !debug
        if(mod(ihst,100000).eq.0)then
          call imhist(lpcm,kend)
          cptm1 = MPI_WTIME()
          write(n6,'(2x,a,f12.3)') 'cpu =', cptm1-cptm0
        endif
        if( kend.eq.1 ) write(n6,*) 'IMSTPBK: Converged!', tstbk
        if( kend.eq.1 ) exit
      enddo                     !tst

      lprf=1
      call imcndv
      ! SY OUT: BKFRAC
      do i = 1, bk_nty
        bk_mag(i) =  totpf*bkfrac(i)
      enddo
      cptm1 = MPI_WTIME()
      cmsg  = "cal."
      write(n6,'(2x,"imbkflx  itim =",i7,2x,a,"  htout =",f8.3,
     >  "  cpu =",f12.3,"  fwcal/fwstp =",2f8.3,"  totpf =",1pe12.4,
     >  "  bk_mag =",1p10e12.4)') itim, cmsg, htout, cptm1-cptm0,
     >  fwcal, fwstp, totpf, bk_mag(1:bk_nty)
!
      write(n6,'(/2x,"*** imbkflx ***  input data of  imp_cal")')
      write(n6,'(2x,5(1x,a4,1x,i7,2x))') (csput(i),npsput(i),i=1,nsput)
      write(n6,'(2x,"lgstk  = ",i2,"  fwstp = ",f8.4,"  ltmax =",i8)')
     >   lgstk, fwstp, ltmax
      write(n6,'(2x,"chwl = ",15(2x,a,1x))') (chwl(ih),ih=1,nhwl)
      write(n6,'(2x,"lwstk= ",15(2x,i2.2,3x))') (lwstk(ih),ih=1,nhwl)

!    SY clear BK particles
      nbks = npext
      nbke = npmax
      wght(nbks:nbke) = 0.0d0
      call imp_pack
!
!::input  (recovery)
      npmaxL(5) = 0
      nsput = sv_nsput
      npsput(1:5) = sv_npsput(1:5)
      lgstk  = sv_lgstk
      fwstp  = sv_fwstp
      ltmax  = sv_ltmax
      lwstk(1:ndwh) = sv_lwstk(1:ndwh)
      itm    = sv_itm
      dtunt  = sv_dtunt
      dtout  = sv_dtout
      tout   = sv_tout
      lpsct  = sv_lpsct
      hmype  = sv_hmype
      htout  = sv_htout
      hstep  = sv_hstep
      hcpum  = sv_hcpum
      hpmax  = sv_hpmax
      itimZ  = sv_itimZ
      timeZ  = sv_timeZ
      npemt  = sv_npemt
      nsrc_spt = sv_nsrc_spt

      return
      end
