!***********************************************************************
      subroutine imp_step( tst, dtm )
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
      use cimcom, only : denz, dtemt, dtout, dtscr, dtsiz, hcpum, hmype
     >    , hpmax, hstep, htout, ichcm, ichmx, ien, igtav, igten, ipcm
     >    , is, ismax, istcm, itimz, lpcm, lpsct, mypcn, myptl, mywgt
     >    , npemt, npmax, npmaxb, npmaxl, ntags, pemt, sptyc, timez
     >    , tout, wght, wstp, nbz
      use cimden, only : nsrc_spt, tdnz
      use cntcom, only : ncmax2, volm
      use cunit,  only : lmype, n6
      use mpi!,    only : mpi_wtime
      implicit none

!::argument
      real(8), intent(in) :: tst, dtm
!

!::local variables
      real(8) :: temt, tnow, tmend, zwgt
      real(8) :: ptz, fav
      integer :: itm, itmmx, ip, ic, iz, npext
      integer :: kcnd
      integer :: ltmx, nv
      real(8) :: cptm0, cptm1
!
      cptm0 = MPI_WTIME()
      ichmx = 4000
!:: time (from IMPV5/imp_step.f)
      tnow = tst
      tmend = tst + dtm
      itm = 0
      itmmx = int(dtm/dtemt + 1.0d-9)
      dtsiz = dtm

!::debug write (from IMPV5)
      write(n6,'(2x,a,2x,a,"  tst/dtm/dtemt =",1p3e11.4,
     > "  itmmx =",i5,"  npmax/npemt =",2i7)') "mdl", sptyc,
     >   tst, dtm, dtemt, itmmx, npmax, npemt

!:: imclear routines are moved from immont.f because we do not use immont
!:: in time-dependent calculation 09.10.2019 SY
!:: imclear("particle") is no longer needed
!:: clear
      call imclear("score")
      call imclear("floss")

!-----------------------------------------------------------------
!::Existing particles
!-----------------------------------------------------------------
      ! SY IF ntags statement assuming no imp_pack routine used
      if( npmaxL(nsrc_spt) > 0 ) then
!::advance time till time = tout
        dtout = dtscr
        tout  = timeZ + dtout
        lpsct = 0
        ltmx = 0      ! maximum of step number (istep)
!::loop ip
        npmaxB = npmax
        do ip = 1, npmax
          if(ntags(ip).ne.nsrc_spt) cycle
          ipcm  = ip
          lpcm  = itm
          ichcm = 0   ! Ar0 <==> Ar+++
          istcm = 0   ! step number in imtrci/imtrcn
          kcnd = 0
 210      continue
          if( ien(ip) .ge. 6 ) goto 310

          ichcm = ichcm + 1  ! change Ar0 <==> Ar+ at imsep
          if( ichcm.ge.ichmx ) goto 920

!::neutral impurity  (loop: istcm)  tnow => tmend
          if( is(ip).eq.0 ) then
            call imtrcn(ip,tmend,kcnd)
            if( kcnd.eq.0 .or. kcnd.ge.6 ) goto 310
          else
!
!::ion impurity  (loop: istcm)
            call imtrci(ip,tmend,kcnd)
            if( kcnd.eq.0 .or. kcnd.ge.6 ) goto 310
          endif
!
!::Ar0 <=> Ar+  Ar0/+ => wall  (loop: icnng)
          if( ien(ip).ge.6 ) goto 310
          call imshft(ip,tmend,kcnd)
          if( ien(ip).ge.6 ) goto 310
          goto 210

!::terminate
 310      continue
          ltmx = max( ltmx, ichcm )  ! istcm => ichcm
        enddo  ! lp(ip)
      endif ! npmax > 0

!-----------------------------------------------------------------
!::emitted partcile  ext : extra particles
!-----------------------------------------------------------------
!::loop (itm)
 100  continue
      itm = itm + 1
      npext = npmax + 1
      temt = tnow
      if( nbz > 0 ) then
        call imp_lnch(temt,tmend)
      endif
      npmaxL(nsrc_spt) = npmax
      tnow = tnow + dtemt

!::advance time till time = tout
      dtout = dtscr
      tout  = timeZ + dtout
      lpsct = 0
      ltmx = 0      ! maximum of step number (istep)

!::loop ip
      do ip = npext, npmax
        ipcm  = ip
        lpcm  = itm
        ichcm = 0   ! Ar0 <==> Ar+++
        istcm = 0   ! step number in imtrci/imtrcn

!::terminate trace
        if( ien(ip).ge.6 ) goto 300

        kcnd = 0
 200    continue
        if( ien(ip).ge.6 ) goto 300
        ichcm = ichcm + 1  ! change Ar0 <==> Ar+ at imsep
        if( ichcm.ge.ichmx ) goto 920

!::neutral impurity  (loop: istcm)  tnow => tmend
        if( is(ip).eq.0 ) then
          call imtrcn(ip,tmend,kcnd)
          if( kcnd.eq.0 .or. kcnd.ge.6 ) goto 300
        else
!
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
 300    continue
        ltmx = max( ltmx, ichcm )  ! istcm => ichcm
      enddo  ! lp(ip)

!::loop (lt)
      if( itm < itmmx )  goto 100

!::particle number mypcn at tmend
      do ip = 1, npmax
        if( ien(ip) == 10 ) cycle  ! cut off of death particle
        zwgt = wght(ip)
        if( wstp(ip) > 0.0d0  ) then
          zwgt = wstp(ip)
          if( wght(ip) /= 0.0d0 ) then
            write(n6,'(2x,"imp_step  ihtst ip,is,wght,wstp =",
     >      i8, i4, 1pe12.4, 1pe12.4)') ip, is(ip), wght(ip), wstp(ip)
            call wexit("imp_step","wght > 0 when wstp > 0")
          endif
        endif

!::Wten = Wsum
        myptl(igten) = myptl(igten) + 1
        mywgt(igten) = mywgt(igten) + zwgt
        mypcn(igten) = mypcn(igten) + pemt(ip)*zwgt
      enddo

!::total particle number calculated from avaraged density
      fav = 1.0d0/dtsiz
      nv = ncmax2    ! incl. vac-region
      ptz = 0.0d0
      do ic = 1, nv
        do iz = 0, ismax
          tdnz(iz,ic) = fav*denz(iz,ic)/volm(ic)
          ptz = ptz + tdnz(iz,ic)*volm(ic)
        enddo
      enddo
      mypcn(igtav) = ptz

!::next step
      itimZ = itimZ + 1
      timeZ = timeZ + dtout
!
!::KSFUJI  hist variables  see New imhist
      cptm1 = MPI_WTIME()
!
!::variables for imp_step
      hmype = lmype
      htout = tout
      hstep = ltmx
      hcpum = cptm1-cptm0
      hpmax = npmax
!
!::scatter test partciles
      return

 920  continue
        write(n6,'(2x,"Stop too many loop (ichncm at imp_step  ",
     >   2i6)') ichcm, ichmx
        call wexit("imp_step","ichcm > ichmx")
      end