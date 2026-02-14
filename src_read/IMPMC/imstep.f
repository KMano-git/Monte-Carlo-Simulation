!***********************************************************************
      subroutine imstep
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
      use cimcns, only : iseed, ndseed
      use cimcom, only : dtimz, dtout, hcpum, hmype, hpmax, hstep, htout
     >    , ien, igi, igtem, ipcm, is, itimz, itg, lpcm, lpsct, ltout
     >    , npmax, sptyi, timez, tout, wght, wpemt, wtemt
      use cimden, only : nwspt
      use csonic, only : lrand
      use cunit,  only : lmype, mype, n6
      use mpi!,    only : mpi_wtime
      implicit none
!
!::local variables
      real(8) xran
      real*8  cptm0, cptm1
      integer ip
      integer ltmx,  jseed, kemt, kcnd
      integer istep, istmx, nerr
      character  cmsg*120
      real*8  tmend
! function
      real(8)    random
!
      cptm0 = MPI_WTIME()
      istmx = 5000
      nerr  = 0
!
!::advance time till t = tout
      dtout = dtimZ*dfloat(ltout)
      tout  = timeZ + dtout
      lpsct = 0
      ltmx = 0
!
!::loop ip
      do ip = 1, npmax
        ipcm  = ip
        istep = 0
!
!::terminate the trace
        if( ien(ip).ge.6 ) goto 400
!
!-----
        if( lrand.eq.1 ) then
          if( itg(ip).gt.ndseed ) call wexit("imstep","itg.gt.ndseed")
          jseed = iseed(itg(ip))
          call srandom(jseed)
          xran = random(0)
        endif
        tmend = tout
!
! generation of neutral or ion impurity
        if( is(ip).lt.0 ) then
          kemt = 1                    ! neutral impurity
          igtem(ip) = 0
          if( sptyi.eq.31 ) then      ! ion impurity
            kemt = 7
            igtem(ip) = -1
          endif
          if( sptyi.eq.41 ) kemt = 8  ! ion impurity source on core edge
          igi(ip) = 1
          call imemit(ip,tmend,kemt)
          wtemt = wtemt + wght(ip)
          wpemt = wpemt + 1.0d0
        endif
        kcnd = 0
!
        do ! istep
          istep = istep + 1
!-----
          if( istep.ge.istmx ) then
            if( mod(istep,istmx).eq.0 ) then
              write(n6,'(2x,
     >         "warning!  too many loop at imstep  istep =",
     >         i7,"  nwspt = ",a,"  npmax =",i6)') 
     >          istep, nwspt, npmax
            endif
            ien(ip) = 9
            nerr = nerr + 1
            write(n6,'(2x,"find error too many loop  ","
     >       lpcm,ipcm =",i10,i8,"  nerr =", i3)') lpcm, ipcm, nerr
            goto 400
          endif
!-----
          kemt = kcnd
          call imemit(ip,tmend,kemt)
          if( ien(ip).ge.6 ) goto 400
!
!::neutral impurity
          if( is(ip).eq.0 ) then
            call imtrcn(ip,tmend,kcnd)
            if( kcnd.eq.0 .or. kcnd.ge.6 ) goto 400
            cycle
          else
!
!::ion impurity
            call imtrci(ip,tmend,kcnd)
            if( kcnd.eq.0 .or. kcnd.ge.6 ) goto 400
            cycle
          endif
        enddo ! istep
!
!::terminate
 400    continue
        ltmx = max( ltmx, istep )
      enddo  ! lp(ip)
!
      itimZ = itimZ + 1
      timeZ = timeZ + dtout
!
!::KSFUJI  hist variables  see New imhist
      cptm1 = MPI_WTIME()
!
!::variables for imstep
      hmype = lmype
      htout = tout
      hstep = ltmx
      hcpum = cptm1-cptm0
      hpmax = npmax
!
      end
