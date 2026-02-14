!***********************************************************************
      subroutine imtrci(ip,tmend,kcnd)
!***********************************************************************
!)
!)   sub. imtrci
!)      Ar+  ==> wall
!)        call imwrfi(krfl)
!)        kcnd <== krfl
!)        ien(ip) = kcnd
!)
!)   krfl = -1  initial set
!)        = 1   No pass for ion    pass at shevron for neutral
!)        = 2   pass gate
!)        = 3   refl/abs at wall when w<wmin
!)        = 4   refl/abs at wall when w>wmin
!)        = 5   stick at wall/gate
!)        = 6   abs
!)
!)   kcnd = -1 (run) : continue calculation
!)        =  0 (end) : reach at the time (tmend)
!)
!)        =  1 (puf) : emit a neutral partcile
!)        =  2 (rec) : recombination of C+  (C+ => C0)
!)        =  3 (ion) : ionization of C0     (C0 => C+)
!)        =  4 (wln) : C0  strikes on wall  (C0 => wall)
!)        =  5 (wli) : Cn+ strike on wall   (C+ => wall)
!)                                     krfl = 2,3,4
!---------
!)    terminate trace
!)      =  6 (stk) : stick at wall  krfl = 5
!)      =  7 (abs) : abs            krfl = 6
!)      =  8 (stk) : juedgment stick at core edge  in imstkce_ion _ntl
!)      =  9 (err) : error
!)
!)         Do not change is(ip) in this routine !
!)         Change is(ip) in imshft (imstep)
!)
!)    trace   if( iswtr > 0 ) call impdt
!)        ==> call lst_trak(ip,"ion_S",ptmx,ptim,dtunt)
!)
!)   kext=1  when ptim2 > ptmx (tmend) => kcnd=0 (reach at the time)
!)
!-----------------------------------------------------------------------
      use cimcns, only : iseed, ndseed
      use cimcom, only : aimas, azmas, dtimz, i6_trac, ien, igi, il
     >    , ipcm, ir, is, iswtr, itg, latom, ldifs, lmrfl, lorbt, lpcm
     >    , lpdif, lscat, lsctmx, pemt, rr, slnv, slw0, wght, tt, zz
      use cimfls, only : lfls
      use cimptl, only : wic0, wis0, wroh0, wrr0, wwt0, wzz0
      use cntcom, only : icwl
      use cntpls, only : dene
      use csize,  only : ndmc
      use csonic, only : lrand
      use cunit,  only : n6
      use mod_shexe, only : impmc_model
      implicit none
!
!::argument
      integer, intent(in)  :: ip
      integer, intent(out) :: kcnd
      real(8), intent(in)  :: tmend
!
!::local variables
      real*8  facdt, cftz
      integer  ltmx, lt, ic, jseed, mscat, lstp, k
      integer  kext, icnd, iw
      integer  kou, icw, lnw, krfl, ker
      real*8   ptim, ptmx, ptim2, dtunt, pdtm
      real*8   tauz, dtstp, tsw, tmw
      character  comt*6
      real*8   xran
      real(8) :: roh
! function
      real(8)    funroh, random
!
      if( impmc_model == 1 ) call lst_trak_head(ip,"ionSta")

      if( ien(ip).ge.6 ) then
        kcnd = ien(ip)
        ptim = tt(ip)
        ptmx  = tmend
        goto 150
      endif
!
!::advance time from ptim to ptmx
      ptim  = tt(ip)
      ptmx  = tmend
      dtunt = dtimZ
      pdtm  = dtunt*1.00001d0
      ltmx  = int((ptmx-ptim)/dtunt + 5.0d0)
      facdt = 1.0d0/20.0d0
      cftz = 1.0d0/dsqrt(aimas)*azmas/(azmas+aimas)
      kcnd = -1
      kext = 0
!
!::debug write
      if( impmc_model == 0 ) then
        if( iswtr.gt.0 ) call impdt(i6_trac,ip,"ion_S",ptmx,ptim,dtunt)
      else
        call lst_trak(ip,"ionL00",ptmx,ptim,dtunt)
      endif
!
!::error
      if( ptim.gt.ptmx ) then
        write(n6,'(2x,"error imtrci  ptim.gt.ptmx ",i10,i8,
     >    1p2e12.3)') lpcm, ipcm, ptim, ptmx
        kcnd = 9
        goto 100
      endif
!
      call imlimt(ip,icnd)
      if( icnd.ne.0 ) then
        write(n6,'(2x,"error imtrci  out of system")')
        kcnd = 9
        goto 100
      endif
!
!-----
      if( lrand.eq.1 ) then
        if( itg(ip).gt.ndseed ) call wexit("imstep","itg.gt.ndseed")
        jseed = iseed(itg(ip))
        call srandom(jseed)
        xran = random(0)
      endif
!-----
!
      do lt = 1, ltmx
        ptim2 = ptim + pdtm    ! avoid too small dtunt (1.0d-15)
!
        if( ptim2.ge.ptmx ) then
          dtunt = ptmx - ptim
          kext  = 1
        endif
        ptim = ptim + dtunt

        if( dtunt .le. 0.0d0 ) then
          write(n6,'(2x,"imtrci dtunt  ptim,ptim2,dtunt =",
     >      1p3e14.6)') ptim, ptim2, dtunt
          kcnd = 0
          ptim = tt(ip)
          ptmx  = tmend
          goto 150
        endif
!
!::step in scattring
        mscat = 1
        ic  = ir(ip)
        if( dene(ic).le.0.0d0 ) then
          lstp = 1
          tauz = 1.0d30
        else
          tauz = 1.0d0/(is(ip)**2*slw0(ic,2)*slnv(ic,2))*cftz
          dtstp = tauz*facdt
          lstp = int(dtunt/dtstp)
          if( lstp.lt.1 ) lstp = 1
          if( lpdif.eq.1 .and. lstp.gt.lsctmx ) mscat = 2
        endif
        dtstp = dtunt/dfloat(lstp)
!
        wrr0 = rr(ip)
        wzz0 = zz(ip)
        wic0 = ir(ip)
        wis0 = is(ip)
        wwt0 = wght(ip)
        wroh0 = funroh(wic0,wrr0,wzz0)
        if((impmc_model == 1).and.(ic.gt.ndmc+10)) then
          write(n6,*) 'faulty ic', ic, 'ip', ip, 'will be deleted'
          kcnd = 5
          krfl = 5
          goto 100
        endif
        call mfndic(wrr0,wzz0,wic0,kou)
!
!----------------------------------------------------------------------
!::hit wall   kou=0 (r,z) in ic  kou<0 hit wall at (r,z) kou=lon
!----------------------------------------------------------------------
        if( kou.lt.0 ) then
          iw = -kou
          ic = icwl(iw)
          rr(ip) = wrr0
          zz(ip) = wzz0
          ir(ip) = ic
          il(ip) = kou
          kcnd = 5
          ien(ip) = kcnd
          goto 100
        endif
!
!----------------------------------------------------------------------
!::various physical process (scatt,orbit,difus,atom)
!----------------------------------------------------------------------
! fforc, tforc
        call imforce(ip,tauz)
        if( mscat.eq.1 ) then
          do k = 1, lstp
            if( lscat.eq.1 ) call imscatt(ip,dtstp)
            if( lorbt.eq.1 ) call imorbit(ip,dtstp)
          enddo  !  == loop(k:local time)
        else
          call impldif(ip,dtunt,tauz)
        endif
!-----
        if( ldifs.eq.1 ) call imdifus(ip,dtunt)
        if( latom.eq.1 ) then
          call imatom(ip,dtunt)
        else 
          call imscor(ip,dtunt)
        endif
!-----
!::Note    hit wall in imatom    kcnd = 5
!
        kcnd = ien(ip)   ! <== return code of imatom/imscor

!::trace till tmend
        if( kcnd.eq.-1 .and. kext.eq.1 ) kcnd = 0
!
!::reflection in the core edge
        if( lmrfl.eq.1 ) call imedrfl(ip,dtunt)
!
!::impurity flux at the core edge
        if( lfls == 1 ) then
          ic = ir(ip)
          roh = funroh(ic,rr(ip),zz(ip))
          if( impmc_model == 0 ) then
            call imflos(ic,is(ip),wght(ip),roh,wroh0,kcnd)
          else
            call imflos(ic,is(ip),pemt(ip)*wght(ip),roh,wroh0,kcnd)
          endif
        endif
!
!::trace
        if( impmc_model == 0 ) then
          if( iswtr.gt.0 ) call impdt(i6_trac,ip,"ion",ptmx,ptim,dtunt)
        else
          call lst_trak(ip,"ionLlt",ptmx,ptim,dtunt)
        endif
!
!::loop (lt:time)
        if( kcnd.ge.0 ) goto 100
      enddo !lt
!
 100  continue
!::reflection at wall
      if( kcnd.eq.5 ) then
        icw = ir(ip)
        lnw = il(ip)
        tsw = tt(ip)  ! dummy
!
!-----
        call imwrfi(ip,icw,lnw,tsw,tmw,krfl,ker)
!-----
        if( krfl.eq.1 ) then
                                 kcnd = 9
        elseif( krfl.eq.2 ) then
                                 kcnd = 5
        elseif( krfl.eq.3 ) then
                                 kcnd = 5
                                 igi(ip) = igi(ip) + 1  ! i-rfl
        elseif( krfl.eq.4 ) then
                                 kcnd = 5
                                 igi(ip) = igi(ip) + 1  ! i-rfl & pas
        elseif( krfl.eq.5 ) then
                                 kcnd = 6
        elseif( krfl.eq.6 ) then
                                 kcnd = 7
        else
                                 kcnd = 9
        endif
      endif
!
      tt(ip) = ptim
      ien(ip) = kcnd
!
 150  continue
      if( impmc_model == 0 .and. iswtr > 0 ) then
        write(comt,'("ionE_",i1)') kcnd
        call impdt(i6_trac,ip,comt,ptmx,ptim,dtunt)
      elseif( impmc_model == 1 ) then
        call lst_trak(ip,"ionEnd",ptmx,ptim,dtunt)
      endif
!
      if( is(ip).eq.0 ) then
        write(n6,'(2x,"error im imtrci  is = 0  lpcm, ipcm =",i10,i8)')
     >      lpcm, ipcm
        call wexit("imtrci","is = 0")
      endif
!
      return
      end