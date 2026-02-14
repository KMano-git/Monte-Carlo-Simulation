!***********************************************************************
      subroutine imatom(ip,dtunt)
!***********************************************************************
!
!    kcnd = -1 (run) : continue calculation
!         =  0 (end) : reach at the time (tmend)
!         =  1 (puf) : emit a neutral partcile
!         =  2 (rec) : recombination of C+  (C+ => C0)
!         =  3 (ion) : ionization of C0     (C0 => C+)
!         =  4 (wln) : C0  strikes on wall  (C0 => wall)
!         =  5 (wli) : Cn+ strike on wall   (C+ => wall)
!         =  6 (wab) : absorption at wall
!         =  9 (err) : error
!
!::spend time
!
!                    tau_atom
!         |<===========================>|
!         |<===>|  dtunt
!  coll.        *
!         C2+   C3+
!  isat   2     3
!  tmat   0.0   t1
!
!         dlt_T(1) = tmat(2) - tmat(1) = dtunt for C2+  Note.
!
!
!                   dtunt
!         |<===========================>|
!         |<============>|<============>|
!   coll.                               *
!        C2+            C2+            C3+
!   isat  2              2              3
!   tmat  t0             t1             t2
!
!        dlt_T(1) = tmat(2) - tmat(1) = t1 - t0  for C2+
!        dlt_T(2) = tmat(3) - tmat(2) = t2 - t1  for C2+  Note.
!
!            dtunt
!      |<===========================>|
!      |<============>|<============>|
!   coll.             *
!      C2+            C3+           C3+
!
!        dlt_T(1) = tmat(2) - tmat(1) = t1 - t0 for C2+
!        dlt_T(2) = tmat(3) - tmat(2) = t2 - t1 for C3+
!
!     The case of latom = 0 is treated in this routine.  2007/12/02
!
!     Error is(ip) = 0  in imtrci
!
!-----------------------------------------------------------------------
      use cimcom, only : cdi, cdr, denz, fforc, friz, ien, il, ionz
     >    , ipcm, ir, irmax, is, itg, latom, lpcm, pemt, recz, rr, temz
     >    , tforc, thfz, tt, vv, vz, vzpz, wght, zz
      use cimptl, only : wic0, wic1, wrr0, wzz0
      use cntcom, only : icwl
      use cunit,  only : lmype, mype, n6
      use mod_shexe, only : impmc_model
      implicit none
!
!::argument
      integer, intent(in) :: ip
      real(8), intent(in) :: dtunt
!
!::local variables for fndway
      integer ndsc; parameter (ndsc=301) !KH160814 101 -> 201
      integer nmv, icmv(ndsc)
      real*8  dtmv(ndsc)
!
!::local variables for scoreing
      integer nat, isat(ndsc)
      real*8  tmat(ndsc)
!
!::local variables (NEW)
      integer nmv1, np
      real*8  tmmv(ndsc)
      real*8  ptm(ndsc), pdt(ndsc)
      integer pis(ndsc), pmv(ndsc)
!
      real*8  sum
      integer  ic_sav
      integer  kcnd
!
!::local variables
!)    dtadv => dtmad  dtscr => dtmsc
      real*8  ztim, zvl2, zsgcx, fcol, dtatm, dt, xx
      real(8)  wrr1, wzz1, dtmad, dtmsc, dt0, fcdtw, zdt0, zwgh
      real*8  xran, xx2, xran2
      integer iz, ic, i, iz0
      integer lon, iw
      integer, parameter :: loop_max = 10000000 ! you can change loop_max value
      real*8, parameter :: fcatm = 1.0d0/20.0d0
! function
      real(8)    random
!
      kcnd = -1
!
      ztim = dtunt
      dt   = 0.0d0
      ic   = ir(ip)
      iz   = is(ip)
      nat  = 1
      isat(nat) = iz
      tmat(nat) = 0.0d0

      zwgh = wght(ip)
      if( impmc_model == 1 ) zwgh = pemt(ip) * zwgh
!
!::neglect atomic process
      if( latom.eq.0 ) then
        nat = 2
        isat(nat) = iz
        tmat(nat) = dtunt
        goto 150 ! skip loop
      endif

!:: atomic process
      do i = 1, loop_max
        if(i.ge.loop_max) call wexit("imatom","too much loop")
!
!::charge exchange recombination
        zvl2  = vv(ip)
        call sgvcxr(ic,iz,zvl2,zsgcx)
!
        if( iz.eq.0 ) then
          fcol = cdi(iz,ic)
        else
          fcol = cdi(iz,ic)+cdr(iz,ic)+zsgcx
        endif
!
        if( ztim == 0.0d0 .or. fcol == 0.0d0 )then
          write(n6,'(2x,"imatom fcol, ztim, zsgcx =",
     >              1p3e12.4)') fcol, ztim, zsgcx
          call wexit("imatom","ztim = 0 or fcol == 0")
        endif
        dtatm = 1.0d0/fcol*fcatm
        dt    = dmin1( ztim, dtatm )
!
!::collision
        xran = random(0)
        xx = xran/(dt*fcol)
!
!-----------------------------------------------------------------------
        iz0 = iz
!
        if( xx.le.1.0d0 ) then
          xran2 = random(0)
          xx2 = xran2/(cdi(iz,ic)/fcol)
!:: Scoring ionization/recombination sources: SY190206
!::ionization
          if( xx2.le.1.0d0 ) then
            ionZ(iz,ic) = ionZ(iz,ic) + zwgh
            iz = iz + 1
!::CXR & recombinaion
          else
            recZ(iz,ic) = recZ(iz,ic) + zwgh
            iz = iz - 1
          endif
        endif
!-----------------------------------------------------------------------
!
!::loop
        ztim = ztim - dt
        if( iz.ne.iz0 .or. ztim.eq.0.0d0 ) then
          nat = nat + 1
          if( nat.gt.ndsc ) call wexit("imatom","nat.gt.ndsc")
          isat(nat) = iz
          tmat(nat) = dtunt - ztim
        endif
        if( (iz.eq.0) .or. (ztim.le.0.0d0) ) exit
      enddo
 150  continue
!
!::particle information
      dtmad = tmat(nat)
      is(ip) = iz
      if( irmax.eq.1 ) goto 200
!
!::trace
      ic_sav = ic !SYamoto
      wrr1 = wrr0 + dtmad*(rr(ip)-wrr0)/dtunt
      wzz1 = wzz0 + dtmad*(zz(ip)-wzz0)/dtunt
!
      call fndway(wrr0,wzz0,wrr1,wzz1,wic0,wic1,ndsc,nmv,icmv,dtmv)
!
!::error
      if( nmv.eq.0 ) then          ! error
        write(n6,'(2x,"error  fndway  nmv =0  mype =",i4,"  ip =",i7,
     >  2x,0p2f11.5," => ",0p2f11.5)') mype, ip, wrr0, wzz0, wrr1, wzz1
        kcnd = 9
        goto 200
      endif
!
!::inside the cell
      if( wic1.gt.0 ) then
        fcdtw = 1.0d0   ! Not hit wall  sumup(dtmv) = 1.0d0
        ic = icmv(nmv)
        lon = 0
!
!::hit the wall   hit point (wrr1,wzz1)
      elseif( wic1.lt.0 ) then
        fcdtw = 0.0d0
        do i = 1, nmv
          fcdtw = fcdtw + dtmv(i)
        enddo
        lon = wic1
        iw  = -lon
        ic  = icwl(iw)
        kcnd = 5       ! <== C+ => wall
!
!::no found the cell
      else
        lon = 0
        kcnd = 9
        write(n6,'(2x,"error at imatom  wic1 = 0 ",2x,"  mype =",i5
     >    "  lpcm =",i10,"  ipcm =",i6)') mype, lpcm, ipcm
        call wexit("imatom","wic1 = 0")
      endif
!
!::final state
      iz = isat(nat)
      tt(ip) = tt(ip) + fcdtw*dtmad
      rr(ip) = wrr1
      zz(ip) = wzz1
      ir(ip) = ic
      is(ip) = iz
      il(ip) = lon
!
!::scoreing (not change charge state)
      if( nat.le.2 ) then
        dtmsc = 0.0d0
        do i = 1, nmv
          dt0   = dtunt*dtmv(i)
          dtmsc = dtmsc + dt0
          ic = icmv(i)
          zdt0  = zwgh * dt0
          denZ(iz,ic) = denZ(iz,ic) + zdt0
          temZ(iz,ic) = temZ(iz,ic) + zdt0 * vv(ip)
          friZ(iz,ic_sav) = friZ(iz,ic_sav) + zdt0 * fforc(ip)
          thfZ(iz,ic_sav) = thfZ(iz,ic_sav) + zdt0 * tforc(ip)
          vzpZ(iz,ic_sav) = vzpZ(iz,ic_sav) + zdt0 * vz(ip)
        enddo
!::scoreing (move to another cell)
!::NEW
      else
        nmv1 = nmv + 1
        if( nmv1.gt.ndsc ) call wexit("imatom","nmv1.gt.ndsc")
        sum = 0.0d0
        do i = 1, nmv
          tmmv(i) = sum
          sum = sum + dtmv(i)*dtmad
        enddo
        tmmv(nmv1) = sum
        icmv(nmv1) = icmv(nmv)
!
!::debug write (NEW)
        call scden(nat,isat,tmat,nmv1,icmv,tmmv
     >                ,np,ptm,pdt,pis,pmv,ndsc)
        if( np.gt.ndsc ) call wexit("imatom","np.gt.ndsc")
!
        dtmsc = 0.0d0
        do i = 1, np-1
          ic = pmv(i)
          iz = pis(i)
          dt0 = pdt(i)
          dtmsc = dtmsc + dt0
          zdt0  = zwgh * dt0
          denZ(iz,ic) = denZ(iz,ic) + zdt0
          temZ(iz,ic) = temZ(iz,ic) + zdt0 * vv(ip)
          friZ(iz,ic_sav) = friZ(iz,ic_sav) + zdt0 * fforc(ip)
          thfZ(iz,ic_sav) = thfZ(iz,ic_sav) + zdt0 * tforc(ip)
          vzpZ(iz,ic_sav) = vzpZ(iz,ic_sav) + zdt0 * vz(ip)
        enddo
      endif
!
 200  continue
      if((kcnd.eq.-1).and.(is(ip).eq.0)) kcnd = 2
!
!::kcnd = 2  recombination  imemit
!::kcnd = 5  hit wall       imemit
!
      if( is(ip).eq.0 ) is(ip) = iz0
      if( is(ip).eq.0 ) is(ip) = 1
      ien(ip) = kcnd
!
!::normal situation  ien = 2(C+=>C0) = 5(C+=>wall)
      if( dtmad < dtunt ) then
        if( ien(ip).ne.2 .and. ien(ip).ne.5 ) then
        write(910,'(2x,"dtmad.lt.dtunt ","lmype =",i5,"  ip =",i7,
     >    "  itg =",i8,"  dtmad,dtunt =",1p2e11.3,"  is,ien =",2i4,
     >    "  normal 2 & 5")')
     >    lmype, ip, itg(ip), dtmad, dtunt, is(ip), ien(ip)
        endif
      endif
!
      return
      end