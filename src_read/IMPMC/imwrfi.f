!***********************************************************************
      subroutine imwrfi(ip,ic,ln,ts,tm,krfl,ker)
!***********************************************************************
!
!        ip   : number of a test partcile
!        ic   : cell number
!        ln   : number of wall segment ln < 0
!        ts   : time  (dummy)
!        tm   : time  (dummy)
!        krfl : 1(pass) 2(gate) 3(refle) 4(refle&pass) 5(stick)
!        ker  : 0 (normal) 1 (error)
!
!    krfl
!     1    pass : Error
!     2    gate : C+ => C0  wght=wt1(pass)   wt1=wt0
!     3,4  refle & wexh     wght=wt1(refle)  wt1=wt0*rflw
!     5    stick            wght=wt0(stick)  wt1=0  (trace stop)
!     6    abs
!
!    Note
!      (vx,vy,vz) => (zvx,zvy,zvz)   vz : common variables of cimcom
!       do not use ip for point number.
!
!       wvaci(mgi,mgo) = 0.0   fin-i = 0.0 at W3g1 and W3g2
!          because no ion can hit the gate.     2012/05/24
!
!VV_5.2  wght > 0.0 when stick => wght = 0.0  wstp = 0.95
!)
!)  No need the statement (if(ien>=6)) in imwrfi/imwrfn
!)
!)     imp_step
!)          |_____if(ien>=6) goto 310
!)          |
!)          |_____imtrci_
!)          |            |___if(ien>=6) goto 150
!)          |            |___imwrfi(kcnd)  set kcnd => ien
!)          |            |___150 continue
!)          |            |___if(kcnd>=6) goto 310
!)          |
!)          |_____imtrcn
!)          |           |____if(ien>=6) goto 550
!)          |           |____imwrfn(kcnd)  set kcnd => ien
!)          |           |____550 continue!
!)          |           |____if(kcnd>=6) goto 310
!)          |
!)          |_____if(ien>=6) goto 310
!)          |_____imshft
!)          |
!)          |_____310 continue
!)
!-----------------------------------------------------------------------
      use cimcom, only : amz, denflxztowall, eneflxztowall, hrfw, ien
     >    , igexh, igpmp, igtno, il, ipcm, ir, is, lpcm, lslf, lwstk
     >    , mypcn, mywgt, pemt, rr, tt, vv, wexh, wght, whit, wrfmin
     >    , wstp, wvaci, wvacn, zz
      use cimntl, only : sint, srn, stb, stm, svx, svy, svz, sx0, sy0
     >    , sz0
      use cntcom, only : chwl, cswl, fnfi, fnth, icwl, igtw, ihgw, ihwl
     >    , nogt, snwl, tcfi, tcth, tsfi, tsth
      use csize,  only : ndmc
      use cunit,  only : mype, n6
      use mod_shexe, only : impmc_model
      implicit none
!
!::argument
      integer, intent(in)    :: ip
      integer, intent(inout) :: ic
      integer, intent(inout) :: ln
      real(8), intent(in)    :: ts
      real(8), intent(out)   :: tm
      integer, intent(out)   :: krfl
      integer, intent(out)   :: ker
!
!::local varables
      integer  ic0, ln0, iz0, iw, ih, igt
      integer  ith, ifi, istk, ien1
      integer  mgi, mgo
      real*8   wt0, wt1, wt2, dwt, rflw, xrn0, xran
      real*8   ptm, wtm, x0, y0, z0
      real*8   v0, vlr1, vlf1, vlz1, csa, sna, csb, snb
      real*8   vlr, vlz, vlf, vlx, vly
      real*8   rw0, zw0, rw1, zw1
      character  cmsg*80
      real*8   v0slf
      real(8)    zfct
! function
      real(8)    random
!
      ic0  = ic
      ln0  = ln
      tm   = ts     !  tm = 0.0 when birth
      krfl = -1
      ker  = 0
      if( impmc_model == 0 ) then
        zfct = 1.0_8
      else
        zfct = pemt(ip)
      endif

!::cond.
      ien1 = 0

      ptm = ts + stb(ip)  ! physical time
      iz0 = is(ip)
      wt0 = wght(ip)
      iw  = -ln0
      ih  = ihwl(iw)
      ic  = icwl(iw)
      igt = igtw(iw)     ! gate number
!
      rflw = hrfw(ih)    ! reflection coefficient
      xrn0 = 0.0d0
!
      if( impmc_model == 0 ) then
        if( ic /= ic0 ) call wexit("imwrfi","ic.ne.ic0")
      else
        if( ic.ne.ic0.and.ic0.lt.ndmc+100 ) then
          write(n6,*) 'ic.ne.ic0', ic, ic0
          call wexit("imwrfi","ic.ne.ic0")
        endif
      endif

!
      istk = lwstk(ih)
      mgi = 0
      if( igtno(ip).gt.0 ) mgi = igtno(ip)
      mgo = ihgw(ih)
!
!::Wg1, Wg2
!::flux into vac-region
      if(ihgw(ih).lt.0 ) then
        mgi = -ihgw(ih)
        igtno(ip) = mgi
        wvaci(mgi,0) = wvaci(mgi,0) + wght(ip)
        wvacn(mgi,0) = wvacn(mgi,0) + wght(ip)  ! Ar2+ => gate ==> Ar0
        mgo = 0                                 !     wvaci   wvacn
      endif
!
!::debug write (A)
!
!::incident flux to wall
      whit(iz0,ih,1) = whit(iz0,ih,1) + wght(ip)   ! n-imp hit
!
! toku, Impurity Flux Onto the Wall --> soldor/pst_wimp.f,  2013.07.04
!                                   --> IMPMC/pst_wimp.f,
! SY modified 2020.03.11: *pemt replaces sflux/emt normalization
      denFlxZtoWall(iz0,iw) = denFlxZtoWall(iz0,iw) + zfct * wght(ip)
      eneFlxZtoWall(iz0,iw) = eneFlxZtoWall(iz0,iw)
     &     + zfct*wght(ip)*0.5d0*amz*vv(ip)
!
!-----------------------------------------------------------------------
!::albedo(-1)
!-----------------------------------------------------------------------
      if( igt.lt.0 ) then
        call wexit("imwrfi","ion hit the crio-panel")
!
!-----------------------------------------------------------------------
!::wall(0)
!-----------------------------------------------------------------------
      elseif( igt.eq.0 ) then
        if( istk.eq.1 ) goto 150        ! stk(5)
        if( lslf.eq.1 ) call slfsput(ic, ip, rflw, v0slf)  ! self-sputtering of C
        if( wght(ip).le.wrfmin ) then
          xrn0 = random(0)
          if( xrn0.le.rflw ) then      ! rfl(3)
            wt1  = wt0
            krfl = 3 ! refle at wall
            goto 130
          else
            goto 160                   ! abs(6)
          endif
        endif
!
        wt1 = wt0*rflw                  ! rfl & abs (4)
        krfl = 4
        goto 130
!
!-----------------------------------------------------------------------
!::gate (1,2)
!-----------------------------------------------------------------------
      elseif( igt.le.nogt ) then
        if( istk.eq.1 ) then           ! stk(5)
! plasma region => Vac region
          rw0 = rr(ip)
          zw0 = zz(ip)
          call tmwgti(rw0,zw0,ic,ln,rw1,zw1,ker)
          if( ker.ne.0 ) goto 940
          rr(ip) = rw1  ! No use
          zz(ip) = zw1  ! NO use
          ir(ip) = ic   ! No use
          il(ip) = ln   ! Use important
          goto 150
        endif
!
        wt1 = wt0                        ! gate(2)
        krfl = 2   ! <== pass gate
        goto 120
!
!-----------------------------------------------------------------------
!::shevron/partition (3,4,5)
!-----------------------------------------------------------------------
      elseif( igt.gt.nogt ) then
        call wexit("imwrfi","ion enter the dome region")
!
!-----------------------------------------------------------------------
!::invalid igt
!-----------------------------------------------------------------------
      else
        goto 920
      endif
!
!-----------------------------------------------------------------------
!::pass  No execution
!-----------------------------------------------------------------------
 110  continue
      call wexit("imwrfi","pass")
!
!-----------------------------------------------------------------------
!::gate
!-----------------------------------------------------------------------
 120  continue
      dwt = wt0 - wt1 ! dwt = 0
      whit(iz0,ih,3) = whit(iz0,ih,3) + wt0
      rw0 = rr(ip)
      zw0 = zz(ip)
      call tmwgti(rw0,zw0,ic,ln,rw1,zw1,ker)
      if( ker.ne.0 ) goto 940
      rr(ip) = rw1  ! No use
      zz(ip) = zw1  ! NO use
      ir(ip) = ic   ! No use
      il(ip) = ln   ! Use important
      goto 200
!
!-----------------------------------------------------------------------
!::reflection
!::refl & ( pass dwt=0.0 or absorp dwt > 0.0 )
!-----------------------------------------------------------------------
 130  continue ! refl at wall abs and refl at wall
      dwt = wt0 - wt1
      whit(iz0,ih,2) = whit(iz0,ih,2) + wt1
      whit(iz0,ih,4) = whit(iz0,ih,4) + dwt   ! 3 => 4
      wexh(iz0,ih) = wexh(iz0,ih) + dwt

! added 1 line integration of ST and TD versionis of IMPMC by kamata 2022/05/14
      if( impmc_model == 1 ) then
!::pump/exhaust
        mywgt(igexh) = mywgt(igexh) + dwt
        mypcn(igexh) = mypcn(igexh) + dwt*pemt(ip)
        if( chwl(ih)(1:1) == "a" ) then
          mywgt(igpmp) = mywgt(igpmp) + dwt
          mypcn(igpmp) = mypcn(igpmp) + dwt*pemt(ip)
        endif
      endif
      wtm = tm    ! <==  Note
      x0 = sx0(ip) + wtm*svx(ip)
      y0 = sy0(ip) + wtm*svy(ip)
      z0 = sz0(ip) + wtm*svz(ip)
      v0 = sqrt(svx(ip)**2+svy(ip)**2+svz(ip)**2)
      if(lslf.eq.1) v0 = v0slf
      iw = -ln
      ic = icwl(iw)
      ith  = int(fnth*random(0)+1.0d0)
      ifi  = int(fnfi*random(0)+1.0d0)
      vlr1 = v0*tsth(ith)*tcfi(ifi)
      vlf1 = v0*tsth(ith)*tsfi(ifi)
      vlz1 = v0*tcth(ith)
      csa  = x0/sqrt(x0**2+y0**2)
      sna =  y0/sqrt(x0**2+y0**2)
      csb = cswl(iw)
      snb = snwl(iw)
      vlr = vlr1*csb - vlz1*snb
      vlz = vlr1*snb + vlz1*csb
      vlf = vlf1
      vlx  = vlr*csa - vlf*sna
      vly  = vlr*sna + vlf*csa
!xx   vlz  = vlz
!
      tt(ip)  = ptm   ! no change
      ir(ip)  = ic
      il(ip)  = ln
!-----
!xx   is(ip)  = 0    ! do not change is in imtrci
      ien(ip) = -1              ! ??? !-----
      sx0(ip) = x0
      sy0(ip) = y0
      sz0(ip) = z0
      svx(ip) = vlx
      svy(ip) = vly
      svz(ip) = vlz
      wght(ip) = wt1
      xran = random(0)
      srn(ip) = -dlog(xran)
      sint(ip) = 0.0d0
      stm(ip) = 0.0d0
      stb(ip) = tt(ip)
      goto 200
!
!-----------------------------------------------------------------------
!::stcik
!-----------------------------------------------------------------------
 150  continue  ! istk = iwstk(ih) = 1
      krfl = 5
      ien1 = 6
      wt1 = 0.0d0
      dwt = wt0 - wt1
      whit(iz0,ih,4) = whit(iz0,ih,4) + dwt  !  3 => 4
      wexh(iz0,ih) = wexh(iz0,ih) + dwt
      goto 200
!
!-----------------------------------------------------------------------
!::abs
!-----------------------------------------------------------------------
 160  continue
      krfl = 6
      wt1 = 0.0d0
      ien1 = 7
      dwt = wt0 - wt1
      whit(iz0,ih,4) = whit(iz0,ih,4) + dwt  !  3 => 4
      wexh(iz0,ih) = wexh(iz0,ih) + dwt
      goto 200
!
!-----------------------------------------------------------------------
!::leave the Vac-region
!-----------------------------------------------------------------------
 200  continue
      if( impmc_model == 1 ) then
        wght(ip) = wt1   ! VV_5.2
        wstp(ip) = 0.0d0
        ien(ip)  = ien1
        dwt = wt0 - wt1
        if( ien1 >= 6 ) then
          wstp(ip) = dwt
        endif
        if( wt1 == 0.0d0 ) then
        write(n6,'(2x,"imwrfi  wt1=0.0",2x,a,2x,i8,i4,1p2e12.4)')
     >    "ip ien wght, wstp =", ip, ien(ip), wght(ip), wstp(ip)
        endif
      endif

!VV_5.2  Particle move out of gate with wt1 (not dwt)
!VV_5.2  Particle is pumped at pump and Vessel dwt

      if( mgo.ge.1 .and. mgo.le.nogt ) igtno(ip) = 0  ! g1/g2
      mgo = ihgw(ih)
      wt2 = dwt
      if( impmc_model == 1 .and. mgo <= nogt ) wt2 = wt1

      if( mgi.gt.0 .and. mgo.gt.0 ) then
        wvaci(mgi,mgo) = wvaci(mgi,mgo) + wt2
      endif

!::debug write (B)
!
      return
!
!::error
 920  continue
      ker = 1
      write(cmsg,'("  mype,lpcm,ipcm =",i5,i11,i6)') mype, lpcm, ipcm
      call wexit("imwrfi","igt invalid number"//trim(cmsg))
!
 930  continue
      ker = 1
      write(cmsg,'("  mype,lpcm,ipcm =",i5,i11,i6)') mype, lpcm, ipcm
      call wexit("imwrfi","next ic"//trim(cmsg))
!
 940  continue
      ker = 1
      write(cmsg,'("  mype,lpcm,ipcm =",i5,i11,i6)') mype, lpcm, ipcm
      call wexit("imwrfi","tmwgtn"//trim(cmsg))
!
      return
      end