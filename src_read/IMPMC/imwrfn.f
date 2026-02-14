!***********************************************************************
      subroutine imwrfn(ip,ic,ln,ts,tm,krfl,ker)
!***********************************************************************
!
!     ipnm : number of a test partcile
!     ic   : cell number
!     ln   : number of segment ln < 0
!     ts   : time
!     tm   : time  tm-ts ~ 1.0d-12
!     krfl : 1(pass) 2(gate) 3(refle) 4(refl&abs) 5(stick) 6(abs)
!     ker  : 0 (normal) 1 (error)
!
!    GREP var.  wt1 (new weit)
!
!    krfl
!     1,2  pass & gate   wght=wt1(pass)   wt1=wt0
!     3,4  refle & wexh  wght=wt1(refle)  wt1=wt0*rflw
!     5    stick         wght=wt0(stick)  wt1=0  (trace stop)
!     6    abs           wght=wt0(abs)    wt1=0  (trace stop)
!
!    difference of krfl = 5 and 6
!      krfl = 5  lwstk(ih) = 1                     ien(ip) = 6
!      krfl = 6  lwstk(ih) = 0 and wght < wrfmin   ien(ip) = 7
!
!)::ien  see sub. IMPMC/imwclas.f
!
!    ien(ip) = -1(run), 0(tmend), 1(puf),  2(rec),  3(ion),
!              4(wln),  5(wli),   6(stk),  7(abs),  9(err)
!     No trace in case of  ien = 6(stk), 7(abs) 9(err)
!     wght(ip) > 0 in those cases  (flux to wall)
!
!)   terminate trace
!)      in case of  ien = 6(stkWL) 7(absMN) 8(stkCE) 9(err)
!)   wght(ip) > 0 in those cases  (need data when cal. flux to wall/CE)
!)Old version   Do not set wght(ip) = 0.0   This is old version.
!)
!)New VV_5.2    introduce wstp  wght(ip) = 0.0
!)
!)  plasma side    |    Vac region   |  plsma side
!)                 |                 |
!)      wght(ip) wt0               wt1
!)               ==>               ==>
!)                 |____||___________|
!)                      \/
!)                 pump  wt2 = dwt = wt0-wt1 (stick ot pump)
!--------------------------------------------------------------------
!
!    impurity flux under the vac region
!       igtw(iw)  =  -1 (pump)  1,2 (gate)  0 (wall) 3,4 (shevron)
!       ihgw(ih)  =  -2,-1 (Wg), 0 (W), 1,2 (Vg), 3 (a), 4 (V), 5 (inV)
!
!     wvacn(mgi,mgo)   mgi : in  mgo : out
!       mgi = igtno(ip)    igt = 1/2 when Vac-in
!       mgo = gate or pump(a,V)  1,2(=nogt),+1(a),+2(V),+3(cal)
!
!       nogt=2:gate  nopm=1:pump ( JT-60SA, slimCS )
!        lgstk => lwstk(ih) 1: 100% stick
!
!       whit(iz,ih,4)  1(inj)/2(rfl)/3(pas)/4(stk/abs)
!
!    igtw(iw) : extended gate number
!        -1               0             1,2            3,4
!       albedo (rfl/stk)  wal(rfl/stk)  gate (pas|stk) shev (rfl/pas)
!
!-----------------------------------------------------------------------
!    Note
!      (vx,vy,vz) => (zvx,zvy,zvz)   vz : common variables of cimcom
!       do not use ip for point number.
!
!    ip krf  chwl    ih     ln   iw  igt     ic   ixiy   wght    vel
!    rr     zz     xrn0    rflw
!
!ip    ih  chwl    ln     iw igt irg  ic     ix   iy   mgi mgo igtno wvacn
!  zwer  zwgh
!       zwer = 1.0e-13/295.0  zwgh = 0.0d0      2012/05/24   K. Shimizu
!
!-----------------------------------------------------------------------
      use cimcom, only : amz, denflxztowall, eneflxztowall, hrfw, hv0w
     >    , ien, igexh, igpmp, igtno, il, ipcm, ir, is, lpcm, lslf
     >    , lwstk, mypcn, mywgt, pemt, rr, tt, v, vv, vz, wexh, wght
     >    , whit, wrfmin, wstp, wvacn, zz
      use cimntl, only : sint, srn, stb, stm, svx, svy, svz, sx0, sy0
     >    , sz0
      use cntcom, only : chwl, cswl, fnfi, fnth, icwl, igtw, ihgw, ihwl
     >    , mknd, mseg, next, nogt, snwl, tcfi, tcth, tsfi, tsth
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
      integer  ic0, ln0, iz0, iw, ih, igt, k, kf, kcn
      integer  ith, ifi, istk, ien1
      integer  mgi, mgo
      real*8   wt0, wt1, wt2, dwt, rflw, xrn0, xran
      real*8   ptm, wtm, x0, y0, z0, r0, zvx, zvy, zvz, x1, y1, z1
      real*8   v0, v02, vlr1, vlf1, vlz1, csa, sna, csb, snb
      real*8   vlr, vlz, vlf, vlx, vly
      character  cmsg*80
      real*8   v0slf
! function
      real(8)    random
!
      ic0  = ic
      ln0  = ln
      tm   = ts     !  tm = 0.0 when birth
      krfl = -1
      ker  = 0
!
      ien1 = 0      !  ==> temporary definition

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
      if( ic.ne.ic0 ) then
        call wexit("imwrfn","ic.ne.ic0")
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
        wvacn(mgi,0) = wvacn(mgi,0) + wght(ip)
        mgo = 0
      endif
!
!::incident flux to wall
      whit(iz0,ih,1) = whit(iz0,ih,1) + wght(ip)   ! n-imp hit
!
! toku, Impurity Flux Onto the Wall --> soldor/pst_wimp.f,  2013.07.04
!                                   --> IMPMC/pst_wimp.f,
! SY modified 2020.03.11: *pemt replaces sflux/emt normalization
      if( impmc_model == 0 ) then
        denFlxZtoWall(iz0,iw) = denFlxZtoWall(iz0,iw) + wght(ip)
        eneFlxZtoWall(iz0,iw) = eneFlxZtoWall(iz0,iw)
     &       + wght(ip)*0.5d0*amz*(svx(ip)**2+svy(ip)**2+svz(ip)**2)
      else
        denFlxZtoWall(iz0,iw) = denFlxZtoWall(iz0,iw)+pemt(ip)*wght(ip)
        eneFlxZtoWall(iz0,iw) = eneFlxZtoWall(iz0,iw)
     &       + pemt(ip)*wght(ip)*0.5d0*amz*vv(ip)
      endif
!
!-----------------------------------------------------------------------
!::albedo(-1) & wall(0)
!-----------------------------------------------------------------------
      if( igt.le.0 ) then
        if( istk.eq.1 ) goto 150       ! stk(5)
        if( lslf.eq.1 ) call slfsput(ic, ip, rflw, v0slf)  ! self-sputtering of C
        if( wght(ip).le.wrfmin ) then
          xrn0 = random(0)
          if( xrn0.le.rflw ) then      ! rfl(3)
            wt1  = wt0
            krfl = 3
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
        if(istk.ne.1) then
          wt1 = wt0                        ! gate(2)
          krfl = 2   ! <== pass gate
          whit(iz0,ih,3) = whit(iz0,ih,3) + wt0
        endif
!::
! plasma region => Vac region
        x0  = sx0(ip)
        y0  = sy0(ip)
        z0  = sz0(ip)
        zvx = svx(ip)
        zvy = svy(ip)
        zvz = svz(ip)
        call tmwgtn(x0,y0,z0,zvx,zvy,zvz,ic,ln,ts,tm,x1,y1,z1,kcn)
        if( impmc_model == 0 ) then
          if( kcn /= 0 ) goto 940
        else
          if( kcn /= 0 ) goto 160 ! absorb
        endif
!
        if( istk.eq.1 ) then
          goto 150  ! stk(5)
        else
          dwt = wt0   ! Note change Vac <=> Plasma
          goto 200
        endif
!-----------------------------------------------------------------------
!::shevron/partition (3,4,5)
!-----------------------------------------------------------------------
      elseif( igt.gt.nogt ) then
        if( istk.eq.1 ) goto 150       ! stk(5)
        xrn0 = random(0)
        if( xrn0.le.rflw ) then
          wt1 = wt0                    ! rfl(3)
          krfl = 3  ! <== reflection
          goto 130
        else
          wt1 = wt0                    ! pas(1)
          krfl = 1  ! <==pass
          if( impmc_model == 1 ) dwt = wt0 - wt1
          goto 110
        endif
!
!-----------------------------------------------------------------------
!::invalid igt
!-----------------------------------------------------------------------
      else
        goto 920
      endif
!
!-----------------------------------------------------------------------
!::pass
!-----------------------------------------------------------------------
 110  continue
      whit(iz0,ih,3) = whit(iz0,ih,3) + wt0
      do k = 1, mseg(ic0)
        kf = k
        if( mknd(ic,k).eq.ln0 ) goto 20
      enddo
      goto 930   !  error shev/part
 20   continue
      ic = next(ic,kf)
      ln = mknd(ic,kf)
      if( ic.le.0 ) goto 930
      dwt = 0.0d0
      goto 200
!
!-----------------------------------------------------------------------
!::reflection
!::refl & ( pass dwt=0.0 or absorp dwt > 0.0 )
!-----------------------------------------------------------------------
 130  continue
      dwt = wt0 - wt1
      whit(iz0,ih,2) = whit(iz0,ih,2) + wt1
      whit(iz0,ih,4) = whit(iz0,ih,4) + dwt   ! 3 => 4
      wexh(iz0,ih) = wexh(iz0,ih) + dwt

      if( impmc_model == 1 ) then
!::pump/exhaust
        mywgt(igexh)  = mywgt(igexh)  + dwt
        mypcn(igexh)  = mypcn(igexh)  + dwt*pemt(ip)
        if( chwl(ih)(1:1) == "a" ) then
          mywgt(igpmp) = mywgt(igpmp) + dwt
          mypcn(igpmp) = mypcn(igpmp) + dwt*pemt(ip)
        endif
      endif

!
      wtm = tm    ! <==  Note
      x0 = sx0(ip) + wtm*svx(ip)
      y0 = sy0(ip) + wtm*svy(ip)
      z0 = sz0(ip) + wtm*svz(ip)
      v0 = hv0w(ih)
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
!
      tt(ip)  = ptm   ! no change
      ir(ip)  = ic
      il(ip)  = ln
      is(ip)  = 0
      ien(ip) = -1
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
!::stick
!-----------------------------------------------------------------------
 150  continue  !  istk = lwstk(ih) = 1
      krfl = 5
      wt1 = 0.0d0
      ien1 = 6
      dwt = wt0 - wt1
      whit(iz0,ih,4) = whit(iz0,ih,4) + dwt  !  3 => 4
      wexh(iz0,ih) = wexh(iz0,ih) + dwt
      goto 200
!
!-----------------------------------------------------------------------
!::abs
!-----------------------------------------------------------------------
 160  continue
      krfl = 6                   ! abs(6)
      ien1 = 7
      wt1 = 0.0d0
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
        wght(ip)= wt1 !SY
        wstp(ip) = 0.0d0
        ien(ip)  = ien1

        if( ien1 >= 6 ) then
          wstp(ip) = dwt  ! VV_5.2
        endif

!::terminate trace  update hit point
        if( ien1 >= 6 ) then
          wtm = tm    ! <==  Note
          x0  = sx0(ip) + wtm*svx(ip)
          y0  = sy0(ip) + wtm*svy(ip)
          z0  = sz0(ip) + wtm*svz(ip)
          r0  = sqrt(x0**2+y0**2)
          v02 = svx(ip)**2+svy(ip)**2+svz(ip)**2
          v0  = sqrt(v02)
          rr(ip) = r0
          zz(ip) = z0
          vv(ip) = v02
          vz(ip) = 0.0d0
          v(ip)  = v0
        endif
      endif

!VV_5.2  Particle move out of gate with wt1 (not dwt)
!VV_5.2  Particle is pumped at pump and Vessel dwt

      if( mgo.ge.1 .and. mgo.le.nogt ) igtno(ip) = 0  ! g1/g2
      mgo = ihgw(ih)
      wt2 = dwt
      if( impmc_model == 1 .and. mgo <= nogt ) wt2 = wt1

      if( mgi.gt.0 .and. mgo.gt.0 ) then
        wvacn(mgi,mgo) = wvacn(mgi,mgo) + wt2
      endif
!
      return
!
!::error
 920  continue
      write(cmsg,'("  mype,lpcm,ipcm =",i5,i11,i6)') mype, lpcm, ipcm
      call wexit("imwrfn","igt invalid number"//trim(cmsg))
!
 930  continue
      write(cmsg,'("  mype,lpcm,ipcm =",i5,i11,i6)') mype, lpcm, ipcm
      call wexit("imwrfn","next ic"//trim(cmsg))
!
 940  continue
 ! it offen occurs by numerical error, so wexit is not called
      write(n6,'(2x,"imwrfn: tmwgtn errr, kcn=",i4)') kcn
      ker = kcn
      write(10000+mype,'(2x,"*** imwrfn: tmwgtn errr ***")')
      write(10000+mype,'(2x,"ip,r0,z0,r1,z1,vx,vy,vz",i5,7f14.8)')
     >  ip, sqrt(x0**2+y0**2), z0, sqrt(x1**2+y1**2), z1, zvx,zvy,zvz

      return
      end