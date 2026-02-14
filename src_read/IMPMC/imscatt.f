!***********************************************************************
      subroutine imscatt(ip,dtunt)
!***********************************************************************
!     simulate coulomb collisions by using the monte-carlo technique
!       m = 2 : collsion with protons  (save cpu-time)
!
!       wci(ic) = 0  in low Te region when using lpdif (new Dif model)
!
!-----------------------------------------------------------------------
      use cimcns, only : gaus0, gaus1
      use cimcom, only : dfm, fxm1, fxm2, fxm3, ir, is, ishvf, lpsct
     >    , ndmp, pemt, slnv, slw0, slw1, v, vr, vv, vvr, vvz, vz, wght
     >    , wsct
      use cntcom, only : fnfi, tsfi
      use cntpls, only : temi, vflw
      use mod_shexe, only : impmc_model
      implicit none
!
!::argument
      integer, intent(in) :: ip
      real(8), intent(in) :: dtunt
!
!::local variables
      integer  ms, ic, ia, iif, ig0, ig1, iph
      real*8   zvf, dt, dts, slx0, slx1, uz, u, fx1, fx2, fx3
      real*8   sl0, sl1, dwz, dww
      real*8   xran1, xran2, xran3
!
! function
      real(8)    random
!
      lpsct = lpsct + 1
!
      ms = 2
!
      ic  = ir(ip)
!
!::neutral impurity and void region
      if( is(ip).eq.0 ) return
      if( temi(ic).le.0.0d0 ) return
!
!::flow
      ia   = 1
      zvf  = vflw(ic,ia)
      if( ishvf.eq.0 ) zvf = 0.0d0
!
!::step size
      dt    = dtunt*is(ip)**2
      dts   = sqrt(dt)
!
      slx0 = dt*slw0(ic,ms)
      slx1 = dts*slw1(ic,ms)
!
!::effect of plasma flow
      uz  = vz(ip)-zvf
      u   = sqrt(vvr(ip)+uz*uz)
!
      iif = int(dfm*slnv(ic,ms)*u+1.0)
      if( iif.ge.1000 ) iif = 1000
      fx1 = fxm1(iif,ms)
      fx2 = fxm2(iif,ms)
      fx3 = fxm3(iif,ms)
!
!--------------------------------------------------------------------
!
      xran1 = random(0)
      xran2 = random(0)
      xran3 = random(0)
      ig0 = int(2000.0d0*xran1+1.0d0)
      sl0 = slx1*fx2*gaus0(ig0)+slx0*fx1
      ig1 = int(1000.0d0*xran2+1.0d0)
      sl1 = slx1*fx3*gaus1(ig1)
      iph = int(fnfi*xran3+1.0d0)
!
      dwz = (sl0*uz   +sl1*tsfi(iph)*vr(ip))/u
      dww = 2.0*sl0*u +sl0*sl0+sl1*sl1 + 2.0*zvf*dwz
!-----
      vz (ip) = vz(ip) + dwz
      vv (ip) = vv(ip) + dww
!-----
      vvz(ip) = vz(ip)*vz(ip)
      vv (ip) = dmax1(vv(ip),vvz(ip))
      vvr(ip) = vv(ip) - vvz(ip)
      vr (ip) = sqrt(vvr(ip))
      v  (ip) = sqrt(vv(ip))
!-----
      if( impmc_model == 0 ) then
        wsct(ic) = wsct(ic) - wght(ip)*dww
      else
        wsct(ic) = wsct(ic) - pemt(ip)*wght(ip)*dww
      endif
!
      end