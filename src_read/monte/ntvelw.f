!***********************************************************************
      subroutine ntvelw(iw_p,vwx,vwy,vwz,vx,vy,vz)
!***********************************************************************
!
!     .vwy      /vwx        nw = 1 : SOL
!       .      /            nw = 3 : PRV
!         .   /
!           ./
!           /A vwz               vfl = vva(j,i,ia)       // b-vector
!          /  *                 vdf = Flux/(na*S^psi)  _|_
!         /     * vny
!        csps(iw_p),snps(iw_p)
!                                vx  = vfl*bx
!     ps : plasma surface        vy  = vfl*by
!                                vz  = vfl*bz
!
!        vwx =  vx*cosw + vy*sinw      SOL surface
!        vwy = -vx*sinw + vy*cosw        vny = -vwy
!        vwz =  vz
!                                      PRV surface
!        vx  = vwx*cosw - vwy*sinw       vny = -vwy
!        vy  = vwx*sinw + vwy*cosw
!        vz  = vwz                     wall: counter clockwise
!
!        cosw = bx/sqrt(bx**2+by**2)
!        sinw = by/sqrt(bx**2+by**2)
!
!        Vdf = Dano/ramd   [m2/s]/[m] = [m/s]
!
!------------------------------------------------------------
!        Vn = V*en = sign(en)*Vf*(-bx*sinw+by*cosw) + Vdif
!               Vdif > 0  -bx*sinw+by*cosw ~ 0
!               sign(en) = - 1 (SOL)  = +1 (PRV)
!
!        Vfwx = Vf*(bx*cosw+by*sinw)
!        Vfwy = Vf*(-bx*sinw+by*cosw)
!        Vfwz = Vf*bz
!
!        vwx = sqrt(T/m)*sqrt(-2*lnY)*cosf + Vfwx
!        vwy = sign(en)*vel_sfn(V*en)
!        vwz = sqrt(T/m)*sqrt(-2*lnY)*sinf + Vfwy
!
!        vx  = vwx*cosw - vwy*sinw
!        vy  = vwx*sinw + vwy*cosw
!        vz  = vwz
!
!----------------------------------------------------------------------
      use cntcom,      only : bvx, bvy, bvz, csbr, flbr, fnfi, icbr
     >    , ikbr, isxp, isyp, iwpbr, snbr, tcfi, tsfi, vion
      use com_phycnst, only : csq2, cspi
      use cplcom,      only : vna, vva
      use cplmet,      only : gare
      implicit none
!
!::argument
      integer, intent(in)  :: iw_p
      real(8), intent(out) :: vx, vy, vz
      real(8), intent(out) :: vwx, vwy, vwz
!
!::local variables
      integer  ik, ia, iw, j, i, ib, ir, ic, ifi
      real*8   sgen, zfi, zar, zna, vfl
      real*8   cosw, sinw, bx, by, bz
      real*8   vfwx, vfwz, vdf, vabs
      real*8   xran
!
!::local variables  (srfsrc_vel)
      real*8   vdr, vth, vn
      real*8   bet, alf2, bet2, x0, gx0, xt, yt
      real*8   xinf, vinf, vmax, xmax, pmax, prb
! function
      real(8)    random
!
      save  sgen, vth, vfwx, vfwz, sinw, cosw
      save  alf2, bet2, vdr, vinf, pmax
!
!::iw   is SOL boundary index(fixed)                                 index of **ps(iw)    see sub. ntwpls.f
!::iw_p is neutral birth SOL boundary index(given by prbr of ntnflx) index of **br(iw_p)  see sub. ntnflx.f
!
      if( iw_p.eq.0 ) goto 150
!
!::ik = 1 (SOL)  = 3 (PRV)
      ik = ikbr(iw_p)
      if( ik.ne.1 .and. ik.ne.3 ) then
        call wexit("ntvelw","ik.ne. 1 or 3")
      endif
      sgen = -1.0d0   ! Note wall circle
!
!::temprally 1-specy
      ia = 1
!
!::parallel velocity & diffusion velocity
      zfi = flbr(iw_p)
      iw  = iwpbr(iw_p)
      j   = isxp(iw)
      i   = isyp(iw)   ! surface
      ib  = i           ! boundary (dummy cell)
      ir  = i + 1       ! real cell
      if( ik.eq.3 ) then
        ib  = i + 1
        ir  = i
      endif
      zar = gare(j,i)
      zna = vna(j,ib,ia)
      vdf = zfi/(zar*zna)
      vfl = vva(j,ir,ia)
!
!::unit vector (wall, en, b-field)
      ic   = icbr(iw_p)
      cosw = csbr(iw_p)
      sinw = snbr(iw_p)
      bx = bvx(ic)
      by = bvy(ic)
      bz = bvz(ic)
!
!::V = Vf*b + Vdf*en  Vn = V*en
      vfwx = vfl*(bx*cosw+by*sinw)
      vfwz = vfl*bz
      vdr  = vdf         ! b*en ~ 0.0
      vth  = vion(ic,ia)
!
!----------------------------------------------------------------------
!::setup for normal velocity
!----------------------------------------------------------------------
!::step 1
      bet  = vth
      bet2 = bet*bet
      x0   = vdr/(csq2*bet)
      gx0  = exp(-x0*x0) + cspi*x0*(1.0d0+erf(x0))
      alf2 = bet2*gx0
!
!::step2
      xt   = 3.0d0
      yt   = exp(-xt*xt) + cspi*x0*(1.0d0-erf(xt))
      yt   = yt/gx0
      xinf = xt
      vinf = vdr + csq2*bet*xinf
!
!::step3
      vmax = 0.5d0*(vdr+sqrt(vdr*vdr+4.0d0*bet2))
      xmax = vmax/vinf
      pmax = vinf**2/alf2*xmax*exp(-(vmax-vdr)**2/(2.0d0*bet2))
!
!----------------------------------------------------------------------
!::normal velocity
!----------------------------------------------------------------------
 150  continue
      xt  = random(0)
      prb = vinf**2/alf2*xt*exp(-(vinf*xt-vdr)**2/(2.0d0*bet2))
      xran = random(0)
      if( xran.gt.prb/pmax ) goto 150
      vn = xt*vinf
!
!----------------------------------------------------------------------
!::coordinate transformation
!----------------------------------------------------------------------
!::vwx,vwy,vwz
      vwy  = sgen*vn
      xran = random(0)
      vabs = vth*sqrt(-2.0d0*dlog(xran))
      ifi  = int(fnfi*random(0) + 1.0d0)
      vwx  = vabs*tcfi(ifi) + vfwx
      vwz  = vabs*tsfi(ifi) + vfwz
!
!::vx,vy,vz
      vx = vwx*cosw - vwy*sinw
      vy = vwx*sinw + vwy*cosw
      vz = vwz
!
      return
      end