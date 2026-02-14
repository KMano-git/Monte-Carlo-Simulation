!***********************************************************************
      subroutine imemit(ip,tmend,kemt)
!***********************************************************************
!
!    kcnd = -1 (run) : continue calculation
!         =  0 (end) : reach at the time (tmend)
!         =  1 (gen) : generation of a neutral impurity (C0 or CD4)
!         =  2 (rec) : recombination of C+  (C+ => C0)
!         =  3 (ion) : ionization of C0     (C0 => C+)
!         =  4 (wln) : C0  strikes on wall  (C0 => wall)
!         =  5 (wli) : Cn+ strikes on wall   (C+ => wall)
!         =  6 (wab) : absorption at wall
!         =  7 (gni) : generation pf ion impurity (Ar+6)
!         =  9 (err) : error
!
!     generation type    see iminit
!     emttyp  "Cchp"  "Cchw"  "Cchd"  "Cphy"  "Arpf"  "Arbk"  "Ipnt"
!     iemtyp   11       12    X 13      14      21      22      31
!
!-----------------------------------------------------------------------
      use cimcns, only : ftom, tomvel
      use cimcom, only : amz, dlmtn, dtimz, icbz, ien, ikbz, i6_trac, il
     >    , ir, is, ismax, iswtr, lslf, rr, sptyc, sptyi, tt, v, v0emtz
     >    , vr, vv, vvr, vvz, vz, wght, zz
      use cimntl, only : sint, srn, stb, stm, svx, svy, svz, sx0, sy0
     >    , sz0
      use cimptl, only : ic_00, is_00, ri_00, zi_00, vv_00, vz_00
      use cntcom, only : cswl, fnfi, fnth, icwl, ipwl, mknd, snwl, tcfi
     >    , tcth, tsfi, tsth, xpnt, ypnt, ncmax
      use cntpls, only : teme
      use cphcns, only : cev, cpi
      use mod_externalgrid
      implicit none
!
!::argument
      integer, intent(in) :: ip, kemt
      real(8), intent(in) :: tmend
!
!::local variables
      integer :: iwem, ic, ik, icx, icy, ln
      integer :: iw, ith, ifi, ipmg, ipmg1, ix, ierr
      real*8  :: x0, y0, z0, v0, vlx, vly, vlz, wgt
      real*8  :: xran, pdtm
      real*8  :: wtm, xi, yi, zi, ri, zvl2, zvlb, tfi  ! tfi  KSFUJI
      real*8  :: vlr1, vlf1, vlz1, vlr, vlf
      real*8  :: xp, yp, zp, viv, vip
      real*8  :: csa, sna, csb, snb, ein, erfl, alf, psr, psz
      real*8  :: xgrid, xgrid_1, ygrid, ygrid_1
      character :: cemt*6
!
      real*8   v0slf, rflw
! function
      integer    imox, imoy
      real(8)    random
!
!----------------------------------------------------------------------
!::generation of neutral impurity
!----------------------------------------------------------------------
      if( kemt.eq.0 ) then
        cemt = "G_emt0"
        goto 100
!
!----------------------------------------------------------------------
!::sputtering at wall
!----------------------------------------------------------------------
      elseif( kemt.eq.1 ) then
        cemt = "G-"//sptyc
        v0 = v0emtZ
!------
        if( v0emtZ.eq.0.0d0 ) then
          xran = random(0)
          ix = int(ftom*xran + 1.0d0)
          v0 = tomvel(ix)
        endif
!-----
        call imstaw(ip,iwem,x0,y0,z0,v0,vlx,vly,vlz)
!
!::cell
        ic  = icbz(iwem)
        ik  = ikbz(iwem)
        icx = imox(ic)
        icy = imoy(ic)
        ln  = mknd(ic,ik)
!
!::check error
        if( .not.use_exdata .and. (icx.le.0 .or. icy.le.0) ) then
          call wexit("imgen_chem"," icx.le.0 .or. icy.le.0")
        endif
!
!::weit
        wgt = 1.0d0
!
!::set common variables
        tt(ip)  = 0.0d0
        ir(ip)  = ic
        il(ip)  = ln
        is(ip)  = 0
        ien(ip) = -1    ! neutral imp.
        sx0(ip) = x0
        sy0(ip) = y0
        sz0(ip) = z0
        svx(ip) = vlx
        svy(ip) = vly
        svz(ip) = vlz
        wght(ip) = wgt
        xran = random(0)
        srn(ip) = -dlog(xran)
        sint(ip) = 0.0d0
        stm(ip) = 0.0d0      ! <== ztim1
        stb(ip) = tt(ip)     ! <== birth time
        goto 100
!
!----------------------------------------------------------------------
!:: ion impurity source on core edge ! toku 2012.10.10
!----------------------------------------------------------------------
      elseif( kemt==8 ) then
        cemt = "G-"//sptyc

!::cell
        call imstace(ip,iwem,x0,y0,z0,vlx,vly,vlz)

        ic  = icbz(iwem)
        icx = imox(ic)
        icy = imoy(ic)

        if( .not.use_exdata .and. (icx.le.0 .or. icy.le.0) ) then
          call wexit("imgen_cedge"," icx.le.0 .or. icy.le.0")
        endif

!::weit
        wgt = 1.0d0
!::set common variables
        tt(ip)  = 0.0d0
        ir(ip)  = ic
        il(ip) = 0
        is(ip)  = ismax
        ien(ip) = 0               ! ion imp.
        sx0(ip) = x0
        sy0(ip) = y0
        sz0(ip) = z0
        svx(ip) = vlx
        svy(ip) = vly
        svz(ip) = vlz
        wght(ip) = wgt
        xran = random(0)
        srn(ip) = -dlog(xran)
        sint(ip) = 0.0d0
        stm(ip) = 0.0d0           ! <== ztim1
        stb(ip) = tt(ip)          ! <== birth time

        ri = sqrt(sx0(ip)**2+sy0(ip)**2)
        zi = sz0(ip)
        zvl2 = svx(ip)**2+svy(ip)**2+svz(ip)**2
        call velion(sx0(ip),sy0(ip),sz0(ip),                              &
     &            svx(ip),svy(ip),svz(ip),zvlb,ierr)
        rr(ip)  = ri
        zz(ip)  = zi
        vv(ip)  = zvl2
        vz(ip)  = zvlb
        vvz(ip) = vz(ip)*vz(ip)
        vv (ip) = dmax1(vv(ip),vvz(ip))
        vvr(ip) = vv(ip) - vvz(ip)
        vr (ip) = sqrt(vvr(ip))
        v  (ip) = sqrt(vv(ip))
        goto 100

!----------------------------------------------------------------------
!::birth of ion impurity for debug
!----------------------------------------------------------------------
      elseif( kemt.eq.7 ) then
        cemt = "G-"//sptyc
        wgt  = 1.0d0
        tt(ip)  = 0.0d0
        ir(ip)  = ic_00
        il(ip)  = 0
        is(ip)  = is_00
        ien(ip) = -1
        rr(ip)  = ri_00
        zz(ip)  = zi_00
        vv(ip)  = vv_00
        vz(ip)  = vz_00
        vvz(ip) = vz(ip)*vz(ip)
        vv (ip) = dmax1( vv(ip), vvz(ip) )
        vvr(ip) = vv(ip) - vvz(ip)
        vr (ip) = sqrt(vvr(ip))
        v (ip)  = sqrt(vv(ip))
        wght(ip) = wgt
        goto 100

!----------------------------------------------------------------------
!::recombination    Ar+ => Ar0
!----------------------------------------------------------------------
      elseif( kemt.eq.2 ) then
      cemt = "gn_rec"
      ri = rr(ip)
      zi = zz(ip)
      tfi = 0.21d0
      viv = sqrt( vv(ip)-vz(ip)**2 )
      vip = vz(ip)
      call velntl(ri,zi,tfi,viv,vip,xp,yp,zp,vlx,vly,vlz)
      is(ip)  = 0
      ien(ip) = -1
      sx0(ip) = xp
      sy0(ip) = yp
      sz0(ip) = zp
      svx(ip) = vlx
      svy(ip) = vly
      svz(ip) = vlz
      xran = random(0)
      srn(ip) = -dlog(xran)
      sint(ip) = 0.0d0
      stm(ip) = 0.0d0
      stb(ip) = tt(ip)
      goto 100
!
!----------------------------------------------------------------------
!::ion    Ar0 => Ar+   (v//, v) are not change
!----------------------------------------------------------------------
      elseif( kemt.eq.3 ) then
      cemt = "gn_ion"
      x0 = sx0(ip)
      y0 = sy0(ip)
      z0 = sz0(ip)
      ic = ir(ip)
      ln = il(ip)
      wtm = tt(ip) - stb(ip)
      xi = sx0(ip) + wtm*svx(ip)
      yi = sy0(ip) + wtm*svy(ip)
      zi = sz0(ip) + wtm*svz(ip)
      ri = sqrt(xi**2+yi**2)
      zvl2 = svx(ip)**2+svy(ip)**2+svz(ip)**2
      call velion(xi,yi,zi,svx(ip),svy(ip),svz(ip),zvlb,ierr)
      !SY dspmtn call for CD dissociation
      if(dlmtn.gt.0.0d0.and.sptyi.ge.11.and.sptyi.le.14) then
         call dspmtn(ip,ic,ln,ri,zi)
      endif
      ir(ip)  = ic
      il(ip)  = ln
      is(ip)  = 1
      ien(ip) = -1
      rr(ip)  = ri
      zz(ip)  = zi
      vv(ip)  = zvl2
      vz(ip)  = zvlb
      vvz(ip) = vz(ip)*vz(ip)
      vv (ip) = dmax1(vv(ip),vvz(ip))
      vvr(ip) = vv(ip) - vvz(ip)
      vr (ip) = sqrt(vvr(ip))
      v  (ip) = sqrt(vv(ip))
      goto 100
!
!----------------------------------------------------------------------
!::reflection at wall and divertor (Ar0 ==> wall)
!       see sub. imtrcn.f    Not passed this part
!----------------------------------------------------------------------
      elseif( kemt.eq.4 ) then
        call wexit("imemit","Not passed this part  Ar0 ==> wall")
!
!----------------------------------------------------------------------
!::reflection at wall and divertor (Ar+ ==> wall)
!----------------------------------------------------------------------
      elseif( kemt.eq.5 ) then
        cemt = "gn_iwal"
        ln = il(ip)
        iw = -ln
        ic = icwl(iw)
        ein = 0.5d0*amz*vv(ip)/cev + 3.0d0*teme(ic)*is(ip)
        erfl = ein*0.8d0
        v0 = 2.0d0*erfl*cev/amz
        v0 = sqrt(v0)
        if( lslf.eq.1 ) then
          call slfsput(ic, ip, rflw, v0slf)  ! self-sputtering of C
          v0 = v0slf
        endif
!
        alf = 2.0d0*cpi*0.0d0    ! no depend on this value
!xx     alf = 2.0d0*cpi*0.721d0  ! confirm  05/02/02
        csa = cos(alf)
        sna = sin(alf)
        xran = random(0)
        ipmg  = ipwl(iw)
        if(use_exdata) then
          ipmg1 = ipwl2(iw)
        else
          ipmg1 = ipwl(iw+1)
        endif
        if(ic.le.ncmax .or. .not.use_exdata) then
          xgrid   = xpnt(ipmg)
          xgrid_1 = xpnt(ipmg1)
          ygrid   = ypnt(ipmg)
          ygrid_1 = ypnt(ipmg1)
        elseif(ic .le. ncmax+vac_ele_size)then
          xgrid   = vac_grid_x(ipmg)
          xgrid_1 = vac_grid_x(ipmg1)
          ygrid   = vac_grid_y(ipmg)
          ygrid_1 = vac_grid_y(ipmg1)
        elseif(ic .le. ncmax+vac_ele_size+pri_ele_size)then
          xgrid   = pri_grid_x(ipmg)
          xgrid_1 = pri_grid_x(ipmg1)
          ygrid   = pri_grid_y(ipmg)
          ygrid_1 = pri_grid_y(ipmg1)
        else
          xgrid   = vgx_EX(ipmg)
          xgrid_1 = vgx_EX(ipmg1)
          ygrid   = vgy_EX(ipmg)
          ygrid_1 = vgy_EX(ipmg1)
        endif
        psr = xgrid + xran*(xgrid_1-xgrid)
        psz = ygrid + xran*(ygrid_1-ygrid)
        x0  = psr*csa
        y0  = psr*sna
        z0  = psz
!
        ith  = int(fnth*random(0)+1.0d0)
        ifi  = int(fnfi*random(0)+1.0d0)
        vlr1 = v0*tsth(ith)*tcfi(ifi)
        vlf1 = v0*tsth(ith)*tsfi(ifi)
        vlz1 = v0*tcth(ith)
!
        csb = cswl(iw)
        snb = snwl(iw)
        vlr = vlr1*csb - vlz1*snb
        vlz = vlr1*snb + vlz1*csb
        vlf = vlf1
!
        vlx  = vlr*csa - vlf*sna
        vly  = vlr*sna + vlf*csa
!
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
        xran = random(0)
        srn(ip) = -dlog(xran)
        sint(ip) = 0.0d0
        stm(ip) = 0.0d0
        stb(ip) = tt(ip)
        goto 100
!
!----------------------------------------------------------------------
!::absorption at wall or error
!----------------------------------------------------------------------
      else
        call wexit("imemit","Not passed this part  wabs or error")
      endif
!
!::debug write
 100  continue
      pdtm = dtimZ
      if( is(ip).eq.0 ) pdtm = stb(ip)
      if( iswtr.eq.1 ) call impdt(i6_trac,ip,cemt,tmend,tt(ip),pdtm)
!
      return
      end
