!***********************************************************************
      subroutine imtrcn(ip,tmend,kcnd)
!***********************************************************************
!
!    initial kcnd     /-1:loop/
!    return           /0:end/,  /3,ion/,  /6,stick/, /9,err/
!    useless          /1:gen/,  /2:rec/,  /4,wal/,
!-----------------------------------------------------------------------
!
!    ip   (i)   number of sample particle
!    tstep(i)   advance time till t = tstep
!    kcnd (o)   condition
!               -1 (loop) / 0(time limit) / 1(ion) / 2(hit) / 9(error)
!      ien(ip)   1 (neut) / 2(time limit) / 0(ion) / 3(hit) / 4(error)
!
!    input data from cimcom
!        ptim         <==  (tt)
!       (x0, y0, z0)  <==  (rr, zz, ff)
!       (vlx,vly,vlz) <==  (vr, vz, vv)
!       ( ic, ln )    <==  (ir, il)
!       (zincx, zint1, ztim1)  <== (vvr, vvz, v)
!
!    output data to cimcom
!        ptim         ==>  (tt)
!       (xi, yi, zi)  ==>  (rr, zz, ff)
!
!         Do not change is(ip) in this routine !
!         Change is(ip) in imemit (imstep)
!
!-----------------------------------------------------------------------
!     subroutine fndips(ipnm,x0,y0,z0,vlx,vly,vlz,ic,ln,
!    >  xi,yi,zi,vi,vib,kcnd, ptim,ptmx)
!-----------------------------------------------------------------------
!
!    ipnm (i)  number of sample particle
!    x0   (i)  x-coordinate of neutral
!    y0   (i)  y-coordinate of neutral
!    z0   (i)  z-coordinate of neutral
!    vlx  (i)  velocity of neutral
!    vly  (i)  velocity of neutral
!    vlz  (i)  velocity of neutral
!    ic   (i/o)  cell number of emitted neutral(i) and ion(o)
!    ln   (i/o)  line number of emitted neutral(i) and ion(o)
!
!    xinz (o)  x-position of ionization
!    yinz (o)  y-position of ionization
!    zinz (o)  z-position of ionization
!    v2inz (o)  velocity**2 of ion
!    vbinz (o)  velocity along b-field
!
!    kcnd (o)  condition 0(ion) / 1(hit) / 2(time limit) / 9(error)
!    ptim (i/o)  current time
!    ptmx (i)    max time
!
!  common  (cimntl)
!    xincx(i)    -log(random)
!    xintg(i)    integ(-freq*dt)
!    xtime(i)    time after neutral particle was born.
!    scic (o)    cell number
!    scdt (o)    scoreing time
!    nsc  (o)    number of passed cell till ionization of hit the wall
!    ndsc (i)    dimension size of scic, scdt
!
!   Note.    (x,y,z) :  Not (r,z,fai)
!
!     ien(ip)
!     kcnd = -1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9
!           run end puf rec ion wln wli stk abs stk err   see imtrci.f
!
!     debug sybcheck list of trace
!     !xx  iswtr > 0 call impdt
!          call lst_trak
!
!-----------------------------------------------------------------------
      use cimcom, only : cdi, cdi0, denz, i6_trac, ien, igi, il, iml
     >    , ionz, ipcm, ir, iswtr, lpcm, ltout, pemt, rf1_ro, rr, temz
     >    , tt, vzpz, wght, zz
      use cimfls, only : lfls
      use cimntl, only : sint, srn, stb, stm, svx, svy, svz, sx0, sy0
     >    , sz0
      use cntcom, only : mrgn, ncmax, mseg
      use csize,  only : ndmc
      use cunit,  only : n6
      use mod_shexe, only : impmc_model
      use mod_externalgrid
      implicit none
!
!::argument
      integer, intent(in)  :: ip
      integer, intent(out) :: kcnd
      real(8), intent(in)  :: tmend
!
!::local variables
      integer  ic, ln, icb, ker, icn, lnn, iz, istz
      integer  icw, lnw, krfl, kmax, ic_tmp
      real*8   ptim, ptim2, ptmx, zincx, zint1, zint2, ztim1, ztim2
      real*8   x0, y0, z0, vlx, vly, vlz, x1, y1, z1
      real*8   ts, tm, trcz, zdlt, zdtw, zvl2, tmw
      real*8   zi, ri
      real*8   roh, rohb
      real(8)    zfct
      character  comt*(6)
      integer i_loop
      integer, parameter :: i_max=1000000 ! you can change this value
      integer ko
      real*8 x_zdlt, y_zdlt, z_zdlt, r_zdlt

! function
      real(8)    funroh
!
      i_loop = 0
!
      if( impmc_model == 0 ) then
        zfct = 1.0_8
      else
        zfct = pemt(ip)
      endif

      if( impmc_model == 1 ) call lst_trak_head(ip,"ntlSta")

      if( ien(ip).ge.6 ) then
        kcnd = ien(ip)
        ptim = tt(ip)
        ptmx = tmend
        goto 550
      endif
!
 100  continue
      ptim = tt(ip)
      ptmx = tmend
      if( ptim.gt.ptmx ) goto 910
!
      x0   = sx0(ip)
      y0   = sy0(ip)
      z0   = sz0(ip)
      vlx  = svx(ip)
      vly  = svy(ip)
      vlz  = svz(ip)
      ic   = ir(ip)
      ln   = il(ip)
      zincx = srn(ip)
      zint1 = sint(ip)
      ztim1 = stm(ip)
!
      zint2 = zint1
      ztim2 = ztim1
!
      rohb = funroh(ic,rr(ip),zz(ip))

!::start
      istz = 0
      icb  = ic
      ker  = 0
      kcnd = -1   !  loop
!
 200  continue
      if(i_loop > i_max) call wexit("imtrcn","too much loop")
      istz   = istz + 1
      tt(ip) = ptim
      ir(ip) = ic
      il(ip) = ln
      if( impmc_model == 0 ) then
        if( iswtr > 0 ) call impdt(i6_trac,ip,"ntl",ptmx,ptim,stb(ip))
      else
        call lst_trak(ip,"ntlLis",ptmx,ptim,stb(ip))  ! TRAC
      endif
!
      icn = ic
      lnn = ln
      ts  = ztim2
! DBG
      if(.not. use_exdata) then ! ordinary mesh
        call tmtrac(x0,y0,z0,vlx,vly,vlz,icn,lnn,ts,tm,x1,y1,z1,ker)
      else ! the mesh made by Gmesh and external tool
        if(icn.le.ncmax) then ! SOL region
          kmax = mseg(icn)
          call tmtrac_ex(x0,y0,z0,vlx,vly,vlz,icn,lnn
     >      ,ts,tm,x1,y1,z1,ker,kmax)
        elseif(icn .le. ncmax+vac_ele_size) then ! vacume region
          ic_tmp = icn-ncmax
          kmax = mseg_vacume(ic_tmp)
          call tmtrac_grid(x0,y0,z0,vlx,vly,vlz,icn,lnn
     >      ,ts,tm,x1,y1,z1,ker
     >      ,vac_element,vac_ele_size,kmax
     >      ,vac_grid_size,vac_grid_x,vac_grid_y,ic_tmp)
        elseif(icn .le. ncmax+vac_ele_size+pri_ele_size) then ! private region
          ic_tmp = icn-(ncmax+vac_ele_size)
          kmax = mseg_pri(ic_tmp)
          call tmtrac_grid(x0,y0,z0,vlx,vly,vlz,icn,lnn
     >      ,ts,tm,x1,y1,z1,ker
     >      ,pri_element,pri_ele_size,kmax
     >      ,pri_grid_size,pri_grid_x,pri_grid_y,ic_tmp)
        else ! sub-diverter
          ic_tmp = icn-(ncmax+vac_ele_size+pri_ele_size)
          kmax = mseg_subdiv(ic_tmp)
          call tmtrac_grid(x0,y0,z0,vlx,vly,vlz,icn,lnn
     >      ,ts,tm,x1,y1,z1,ker
     >      ,subdiv_cell,cell_size,kmax
     >      ,grid_size,vgx_EX,vgy_EX,ic_tmp)
        endif
      endif
      if( ker.ne.0 ) goto 920
!
      trcz  = cdi(0,ic)    ! <== C, Ar
      if( iml(ip).eq.1 ) then
        trcz  = cdi0(ic)     ! <== CD4
      endif
      ztim1 = tm
      zdlt  = ztim1 - ztim2
!-----
!::end
      ptim2 = ptim + zdlt
      if( ptim2.ge.ptmx ) then
        zdlt  = ptmx - ptim
        kcnd = 0  ! <== end
        icn  = ic
        lnn  = 0
      endif
!-----
!::ion
      zint1 = zint1 + zdlt*trcz
      if( zint1.ge.zincx ) then
        zdlt = (zincx-zint2)/trcz
        kcnd = 3   ! <== ion
        iml(ip) = 0
        ionZ(0,ic) = ionZ(0,ic) + zfct * wght(ip)
        ! check particle is in the new cell or not
        if(is_compareEX .or. use_exdata) then
          x_zdlt = x0+vlx*(ztim2+zdlt)
          y_zdlt = y0+vly*(ztim2+zdlt)
          z_zdlt = z0+vlz*(ztim2+zdlt)
          r_zdlt = sqrt(x_zdlt*x_zdlt+y_zdlt*y_zdlt)
          call mchkin(r_zdlt,z_zdlt,ic,ko)
          if(ko.eq.0) then
            icn = ic
            lnn = 0
          endif
        endif
!:
      endif
!
      ztim1 = ztim2 + zdlt
      zint1 = zint2 + zdlt*trcz
      ptim  = ptim  + zdlt
!
!::scoreing
      iz   = 0
      zdtw = zdlt*wght(ip)
      zvl2 = svx(ip)**2+svy(ip)**2+svz(ip)**2
      denZ(iz,ic) = denZ(iz,ic) + zfct * zdtw
      temZ(iz,ic) = temZ(iz,ic) + zfct * zdtw * zvl2
      if( impmc_model == 0 ) then
        vzpZ(iz,ic) = vzpZ(iz,ic) + zdtw      !debug
      else
        vzpZ(iz,ic) = vzpZ(iz,ic) + 1.d0 !debug
      endif
!
!::next cell
      ztim2 = ztim1
      zint2 = zint1
      icb = ic
      ic  = icn
      ln  = lnn
!
      if( 0 < ic .and. ic <= ndmc ) then
        if( mrgn(ic).eq.6 .or. mrgn(ic).eq.7 ) then
          ri = sqrt(x1**2+y1**2)
          zi = z1
          roh = funroh(ic,ri,zi)
!::impurity flux at the core edge
          if( lfls == 1 ) then
            if( impmc_model == 0 ) then
              call imflos(ic,iz,wght(ip),roh,rohb,kcnd)
            else
              call imflos(ic,iz,pemt(ip)*wght(ip),roh,rohb,kcnd)
            endif
          endif
!::forced ionization due to deep core
          if( roh.lt.rf1_ro ) then
            kcnd = 3   ! <== ion
            iml(ip) = 0
            ionZ(0,ic) = ionZ(0,ic) + zfct * wght(ip)
          endif
          rohb = roh
        endif
      else
        if(ln>0) call wexit("imtrcn","trace fails")
      endif
!
!----------------------------------------------------------------------
!::wall    under development   2009/12/15   2010/05/08
!----------------------------------------------------------------------
      if( ln.lt.0 ) then
        if( impmc_model == 0 ) then
          if( iswtr > 0 )
     >      call impdt(i6_trac,ip,"n_hitw",ptmx,ptim,stb(ip))
        else
          call lst_trak(ip,"n_hitw",ptmx,ptim,stb(ip))  ! TRAC
        endif
        icw = icb
        lnw = ln
        call imwrfn(ip,icw,lnw,ztim2,tmw,krfl,ker)
        if( ker.ne.0 ) goto 920
        ic = icw
        ln = lnw
!
!::pass & gate
        if( krfl.eq.1 .or. krfl.eq.2 ) then
          ztim2 = tmw   ! <== important
          ir(ip) = ic
          il(ip) = ln
          if( impmc_model == 0 ) then
            if( iswtr > 0 )
     >        call impdt(i6_trac,ip,"n_psgt",ptmx,ptim,stb(ip))
          else
            call lst_trak(ip,"n_psgt",ptmx,ptim,stb(ip))  ! TRAC
          endif
          i_loop = i_loop+1
          goto 200
!
!::refle & exh
        elseif( krfl.eq.3 .or. krfl.eq.4 ) then
          igi(ip) = igi(ip) + 1      ! <== n-rfl & pas
          if( impmc_model == 0 ) then
            if( iswtr > 0 )
     >        call impdt(i6_trac,ip,"n_wrfl",ptmx,ptim,stb(ip))
          else
            call lst_trak(ip,"n_wrfl",ptmx,ptim,stb(ip))  ! TRAC
          endif
          i_loop = i_loop+1
          goto 100
!
!::stick
        elseif( krfl.eq.5 ) then
          kcnd = 6
          if( impmc_model == 0 ) then
            if( iswtr > 0 )
     >        call impdt(i6_trac,ip,"n_wstk",ptmx,ptim,stb(ip))
          else
            call lst_trak(ip,"n_wstk",ptmx,ptim,stb(ip))  ! TRAC
          endif
          goto 500

!::abs
        elseif( krfl.eq.6 ) then
          kcnd = 7
          if( impmc_model == 0 ) then
            if( iswtr > 0 )
     >        call impdt(i6_trac,ip,"n_wstk",ptmx,ptim,stb(ip))
          else
            call lst_trak(ip,"n_wabs",ptmx,ptim,stb(ip))  ! TRAC
          endif
          goto 500
!
!::others
        else
          call wexit("imtrcn","krfl invalid number")
        endif
      endif
!
!::loop
      if( kcnd.eq.-1 ) then
        if(istz > ltout) then !.. prevent infinity loop
        !... debug k-itoh
        ! call wexit("imtrcn","kcnd.eq.-1 and istz > ltout")
          kcnd = 7 !.. absorption in this case
          if(iswtr.gt.0) then
            if( impmc_model == 0 ) then
              call impdt(i6_trac,ip,"n_wstk",ptmx,ptim,stb(ip))
            else
              call lst_trak(ip,"n_wstk",ptmx,ptim,stb(ip))  ! TRAC
            endif
          endif
          write(n6,'(2x,"error in imtrcn:  istz > ltout")')
          goto 500
        endif
!
        i_loop = i_loop+1
        goto 200
      endif
!
!::termination of trace
 500  continue
      tt(ip)  = ptim
      ir(ip)  = ic     ! gate in vac
      il(ip)  = ln
      ien(ip) = kcnd   ! 0(ion)/2(tim)/3(hit)/4(err)
      sint(ip) = zint1
      stm(ip) = ztim1  ! time of before step
!
 550  continue
      if( iswtr.gt.0 ) then
        if( impmc_model == 0 ) then
          write(comt,'("ntlE_",i1)') kcnd
          call impdt(i6_trac,ip,comt,ptmx,ptim,stb(ip))
        else
          call lst_trak(ip,"ntlEnd",ptmx,ptim,stb(ip))  ! TRAC
        endif
      endif
      if( kcnd.eq.3 ) call imout_ipos(ip,ptmx,ptim)
      return
!
!::error
 910  continue
      kcnd = 9
      ien(ip) = kcnd
      write(n6,'(2x,"error in imtrcn  ptim.gt.ptmx  ",i10,i8,
     >  "  ptim,ptmx =",1p2e12.3)') lpcm, ipcm, ptim, ptmx
      return
!
 920  continue
      kcnd = 9
      ien(ip) = kcnd
      write(n6,'(2x,"error in imtrcn  ker.ne.0  ",i10,i8,
     >  2x,"  ker =",i3)') lpcm, ipcm, ker
      return
      end
