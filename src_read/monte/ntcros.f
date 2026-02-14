!**********************************************************************
      subroutine ntcros
!**********************************************************************
!
!    calculate cross section
!
!    sub. ntcrs_XXX(ic,pvel)
!
!    input
!         XXX  :  = MH0 (D0), = MH2 (D2), = MHP (D2+), AHP = (D+)
!         ic   :  cell number
!         pvel :  neutral velocity         [m/s]
!
!    output (/cntcom/)
!         mrct       :  = 1 (D0), = 2 (D2), = 3 (D2+), = 4 (D+)
!         nrct       :  number of reaction
!         dmfp       :  mean free path           [m]
!         trct       :  sum_j{den_j*<sigv>_j}    [1/s]
!         trcx       :  sum_j{den_j*<sigv>_j} except for inonization
!         wrct       :  after collission weit = weit*wrct
!         frct(nrct) :  cumulative probability   []
!         srct(nrct) :  den_j*<sigv>_j           [1/s]
!         elrc(nrct) :  energy loss of electron  [eV]   when DS/EI/RC
!         derc(nrct) :  dissociation energy      [eV]   when DS
!         irct(nrct) :  index of reaction
!         igrc(nrct) :  index of plasma ion = -1 (electron)
!
!    Note.
!     wrct = sum_cx(freq) / sum_all(freq) = trcx/trct
!                                   (suppresion of absorption)
!     elrc : positive means loss
!     label of reaction   lbrc(irct(i))
!     reaction            ctrc(irct(i))
!
!    reaction
!      mrct  irct  monte  label    reacion
!                  -part.
!       1            D0
!             1           "EI"    [ e  + H0  =>  H+  + 2e ]
!             2           "CX"    [ H+ + H0  =>  H0  + H+ ]
!             3           "EL1"   [ H+ + H0  =>  H+  + H0 ]
!       2            D2
!             4           "DS1"   [ e  + H2  =>  H0  + H0 + e ]
!             5           "DS2"   [ e  + H2  =>  H2+ + 2e ]
!             6           "DS3"   [ e  + H2  =>  H+  + H0 + e ]
!             7           "EL2"   [ H+ + H2  =>  H+ + H2 ]
!       3            D2+
!             8           "DS4"   [ e  + H2+ =>  H+  + H0 + e ]
!             9           "DS5"   [ e  + H2+ =>  H+  + H+ + 2e ]
!            10           "DS6"   [ e  + H2+ =>  H0  + H0 ]
!       4            D+
!            11           "RC"    [ e  + H+  =>  H0 ]
!
! 2015.10.02, n-n collision implementation by toku
!                  -part.
!       1            D0
!            12           "EL3"   [ H0 + H0  =>  H0  + H0 ] kel=3
!            13           "EL4"   [ H0 + H2  =>  H0  + H2 ] kel=5
!       2            D2
!            14           "EL5"   [ H2 + H0  =>  H2  + H0 ] kel=4
!            15           "EL6"   [ H2 + H2  =>  H2  + H2 ] kel=6
!
!      c.x reaction   Ni*<sig*v>cx          -> ->
!              vrel**2 = 8*Ti*eV/(pi*Mi) + (v0-Vf)**2     2003/2/13
!
!----------------------------------------------------------------------
      use cntcom, only : ctrc, lbrc, lelmd, lnnel, ndrmx, nrcmx, rcfc
     >    , scfc
      use cunit,  only : n6
      implicit none
!
      integer lopcrs
      parameter (lopcrs=1)
!
!::local variables
      integer  nsiz, i
! function
      integer  lenx
!
! <H0>
      lbrc( 1) = "EI";  ctrc( 1) = "[ e  + H0  =>  H+  + 2e ]"
      lbrc( 2) = "CX";  ctrc( 2) = "[ H+ + H0  =>  H0  + H+ ]"
      lbrc( 3) = "EL1"; ctrc( 3) = "[ H+ + H0  =>  H+  + H0 ]"
! <H2>
      lbrc( 4) = "DS1"; ctrc( 4) = "[ e  + H2  =>  H0  + H0 + e ]"
      lbrc( 5) = "DS2"; ctrc( 5) = "[ e  + H2  =>  H2+ + 2e ]"
      lbrc( 6) = "DS3"; ctrc( 6) = "[ e  + H2  =>  H+  + H0 + e ]"
      lbrc( 7) = "EL2"; ctrc( 7) = "[ H+ + H2  =>  H+ + H2 ]"
! <H2+>
      lbrc( 8) = "DS4"; ctrc( 8) = "[ e  + H2+ =>  H+  + H0 + 2e ]"
      lbrc( 9) = "DS5"; ctrc( 9) = "[ e  + H2+ =>  H+  + H+ + e ]"
      lbrc(10) = "DS6"; ctrc(10) = "[ e  + H2+ =>  H0  + H0 + e ]"
! <H+>
      lbrc(11) = "RC";  ctrc(11) = "[ e  + H+  =>  H0 ]"
! <D0> n-n
      lbrc(12) = "EL3"; ctrc(12) = "[ D0 + D0  =>  D0  + D0 ]"
      lbrc(13) = "EL4"; ctrc(13) = "[ D0 + D2  =>  D0  + D2 ]"
! <D2> n-n
      lbrc(14) = "EL5"; ctrc(14) = "[ D2 + D0  =>  D2  + D0 ]"
      lbrc(15) = "EL6"; ctrc(15) = "[ H2 + H2  =>  H2  + H2 ]"
!
      nrcmx = 15
!
      write(n6,'(/2x,"*** ntcros ***   cross section  New vrel2")')
      write(n6,'(4x,"cx reaction   vrel2 = 8Ti/(pi*Mi) + (v0-Vf)**2")')
!
!----------------------------------------------------------------------
      if( lopcrs.eq.0 ) then
      write(n6,'(4x,"set rcsc & scfc ")')
      nsiz = ndrmx
      call setd( rcfc, nsiz, 0.0d0 )
      call setd( scfc, nsiz, 0.0d0 )
!
      do i = 1, nrcmx
      rcfc(i) = 1.0d0
      scfc(i) = 1.0d0
      enddo
      endif
!
!::option   elastic collision
      write(n6,'(2x,"lelmd =",i3)') lelmd
      write(n6,'(2x,"when lelmd = 0, no effect elastic collifion")')
      if( lelmd.eq.0 ) then
        rcfc(3) = 0.0d0; rcfc(7) = 0.0d0
        scfc(3) = 0.0d0; scfc(7) = 0.0d0
      endif
!
      if( lnnel==0)then
            rcfc(12:15) = 0.0d0
            scfc(12:15) = 0.0d0
            write(n6,'(2x,"lnnel=0. scfc set to zero.")')
      endif
!
!::debug write
      do i = 1, nrcmx
      write(n6,'(2x,i3,2x,a,2f6.1,4x,a)') i,lbrc(i),rcfc(i),scfc(i)
     >  ,ctrc(i)(1:lenx(ctrc(i)))
      enddo
!
      return
      end
!
!**********************************************************************
      subroutine ntcrs_AH0(ic,pvel)
!**********************************************************************
      use BGProf, only : den0BG, dengBG, eng0BG, enggBG, nfl0BGx
     >    , nfl0BGy, nfl0BGz, nflgBGx, nflgBGy, nflgBGz
      use celcom, only : el_atm, el_cros
      use cntcom, only : bvx, bvy, bvz, derc, dmfp, elrc, frct, igrc
     >    , irct, lnnel, mrct, ngas, nrct, rcfc, rmas, srct, tbel_ei
     >    , tbsg_ei, trct, trcx, velx, vely, velz, wrct
      use cntpls, only : dene, deni, temi, vflw
      use cphcns, only : cev, cmp, cpi
      use csize,  only : ndgs
      use KrsticSchultz3, only : sigma_DxD2_t, sigma_DxD_t
      implicit none
!
!::argument
      integer, intent(in) :: ic
      real(8), intent(in) :: pvel
!
!::local variables
      integer ii, irc, ig, i
      real*8  zne, zti, zngbg, zn0bg
      real*8  zvntl2, zvrel2, zerel, zvrel, zscx, zscxv, zsum
      real*8  tb_erel(ndgs),tb_vrel(ndgs)
      real(8) :: v00rel2, v00rel, E00rel, v0grel2, v0grel, E0grel,
     >     sig00, sigv00, sig0g, sigv0g, lnE
! function
      real(8)    xdt_eval1
!
      mrct = 1
      ii = 0
!
!--plasma parameter
      zne = dene(ic)
      zti = temi(ic)
      zngbg= dengBG(ic,1)
      zn0bg= den0BG(ic,1)
      if(lnnel==1) then
        if( zne+zn0bg+zngbg.le.0.0d0 ) goto 100
      else
        if( zne.le.0.0d0 ) goto 100
      endif
!
!--EI [ e  + H0  =>  H+  + 2e ]
      irc = 1
      ii  = ii + 1
      irct(ii) = irc
      srct(ii) = tbsg_ei(ic)
!xx   elrc(ii) = elhion(zte,zne)
      elrc(ii) = tbel_ei(ic)
      derc(ii) = 0.0d0
      igrc(ii) = -1
!
!--CX [ H+ + H0  =>  H0  + H+ ]
      irc = 2
      do ig = 1, ngas
!xx   zvntl2 = pvel*pvel  ! 2003/02/13
      zvntl2 = (velx-vflw(ic,ig)*bvx(ic))**2
     >        +(vely-vflw(ic,ig)*bvy(ic))**2
     >        +(velz-vflw(ic,ig)*bvz(ic))**2
      zvrel2 = 8.0d0*zti*cev/(rmas(ig)*cpi) + zvntl2
      zerel  = 0.5d0*cmp*zvrel2/cev
      zvrel  = dsqrt(zvrel2)
      zscx   = 0.6937d-18*(1.0-0.155*dlog10(zerel))**2
     >         /(1.0+0.1112d-14*zerel**3.3)
      zscxv  = zscx*zvrel
      tb_erel(ig) = zerel
      tb_vrel(ig) = zvrel
      ii = ii + 1
      irct(ii) = irc
      srct(ii) = deni(ic,ig)*zscxv
      elrc(ii) = 0.0d0
      derc(ii) = 0.0d0
      igrc(ii) = ig
      enddo
!
!--EL1  [ H+ + H0  =>  H+  + H0 ]
      if( rcfc(3).gt.0.0d0 ) then
      irc = 3
      do ig = 1, ngas
      zerel = tb_erel(ig)
      zvrel = tb_vrel(ig)
      zscx  = xdt_eval1(el_cros(el_atm),zerel)
      zscx  = zscx*1.0d-4    ! Note.  erel [eV/amu], scx [cm2]
      zscxv = zscx*zvrel
      ii = ii + 1
      irct(ii) = irc
      srct(ii) = deni(ic,ig)*zscxv
      elrc(ii) = 0.0d0
      derc(ii) = 0.0d0
      igrc(ii) = ig
      enddo
      endif
!
      ! ------ NNC by toku 160625 ----------------------------------------
!--EL3  [ H0 + H0  =>  H0  + H0 ]
      if( rcfc(12).gt.0.0d0 ) then
      irc = 12
      do ig = 1, ngas
      v00rel2 = (velx-nfl0BGx(ic,ig))**2
     >         +(vely-nfl0BGy(ic,ig))**2
     >         +(velz-nfl0BGz(ic,ig))**2
      v00rel2 = 8.0d0*eng0BG(ic,ig)*cev/(rmas(ig)*cpi) + v00rel2
      E00rel  = 0.5d0*cmp*v00rel2/cev
      v00rel  = dsqrt(v00rel2)
      lnE=log(E00rel)
      call sigma_DxD_t(lnE,sig00)
      sigv00  = dble(lnnel)*1.d-2*sig00*v00rel
      ii = ii + 1 ! number of collision type X combination of species
      irct(ii) = irc
      srct(ii) = den0BG(ic,ig)*sigv00
      elrc(ii) = 0.0d0
      derc(ii) = 0.0d0
      igrc(ii) = ig
      enddo
      endif
!
!--EL4  [ H0 + H2  =>  H0  + H2 ]
      if( rcfc(13).gt.0.0d0 ) then
      irc = 13
      do ig = 1, ngas
      v0grel2 = (velx-nflgBGx(ic,ig))**2
     >         +(vely-nflgBGy(ic,ig))**2
     >         +(velz-nflgBGz(ic,ig))**2
      v0grel2 = 8.0d0*enggBG(ic,ig)*cev/(2.0d0*rmas(ig)*cpi) + v0grel2
      E0grel  = 0.5d0*cmp*v0grel2/cev
      v0grel  = dsqrt(v0grel2)
      lnE=log(E0grel)
      call sigma_DxD2_t(lnE,sig0g)
      sigv0g  = dble(lnnel)*1.d-2*sig0g*v0grel ! tentative
      ii = ii + 1 ! number of collision type X combination of species
      irct(ii) = irc
      srct(ii) = dengBG(ic,ig)*sigv0g
      elrc(ii) = 0.0d0 ! ?
      derc(ii) = 0.0d0 ! ?
      igrc(ii) = ig
      enddo
      endif
! ------------------------------------------------------------------
!
      nrct = ii
!
!::probability (except for ionization)
      zsum = 0.0d0
      frct(1) = zsum
      do i = 2, nrct
      zsum = zsum + srct(i)
      frct(i) = zsum
      enddo
      trcx = zsum
      do i = 1, nrct
      frct(i) = frct(i)/zsum
      enddo
      frct(nrct) = 1.0001d0
!
!::total collision frequency
      trct = srct(1) + trcx
      wrct = trcx/trct
      dmfp = pvel/trct
      return
!
!::vacume
 100  continue
      nrct = 0;  trct = 0.0d0;  wrct = 1.0d0;  dmfp = 1.0d30
      return
!
      end
!
!**********************************************************************
      subroutine ntcrs_MH0(ic,pvel)
!**********************************************************************
      use BGProf, only : den0BG, dengBG, eng0BG, enggBG, nfl0BGx
     >    , nfl0BGy, nfl0BGz, nflgBGx, nflgBGy, nflgBGz
      use celcom, only : el_cros, el_mol
      use cntcom, only : bvx, bvy, bvz, derc, dmfp, elrc, frct, igrc
     >    , irct, lnnel, mrct, ngas, nrct, rcfc, rmas, srct, tbde_ds
     >    , tbel_ds, tbsg_ds, trct, velx, vely, velz, wrct
      use cntpls, only : dene, deni, temi, vflw
      use cphcns, only : cev, cmp, cpi
      use KrsticSchultz3, only : sigma_DxD2_t
      use Phelps, only : sigma_H2xH2_t
      implicit none
!
!::argument
      integer, intent(in) :: ic
      real(8), intent(in) :: pvel
!
!::local variables
      integer  ii, ir, irc, ig, i
      real*8   zne, zti, zn0bg, zngbg
      real*8   zvntl2, zvrel2, zerel, zvrel, zscx, zscxv, zsum
! function
      real(8)    xdt_eval1
      real(8) :: vg0rel2, vg0rel, Eg0rel, vggrel2, vggrel, Eggrel,
     >     sigg0, sigvg0, siggg, sigvgg, lnE
!
      mrct = 2
      ii = 0
!
!--plasma parameter
      zne = dene(ic)
      zti = temi(ic)
      zn0bg= den0BG(ic,1)
      zngbg= dengBG(ic,1)
      if(lnnel==1) then
        if( zne+zn0bg+zngbg.le.0.0d0 ) goto 100
      else
        if( zne.le.0.0d0 ) goto 100
      endif
!
!--cross
!xx   call sgvmhy( zte, asgv, aels, aeng )
!xx   set tbsg_ds, tbel_ds, tbde_ds sub.  ntcrds
!
!--DS1  [ e  + H2  =>  H0  + H0 + e ]
!--DS2  [ e  + H2  =>  H2+ + 2e ]
!--DS3  [ e  + H2  =>  H+  + H0 + e ]
      ir = 0
      do irc = 4, 6
      ir = ir + 1
      ii = ii + 1
      irct(ii) = irc
!x    srct(ii) = zne*asgv(ir)
!x    elrc(ii) = aels(ir)
!x    derc(ii) = aeng(ir)
      srct(ii) = tbsg_ds(ic,ir)
      elrc(ii) = tbel_ds(ic,ir)
      derc(ii) = tbde_ds(ic,ir)
      igrc(ii) = -1
      enddo
!
!--EL2  [ H+ + H2  =>  H+ + H2 ]
      if( rcfc(7).gt.0.0d0 ) then
      irc = 7
      do ig = 1, ngas
!xx   zvntl2 = pvel*pvel                      !  2003/02/13
      zvntl2 = (velx-vflw(ic,ig)*bvx(ic))**2
     >        +(vely-vflw(ic,ig)*bvy(ic))**2
     >        +(velz-vflw(ic,ig)*bvz(ic))**2
      zvrel2 = 8.0d0*zti*cev/(rmas(ig)*cpi) + zvntl2
      zerel  = 0.5d0*cmp*zvrel2/cev
      zvrel  = dsqrt(zvrel2)
      zscx  = xdt_eval1(el_cros(el_mol),zerel)
      zscx  = zscx*1.0d-4   ! Note.  erel [eV/amu], scx [cm2]
      zscxv = zscx*zvrel
      ii = ii + 1
      irct(ii) = irc
      srct(ii) = deni(ic,ig)*zscxv
      elrc(ii) = 0.0d0
      derc(ii) = 0.0d0
      igrc(ii) = ig
      enddo
      endif
!
!--EL5  [ D2 + D0  =>  D2  + D0 ]
      if( rcfc(14).gt.0.0d0 ) then
      irc = 14
      do ig = 1, ngas
      vg0rel2 = (velx-nfl0BGx(ic,ig))**2
     >         +(vely-nfl0BGy(ic,ig))**2
     >         +(velz-nfl0BGz(ic,ig))**2
      vg0rel2 = 8.0d0*eng0BG(ic,ig)*cev/(rmas(ig)*cpi) + vg0rel2 ! ???
      Eg0rel  = 0.5d0*cmp*vg0rel2/cev
      vg0rel  = dsqrt(vg0rel2)
      lnE=log(Eg0rel)
      call sigma_DxD2_t(lnE,sigg0)
      sigvg0  = dble(lnnel)*1.d-2*sigg0*vg0rel
      ii = ii + 1 ! number of collision type X combination of species
      irct(ii) = irc
      srct(ii) = den0BG(ic,ig)*sigvg0
      elrc(ii) = 0.0d0
      derc(ii) = 0.0d0
      igrc(ii) = ig
      enddo
      endif
!
!--EL6  [ H2 + H2  =>  H2  + H2 ]
      if( rcfc(15).gt.0.0d0 ) then
      irc = 15
      do ig = 1, ngas
      vggrel2 = (velx-nflgBGx(ic,ig))**2
     >         +(vely-nflgBGy(ic,ig))**2
     >         +(velz-nflgBGz(ic,ig))**2
      vggrel2 = 8.0d0*enggBG(ic,ig)*cev/(2.0d0*rmas(ig)*cpi) + vggrel2 ! ???
      Eggrel  = 0.5d0*cmp*vggrel2/cev
      vggrel  = dsqrt(vggrel2)
      lnE=log(Eggrel)
      call sigma_H2xH2_t(lnE,siggg)
!      siggg=0.d0
      sigvgg  = dble(lnnel)*1.d-2*siggg*vggrel
      ii = ii + 1 ! number of collision type X combination of species
      irct(ii) = irc
      srct(ii) = dengBG(ic,ig)*sigvgg
      elrc(ii) = 0.0d0
      derc(ii) = 0.0d0
      igrc(ii) = ig
      enddo
      endif
!
      nrct = ii
!
!::probability
      zsum = 0.0d0
      do i = 1, nrct
      zsum = zsum + srct(i)
      frct(i) = zsum
      enddo
      do i = 1, nrct
      frct(i) = frct(i)/zsum
      enddo
      frct(nrct) = 1.0001d0
!
!::total collision frequency
      trct = zsum
      wrct = 1.0d0
      dmfp = pvel/trct
      return
!
!::vacume
 100  continue
      nrct = 0;  trct = 0.0d0;  wrct = 1.0d0;  dmfp = 1.0d30
      return
      end
!
!**********************************************************************
      subroutine ntcrs_MHP(ic,pvel)
!**********************************************************************
      use cntcom, only : derc, dmfp, elrc, frct, igrc, irct, mrct, nrct
     >    , srct, tbde_ds, tbel_ds, tbsg_ds, trct, wrct
      use cntpls, only : dene
      implicit none
!
!::argument
      integer, intent(in) :: ic
      real(8), intent(in) :: pvel
!
!::local variables
      integer ii, ir, irc, i
      real*8  zne, zsum
!
      mrct = 3
      ii = 0
!
!--plasma parameter
      zne = dene(ic)
      if( zne.le.0.0d0 ) goto 100
!
!--cross
!--DS4  [ e  + H2+ =>  H+  + H0 + 2e ]
!--DS5  [ e  + H2+ =>  H+  + H+ + e ]
!--DS6  [ e  + H2+ =>  H0  + H0 + e" ]
      ir = 3
      do irc = 8, 10
      ir = ir + 1
      ii = ii + 1
      irct(ii) = irc
!x    srct(ii) = zne*asgv(ir)
!x    elrc(ii) = aels(ir)
!x    derc(ii) = aeng(ir)
      srct(ii) = tbsg_ds(ic,ir)
      elrc(ii) = tbel_ds(ic,ir)
      derc(ii) = tbde_ds(ic,ir)
      igrc(ii) = -1
      enddo
!
      nrct = ii
!
!::probability
      zsum = 0.0d0
      do i = 1, nrct
      zsum = zsum + srct(i)
      frct(i) = zsum
      enddo
      do i = 1, nrct
      frct(i) = frct(i)/zsum
      enddo
      frct(nrct) = 1.0001d0
!
!::total collision frequency
      trct = zsum
      wrct = 1.0d0
      dmfp = pvel/trct
      return
!
!::vacume
 100  continue
      nrct = 0;  trct = 0.0d0;  wrct = 1.0d0;  dmfp = 1.0d30
      return
!
      end
!
!
!**********************************************************************
      subroutine ntcrtb
!**********************************************************************
!
!      cross-section of el-ionization
!                       el-dissociation
!
!----------------------------------------------------------------------
      use cntcom, only : ncmax, tbde_ds, tbel_ds, tbel_ei, tbsg_ds
     >    , tbsg_ei, temin_ion
      use cntpls, only : dene, teme
      use csize,  only : ndmc
      implicit none
!
!::local variables
      integer  ir, nsiz, ic
      real*8   zne, zte
      real*8   asgv(6),aels(6),aeng(6)
! function
      real(8)    elhion, svhion
!
!::clear
      nsiz = ndmc+1
      call setd( tbsg_ei, nsiz,  0.0d0 )
      call setd( tbel_ei, nsiz,  0.0d0 )
      nsiz = (ndmc+1)*6
      call setd( tbsg_ds, nsiz,  0.0d0 )
      call setd( tbel_ds, nsiz,  0.0d0 )
      call setd( tbde_ds, nsiz,  0.0d0 )
!
!::ne*<sig*v>,Eion
      do ic = 1, ncmax
        if( teme(ic).gt.0.0d0 .and. dene(ic).gt.0.0d0) then
          zne = dene(ic)
          zte = dmax1(teme(ic),temin_ion)   !---today=09.24
          tbsg_ei(ic) = dmax1(zne*svhion(zte,zne),1.0d-99)
          tbel_ei(ic) = elhion(zte,zne)
          call sgvmhy( zte, zne, asgv, aels, aeng )
          do ir = 1, 6
            tbsg_ds(ic,ir) = dmax1(zne*asgv(ir),1.0d-99)
            tbel_ds(ic,ir) = aels(ir)
            tbde_ds(ic,ir) = aeng(ir)
          enddo
        endif
      enddo
!
      return
      end
