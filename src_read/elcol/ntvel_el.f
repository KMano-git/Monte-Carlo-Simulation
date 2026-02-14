!**********************************************************************
      subroutine ntvel_el(irc,ipl,xpos1,ypos1,wpox1,wpoy1)
!**********************************************************************
!
!     new velocity after elastic collision
!
!       irc :  index of reaction    see sub. ntcros
!                 = 3  (D0 + D+)  = 7  (D2 + D+)
!                 =12  (D0 + D0)  =13  (D0 + D2)
!                 =14  (D2 + D0)  =15  (D2 + D2)
!       ipl :  plasma ion
!
!
!   Procedure to calculate the velocity after collision
!
!        A : the original co-rodinate
!        B : the frame moving with the flow velocity
!
!                          test      ion      mass of center
!-----------------------------------------------------------
!   before col.  in <A>    vtx0      vix0
!                in <B>    vtx       vix      vrx
!-----------------------------------------------------------
!   after  col.  in <B>    wtx       wix      wrx
!                in <A>    wtx0      wix0
!
!-----------------------------------------------------------------------
      use BGProf, only : vlp0BG, vlpgBG, eng0BG, enggBG, nfl0BGx
     >    , nfl0BGy, nfl0BGz, nflgBGx, nflgBGy, nflgBGz
      use celcom, only : agmi, dcl_eng, dcl_ipl, dcl_mvp, dcl_mvx
     >    , dcl_mvy, dcl_mvz, el_angl, el_cros, pcl_eng, pcl_mvp
     >    , pcl_mvx, pcl_mvy, pcl_mvz, pcl_sgv, sgmx
      use cntcom, only : bvx, bvy, bvz, evel, fnor2, icpo, rmas, tmas
     >    , tnor2, vel, vel2, velp, velx, vely, velz, vion
      use cntpls, only : vflw
      use cphcns, only : cev, cmp, cpi
      use cunit,  only : mype, n6
      use KrsticSchultz3, only : sigma_DxD2_t, sigma_DxD_t
      use Phelps, only : sigma_H2xH2_t
      implicit none
!
!::arguments
      integer, intent(in) :: irc, ipl
      real(8), intent(in) :: xpos1, ypos1, wpox1, wpoy1 ! dummy, for nnc examination
!
!::local variables
      real*8  vfx0,vfy0,vfz0
      real*8  vtx0,vty0,vtz0, vix0,viy0,viz0
      real*8  vtx, vty, vtz,  vix, viy, viz,  vrx,vry,vrz
      real*8  wtx, wty, wtz,  wix, wiy, wiz,  wrx,wry,wrz
      real*8  wtx0,wty0,wtz0, wix0,wiy0,wiz0
      real*8  vtp0,vgx,vgy,vgz,vrv2,vrv,vr2,vr
      real*8  ami,zvi,fvt,fvi,fwt,fwi,fai,sth,cth,sfi,cfi
      real*8  amn
      integer ivx,ivy,ivz, kel, ixs_angl, ixs_cros
      real*8  sgvmax, angmin, sgma, sgvr, xran, prb
      real*8  zvr,zwr,zct
      real*8  zimvx,zimvy,zimvz,zimvp,zieng
      integer :: ierr = 0
      save ierr
! function
      integer    imox, imoy
      real(8)    random, xdt_eval1, xdt_eval2
!
!::temp common for test program
      real*8  erl, theta, lnErel 
      integer lp; data lp/0/
      integer lrj, lrjmax
      character*18 :: cErmsg
!
!----------------------------------------------------------------------
!
      lrj    = 0
      lrjmax = 1000
!
!::set ixs_angl
      if( irc.eq.3 ) then
        kel = 1
      elseif( irc.eq.7 ) then
        kel = 2
      elseif( irc.eq.12) then
        kel = 3 ! D+D
        amn=rmas(ipl)
      elseif( irc.eq.14) then
        kel = 4 ! D2+D  !!! Sorted by Background Species!!
        amn=rmas(ipl)
      elseif( irc.eq.13) then
        kel = 5 ! D+D2
        amn=2.*rmas(ipl)
      elseif( irc.eq.15) then
        kel = 6 ! D2+D2
        amn=2.*rmas(ipl)
      else
        write(cErmsg,'(a14,i2)') 'incorrect irc=',irc
        call wexit("ntvel_el",cErmsg)
      endif
!
!::test particle
      vtx0 = velx
      vty0 = vely
      vtz0 = velz
      vtp0 = velp

!::set ixs, sigv_max, angle_min
      if(kel<3)then
        ixs_angl = el_angl(kel)
        ixs_cros = el_cros(kel)
        sgvmax   = sgmx(kel)
        angmin   = agmi(kel)
!
!::plasma flow
        vfx0 = vflw(icpo,ipl)*bvx(icpo)
        vfy0 = vflw(icpo,ipl)*bvy(icpo)
        vfz0 = vflw(icpo,ipl)*bvz(icpo)
        ami  = rmas(ipl)
        zvi  = vion(icpo,ipl)
!
!----------------------------------------
      else ! NNC by toku
        ixs_angl = el_angl(1) ! tentative
        ixs_cros = el_cros(1)  ! tentative
        sgvmax   = 1.d-1*sgmx(1) ! tentative
        angmin   = agmi(1) ! tentative
!
!::Neutral BG flow
        if(kel<5)then
          vfx0 = nfl0BGx(icpo,ipl)
          vfy0 = nfl0BGy(icpo,ipl)
          vfz0 = nfl0BGz(icpo,ipl)
          ami  = rmas(ipl)    ! D
          zvi  = sqrt(eng0BG(icpo,ipl)*cev/ami)
        elseif(kel<7)then
          vfx0 = nflgBGx(icpo,ipl)
          vfy0 = nflgBGy(icpo,ipl)
          vfz0 = nflgBGz(icpo,ipl)
          ami  = rmas(ipl)*2.d0 ! D2
          zvi  = sqrt(enggBG(icpo,ipl)*cev/(2.*ami))
        endif
!----------------------------------------
      endif
!
!--frame moving with the flow velocity
!
!::test velocity
      vtx = vtx0 - vfx0
      vty = vty0 - vfy0
      vtz = vtz0 - vfz0
!
!::ion velocity   ! rare case  (vty = viy, vtz = viz)
 100  continue
      ivx = int(fnor2*random(0) + 1.0d0)
      ivy = int(fnor2*random(0) + 1.0d0)
      ivz = int(fnor2*random(0) + 1.0d0)
      vix = zvi*tnor2(ivx)  ;  vix0 = vix + vfx0
      viy = zvi*tnor2(ivy)  ;  viy0 = viy + vfy0
      viz = zvi*tnor2(ivz)  ;  viz0 = viz + vfz0
!
!::center of mass velocity
      fvt = tmas/(tmas+ami)
      fvi = ami/(tmas+ami)
      vgx = fvt*vtx+fvi*vix
      vgy = fvt*vty+fvi*viy
      vgz = fvt*vtz+fvi*viz
!
!::relative velocity
      vrx = vtx - vix
      vry = vty - viy
      vrz = vtz - viz
      vrv2= vry**2+vrz**2
!-----
      if( vrv2.le.1.0d-30 ) then
        ierr = ierr + 1
        if( ierr.gt.100 ) call wexit("ntvel_el","vrv2~0")
        write(n6,'(2x,"Warning at ntvel_el vrv2~0 ",1pe12.4,
     >  "  vtx =",1p3e12.4,"  vix =",1p3e12.4,
     >  "  vtx0 =",1p3e12.4,"  vfx0 =",1p3e12.4,"  zvi,icpo,ipl =",
     >  1p2e12.4,i7,i3,2i5)') vrv2,
     >   vtx,vty,vtz, vix,viy,viz,  vtx0,vty0,vtz0, vfx0,vfy0,vfz0,
     >   zvi,vflw(icpo,ipl),icpo,ipl, imox(icpo),imoy(icpo)
        vrv2= 1.0d-30
      endif
!-----
      vrv = dsqrt(vrv2)
      vr2 = vrx**2+vrv2
      vr  = dsqrt(vr2)
      erl = 0.5d0*cmp*vr2/cev
!
!::sgvr(vrel)
      if(kel<3)then
        sgma = xdt_eval1( ixs_cros, erl )
        sgma = sgma*1.0d-4  !  cgs => MKS unit
!
! ---- NNC -----------------------------------------------
      else
!         erl = erl*amn/(amn+amt) ! E_lab -> E_CoM
         lnErel=log(erl)
         if(kel==3)then
            call sigma_DxD_t(lnErel,sgma)
         elseif(kel==4)then
            call sigma_DxD2_t(lnErel,sgma)
         elseif(kel==5)then
            call sigma_DxD2_t(lnErel,sgma)
         elseif(kel==6)then
            lnErel=log10(erl)
            call sigma_H2xH2_t(lnErel,sgma)
         endif
      endif
! --------------------------------------------------------
      sgvr = sgma*vr
!
!::rejection technique
!-----
      lrj  = lrj + 1
      if( lrj.gt.lrjmax ) then
        write(n6,'(2x,"--- Warning ntvel_el  lrj.gt.lrjmax ",2i6)')
     >    lrj,lrjmax
        goto 150
      endif
!-----
      xran = random(0)
      if( sgvr.lt.xran*sgvmax ) goto 100
 150  continue
!================================================================================
!
!::scattering angle & polar angle
      prb = random(0)
      if(kel<3)then
        theta = xdt_eval2(ixs_angl,prb,erl) ! Note unit of erl [eV/amu]
      else
! ---- NNC (tentative) -------------------------------------------------
        theta = xdt_eval2(ixs_angl,prb,erl) ! Note unit of erl [eV/amu]
! ----------------------------------------------------------------------
      endif
      if( theta.lt.angmin ) theta = angmin
      fai = 2.0d0*cpi*random(0)
!
!::trigonometric function
      sth = sin(theta)
      cth = cos(theta)
      sfi = sin(fai)
      cfi = cos(fai)
!
!::relative velocity after collision
      wrx = vrx*cth-sth*sfi*vrv
      wry = vry*cth+sth*(vrx*vry*sfi-vr*vrz*cfi)/vrv
      wrz = vrz*cth+sth*(vrx*vrz*sfi+vr*vry*cfi)/vrv
!
!::velocity
      fwt = ami/(tmas+ami) ;  fwi = tmas/(tmas+ami)
      wtx = vgx + fwt*wrx  ;  wix = vgx - fwi*wrx
      wty = vgy + fwt*wry  ;  wiy = vgy - fwi*wry
      wtz = vgz + fwt*wrz  ;  wiz = vgz - fwi*wrz
!
!--the original co-ordinate
!
!::add flow velocity
      wtx0 = wtx + vfx0    ;  wix0 = wix + vfx0
      wty0 = wty + vfy0    ;  wiy0 = wiy + vfy0
      wtz0 = wtz + vfz0    ;  wiz0 = wiz + vfz0
!
!::set velocity
      velx = wtx0
      vely = wty0
      velz = wtz0
      velp = velx*bvx(icpo)+vely*bvy(icpo)+velz*bvz(icpo)
      vel2 = velx*velx + vely*vely + velz*velz
      evel = 0.5d0*tmas*vel2/cev
      vel  = dsqrt(vel2)
!
!::   collision estimator (momentum & energy change of test particle)
      dcl_ipl = ipl
!-----
      pcl_mvx = tmas*(wtx0-vtx0)
      pcl_mvy = tmas*(wty0-vty0)
      pcl_mvz = tmas*(wtz0-vtz0)
      pcl_mvp = tmas*(velp-vtp0)
      pcl_eng = (evel*cev - 0.5d0*tmas*(vtx0**2+vty0**2+vtz0**2))
!-----
      dcl_mvx = pcl_mvx*sgvr
      dcl_mvy = pcl_mvy*sgvr
      dcl_mvz = pcl_mvz*sgvr
      dcl_mvp = pcl_mvp*sgvr
      dcl_eng = pcl_eng*sgvr
!-----
!
!::debug check
      pcl_sgv = sgvr
      return
!
!-----------------------------------------------------------------------
!::debug write  <relative veocity>
!-----------------------------------------------------------------------
!       if( mype.ne.0 ) return
! !
!       if( lp.eq.1 ) then
!       write(36,'(/2x,"*** ntvel_el ***  check  |vr|=|wr|, cos(theta)")')
!       write(36,'(4x,"lp",4x,"erl",9x,"th",10x,"|vr|",8x,"|vp|",8x,
!      >    "ang_v",7x,"th",10x,"er_|v|",6x,"er_cth")')
!       endif
!       zvr = dsqrt(vrx*vrx+vry*vry+vrz*vrz)
!       zwr = dsqrt(wrx*wrx+wry*wry+wrz*wrz)
!       zct = (vrx*wrx+vry*wry+vrz*wrz)/(zvr*zwr)
!       write(36,'(2x,i6,1p8e12.3)')
!      > lp,erl,theta,zvr,zwr,zct,cth,zvr-zwr,zct-cth
! !
! !-----------------------------------------------------------------------
! !::debug write <conservation>
! !-----------------------------------------------------------------------
!       zimvx = ami*(wix0-vix0)*sgvr
!       zimvy = ami*(wiy0-viy0)*sgvr
!       zimvz = ami*(wiz0-viz0)*sgvr
!       zimvp = zimvx*bvx(icpo)+zimvy*bvy(icpo)+zimvz*bvz(icpo)
!       zieng =(0.5d0*ami*(wix0**2+wiy0**2+wiz0**2)
!      >       -0.5d0*ami*(vix0**2+viy0**2+viz0**2))*sgvr
!       if( lp.eq.1 ) then
!       write(37,'(/2x,"*** ntvel_el ***  check  conservation")')
!       write(37,'(4x,"lp",4x,"tmvx",8x,"dmvx",8x,"tmvy",8x,"dmvy",8x
!      >  ,"tmvz",8x,"dmvz",8x,"tmvp",8x,"dmvp",8x,"teng",8x,"deng")')
!       endif
!       write(37,'(2x,i6,1p10e12.2)') lp
!      > ,dcl_mvx,(dcl_mvx+zimvx)/dcl_mvx, dcl_mvy,(dcl_mvy+zimvy)/dcl_mvy
!      > ,dcl_mvz,(dcl_mvz+zimvz)/dcl_mvz, dcl_mvp,(dcl_mvp+zimvp)/dcl_mvp
!      > ,dcl_eng,(dcl_eng+zieng)/dcl_eng
! !
!       return
      end
