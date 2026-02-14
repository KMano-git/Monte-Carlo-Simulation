!***********************************************************************
      subroutine imclear(ctyp)
!***********************************************************************
!
!    ir : cell number   is : charge state   il : cell surface
!    ivs: shallow       ivd: deep           itg: ???
!    igs: sputtering     igi: first ionization
!    iml: CD4           ien: condition
!    igtno: gate No     igtem: passing gate
!
!    tt : time          rr : position (R)   zz : position (Z)
!    fi : angle         vr : velocity (Vv)  vz : velocity (Vb)
!    vv : energy (V2)  wght: weight        wspt: weight at stick/error
!    vvr: Vr**2        vvz : Vz**2          v  : velocity
!
!       ctyp = "particle"  "score"  "floss"
!-----------------------------------------------------------------------
      use cimcom, only : denz, fforc, fi, friz, hcnpt, hcnwt, hcpum
     >    , hgrpn, hgrpp, hgrpw, hmype, hpmax, hregn, hregp, hregw
     >    , hrgpa, hrgpc, hstep, htout, ien, igi, igs, igtem, igtno
     >    , il, iml, ionz, ir, is, itg, ivd, ivs, mypcn, myptl, mywgt
     >    , ndext, ndmp, ntags, recz, rr, sdmy1, sdmy2, sdmy3, sdtmz
     >    , sflux, sitmp, sptcl, stimp, stimz, swtot, temz, tforc, thfz
     >    , tt, v, vr, vv, vvr, vvz, vz, vzpz, wexh, wght, whit, wpemt
     >    , wsct, wstp, wtemt, wvaci, wvacn, zz
! added hgrpn, hgrpp, hgrpw, hregn, hregp, hregw, mypcn, myptl, mywgt, ntags, wstp from IMPMC_TD by kamata 2022/05/14
      use cimctl, only : icalz
      use cimfls, only : fls_flxi, fls_flxo
      use cimntl, only : ndmp2, sint, srn, stb, stm, svx, svy, svz, sx0
     >    , sy0, sz0
      use csize,  only : ndgt, ndmc, ndmis
      use csonic, only : itim
      use cunit,  only : n6
      use mod_shexe, only : impmc_model
      implicit none

!::argument
      character, intent(in) :: ctyp*(*)

      if( ctyp.eq."score" ) then
      write(n6,'(2x,"*** imclear: score ***  icalZ/itim =",i5,i7)')
     >   icalZ, itim

!::/cimcom_11/
      if( impmc_model == 0 ) then
      sflux = 0.0d0
      swtot = 0.0d0
      endif
      sptcl = 0.0d0
      stimp = 0.0d0
      sitmp = 0.0d0
      stimz = 0.0d0
      sdtmz = 0.0d0
      sdmy1 = 0.0d0
      sdmy2 = 0.0d0
      sdmy3 = 0.0d0
      denZ(0:ndmis,0:ndmc) = 0.0d0
      temZ(0:ndmis,0:ndmc) = 0.0d0
      wsct(0:ndmc) = 0.0d0
      friZ(0:ndmis,0:ndmc) = 0.0d0
      thfZ(0:ndmis,0:ndmc) = 0.0d0
      vzpZ(0:ndmis,0:ndmc) = 0.0d0
      ionZ(0:ndmis,0:ndmc) = 0.0d0
      recZ(0:ndmis,0:ndmc) = 0.0d0

!::/cimcpm_12/
      wtemt = 0.0d0
      wpemt = 0.0d0
      whit(0:ndmis,1:30,1:4) = 0.0d0  ! 4 : inj/rfl/pas/abs
      wexh(0:ndmis,1:30) = 0.0d0
      myptl(1:15) = 0.0d0
      mywgt(1:15) = 0.0d0
      mypcn(1:15) = 0.0d0

!::/cimcom_12b/
      hmype = 0.0d0
      htout = 0.0d0
      hstep = 0.0d0
      hcpum = 0.0d0
      hpmax = 0.0d0
      hrgpa(0:10) = 0.0d0
      hrgpc(0:10) = 0.0d0
      hcnpt(0:10) = 0.0d0
      hcnwt(0:10) = 0.0d0
      hgrpp(1:15) = 0.0d0
      hgrpw(1:15) = 0.0d0
      hgrpn(1:15) = 0.0d0
      hregp(1:10) = 0.0d0
      hregw(1:10) = 0.0d0
      hregn(1:10) = 0.0d0

!::/cimcom_12c/ :: flow under dome
      wvacn(0:ndgt,0:ndext) = 0.0d0
      wvaci(0:ndgt,0:ndext) = 0.0d0
      endif

!::particle variables
      if( ctyp.eq."particle") then
      write(n6,'(2x,"*** imclear: partc ***  icalZ/itim =",i5,i7)')
     >   icalZ, itim
      ir   (1:ndmp)  = 0
      is   (1:ndmp)  = 0
      il   (1:ndmp)  = 0
      ivs  (1:ndmp)  = 0
      ivd  (1:ndmp)  = 0
      itg  (1:ndmp)  = 0
      igs  (1:ndmp)  = 0
      igi  (1:ndmp)  = 0
      iml  (1:ndmp)  = 0
      ien  (1:ndmp)  = 0
      igtno(1:ndmp)  = 0
      igtem(1:ndmp)  = -9
      ntags(1:ndmp)  = 0

      tt   (1:ndmp)  = 0.0d0
      rr   (1:ndmp)  = 0.0d0
      zz   (1:ndmp)  = 0.0d0
      fi   (1:ndmp)  = 0.0d0
      vr   (1:ndmp)  = 0.0d0
      vz   (1:ndmp)  = 0.0d0
      vv   (1:ndmp)  = 0.0d0
      wght (1:ndmp)  = 0.0d0
      wstp (1:ndmp)  = 0.0d0
      vvr  (1:ndmp)  = 0.0d0
      vvz  (1:ndmp)  = 0.0d0
      v    (1:ndmp)  = 0.0d0
      fforc(1:ndmp)  = 0.0d0
      tforc(1:ndmp)  = 0.0d0

      sx0 (1:ndmp2) = 0.0d0
      sy0 (1:ndmp2) = 0.0d0
      sz0 (1:ndmp2) = 0.0d0
      svx (1:ndmp2) = 0.0d0
      svy (1:ndmp2) = 0.0d0
      svz (1:ndmp2) = 0.0d0
      srn (1:ndmp2) = 0.0d0
      sint(1:ndmp2) = 0.0d0
      stm (1:ndmp2) = 0.0d0
      stb (1:ndmp2) = 0.0d0
      endif

!::flux at core edge
      if( ctyp.eq."floss" ) then
      write(n6,'(2x,"*** imclear: floss ***  icalZ/itim =",i5,i7)')
     >   icalZ, itim
      fls_flxI(0:ndmis,1:5) = 0.0d0  ! flux at core edge
      fls_flxO(0:ndmis,1:5) = 0.0d0  ! flux at core edge
!      fls_recI(1:ndmp,1:5) = 0     ! recode of flux to core
!      fls_recO(1:ndmp,1:5) = 0     ! recode of flux from core
      endif

      return
      end

!***********************************************************************
      subroutine chk_wpemt(cmsg)
!***********************************************************************
      use cimcom, only : hcnpt, idemt, wpemt
      use cunit,  only : n6
      implicit none

      character, intent(in) :: cmsg*(*)

      write(n6,'(2x,"===",a,"===   wpemt,stpt =",1p2e12.3)')
     >  trim(cmsg), wpemt, hcnpt(idemt)

      return
      end
