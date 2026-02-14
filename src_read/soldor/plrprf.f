!**********************************************************************
! modified 1/1 lines integrated calculation with TOPICS by kamata 2020/11/16
!ik   subroutine plrprf(ldummy)
      subroutine plrprf_t(ldummy)
!**********************************************************************
!
!     pointer of ic
!         icaxs : plasma axis
!         icmpe : plasma edge
!         icmps : plasma surface
!         icspx : flux tube close to separatrix
!         icwl1 : sol wall
!
!    radial profile in sol
!     fac(ic) = (dmid(ic)-dmid(icspx))/(dmid(icwl1)-dmid(icspx))
!      na(ic) = nasx - (nasx-naws)*fac(ic)
!      Ti(ic) = Tisx - (Tisx-Tiws)*fac(ic)
!      Te(ic) = Tesx - (Tesx-Tews)*fac(ic)
!
!       where Tisx = Ti(icspx), Tiws = Ti(icwl1)
!
!       fprna ==> fcna    2005/7/10
!
!     1D radial profile
!                         main  anamp, atimp, atemp
!                         SOL   anasl, atisl, atesl
!                         prv   anapv, atipv, atepv
!     1D radial diffusion
!                         main  adamp, aetmp, aximp, axemp
!                         SOL   adasl, aetsl, axisl, axesl
!                         prv   adapv, aetpv, axepv, axipv
!     flux at core edge   flxni, flxqi, flxqe
!     total rad. power    trad_min, trad_max
!
!-----------------------------------------------------------------------
      use cplcom, only : anamp, anapv, anasl, anemp, anepv, anesl, animp
     >    , anipv, anisl, atemp, atepv, atesl, atimp, atipv, atisl, cfti
     >    , fcna, mdl_wrd, nfti, nion, qtim, trad_max, trad_min
      use cplmet, only : icaxs, icmpe, icmps, icspx, icwl1, icwl2
      use cpmpls, only : clsav, clsfc, clsni, clste, clsti, clsxp, dn0c
     >    , dn0s, fflna, flxni, flxqe, flxqi, fprna, lnedg, prfni, prfte
     >    , prfti, pufs, rmdc, rmds, taup, wlose, wlosi, ltedg
     >    , nty_ctl, rad_const_step, rad_const_fac
      use csize,  only : ndsp, ndy
      use csonic, only : time
      use cunit,  only : cdgr, mygrp, n5, n6
      use topics_mod, only : dtcal, dtduc, gtim
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  ldummy
      integer, intent(in) :: ldummy ! dummy
!
!::Use
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  ia, iy, icprv, kcon
      integer  ia, iy, icprv
      real*8   facspx, facwl1, facprv, facwl2
      real*8   x, fac
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   real*8   exd_tim0
!
!::NO use local variables
      real*8  bnip, btip, btep, bnix, btix, btex
      real*8  dtstep
      integer ktmax, ntrmax
!
!::input data
      namelist /umpinp/
     >  flxni, flxqi, flxqe, fflna,
     >  trad_min, trad_max,
     >  lnedg, ltedg, nfti, cfti,
!--------------------------------------------------------------
     >  facspx, facwl1, facprv, facwl2,   ! NEW
!--------------------------------------------------------------
     >  prfni, prfti, prfte, fprna,
     >  bnip, btip, btep, bnix, btix, btex,
     >  ktmax, ntrmax, dtstep, rmdc, dn0c, taup, wlosi, wlose,
     >  rmds, dn0s, pufs,
     >  clsav, clsxp, clsfc, clsni, clsti, clste
!::radiation power constraint
     >  , nty_ctl, rad_const_step, rad_const_fac
!
      if( mygrp.ne.cdgr(1) ) return
!
! modified 1/1 lines integrated calculation with TOPICS by kamata 2020/11/16
!ik   write(n6,'(/2x,"***  plrprf (2010/08/29) ***")')
      write(n6,'(/2x,"***  plrprf_t ***")')
!
!-----------------------------------------------------------------------
!::input data
!-----------------------------------------------------------------------
      facspx = 0.95d0
      facwl1 = 0.50d0
      facprv = 0.90d0
      facwl2 = 0.40d0
!
      read(n5,umpinp)
!xx   write(n6,umpinp)
!
!::lscal (=0:read file/=1:calculation)
!::lscal ==> dummy
!
!::check input data
      if( nion.eq.1 ) fflna(1) = 1.0d0
!
!::clear
      anamp(1:ndy,1:ndsp) = 0.0d0
      animp(1:ndy) = 0.0d0
      anemp(1:ndy) = 0.0d0
      atimp(1:ndy) = 0.0d0
      atemp(1:ndy) = 0.0d0
!
      anasl(1:ndy,1:ndsp) = 0.0d0
      anisl(1:ndy) = 0.0d0
      anesl(1:ndy) = 0.0d0
      atisl(1:ndy) = 0.0d0
      atesl(1:ndy) = 0.0d0
!
      anapv(1:ndy,1:ndsp) = 0.0d0
      anipv(1:ndy) = 0.0d0
      anepv(1:ndy) = 0.0d0
      atipv(1:ndy) = 0.0d0
      atepv(1:ndy) = 0.0d0

!
!::Main plasma
! modified 3/5 lines integrated calculation with TOPICS by kamata 2020/11/16
!ik   call exd_tokrcv(kcon)
!ik   write(n6,'(2x,"TOP  gtim =",1pe14.6,"  gstp =",1pe14.6)')
!ik  >   gtim, gstp
      call rcvftop( 2 )
      call set_top_sr( 3 )
      dtduc = dtcal
      write(n6,'(2x,"TOP  gtim =",es14.6,"  dtcal =",es14.6)')
     >   gtim, dtcal
!
! deleted 6 lines integrated calculation with TOPICS by kamata 2020/11/16
!ik   if( gstp.gt.0.0d0 ) then
!ik   exd_time = gtim
!ik   exd_dtim = gstp
!ik   time = exd_time
!ik   exd_tim0 = exd_time
!ik   exd_time = exd_time + exd_dtim
!
! deleted 4 lines integrated calculation with TOPICS by kamata 2020/11/16
!ik   write(n6,'(2x,"TOP  time = ",1pe14.6,"  from topics")') time
!ik   write(n6,'(2x,"TOP  time,exd_tim0 =",1p2e14.6,"  exd_dtim =",
!ik  >  1pe14.6)') exd_time, exd_tim0, exd_dtim
!ik   endif
!
!::plot
! modified 1/1 lines integrated calculation with TOPICS by kamata 2020/11/16
!ik   if( gstp.eq.0.0d0 ) then
      if( dtcal == 0.0_8 ) then
      time = qtim ! pldisk
      write(n6,'(2x,"No data of exd_tokrcv for plotf-process")')
      write(n6,'(2x,"time = exd_time  ",1pe14.6)') time
      write(n6,'(2x,"set value of anamp,animp,anemp,atimp,atemp")')
      do ia = 1, nion
      anamp(1:ndy,ia) = 1.0d19
      enddo
      animp(1:ndy) = 1.0d19
      anemp(1:ndy) = 1.0d19
      atimp(1:ndy) = 180.0d0
      atemp(1:ndy) = 80.0d0
      endif
!
!::SOL plasma
      do iy = icwl1, icspx
        x = dfloat(iy-icwl1)/dfloat(icspx-icwl1)
        fac = facwl1 + x*(facspx-facwl1)
        anisl(iy) = animp(icmps)*fac
        anesl(iy) = anemp(icmps)*fac
        atisl(iy) = atimp(icmps)*fac
        atesl(iy) = atemp(icmps)*fac
        do ia = 1, nion
        anasl(iy,ia) = anisl(iy)*fcna(ia)
        enddo
      enddo
!
!::PRV plasma
      icprv = icspx + 1
      do iy = icprv, icwl2
      x = dfloat(iy-icwl2)/dfloat(icprv-icwl2)
      fac = facwl2 + x*(facprv-facwl2)
        anipv(iy) = animp(icmps)*fac
        anepv(iy) = anemp(icmps)*fac
        atipv(iy) = atimp(icmps)*fac
        atepv(iy) = atemp(icmps)*fac
        do ia = 1, nion
        anapv(iy,ia) = anipv(iy)*fcna(ia)
        enddo
      enddo
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
!
!::flux from the main plasma
      write(n6,'(5x,"lnedg =",i2)') lnedg
      write(n6,'(5x,"fflna =",10f7.3)') (fflna(ia),ia=1,nion)
      write(n6,'(5x,"flxni =",1pe11.3,"  flxqi =",1pe11.3,"  flxqe =",
     >  1pe11.3)') flxni,flxqi,flxqe
!
!::radiation
      write(n6,'(5x,"mdl_wrd =",i2,"  trad_min =",1pe12.3,
     > "  trad_max =",1pe12.3)') mdl_wrd, trad_min, trad_max
!
!::radial profile
      write(n6,'(2x)')
      write(n6,'(2x,"icwl1 =",i3,f8.4,"  icspx =",i3,f8.4)')
     >   icwl1, facwl1, icspx, facspx
      write(n6,'(2x,"icprv =",i3,f8.4,"  icwl2 =",i3,f8.4)')
     >   icprv, facprv, icwl2, facwl2
      write(n6,'(2x,"icmps =",i3,"  icmpe =",i3,"  icaxs =",i3)')
     >   icmps, icmpe, icaxs
!
      write(n6,'(3x,"iy",3x,"anasl",7x,"atisl",7x,"atesl",7x,
     >   "anapv",7x,"atipv",7x,"atepv",7x,"anamp",7x,"atimp",7x
     >   "atemp")')
      ia = 1
      do iy = icwl1, icaxs
      if( iy.le.icmpe .or. mod(iy,10).eq.0 .or. iy.eq.icaxs ) then
      write(n6,'(2x,i3,1p9e12.3)')
     >  iy, anasl(iy,ia), atisl(iy), atesl(iy),
     >      anapv(iy,ia), atipv(iy), atepv(iy),
     >      anamp(iy,ia), atimp(iy), atemp(iy)
      endif
      enddo
!
      return
      end
