!***********************************************************************
      subroutine plsrcn_man
!***********************************************************************
      use cntmnt, only : sn0, ssn, swe, swi
      use cphcns, only : cev
      use cplcom, only : ama, nion, vcs, vne, vni, vte, vti
      use cplmet, only : hvol, icaxs, icmpe, icmps, icspx, icwl1
      use cpmpls, only : alsfr, an0mp, arhmp, armp, asnmp, clsav, clsfc
     >    , clsni,  clste, clsti, clsxp, dvmp, jmd1, mdp_ni, mdp_te
     >    , sfmp, taup, vlmp, wlmp, wlose, wlosi
      use csize,  only : ndy
      use cunit,  only : n6
      implicit none
!
!::local variable
      integer ix, iy, ia
      real*8  we0, wn0, wv0, ws0, wf0, werad, wirad
      real*8  zwi, zwe, gsi, znrm, znrms, tni
      real*8  zfa, zfb, zdl, zne, zte, zal, zal2, zsgv, zrmd, zsum, zvl
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8  zsi, zn0, zsg, tsi, tfd, zln, zvp, ztp, zfr, zti, avni
!ik   real*8  tno, tvl, twi, twe, tsi2, twi2, twe2, zeizi, zeize, ztau
      real*8  zsi, zn0, tsi, tfd, zln, zvp, ztp, zfr, avni
      real*8  tvl, twi, twe, tsi2, twi2, twe2, zeizi, zeize, ztau
      real*8  xn0(ndy), xsi(ndy)
      real*8  log_huge
!
      if( nion.gt.1 ) then
      call wexit("plsrn_man","nion > 1")
      endif

      log_huge = log(huge(zsum)) - 10.0
!
!
!----------------------------------------------------------------------
!::input parameter
!----------------------------------------------------------------------
      write(n6,'(/2x,"*** plsrcn_man ***")')
      we0 = 30.0d0  ! including charge exchange
      wn0 = 1.0d16  ! temprary
      wv0 = sqrt(we0*cev/ama(1))
      ws0 = sfmp(icspx)
      wf0 = 1.0d0/4.0d0*ws0*wn0*wv0
      wirad = 10.0d0
      werad = 40.0d0
      write(n6,'(2x,"e0 =",1pe11.3,"  n0 =",1pe11.3,"  sf =",1pe11.3,
     >  "  flx0 =",1pe11.3,"  wirad =",1pe10.3,"  werad =",1pe10.3)')
     >   we0, wn0, ws0, wf0, wirad, werad
!
!----------------------------------------------------------------------
!::(N0,Si) in the main plasma
!----------------------------------------------------------------------
      write(n6,'(2x,2x,"iy",3x,"Ne",10x,"Te",10x,"Fa",10x,"Fb",10x,
     >   "dl",10x,"dvl",9x,"N0",10x,"Si",10x,"tSi")')
      ix = jmd1
      zfb = wf0
      zsum = 0.0d0
      tsi = 0.0d0
      tni = 0.0d0
      tvl = 0.0d0
      tsi2 = 0.0d0
      do iy = icmps, icaxs
      zdl = armp(icspx) - arhmp(iy)
      zne = mdp_ni(iy)
      zte = mdp_te(iy)
      zal = zte/10.0d0
      zal2 = zal**2
      zsgv = 3.0d-14*zal2/(3.0d0+zal2)
      zrmd = wv0/(zne*zsgv)
      zsum = zsum + zdl/zrmd
      zfa  = zfb
      zvl  = dvmp(iy)
      if(zsum < log_huge) then
        zn0  = wn0*exp(-zsum)
      else
        zn0 = 0.0_8
      endif
      zsi  = zne*zn0*zsgv
      zfb  = 1.0d0/4.0d0*sfmp(iy)*zn0*wv0
      tvl  = tvl + zvl
      tsi  = tsi + zsi*zvl
      tni  = tni + zne*zvl
      if( iy.lt.icmpe ) then
      tsi2 = tsi2 + zsi*zvl
      endif
      xn0(iy) = zn0
      xsi(iy) = zsi
      write(n6,'(2x,i4,1p9e12.3)')
     >   iy, zne, zte, zfa, zfb, zdl, zvl, zn0, zsi, tsi
      enddo
!
      avni = tni/tvl
      ztau = tni/tsi
      znrm = ztau/taup
      tsi2 = tsi2*znrm
      zeizi = wlosi/(tsi2*cev)
      zeize = wlose/(tsi2*cev)
!
      write(n6,'(2x,"N0(spx) =",1pe10.3,"  Vol =",1p2e10.3,"  <Ne> =",
     >  1pe10.3,"  Si =",1pe10.3,"  Taup =",1p2e10.3,"  Norm =",
     >  1pe10.3)') wn0, vlmp(icspx), tvl, avni, tsi, ztau, taup, znrm
      write(n6,'(2x,"wlosi =",1pe10.3,"  Eioni =",1pe10.3,"  wlose =",
     >  1pe10.3,"  Eione =",1pe10.3)') wlosi,zeizi,wlose,zeize
!
!::absolute value
      ia = 1
      tni = 0.0d0
      tsi = 0.0d0
      twi = 0.0d0
      twe = 0.0d0
      tsi2 = 0.0d0
      twi2 = 0.0d0
      twe2 = 0.0d0
      do iy = icmps, icaxs
      zne = mdp_ni(iy)
      zn0 = znrm*xn0(iy)
      zsi = znrm*xsi(iy)
      zwi = -zsi*zeizi*cev
      zwe = -zsi*zeize*cev
      zvl = dvmp(iy)
      tni = tni + zne*zvl
      tsi = tsi + zsi*zvl
      twi = twi + zwi*zvl
      twe = twe + zwe*zvl
      an0mp(iy) = zn0
      asnmp(iy) = zsi
      if( iy.le.icmpe ) then
      sn0(ix,iy,ia) = zn0
      ssn(ix,iy,ia) = zsi
      swi(ix,iy)    = zwi
      swe(ix,iy)    = zwe
      if( iy.eq.icmpe ) cycle
      tsi2 = tsi2 + zsi*hvol(ix,iy)
      twi2 = twi2 + zwi*hvol(ix,iy)
      twe2 = twe2 + zwe*hvol(ix,iy)
      endif
      enddo
!
      write(n6,'(2x,"Main  Si =",1pe10.3,"  Wi =",1pe10.3,"  We =",
     >  1pe10.3,"  taup =",1pe10.3)') tsi, twi, twe, tni/tsi
      write(n6,'(2x,"Edge  Si =",1pe10.3,"  Wi =",1pe10.3,"  We =",
     >  1pe10.3)') tsi2, twi2, twe2
!
!----------------------------------------------------------------------
!::(N0,Si) in SOL
!----------------------------------------------------------------------
      write(n6,'(/2x,2x,"iy",4x,"Ne",10x,"Te",10x,"Tp",10x,"Vl",10x,
     >  "N0",10x,"Si",10x,"tSi",9x,"tFd")')
      ix  = jmd1
      tvl = 0.0d0
      tsi = 0.0d0
      tni = 0.0d0
      tfd = 0.0d0
      do iy = icwl1+1, icspx
      zne = mdp_ni(iy)
      zte = mdp_te(iy)
      zal = zte/10.0d0
      zal2 = zal**2
      zsgv = 3.0d-14*zal2/(3.0d0+zal2)
      zvl  = dvmp(iy)
      zn0  = wn0
      zsi  = zne*zn0*zsgv
      tvl  = tvl + zvl
      tsi  = tsi + zsi*zvl
      tni  = tni + zne*zvl
      xn0(iy) = zn0
      xsi(iy) = zsi
!
      zln = 0.5d0*wlmp(iy)
      zvp = clsav*vcs(ix,iy)
      ztp = zln/zvp
      zfr = clsxp/ztp*zvl
      alsfr(iy) = zfr
      tfd = tfd + vni(ix,iy)*zfr*clsni
      write(n6,'(2x,i4,1p8e12.3)')
     >   iy, zne, zte, ztp, zvl, zn0, zsi, tsi, tfd
      enddo
!
      gsi = tfd*clsfc
      znrms = gsi/tsi
      write(n6,'(2x,"N0(spx) =",1pe10.3,"  Si =",1pe10.3,"  Fd =",
     >  1pe10.3,"  Norm =",1pe10.3,"  gSi =",1pe10.3)')
     >   wn0, tsi, tfd, znrms, gsi
!
!::absolute value
      ia = 1
      tni = 0.0d0
      tsi = 0.0d0
      twi = 0.0d0
      twe = 0.0d0
      tsi2 = 0.0d0
      twi2 = 0.0d0
      twe2 = 0.0d0
      do iy = icwl1+1, icspx
      zne = mdp_ni(iy)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   zti = mdp_ti(iy)
      zte = mdp_te(iy)
      zn0 = znrms*xn0(iy)
      zsi = znrms*xsi(iy)
      zwi = -zsi*wirad*cev
      zwe = -zsi*werad*cev
      zvl = dvmp(iy)
      tsi = tsi + zsi*zvl
      twi = twi + zwi*zvl
      twe = twe + zwe*zvl
      an0mp(iy) = zn0
      asnmp(iy) = zsi
      sn0(ix,iy,ia) = zn0
      ssn(ix,iy,ia) = zsi
      swi(ix,iy)    = zwi
      swe(ix,iy)    = zwe
      zfr = alsfr(iy)
      tsi2 = tsi2 - vni(ix,iy)*zfr*clsni
      twi2 = twi2 - 1.5d0*vni(ix,iy)*vti(ix,iy)*cev*zfr*clsti
      twe2 = twe2 - 1.5d0*vne(ix,iy)*vte(ix,iy)*cev*zfr*clste
      enddo
!
      write(n6,'(2x,"tsi, tsi2 =",1p2e11.3,"  twi, twi2 =",1p2e11.3,
     >  "  twe, twe2 =",1p2e11.3)') tsi,tsi2,twi,twi2,twe,twe2
!
      return
      end
