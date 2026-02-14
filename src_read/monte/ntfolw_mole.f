!**********************************************************************
      subroutine hitwall_mole(iend,zincx,zint1,zint2,ztim1,ztim2
     >       ,wpox,wpoy,mtyp)
!**********************************************************************
      use cntcom, only : csem, cswl, e0em, e0wl
     >    , evel, icpo, icwl, igas
     >    , igtw, ihwl, ilon, irfw, istp, ixpo
     >    , iypo, migx, migy, mrgn, nhwl, rcwl
     >    , rmas, snem, snwl, tmas
     >    , velx, vely
     >    , tywl, weit, xpos, ypos
      use cntwcn, only : wehtm, wfhtm,whtm, wwal
      use cunit,  only : n6
      implicit none

!::arguments
      integer, intent(out) :: iend,mtyp
      real*8, intent(inout) :: zincx(0:1),zint1(0:1),zint2(0:1)
     >  ,ztim1(0:1),ztim2(0:1)
      real*8, intent(in) :: wpox,wpoy

!::local variables
      integer ker, ic, ir, ih, iw
      real*8  zrcy, zrfl, weit0, weitd, wt
!
! function
      real(8)    random

!----------------------------------------------------------------------
!::hit the wall
!----------------------------------------------------------------------
      iw = -ilon
!
!--recycling & wall temperature
      zrcy = rcwl(iw)
      e0em = e0wl(iw)
!
!::gate
      if( igtw(iw).gt.0 ) then
        zrfl = rcwl(iw)
        zrcy = 1.0d0
        if( zrfl.gt.0.0d0 ) then
          if( random(0).le.zrfl ) goto 320
        endif
!
!--pass the gate
        ih = ihwl(iw)
        whtm(ih,1) = whtm(ih,1) + weit*2.0d0
        whtm(ih,3) = whtm(ih,3) + weit*2.0d0
        wt = ztim1(1)
        call ntwarp(xpos,ypos,velx,vely,icpo,ilon,wt,wpox,wpoy,ker)
        ztim1(1) = wt
!
        if( ker.ne.0 ) then
          ih = nhwl + 1
          ic = icwl(iw)
          ir = mrgn(ic)
          whtm(ih,1) = whtm(ih,1) + weit*2.0d0
          write(n6,
     >     '(2x,"error occured at ntwarp(Mol)  ic,ix,iy,ir,ih =",
     >      i7,2i5,i4,i6)') ic,migx(ic),migy(ic),ir,ih
          call terminate_mole(iend,5)
          return
        endif
!
!::no effect
        if( irfw(ih).ne.1 ) then
          return
        endif
!
!::Mole ==> Mole with new evel at the gate
        iw   = -ilon
        ixpo = migx(icpo)
        iypo = migy(icpo)
        xpos = wpox
        ypos = wpoy
        csem = cswl(iw)
        snem = snwl(iw)
        e0em = e0wl(iw)
        call ntrc680( istp,xpos,ypos,weit,icpo,"ScatwM-"//tywl(iw) )
        call set_velocity(iw)
        call initialize_rand(1,zincx,zint1,zint2,ztim1,ztim2)
        return
      endif     ! <== gate
!
!--hit point
 320  continue
      icpo = icwl(iw)
      ixpo = migx(icpo)
      iypo = migy(icpo)
      xpos = xpos + velx*ztim1(1)
      ypos = ypos + vely*ztim1(1)
      csem = cswl(iw)
      snem = snwl(iw)
!
      call ntrc680( istp,xpos,ypos,weit,icpo,"hitD2-"//tywl(iw) )
!
!::neutral flux onto the wall (lopht)
      wfhtm(iw) = wfhtm(iw) + weit
      wehtm(iw) = wehtm(iw) + weit*evel
      call ntwedf(iw, weit, weit*evel, 2)
!
!--new weit
      weit0 = weit
      weit  = weit*zrcy
      weitd = weit0 - weit
      wwal(iw) = wwal(iw) + weitd*2.0d0
!
      ih = ihwl(iw)
      whtm(ih,1) = whtm(ih,1) + weit0*2.0d0
      whtm(ih,2) = whtm(ih,2) + weit*2.0d0
      whtm(ih,3) = whtm(ih,3) + weitd*2.0d0
!
!--absorpotion
      if( weit.le.0.0d0 ) then
        call terminate_mole(iend,3)
        return
      endif
!
      !::energy & angle of wall segment
      mtyp = 1
      tmas = rmas(igas)*2.0d0
      call set_velocity(iw)
      call initialize_rand(1,zincx,zint1,zint2,ztim1,ztim2)
      end subroutine hitwall_mole

!**********************************************************************
      subroutine collision_mole(iend,zincx,zint1,zint2,ztim1,ztim2
     >  ,mtyp)
!**********************************************************************
      use cntcom, only : bvx, bvy, bvz, derc
     >    , evel, evelb, fnfi, frct, fvnth, icpo, igas
     >    , igasb, igrc, ilon, irct, irfw, istp
     >    , lbrc, migx, migy, nrct
     >    , rmas, tcfi, tmas, tmasb, trct, tsfi
     >    , tvcth, tvsth, vel, vel2, velp, velpb, velx, vely, velz
     >    , weit, weitb, xpos, ypos
      use cntwcn, only : wgden, wgeng, wtion
      use cphcns, only : cev
      use cunit,  only : n6
      implicit none

!:: arguments
      integer,intent(out) :: iend
      integer,intent(inout) :: mtyp
      real*8,intent(inout) :: zincx(0:1),zint1(0:1),zint2(0:1)
     > ,ztim1(0:1),ztim2(0:1)

!::local variables
      integer i, isel, irc, ipl
      integer ith, ifi
      real*8  xran, zdlt, zdn, zde, edis
      character  clty*4
      real*8  wpox, wpoy ! dummy

! function
      real(8)    random

!----------------------------------------------------------------------
!::collision occured
!----------------------------------------------------------------------
!:: check trct
      if( trct.eq.0.0d0 ) then
         write(n6,'(2x,"*** ntmole   col. occured in void  ",i5)') icpo
         call terminate_mole(iend,5)
         return
      endif   
      zdlt  = (zincx(1)-zint2(1))/trct
      ztim1(1) = ztim2(1) + zdlt
!
!--save particle information before collision
      velpb=velp; evelb=evel; tmasb=tmas; weitb=weit; igasb=igas
!
!--new start position
      xpos = xpos + velx*ztim1(1)
      ypos = ypos + vely*ztim1(1)
!
!
!--scoreing density and temperature
      zdn = zdlt*weit
      zde = zdn*vel2
      wgden(icpo,mtyp) = wgden(icpo,mtyp) + zdn
      wgeng(icpo,mtyp) = wgeng(icpo,mtyp) + zde
!
!--track estimator  before H2 collide
      call ntscor_tr(zdlt)
!
!--what type of collision
 410  continue
      xran = random(0)
      do i = 1, nrct
        isel = i
        if( xran.le.frct(i) ) goto 430
      enddo
      isel = nrct
 430  continue
      edis = derc(isel)
      irc  = irct(isel)
      ipl  = igrc(isel)
      clty = lbrc(irc)
!
      call ntrc680( istp,xpos,ypos,weit,icpo,"collD2"//clty )

!--irc = 4    h2dis      [e + H2  => e + H + H]
      if( clty.eq."DS1 " ) then
        tmas = rmas(igas)
        call ntscor_cl(irc)    !  DS1
        evel = 0.5d0*tmas*vel2/cev + edis
        weit = 2.0d0*weit
        goto 440
!
!--irc = 5    h2ion      [e + H2  => e + H2+ + e]
      elseif( clty.eq."DS2 " ) then
        tmas = rmas(igas)*2.0d0
        call ntscor_cl(irc)    !  DS2
!
!  assumed that the collision occured immediately
!  the displcement of H2+ ion is assumed to be small.
        mtyp = 2
        call ntcrs_MHP(icpo,vel)
!
        if( trct.eq.0.0d0 ) then
          write(n6,
     >      '(2x,"*** ntmole   col.(H20 => H2+) occured in void  "
     >      ,3i6)') icpo,migx(icpo),migy(icpo)
          call terminate_mole(iend,5)
          return
        endif
!
!::scorering
        velp = velx*bvx(icpo)+vely*bvy(icpo)+velz*bvz(icpo)
        zdlt = 1.0d0/trct
        zdn = zdlt*weit
        zde = zdn*vel2
        wgden(icpo,mtyp) = wgden(icpo,mtyp) + zdn
        wgeng(icpo,mtyp) = wgeng(icpo,mtyp) + zde
!
!--track estimator   before H2+ collide
        call ntscor_tr(zdlt)
!
        goto 410
!
!--irc = 6    h2dision   [e + H2  => e + H+ + H + e]
      elseif( clty.eq."DS3 " ) then
        tmas = rmas(igas)
        call ntscor_cl(irc)    !  DS3
        evel = 0.5d0*tmas*vel2/cev + edis
        wtion(icpo,igas) = wtion(icpo,igas) + weit
        goto 440
!
!--irc = 8    h2pdis     [e + H2+ => e + H+ + H]
      elseif( clty.eq."DS4 " ) then
        tmas = rmas(igas)
        call ntscor_cl(irc)    !  DS4
        evel = 0.5d0*tmas*vel2/cev + edis
        wtion(icpo,igas) = wtion(icpo,igas) + weit
        goto 440
!
!--irc = 9    h2pdision  [e + H2+ => e + H+ + H+ + e]
      elseif( clty.eq."DS5 " ) then
        tmas = rmas(igas)
        call ntscor_cl(irc)    !  DS5
        wtion(icpo,igas) = wtion(icpo,igas) + 2.0d0*weit
        weit = 0.0d0
        call terminate_mole(iend,4)
        return
!
!--irc =10    h2pdisrec  [e + H2+ => H + H *]
      elseif( clty.eq."DS6 " ) then
        tmas = rmas(igas)
        call ntscor_cl(irc)    !  DS6
        evel = 0.5d0*tmas*vel2/cev + edis
        weit = 2.0d0*weit
        goto 440
!
!--irc = 7    elastic    [H+ + H2 => H+ + H2]
      elseif( clty.eq."EL2 " ) then
        call ntvel_el(irc,ipl,xpos,ypos,wpox,wpoy)
        call ntscor_cl(irc)    !  EL2
        ilon = 0
        call initialize_rand(1,zincx,zint1,zint2,ztim1,ztim2)
        return
!
!--irc = 14    elastic    [D2 + D => D2 + D]
      elseif( clty.eq."EL5 " ) then
        call ntvel_el(irc,ipl,xpos,ypos,wpox,wpoy)
        ilon = 0
        call initialize_rand(1,zincx,zint1,zint2,ztim1,ztim2)
        return
!
!--irc = 15    elastic    [H2 + H2 => H2 + H2]
      elseif( clty.eq."EL6 " ) then
        call ntvel_el(irc,ipl,xpos,ypos,wpox,wpoy)
        ilon = 0
        call initialize_rand(1,zincx,zint1,zint2,ztim1,ztim2)
        return
!
!--irc = ?
      else
        call wexit("ntmole","wrong selection of reaction")
      endif
!
!--new velocity of atom
 440  continue
      vel2 = 2.0d0*evel*cev/tmas
      vel  = dsqrt(vel2)
      ith  = int(fvnth*random(0) + 1.0d0)
      ifi  = int(fnfi *random(0) + 1.0d0)
      velx = vel*tvsth(ith)*tcfi(ifi)
      vely = vel*tvsth(ith)*tsfi(ifi)
      velz = vel*tvcth(ith)
      vel2 = velx*velx+vely*vely+velz*velz
      evel = 0.5d0*tmas*vel2/cev
      vel  = dsqrt(vel2)
      ilon = 0
      call terminate_mole(iend,1)
      return
      end subroutine collision_mole
!
!**********************************************************************
      subroutine terminate_mole(iend,iend_in)
!**********************************************************************
      use cntcom, only : weit,ievt,ixpo,iypo,iptl
     > , icpo, migx, migy, istp
      use cntwcn, only : wend
      use csonic, only : itim
      use cunit,  only : mype, n6
      implicit none
! arguments
      integer, intent(inout) :: iend
      integer, intent(in) :: iend_in
!
      iend = iend_in
!
      if(iend == 1) then  
        continue
!--too many step
      elseif(iend == 2) then  
        wend(4) = wend(4) + weit*2.0d0
        weit  = 0.0d0
        write(n6,'(2x,"too many step  ",2i8,2i4)') istp,ievt,ixpo,iypo
!--weit = 0
      elseif(iend == 3) then  
        continue
!--weit = 0
      elseif(iend == 4) then        
        continue
!--error
      elseif(iend == 5) then
        wend(5) = wend(5) + weit*2.0d0
        write(n6,'(2x,"Error ntmole  mype =",i5,"  itim =",i6,
     >    "  iptl =",i6,"  ic =",3i5,"  weit,werr =",1p2e12.3)')
     >    mype,itim,iptl,icpo,migx(icpo),migy(icpo),weit*2.0d0,wend(5)
        weit = 0.0d0
      endif
      end subroutine terminate_mole

!**********************************************************************
      subroutine set_velocity(iw)
!**********************************************************************
      use cntcom, only : csem, e0em
     >    , evel, fnfi, fnth, icpo
     >    , istp, snem, tcfi, tcth, tmas, tsfi, tsth
     >    , vel, vel2, velx, vely, velz
     >    , weit, xpos, ypos
      use cntwcn, only : weemm
      use cphcns, only : cev
      implicit none
!
!::arguments
      integer, intent(in) :: iw
!
!::local variables
      real*8  cosw, sinw, ctx, stx, cfx, sfx
      real*8  vxp, vyp, vzp
      integer ith, ifi
!
! function
      real(8)    random

!::energy & angle of wall segment
      evel = e0em
      cosw = csem
      sinw = snem
      call ntrc680( istp,xpos,ypos,weit,icpo,"mole_emit" )
!
!::velocity (cosine distribution)
      vel2 = 2.0d0*evel*cev/tmas
      vel  = dsqrt(vel2)
      ith  = int(fnth*random(0) + 1.0d0)
      ctx  = tcth(ith)
      stx  = tsth(ith)
      ifi  = int(fnfi*random(0) + 1.0d0)
      cfx  = tcfi(ifi)
      sfx  = tsfi(ifi)
      vxp  = vel*stx*sfx
      vyp  = vel*ctx
      vzp  = vel*stx*cfx
!
      velx = vxp*cosw - vyp*sinw
      vely = vxp*sinw + vyp*cosw
      velz = vzp
      vel2 = velx*velx+vely*vely+velz*velz
      vel  = dsqrt(vel2)
      evel = 0.5d0*tmas*vel2/cev
!
      if(iw.ge.0) weemm(iw) = weemm(iw) + evel*weit
      end subroutine set_velocity
