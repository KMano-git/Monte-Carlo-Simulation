!**********************************************************************
      subroutine hitwall_atom(iend,zincx,zint1,zint2,ztim1,ztim2
     >       ,wpox,wpoy)
!**********************************************************************
      use cntcom, only : csem, cswl, e0em, e0wl, evel
     >    , icpo, icwl, igtw
     >    , ihwl, ilon, irfw, istp, ixpo, iypo, lemd2
     >    , migx, migy, mrgn, nhwl, rcwl, snem, snwl
     >    , tywl
     >    , velx, vely, velz, weit, xpos, ypos
     >    , is_atom_or_mole, ntmont_flg
      use cntwcn, only : wehta, wfhta, whta, wwal
      use cunit,  only : n6
      implicit none
!::arguments
      integer, intent(out) :: iend
      real*8, intent(in) :: wpox,wpoy
      real*8, intent(inout) :: zincx(0:1),zint1(0:1),zint2(0:1)
     >  ,ztim1(0:1),ztim2(0:1)
!
!::local variables
      integer ker, ic, ir, ih, iw
      real*8  zrcy, zrfl, weit0, weitd, ein, wt
!
! function
      real(8)    random

!----------------------------------------------------------------------
!::hit the wall
!----------------------------------------------------------------------
      iw = -ilon
!--recycling & wall temperature
      zrcy = rcwl(iw)
      e0em = e0wl(iw)
!
!::igtw : gate(1,2)/shevron(3,4)/partition(5)
!::rcwl :   0.0    /  0.5       / 1.0    reflection
!::zrcy :   1.0       1.0         1.0    new weight
!
      if( igtw(iw).gt.0 ) then  !<== gate
        zrfl = rcwl(iw)
        zrcy = 1.0d0
        if( zrfl.gt.0.0d0 ) then
          if( random(0).le.zrfl ) goto 320
        endif
!
!--pass the gate/shevron
        ih = ihwl(iw)
        whta(ih,1) = whta(ih,1) + weit  ! (flux-in)
        whta(ih,3) = whta(ih,3) + weit  ! (flux-pass)
        wt = ztim1(0)
        call ntwarp(xpos,ypos,velx,vely,icpo,ilon,wt,wpox,wpoy,ker)
        ztim1(0) = wt
        if( ker.ne.0 ) then
          ih = nhwl + 1
          whta(ih,1) = whta(ih,1) + weit
          ic = icwl(iw)
          ir = mrgn(ic)
          write(n6,'(2x,"error occured at ntwarp(Hyd)
     >      ic,ix,iy,ir,ih =",
     >      i7,2i5,i4,i6)') ic,migx(ic),migy(ic),ir,ih
          call terminate_atom(iend,5)
          return
        endif
!
!--no effect
        if( irfw(ih).ne.1 ) then
          return
        endif
!
!--Atom ==> Mole at the gate
        iw   = -ilon
        ixpo = migx(icpo)
        iypo = migy(icpo)
        xpos = wpox
        ypos = wpoy
        csem = cswl(iw)
        snem = snwl(iw)
        e0em = e0wl(iw)
        call ntrc680( istp,xpos,ypos,weit,icpo,"ScatwA-"//tywl(iw) )
        ir  = mrgn(icpo)
        ein = 0.0d0
        call ntrefl( iw, ein,csem,snem )
        if(is_atom_or_mole == 1) then
          ntmont_flg = 2
          return
        endif
        if( weit.le.0.0d0 ) then
          call terminate_atom(iend,6)
          return
        endif
        call initialize_rand(0,zincx,zint1,zint2,ztim1,ztim2)
        return
      endif  ! <== gate
!
!--hit point
 320  continue
      icpo = icwl(iw)
      ixpo = migx(icpo)
      iypo = migy(icpo)
      xpos = xpos + velx*ztim1(0)
      ypos = ypos + vely*ztim1(0)
      csem = cswl(iw)
      snem = snwl(iw)
!
      call ntrc680( istp,xpos,ypos,weit,icpo,"hitw-"//tywl(iw) )
!
!::neutral flux onto the wall (lopht)
      wfhta(iw) = wfhta(iw) + weit
      wehta(iw) = wehta(iw) + weit*evel
      call ntwedf(iw, weit, weit*evel, 1)
!
!--new weit
      weit0 = weit
      weit  = weit*zrcy
      weitd = weit0 - weit
      wwal(iw) = wwal(iw) + weitd
!
      ih = ihwl(iw)
      whta(ih,1) = whta(ih,1) + weit0   ! (flux-in)
      whta(ih,2) = whta(ih,2) + weit    ! (flux-refl)
      whta(ih,3) = whta(ih,3) + weitd   ! (flux-pass/absorp)
!
!--absorpotion
      if( weit.le.0.0d0 ) then
        call terminate_atom(iend,6)
        return
      endif
!
!--reflection
!--new velocity of mole
      ein = evel
      if( irfw(ih).eq.1 ) ein = 0.0d0
      if( lemd2.eq.1 ) ein = 0.0d0
      call ntrefl( iw, ein,csem,snem )
      if(is_atom_or_mole == 1) then
        ntmont_flg = 2
        return
      endif
      if( weit.le.0.0d0 ) then
        call terminate_atom(iend,6)
        return
      endif
!
!--reflection/absorption
      call initialize_rand(0,zincx,zint1,zint2,ztim1,ztim2)
      end subroutine hitwall_atom

!**********************************************************************
      subroutine collision_atom(iend,zincx,zint1,zint2,ztim1,ztim2)
!**********************************************************************
      use cntcom, only : bvx, bvy, bvz, evel
     >    , evelb, fnor, frct, icpo, igas, igasb, igrc
     >    , ilon, irct, istp, nrct, rmas
     >    , tmas, tmasb, tnor, trct, trcx, vel, vel2, velp, velpb
     >    , velx, vely, velz, vion, weit, weitb, wtmin, xpos, ypos
      use cntpls, only : vflw
      use cntwcn, only : wden, weng, wnfl0x
     >    , wnfl0y, wnfl0z, wtion, wvlp
      use cphcns, only : cev
      use cunit, only:n6
      implicit none

!:: arguments
      integer,intent(inout) :: iend
      real*8, intent(inout) :: zincx(0:1),zint1(0:1),zint2(0:1)
     >  ,ztim1(0:1),ztim2(0:1)

!::local variables
      integer i, isel, irc, ipl, ivx, ivy, ivz
      real*8  zdlt, zdn, zde, zdv, weitd, zeps, zvi
      real*8  wpox, wpoy ! dummy

! function
      real(8)    random

!----------------------------------------------------------------------
!::collision occured
!----------------------------------------------------------------------
!::check trct
      if( trct.eq.0.0d0 ) then
        write(n6,'(2x,"*** ntfolw   col. occured in void  "
     >       ,i5)') icpo
        call terminate_atom(iend,5)
        return
      endif
!
      zdlt  = (zincx(0)-zint2(0))/trct
      ztim1(0) = ztim2(0) + zdlt
!
!--new start position
      xpos = xpos + velx*ztim1(0)
      ypos = ypos + vely*ztim1(0)
!
!--scoreing density and temperature
      velp = velx*bvx(icpo)+vely*bvy(icpo)+velz*bvz(icpo)
      zdn = zdlt*weit
      zde = zdn*vel2
      zdv = zdn*velp
      wden(icpo,igas) = wden(icpo,igas) + zdn
      weng(icpo,igas) = weng(icpo,igas) + zde
      wvlp(icpo,igas) = wvlp(icpo,igas) + zdv
! count neutral flow, 2016.06.23 toku
      wnfl0x(icpo,igas)=wnfl0x(icpo,igas)+zdn*velx
      wnfl0y(icpo,igas)=wnfl0y(icpo,igas)+zdn*vely
      wnfl0z(icpo,igas)=wnfl0z(icpo,igas)+zdn*velz

!--track estimator
      call ntscor_tr(zdlt)
!
!--save particle information before collision
      velpb=velp; evelb=evel; tmasb=tmas; weitb=weit; igasb=igas
!
!--New weit
!--weit is too small
      if( weit.lt.wtmin ) then
        zeps = random(0)*trct
        if( zeps.gt.trcx ) then
          weit = 0.0d0  ! (ei)
        endif
      else
        weit  = weit*trcx/trct  ! (ei/cx)
      endif
!
!--collision estimator
      call ntscor_cl(1)      !  EI
!
!--ionization
      weitd = weitb - weit
      wtion(icpo,igas) = wtion(icpo,igas) + weitd
!
      if( weit.le.0.0d0 ) then
        call ntrc680( istp,xpos,ypos,weit,icpo,"EI" )
        call terminate_atom(iend,1)
        return
      endif
!
!----------------------------------------------------------------------
!::choose collision type CX/EL and ion species
!----------------------------------------------------------------------
      zeps = random(0)
      do i = 2, nrct    !  excluding ionziation process (EI:i=1)
        isel = i
        if( zeps.le.frct(i) ) goto 470
      enddo
      isel = nrct
 470  continue
      irc = irct(isel)
      ipl = igrc(isel)
!
!----------------------------------------------------------------------
!::charge exchange [2]
!----------------------------------------------------------------------
      if( irc.eq.2 ) then
        call ntrc680( istp,xpos,ypos,weit,icpo,"CX" )
!
        igas = ipl
        tmas = rmas(igas)
!
!--new velocity
        ivx  = int(fnor*random(0) + 1.0d0)
        ivy  = int(fnor*random(0) + 1.0d0)
        ivz  = int(fnor*random(0) + 1.0d0)
        zvi  = vion(icpo,igas)
        velx = zvi*tnor(ivx) + vflw(icpo,igas)*bvx(icpo)
        vely = zvi*tnor(ivy) + vflw(icpo,igas)*bvy(icpo)
        velz = zvi*tnor(ivz) + vflw(icpo,igas)*bvz(icpo)
        velp = velx*bvx(icpo)+vely*bvy(icpo)+velz*bvz(icpo)
        vel2 = velx*velx + vely*vely + velz*velz
        evel = 0.5d0*tmas*vel2/cev
        vel  = dsqrt(vel2)
!
!--collision estimator
        call ntscor_cl(irc)    !  CX
!
!----------------------------------------------------------------------
!::elastic collision [3]
!----------------------------------------------------------------------
      elseif( irc.eq.3 ) then
        call ntrc680( istp,xpos,ypos,weit,icpo,"EL" )
!
!--new velocity
        call ntvel_el(irc,ipl,xpos,ypos,wpox,wpoy)
!
!--collision estimator
        call ntscor_cl(irc)    !  EL1
!
!----------------------------------------------------------------------
!::n-n elastic collision [12] EL3, D+D
! COLLISION ESTIMATOR HAS NOT BEEN INPLEMENTED!!
      elseif( irc.eq.12 ) then
        call ntrc680( istp,xpos,ypos,weit,icpo,"EL3" )
        call ntvel_el(irc,ipl,xpos,ypos,wpox,wpoy)
!::n-n elastic collision [13] EL4, D+D2
      elseif( irc.eq.13 ) then
        call ntrc680( istp,xpos,ypos,weit,icpo,"EL4" )
        call ntvel_el(irc,ipl,xpos,ypos,wpox,wpoy)
!::error collision [?]
      else
        call wexit("ntfolw","wrong selection of reaction")
      endif
!
!--go back to start point
      ilon = 0
      call initialize_rand(0,zincx,zint1,zint2,ztim1,ztim2)
      end subroutine collision_atom

!**********************************************************************
      subroutine terminate_atom(iend,iend_in)
!**********************************************************************
      use cntcom, only : weitb,weit,ievt,lsmax,ixpo,iypo,iptl
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
!--ionization
      if(iend == 1) then      
        wend(1) = wend(1) + weitb
!--too many step
      elseif(iend == 4) then        
        wend(4) = wend(4) + weit
        write(n6,'(2x,"too many step  istp,ievt > lsmax ",3i6,
     >  "  ixpo,iypo =",2i6)') istp,ievt,lsmax,ixpo,iypo
!--error
      elseif(iend == 5) then
        wend(5) = wend(5) + weit
        write(n6,'(2x,"Error ntfolw  mype =",i5,
     >    "  itim =",i6,"  iptl ="
     >    ,i6,"  ic =",3i5,"  weit,werr =",1p2e12.3)')
     >    mype,itim,iptl,icpo,migx(icpo),migy(icpo),weit, wend(5)
!--normal end (weit = 0.0)
      else
        iend = 6
      endif
      end subroutine terminate_atom