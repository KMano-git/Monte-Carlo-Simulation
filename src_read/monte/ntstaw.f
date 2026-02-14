!**********************************************************************
      subroutine ntstaw
!**********************************************************************
!
!      emit a sample particle from plate/wall/puff
!
!       need define ifgwl, iwbr(iw_p)
!         = 0 : wl D0 3eV inward
!         = 1 : wl D0 Ti  outward
!         = 3 : dp
!         = 4 : pf D0 3eV inward
!         = 4 : pf D2 0.038 inward
!
!         icpo : position of particle   in neut2d
!         ixpo : x-position of particle in soldor
!         iypo : y-position of particle in soldor
!         izpo : z-position of particle   (dummy)
!         ilon : type of cell boundary  in neut2d
!
!----------------------------------------------------------------------
      use cntcom, only : csbr, csem, e0br, e0em, e0pfa, e0pfm, eibr
     >    , evel, icbr, icpo, iebr, ievt, igas, ikbr, ikon, ilon, isbr
     >    , istp, iwbr, iwpbr, ixpo, iypo, migx, migy, mknd, nbr, next
     >    , ntrc, prbr, rmas, snbr, snem, tmas, vel, vel2, velx, vely
     >    , velz, weit, xpnt, xpos, ypnt, ypos, zpos
     >    , is_atom_or_mole, ntmont_flg
     >    , iwpbr_wall,ipwl
      use cntwcn, only : wees, wsbr, wtot
      use cphcns, only : cev
      use mod_externalgrid, only : use_exdata
      implicit none
!
!::model
      integer  lp
      real*8   vwx, vwy, vwz
!
!::local variables
      real*8  xran, ein, cosw, sinw
      integer i, iw_p, ik, ip, ip1, ic_new, ker
      character  clin*120
      integer  nxon(4)
      data nxon/3,4,1,2/
! function
      real(8)    random
!
!::check
      ntrc = 0
      call ntrc680( 0, 0.0d0,0.0d0,0.0d0, 0, "[ntstaw]" )
!
!----------------------------------------------------------------------
!::wall number
!----------------------------------------------------------------------
      xran = random(0)
      do i = 1, nbr
        iw_p = i
        if( prbr(i).ge.xran ) goto 115
      enddo
      iw_p = nbr
 115  continue
!
!----------------------------------------------------------------------
!::cell number
!----------------------------------------------------------------------
      icpo = icbr(iw_p)
      ixpo = migx(icpo)
      iypo = migy(icpo)
!
      ik   = ikbr(iw_p)
      ikon = ik
      ilon = mknd(icpo,ik)
!
      if( ixpo.le.0 .or. iypo.le.0 ) then
        write(clin,'("-- ntstaw --  iw_p,ixpo,iypo =",3i5)')
     >    iw_p,ixpo,iypo
        call wexit("ntstaw",clin)
      endif
!
!----------------------------------------------------------------------
!::start position
!----------------------------------------------------------------------
      istp = 0
      ievt = 0
      xran = random(0)
      ip   = isbr(iw_p)
      ip1  = iebr(iw_p)
      xpos = xpnt(ip) + xran*(xpnt(ip1)-xpnt(ip))
      ypos = ypnt(ip) + xran*(ypnt(ip1)-ypnt(ip))
      zpos = 0.0d0
!
!----------------------------------------------------------------------
!::particle type
!----------------------------------------------------------------------
      igas = 1
      tmas = rmas(igas)
      evel = 0.0d0        ! temporary
!
!----------------------------------------------------------------------
!::weit
!----------------------------------------------------------------------
      weit = 1.0d0
      wtot = wtot + weit
      wsbr(icpo,igas) = wsbr(icpo,igas) + weit
!
!----------------------------------------------------------------------
!::new velocity
!----------------------------------------------------------------------
      csem = csbr(iw_p)
      snem = snbr(iw_p)
      e0em = e0br(iw_p)
!
!::wall (D0 3eV inward) OLD
      if( iwbr(iw_p).eq.1 ) then
        call ntrc680( 0,xpos,ypos,weit,icpo,"staw=wli" )
        evel = e0em
        vel2 = 2.0d0*evel*cev/tmas
        vel  = dsqrt(vel2)
        cosw = csem
        sinw = snem
        velx = - vel*sinw
        vely = + vel*cosw
        velz = 0.0d0
        vel2 = velx*velx+vely*vely+velz*velz
        vel  = dsqrt(vel2)
        evel = 0.5d0*tmas*vel2/cev
        wees(iwpbr(iw_p)) = wees(iwpbr(iw_p)) + evel
!
!::wall (D0 Ti outward) EIRENE type (NEW)
      elseif( iwbr(iw_p).eq.2 ) then
        call ntrc680( 0,xpos,ypos,weit,icpo,"staw=wlo" )
        call ntvelw(iw_p,vwx,vwy,vwz,velx,vely,velz)
        vel2 = velx*velx+vely*vely+velz*velz
        vel  = dsqrt(vel2)
        evel = 0.5d0*tmas*vel2/cev
!
        ic_new = next(icpo,ik)
        if(use_exdata) then
          if(ic_new.eq.0) then
            call crossBorder(xpos,ypos,icpo,ic_new,ker)
          endif
        else
          ik = nxon(ik)
        endif
        icpo = ic_new
        ikon = ik
        ilon = mknd(icpo,ik)
        call ntrc680( 0,xpos,ypos,weit,icpo,"staw_nxt" )
        wees(iwpbr(iw_p)) = wees(iwpbr(iw_p)) + evel
!
!::d-plate
      elseif( iwbr(iw_p).eq.3 ) then
        call ntrc680( 0,xpos,ypos,weit,icpo,"staw=dpl" )
        ein = eibr(iw_p) ! ein=0 for lemd2=1, see ntnflx
        call ntrefl(iwpbr_wall(iw_p), ein,csem,snem)
        if(is_atom_or_mole==1) ntmont_flg = 1
!
!::puff (D0 3eV inward)
      elseif( iwbr(iw_p).eq.4 ) then
        call ntrc680( 0,xpos,ypos,weit,icpo,"staw=pfA" )
        e0em = e0pfa
        evel = e0em
        vel2 = 2.0d0*evel*cev/tmas
        vel  = dsqrt(vel2)
        cosw = csem
        sinw = snem
        velx = - vel*sinw
        vely = + vel*cosw
        velz = 0.0d0
        vel2 = velx*velx+vely*vely+velz*velz
        vel  = dsqrt(vel2)
        evel = 0.5d0*tmas*vel2/cev
!
!::puff (D2 e0pf inward)
      elseif( iwbr(iw_p).eq.5 ) then
        call ntrc680( 0,xpos,ypos,weit,icpo,"staw=pfM" )
        e0em = e0pfm
        ! ntmole change
        weit = weit*0.5d0
        tmas = rmas(igas)*2.0d0
        call set_velocity(-1)
        is_atom_or_mole = 1 ! molecular
        ntmont_flg = 1
!
!::error
      else
        call wexit("ntstaw","invalid iwbr")
      endif
!
!----------------------------------------------------------------------
!::return
!----------------------------------------------------------------------
!::debug
      if(ntmont_flg.eq.0) then
        call ntrc680( 1,xpos,ypos,weit,icpo,"staw=end" )
      endif
      return
      end
