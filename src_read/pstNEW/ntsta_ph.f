!**********************************************************************
      subroutine ntsta_ph
!**********************************************************************
!
!      generate a sample particle due to photon
!
!         icpo : position of particle   in neut2d
!         ixpo : x-position of particle in soldor
!         iypo : y-position of particle in soldor
!         izpo : z-position of particle   (dummy)
!         ilon : type of cell boundary  in neut2d
!
!         ixbc : x-position of boundary in soldor
!         iybc : y-position of boundary in soldor
!
!----------------------------------------------------------------------
      use cntcom, only : evel, icbr, icpo, ievt, igas, ilon, iptl, istp
     >    , ixpo, iypo, ltrc, migx, migy, nbr, prbr, rmas, tmas, weit
     >    , vel, vel2, velx, vely, velz, xpos, ypos, zpos
      use cntwcn, only : wsbr, wtot
      use cphcns, only : cpi
      use cunit,  only : mype, n6
      implicit none
!
!::local variables
      real*8  xran
      real*8  ss, tt, zth, zfi, cth, sth, cfi, sfi
      real*8  wtim, wpox, wpoy
      integer i, ii, iw, ier
      character  clin*128
! function
      real(8)    random
!
!----------------------------------------------------------------------
!::wall number (cell number)
!----------------------------------------------------------------------
      xran = random(0)
      do i = 1, nbr
      ii = i
      if( prbr(i).ge.xran ) goto 115
      enddo
      ii = nbr
 115  continue
      iw = ii
!
!----------------------------------------------------------------------
!::cell number
!----------------------------------------------------------------------
      icpo = icbr(iw)
      ixpo = migx(icpo)
      iypo = migy(icpo)
!
      ilon = 0
!
      if( ixpo.le.0 .or. iypo.le.0 ) then
      write(clin,'("-- ntstav --  ii,icpo,ixpo,iypo =",4i6)')
     >   ii,icpo,ixpo,iypo
      call wexit("ntsta_ph",clin)
      endif
!
!----------------------------------------------------------------------
!::start position
!----------------------------------------------------------------------
      istp = 0
      ievt = 0
      ss   = random(0)
      tt   = random(0)
      call ctorsp(icpo,ss,tt,xpos,ypos,ier)
      zpos = 0.0d0
!
!----------------------------------------------------------------------
!::particle type
!----------------------------------------------------------------------
      igas = 1
      tmas = rmas(igas)
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
      zth = cpi*random(0)
      zfi = 2.0d0*cpi*random(0)
      cth = cos(zth)
      sth = sin(zth)
      cfi = cos(zfi)
      sfi = sin(zfi)
      velx = sth*cfi
      vely = sth*sfi
      velz = cth
      vel2 = velx*velx + vely*vely + velz*velz
      evel = vel2
      vel  = dsqrt(vel2)
!
!----------------------------------------------------------------------
!::return
!----------------------------------------------------------------------
      wtim = 0.0d0
      wpox = xpos
      wpoy = ypos
!
      if( ltrc.eq.1 .and. iptl.le.20 ) then
      write(n6,602) istp
     > ,wtim,wpox,wpoy,velx/vel,vely/vel,velz/vel,evel,icpo,ixpo,iypo,
     >  ilon,iptl,mype
      endif
!
 602  format(2x,i5,1p7e12.4,6i6,2x,a)
!
      return
      end
