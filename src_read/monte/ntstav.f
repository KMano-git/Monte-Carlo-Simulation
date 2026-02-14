!**********************************************************************
      subroutine ntstav
!**********************************************************************
!
!      generate a sample particle due to recombination
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
!::serious error    2010/03/06
!xx   include 'cplcom'  contains ss
!
!----------------------------------------------------------------------
      use cntcom, only : bvx, bvy, bvz, evel, fnor, icbr, icpo, ievt
     >    , igas, ikon, ilon, istp, ixpo, iypo, migx, migy, nbr, ntrc
     >    , prbr, rmas, tmas, tnor, vel, vel2, velx, vely, velz, vion
     >    , weit, xpos, ypos, zpos
      use cntpls, only : vflw
      use cntwcn, only : wsbr, wtot
      use cphcns, only : cev
      use mod_externalgrid, only : use_exdata
      implicit none
!
!::local variables
      integer i, ii, ier, ivx, ivy, ivz
      real*8  xran, zvi
      real*8  ss, tt
! function
      real(8)    random
!
      character clin*120
!
!::check
      ntrc = 0
      call ntrc680( 0, 0.0d0,0.0d0,0.0d0, 0, "ntstav" )
!
!----------------------------------------------------------------------
!::wall number (cell number)
!----------------------------------------------------------------------
      xran = random(0)
      do 110 i = 1, nbr
      ii = i
      if( prbr(i).ge.xran ) goto 115
 110  continue
      ii = nbr
 115  continue
!
!----------------------------------------------------------------------
!::cell number
!----------------------------------------------------------------------
      icpo = icbr(ii)
      ixpo = migx(icpo)
      iypo = migy(icpo)
!
      ilon = 0
      ikon = 0
!
      if( .not.use_exdata .and. (ixpo.le.0 .or. iypo.le.0) ) then
        write(clin,'("-- ntstav --  ii,icpo,ixpo,iypo =",4i6)')
     >   ii,icpo,ixpo,iypo
        call wexit("ntstav",clin)
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
      ivx  = int(fnor*random(0) + 1.0d0)
      ivy  = int(fnor*random(0) + 1.0d0)
      ivz  = int(fnor*random(0) + 1.0d0)
      zvi  = vion(icpo,igas)
      velx = zvi*tnor(ivx) + vflw(icpo,igas)*bvx(icpo)
      vely = zvi*tnor(ivy) + vflw(icpo,igas)*bvy(icpo)
      velz = zvi*tnor(ivz) + vflw(icpo,igas)*bvz(icpo)
      vel2 = velx*velx + vely*vely + velz*velz
      evel = 0.5d0*tmas*vel2/cev
      vel  = dsqrt(vel2)
!
!----------------------------------------------------------------------
!::return
!----------------------------------------------------------------------
!::debug
      call ntrc680( 1,xpos,ypos,weit,icpo,"stav-end" )
 602  format(2x,i5,1p7e14.4,6i6,2x,a)
!
      return
      end
