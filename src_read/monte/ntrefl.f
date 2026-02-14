!**********************************************************************
      subroutine ntrefl( iw, ein,cosw,sinw )
!**********************************************************************
      use cntcom, only : evel, fnfi, fnth, icpo, ilon, iptl, istp, ixpo
     >    , iypo, ntrc, tcfi, tcth, tmas, tsfi, tsth, vel, vel2, velx
     >    , vely, velz, weit, xpos, ypos
     >    , rmas, igas, is_atom_or_mole
      use cntwcn, only : weema
      use cphcns, only : cev
      use cunit,  only : n6
      implicit none
!
!::argument
      integer, intent(in) :: iw
      real(8), intent(in) :: ein, cosw, sinw
!
!::local variables
      real*8  erfl, rn, re, xran, ctx, stx, cfx, sfx
      real*8  vxp, vyp, vzp, ztim1
      integer  ith, ifi
! function
      real(8)    random
!
!----------------------------------------------------------------------
!::emission H2 / reflection
!----------------------------------------------------------------------
      if( ein.le.0.0d0 ) goto 100
!
!----------------------------------------------------------------------
!::reflection coefficients
!----------------------------------------------------------------------
      call refcof( ein, erfl, rn, re )
!
!----------------------------------------------------------------------
!::reflection
!----------------------------------------------------------------------
      xran = random(0)
      if( xran.le.rn ) then
        evel = erfl
!
!--cosine distribution
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
        weema(iw) = weema(iw) + weit*evel
        goto 200
      endif
!
!----------------------------------------------------------------------
!::emission of H2
!----------------------------------------------------------------------
 100  continue
      ! ntmole change
      weit = weit*0.5d0
      tmas = rmas(igas)*2.0d0
      call set_velocity(iw)
      is_atom_or_mole = 1
!
!----------------------------------------------------------------------
!::return
!----------------------------------------------------------------------
 200  continue
!
!::debug write
      if( ntrc.eq.1 .and. vel.ne.0.0d0 ) then
        ztim1 = 0.0d0
        write(n6,'(2x,"new velocity after ntrefl")')
        write(n6,602) istp
     >   ,ztim1,xpos,ypos,velx/vel,vely/vel,vel,evel,icpo,ixpo,iypo
     >   ,ilon,iptl
      endif
 602  format(2x,i5,1pe17.9,0p4f22.18,1p2e17.9,i6,2i5,i6,i8)
!
      return
      end
