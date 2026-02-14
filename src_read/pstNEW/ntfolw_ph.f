!**********************************************************************
      subroutine ntfolw_ph(xst,yst,xen,yen,dlen,icon)
!**********************************************************************
!
!      folow sample particle
!
!         istp : event of a particle
!         ievt : event in a cell                 !  2003/01/20
!
!::need many steps to hit wall in vacume region
!::need more steps to include EL-collision
!
!   if( weit.le.0.0d0 ) goto 600  ! Fatal error  find 03/07/30
!
!::poloidal plane: x-y, toroidal direction: z
!::velocity (vx,vy,vz)
!::position (x,y,z)
!        x (r) = x0 + velx*time
!        y (z) = y0 + vely*time
!
!----------------------------------------------------------------------
      use cntcom, only : evel, icpo, icwl, ievt, ihwl, ilon, iptl, istp
     >    , ixpo, iypo, lsmax, ltrc, migx, migy, vel, velx, vely, velz
     >    , weit, weitb, xpos, ypos
      use cntwcn, only : wehta, wend, wfhta, whta, wwal
      use csonic, only : itim
      use cunit,  only : mype, n6
      implicit none
!
!::argument
      real(8), intent(out) :: xst, yst, xen, yen, dlen
      integer, intent(out) :: icon
!
!::local variables
      real*8  zint1, ztim1, ztim2, wtim, wpox, wpoy, zdlt
      real*8  zrcy, weit0, weitd, zvx1, zvy1, zvx2, zvy2
      integer  icpn, ker, iw, ih
      real*8  zvlp
!
      if( weit.le.0.0d0 ) goto 570
!
      ievt = 0
!
!::poloidal velocity
      zvlp = sqrt(velx**2+vely**2)/sqrt(velx**2+vely**2+velz**2)
      xst = xpos
      yst = ypos
      icon = 99
!
!----------------------------------------------------------------------
!::start after lounch or charge exchange event
!----------------------------------------------------------------------
 100  continue
      zint1 = 0.0d0
      ztim1 = 0.0d0
      ztim2 = 0.0d0
!
!----------------------------------------------------------------------
!::loop of track
!----------------------------------------------------------------------
 200  continue
      istp = istp + 1
      ievt = ievt + 1
!
!--debug
      wtim = ztim1
      wpox = xpos + velx*wtim
      wpoy = ypos + vely*wtim
      ixpo = migx(icpo)
      iypo = migy(icpo)
!
      if( ltrc.eq.1 ) then
      write(n6,602) istp
     > ,wtim,wpox,wpoy,velx/vel,vely/vel,velz/vel,evel,icpo,ixpo,iypo,
     >  ilon,iptl,mype,"ntfolw_ph"
      endif
 602  format(2x,i5,1p7e12.4,6i6,2x,a)
!
!--spend time in cell
      icpn = icpo
      call sptim(xpos,ypos,velx,vely,icpn,ilon,wtim,wpox,wpoy,ker)
      if( ker.ne.0 ) goto 560
!
      ztim1 = wtim
      zdlt  = ztim1 - ztim2
      zint1 = zint1 + zdlt*zvlp
!
!--scoreing density & temperature
!--track estimator
!
!--next cell
      ztim2 = ztim1
      icpo  = icpn
!
!--hit the wall
      if( ilon.lt.0 ) goto 300
!
!--too many loop
      if( ievt.gt.lsmax ) goto 550
      if( istp*ltrc.gt.lsmax ) goto 550
!
      ievt = 0
      goto 200
!
!----------------------------------------------------------------------
!::hit the wall
!----------------------------------------------------------------------
 300  continue
      if( ltrc.eq.1 ) write(n6,'(2x,"hit the wall")')
      iw = -ilon
!
!--recycling & wall temperature
!::gate
!
!--hit point
 320  continue
      icpo = icwl(iw)
      ixpo = migx(icpo)
      iypo = migy(icpo)
      xpos = xpos + velx*ztim1
      ypos = ypos + vely*ztim1
!
      xen = xpos
      yen = ypos
!
!::debug
      if( ltrc.eq.1 )then
      write(n6,602) istp
     > ,wtim,wpox,wpoy,velx/vel,vely/vel,velz/vel,evel,icpo,ixpo,iypo,
     >  ilon,iptl,mype,"hit wall"
      endif
!
!::neutral flux onto the wall (lopht)
      wfhta(iw) = wfhta(iw) + weit
      wehta(iw) = wehta(iw) + weit*evel
!
!--new weit
      zrcy  = 0.0d0
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
      if( weit.le.0.0d0 ) goto 570
!
!--reflection
!--new velocity of mole
!--reflection/absorption
      goto 100
!
!----------------------------------------------------------------------
!::collision occured
!----------------------------------------------------------------------
 400  continue
      call wexit("ntfolw_ph","collision occured")
!
!----------------------------------------------------------------------
!::terminate
!----------------------------------------------------------------------
!
!--ionization
 520  continue
      wend(1) = wend(1) + weitb
      if( ltrc.eq.1 ) write(n6,'(2x,"ionize")')
      goto 600
!
!--absorption
 530  continue
      wend(2) = wend(2) + weit
      if( ltrc.eq.1 ) write(n6,'(2x,"absorption")')
      goto 600
!
!--too many hit wall
 540  continue
      wend(3) = wend(3) + weit
      if( ltrc.eq.1) write(n6,'(2x,"too many hit wall")')
      goto 600
!
!--too many step
 550  continue
      wend(4) = wend(4) + weit
      if( ltrc.eq.1 ) write(n6,'(2x,"too many step  ",2i8,2i4)')
     >   istp,ievt,ixpo,iypo
      goto 600
!
!--error
 560  continue
      wend(5) = wend(5) + weit
!
      write(n6,'(2x,"Error ntfolw  mype =",i5,"  itim =",i6,"  iptl ="
     > ,i6,"  ic =",3i5,"  weit,werr =",1p2e12.3)')
     >      mype,itim,iptl,icpo,migx(icpo),migy(icpo),weit, wend(5)
      goto 600
!
!--normal end (weit = 0.0)
 570  continue
      icon = 0
!
!----------------------------------------------------------------------
!::return
!----------------------------------------------------------------------
 600  continue
      if( ltrc.eq.1 ) then
      write(n6,'(2x,"ztim1 =",1pe12.3,"  zint1 =",1pe12.3)')
     >  ztim1, zint1
      endif
!
      dlen = sqrt((xen-xst)**2+(yen-yst)**2)
      zvx1 = velx/zvlp
      zvy1 = vely/zvlp
      zvx2 = (xen-xst)/dlen
      zvy2 = (yen-yst)/dlen
!
      if( icon.ne.0 .or. mod(iptl,1000).eq.0 ) then
      write(n6,'(2x,i6,i4,1p10e12.3,2x,i7,2i5,i7)') iptl, icon,
     >  xst, yst, xen, yen, dlen, zint1,
     >  zvx1, zvx2, zvy1, zvy2, icpo, ixpo, iypo, ilon
      endif
!
      return
      end
