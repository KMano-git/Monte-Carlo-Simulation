!**********************************************************************
      subroutine ntfolw_teff(rst,zst,ren,zen,icst,icen,iwen,dlen)
!**********************************************************************
!
!      folow sample particle
!
!         istp : event of a particle
!         ievt : event in a cell                 !  2003/01/20
!
!::need many steps to hit wall in vacume region
!::need more steps to include EL-collision
!x    data lsmax/6000/  !   <=== lsmax/2000/
!x    when trace, istp > lsmax  goto 550
!
!   if( weit.le.0.0d0 ) goto 600  ! Fatal error  find 03/07/30
!
!::birth near O-dplate and hit O-dplate   ~1/2  O.K.  2011/01/18
!
!----------------------------------------------------------------------
      use cimcom, only : isw2tr, iswtr, lpcm
      use cntcom, only : evel, icpo, icwl, ievt, ihwl, ilon, iptl, istp
     >    , ixpo, iypo, lsmax, ltrc, migx, migy, weit, vel, velx, vely
     >    , velz, weitb, xpos, ypos, zpos, ncmax, mseg
      use cntwcn, only : wehta, wend, wfhta, whta, wwal
      use csonic, only : itim
      use cunit,  only : mype, n6
      use mod_externalgrid
      implicit none
!
!::argument
      real(8), intent(out) :: rst, zst, ren, zen, dlen
      integer, intent(out) :: icst, icen, iwen
!
!::local variables
      real*8  ztim1, ztim2, wtim, wpox, wpoy
      real*8  zrcy, weit0, weitd
      real*8  xst, yst, xen, yen, wpor, wpoz, xp1, yp1, zp1
      integer  icpn, ker, iw, ih, kmax, ic_tmp
      integer, parameter :: loop_max = 1000000
      integer  i, j

      xen = 0.0
      xst = 0.0
      yen = 0.0
      yst = 0.0
      zen = 0.0
      zst = 0.0

      if( weit.le.0.0d0 ) goto 570
!
      ievt = 0
!
      vel = sqrt(velx**2+vely**2+velz**2)
      xst = xpos
      yst = ypos
      rst = sqrt(xpos**2+ypos**2)
      zst = zpos
      icst = icpo
!
      lpcm = 0    ! loop  lpcm = lt < ltmax  2012/01/19
      iswtr = 0   ! debug
      isw2tr = 0
!
!----------------------------------------------------------------------
!::start after lounch or charge exchange event
!----------------------------------------------------------------------
      do j = 1, loop_max
      ztim1 = 0.0d0
      ztim2 = 0.0d0
!
!----------------------------------------------------------------------
!::loop of track
!----------------------------------------------------------------------
      do i = 1, loop_max
      istp = istp + 1
      ievt = ievt + 1
!
!--debug
      wtim = ztim1
      wpox = xpos + velx*wtim
      wpoy = ypos + vely*wtim
      wpoz = zpos + velz*wtim
      wpor = sqrt(wpox**2+wpoy**2)
      ixpo = migx(icpo)
      iypo = migy(icpo)
!
      if( ltrc.eq.1 ) then
      write(300,602) istp
     > ,wtim,wpor,wpoz,velx/vel,vely/vel,velz/vel,evel,icpo,ixpo,iypo,
     >  ilon,iptl,mype,"ntfolw_teff"
      endif
 602  format(2x,i5,1p7e12.4,6i6,2x,a)
!
!--spend time in cell
      icpn = icpo
      if(.not. use_exdata) then ! ordinary mesh
        call tmtrac(xpos,ypos,zpos,velx,vely,velz,icpn,ilon,
     >    ztim2,wtim,xp1,yp1,zp1,ker)
      else ! the mesh made by Gmesh and external tool
        if(icpn.le.ncmax) then ! SOL region
          kmax = mseg(icpn)
          call tmtrac_ex(xpos,ypos,zpos,velx,vely,velz,icpn,ilon,
     >      ztim2,wtim,xp1,yp1,zp1,ker,kmax)
        elseif(icpn .le. ncmax+vac_ele_size) then ! vacume region
          ic_tmp = icpn-ncmax
          kmax = mseg_vacume(ic_tmp)
          call tmtrac_grid(xpos,ypos,zpos,velx,vely,velz,icpn,ilon
     >      ,ztim2,wtim,xp1,yp1,zp1,ker
     >      ,vac_element,vac_ele_size,kmax
     >      ,vac_grid_size,vac_grid_x,vac_grid_y,ic_tmp)
        elseif(icpn .le. ncmax+vac_ele_size+pri_ele_size) then ! private region
          ic_tmp = icpn-(ncmax+vac_ele_size)
          kmax = mseg_pri(ic_tmp)
          call tmtrac_grid(xpos,ypos,zpos,velx,vely,velz,icpn,ilon
     >      ,ztim2,wtim,xp1,yp1,zp1,ker
     >      ,pri_element,pri_ele_size,kmax
     >      ,pri_grid_size,pri_grid_x,pri_grid_y,ic_tmp)
        else ! sub-diverter
          ic_tmp = icpn-(ncmax+vac_ele_size+pri_ele_size)
          kmax = mseg_subdiv(ic_tmp)
          call tmtrac_grid(xpos,ypos,zpos,velx,vely,velz,icpn,ilon
     >      ,ztim2,wtim,xp1,yp1,zp1,ker
     >      ,subdiv_cell,cell_size,kmax
     >      ,grid_size,vgx_EX,vgy_EX,ic_tmp)
        endif
      endif
      if( ker.ne.0 ) goto 560
!
      ztim1 = wtim
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
!
      ievt = 0
      enddo
      call wexit("ntfolw_teff","too many loop")
!
!----------------------------------------------------------------------
!::hit the wall
!----------------------------------------------------------------------
 300  continue
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
      zpos = zpos + velz*ztim1
!
      xen = xpos
      yen = ypos
      ren = sqrt(xpos**2+ypos**2)
      zen = zpos
      icen = icpo
      iwen = iw
!
!::debug
      if( ltrc.eq.1 )then
      write(300,602) istp
     > ,wtim,ren,zen,velx/vel,vely/vel,velz/vel,evel,icpo,ixpo,iypo,
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
      enddo
      call wexit("ntfolw_teff","too many loop refrection")
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
!--too many step
 550  continue
      wend(4) = wend(4) + weit
      if( ltrc.eq.1 ) write(n6,'(2x,"too many step  ",2i8,2i4)')
     >   istp,ievt,ixpo,iypo
      call wexit("ntfolw_teff","too many step")
!
!--error
 560  continue
      wend(5) = wend(5) + weit
!
      write(n6,'(2x,"Error ntfolw  mype =",i5,"  itim =",i6,"  iptl ="
     > ,i6,"  ic =",3i5,"  weit,werr =",1p2e12.3)')
     >      mype,itim,iptl,icpo,migx(icpo),migy(icpo),weit, wend(5)
      return
!
!--normal end (weit = 0.0)
 570  continue
!
!----------------------------------------------------------------------
!::return
!----------------------------------------------------------------------
 600  continue
      dlen = sqrt((xen-xst)**2+(yen-yst)**2+(zen-zst)**2)
!
      return
      end