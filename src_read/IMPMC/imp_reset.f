!**********************************************************************
      subroutine imp_reset
!**********************************************************************
!)
!) ion
!)  tt = 0.5d0 => 0.0d0
!)
!) neutral
!)  tt  = 0.5d0 => 0.0d0
!)  sx0 = sx0 + svx*stm
!)  sy0 = sy0 + svy*stm
!)  sz0 = sz0 + svz*stm
!)  (svx,svy,svz)  steady velocity
!)  srn = -dlog(random)
!)  sint = 0.0
!)  stm = 0.0
!)  stb = 0.0
!)              see IMPV5/imemt_wall.f
!)
!)  monte/ntfolw.f
!::collision occured
!) 400 continue
!)  zdlt  = (zincx-zint2)/trct
!)  ztim1 = ztim2 + zdlt
!)  new start position
!)  xpos = xpos + velx*ztim1
!)  ypos = ypos + vely*ztim1
!)
!)    stm(ip) = ztim1
!)
!)--------------------------------------------------------------------
      use cimcom, only : il, ir, is, npmax, tt
      use cimntl, only : sint, srn, stb, stm, svx, svy, svz, sx0, sy0
     >    , sz0
      use cunit,  only : n6
      implicit none
!
      integer :: ip, ic, ko, ic1, ko1, ier
      real(8) :: x0, y0, z0, r0, vlx, vly, vlz, ztm
      real(8) :: x1, y1, z1, r1, xran
! function
      integer    imox, imoy
      real(8)    random

      write(n6,'("imp_reset  npmax =",i5)') npmax
      write(n6,'(4x,a)') "tt = 0.5 => 0.0  New position of neutral"//
     >     " when ko1 < 0  particle on a wall-segment"

!::start point for neutral particle
!:: for ions, unchanged excep. tt
      ier = 0
      do ip = 1, npmax
      if( is(ip) == 0 ) then
        x0  = sx0(ip)
        y0  = sy0(ip)
        z0  = sz0(ip)
        r0  = sqrt(x0**2+y0**2)
        vlx = svx(ip)
        vly = svy(ip)
        vlz = svz(ip)
        ztm = stm(ip)
        ztm = 0.9999d0*ztm  ! wall point => interior point
        x1  = x0 + ztm*vlx
        y1  = y0 + ztm*vly
        z1  = z0 + ztm*vlz
        r1  = sqrt(x1**2+y1**2)

!::ic number to contain (r1,z1) point
        ic  = ir(ip)
        call mchkin(r1,z1,ic,ko)
        if( ko /= 0 ) then
          ier = ier + 1
          ic1 = ic
          if(ic.eq.0) write(n6,*) 'Err ic=0, ip=', ip
          call mfndic(r1,z1,ic1,ko1)
          write(n6,'(6x,i5,4f12.5,i7,i6,i5,i4," => ",i7,i6,i5,i6)')
     >      ip, r0, z0, r1, z1, ic, imox(ic), imoy(ic), ko,
     >                          ic1, imox(ic1), imoy(ic1), ko1
          ic = ic1
          ko = ko1
        endif

!::set common variables
        tt(ip)  = 0.0d0
        ir(ip)  = ic
        il(ip)  = 0
        is(ip)  = 0
        sx0(ip) = x1
        sy0(ip) = y1
        sz0(ip) = z1
        svx(ip) = vlx
        svy(ip) = vly
        svz(ip) = vlz

!::unchange
!::  pemt(ip), wght(ip), wstp(ip)

        xran = random(0)
        srn(ip) = -dlog(xran)
        sint(ip) = 0.0d0
        stm(ip) = 0.0d0      ! <== ztim1
        stb(ip) = 0.0d0      ! <== birth time
      else
        tt(ip) = 0.0d0
      endif
      enddo

      if( ier > 0 ) then
      write(n6,'(7x,"ip",5x,"R0",10x,"Z0",10x,"R1",10x,"Z1",10x,
     >   "ic",4x,"ix",3x,"iy",2x,"ko",7x,
     >   "ic1",4x,"ix1",2x,"iy1",3x,"ko1")')
      write(n6,'(4x,"ier = ",i5,"  invalid ic")') ier
      endif
!      write(n6,'(4x,"ier = ",i5,"  invalid ic")') ier

      return
      end

!**********************************************************************
      subroutine debg_regnw(cmsg)
!**********************************************************************
      use cimcom, only : ien, ir, npmax, wght
      use cntcom, only : mrgn
      use cunit,  only : n6
      implicit none

      character(*), intent(in) :: cmsg
      real(8) :: zregp(10), zregw(10)
      integer :: ip, npcal, ic, irg

      write(n6,'(2x,"debg_regn ---",2x,a)') trim(cmsg)

      npcal = 0
      zregp(1:10) = 0.0d0
      zregw(1:10) = 0.0d0

      do ip = 1, npmax
        if( ien(ip) /= 0 ) cycle
        ic = ir(ip)
        irg = mrgn(ic)
        npcal = npcal + 1
        zregp(irg) = zregp(irg) + 1.0d0
        zregw(irg) = zregw(irg) + wght(ip)
      enddo

      write(n6,'(2x,"debg_regn  npmax,npcal =",2i6)') npmax, npcal
      write(n6,'(2x,"debg regn  zregp =",1p9e12.4)') zregp(1:8)
      write(n6,'(2x,"debg regn  zregw =",1p9e12.4)') zregw(1:8)

      return
      end
