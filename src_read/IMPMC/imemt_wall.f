!***********************************************************************
      subroutine imemt_wall(ip)
!***********************************************************************
      use cimcns, only : ftom, tomvel
      use cimcom, only : emptl, icbz, ien, ikbz, il, ir, is, pemt, tt
     >    , v0emtz, wght, wstp
      use cimntl, only : sint, srn, stb, stm, svx, svy, svz, sx0, sy0
     >    , sz0
      use cntcom, only : mknd
      use cunit, only : n6
      use mod_externalgrid, only : use_exdata
      implicit none

!::arguments
      integer, intent(in) :: ip

!::local variables
      integer :: iwem
      real(8) :: xran
      integer :: ix
      real(8) :: x0, y0, z0, v0, vlx, vly, vlz, wgt
      integer :: ic, icx, icy, ik, ln
! function
      real(8)    random
      integer    imox, imoy

!----------------------------------------------------------------------
!::[1] sputtering at wall    imemit.f => imshft.f
!----------------------------------------------------------------------
!VV_5.1  v0emtZ => v0emt
      v0 = v0emtZ
      if( v0emtZ.eq.0.0d0 ) then
        xran = random(0)
        ix = int(ftom*xran + 1.0d0)
        v0 = tomvel(ix)
      endif

      call imstaw(ip,iwem,x0,y0,z0,v0,vlx,vly,vlz)

!::cell
      ic  = icbz(iwem)
      ik  = ikbz(iwem)
      ln  = mknd(ic,ik)

!:: debug check
      if(.not. use_exdata) then
        icx = imox(ic)
        icy = imoy(ic)
        if( icx.le.0 .or. icy.le.0 .or. ic.eq.0) then
          write(n6,'(2x,"ic,icx,icy=",3i6)') ic,icx,icy
          call wexit("imemt_wall",
     >     "icx.le.0 .or. icy.le.0 .or. ic.eq.0")
        endif
      endif

!::weit
      wgt = 1.0d0

!::set common variables
      tt(ip)  = 0.0d0
      ir(ip)  = ic
      il(ip)  = ln
      is(ip)  = 0
      sx0(ip) = x0
      sy0(ip) = y0
      sz0(ip) = z0
      svx(ip) = vlx
      svy(ip) = vly
      svz(ip) = vlz
      wght(ip) = wgt
      wstp(ip) = 0.0d0
      ien(ip)  = 0
      pemt(ip) = emptl

      xran = random(0)
      srn(ip) = -dlog(xran)
      sint(ip) = 0.0d0
      stm(ip) = 0.0d0      ! <== ztim1
      stb(ip) = tt(ip)     ! <== birth time
      return
      end
