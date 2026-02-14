!***********************************************************************
      subroutine imstaw(ipnm,iwem,x0,y0,z0,v0,vlx,vly,vlz)
!***********************************************************************
!
!      ipnm (i)  number of sample particle
!      iwem (o)  number of wall segment  (flux when iwem = 0 )
!      x0   (o)  x-coordinate of emitted particle from wall
!      y0   (o)  y-coordinate of emitted particle
!      z0   (o)  z-coordinate of emitted particle
!      v0   (i)  velocity
!      vx   (o)  x-velocity
!      vy   (o)  y-velocity
!      vz   (o)  z-velocity
!
!-----------------------------------------------------------------------
      use cimcom, only : csbz, iebz, isbz, nbz, prbz, snbz, icbz
     > , ikbz 
      use cntcom, only : fnfi, fnth, tcfi, tcth, tsfi, tsth, xpnt, ypnt
     > , ncmax, mseg, mgrd, mknd, igtw
      use cphcns, only : cpi
      use mod_externalgrid
      implicit none
!
!::argument
      integer, intent(in)  :: ipnm ! dummy
      integer, intent(out) :: iwem
      real(8), intent(in)  :: v0
      real(8), intent(out) :: x0, y0, z0, vlx, vly, vlz
!
!::local variables
      real*8  xran
      integer i, ith, ifi
      integer ipmg, ipmg1, ic
      real*8  alf, csa, sna, psr, psz
      real*8  vlr1, vlf1, vlz1, csb, snb, vlr, vlf
      real*8  xgrid, xgrid_1, ygrid, ygrid_1
! function
      real(8)    random
!
!::wall segment
      xran = random(0)
      do i = 1, nbz
        iwem = i
        if( prbz(i).ge.xran ) goto 10
      enddo
      iwem = nbz
 10   continue
!
!::start position
      alf = 2.0d0*cpi*0.0d0    ! no depend on this value
!xx   alf = 2.0d0*cpi*0.721d0  ! confirm  05/02/02
      csa = cos(alf)
      sna = sin(alf)
      xran = random(0)
      ipmg  = isbz(iwem)
      ipmg1 = iebz(iwem)
      ic = icbz(iwem)
      if(ic.le.ncmax .or. .not.use_exdata) then
        xgrid   = xpnt(ipmg)
        xgrid_1 = xpnt(ipmg1)
        ygrid   = ypnt(ipmg)
        ygrid_1 = ypnt(ipmg1)
      elseif(ic .le. ncmax+vac_ele_size)then
        xgrid   = vac_grid_x(ipmg)
        xgrid_1 = vac_grid_x(ipmg1)
        ygrid   = vac_grid_y(ipmg)
        ygrid_1 = vac_grid_y(ipmg1)
      elseif(ic .le. ncmax+vac_ele_size+pri_ele_size)then
        xgrid   = pri_grid_x(ipmg)
        xgrid_1 = pri_grid_x(ipmg1)
        ygrid   = pri_grid_y(ipmg)
        ygrid_1 = pri_grid_y(ipmg1)
      else
        call wexit("imstaw","invalid range for icbz")
      endif
!
      psr = xgrid + xran*(xgrid_1-xgrid)
      psz = ygrid + xran*(ygrid_1-ygrid)
      x0  = psr*csa
      y0  = psr*sna
      z0  = psz
!
!::new velocity (ez1:normal direction of wall)
      ith  = int(fnth*random(0)+1.0d0)
      ifi  = int(fnfi*random(0)+1.0d0)
!:: tsth:0~pi/2, tsfi:0~2pi
      vlr1 = v0*tsth(ith)*tcfi(ifi) ! x
      vlf1 = v0*tsth(ith)*tsfi(ifi) ! y
      vlz1 = v0*tcth(ith)           ! z (> 0)
!
      csb = csbz(iwem)
      snb = snbz(iwem)
      vlr = vlr1*csb - vlz1*snb
      vlz = vlr1*snb + vlz1*csb
      vlf = vlf1
!
      vlx  = vlr*csa - vlf*sna
      vly  = vlr*sna + vlf*csa
!
      return
      end
