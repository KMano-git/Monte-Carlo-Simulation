!***********************************************************************
      subroutine imstace(ipnm,iwem,x0,y0,z0,vlx,vly,vlz)
!***********************************************************************
!
!      ipnm (i)  number of sample particle
!      iwem (o)  number of wall segment  (flux when iwem = 0 )
!      x0   (o)  x-coordinate of emitted particle from wall
!      y0   (o)  y-coordinate of emitted particle
!      z0   (o)  z-coordinate of emitted particle
!      vx   (o)  x-velocity
!      vy   (o)  y-velocity
!      vz   (o)  z-velocity
!
!-----------------------------------------------------------------------
      use cimcom, only : amz, icbz, nbz, prbz, xpbz, ypbz
      use cntcom, only : bvx, bvy, bvz, fnor, tnor
      use cntpls, only : temi, vflw
      use cphcns, only : cev, cpi
      implicit none
!
!::argument
      integer, intent(in)  :: ipnm ! dummy
      integer, intent(out) :: iwem
      real(8), intent(out) :: x0, y0, z0, vlx, vly, vlz
!
!::local variables
      real*8  xran
      integer i, ic
      real*8  alf, csa, sna, psr, psz
      real*8  vlf
      integer ivx, ivy, ivz
      real*8  zti, zvz
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
!      alf = 2.0d0*cpi*0.721d0   ! confirm  05/02/02
      csa = cos(alf)
      sna = sin(alf)
      xran = random(0)
      psr = xpbz(iwem) + xran*(xpbz(iwem+1)-xpbz(iwem))
      psz = ypbz(iwem) + xran*(ypbz(iwem+1)-ypbz(iwem))
      x0  = psr*csa
      y0  = psr*sna
      z0  = psz
!
!::new velocity (ez1:normal direction of wall)
!-----------------------------------------------
! Gaussian distribution
      ic  = icbz(iwem)
      zti  = temi(ic)
      vlf  = vflw(ic,1)
      zvz = sqrt(zti*cev/amz)
      ivx = int(fnor*random(0) + 1.0d0)
      ivy = int(fnor*random(0) + 1.0d0)
      ivz = int(fnor*random(0) + 1.0d0)
      vlx = zvz*tnor(ivx) + vlf*bvx(ic)
      vly = zvz*tnor(ivy) + vlf*bvy(ic)
      vlz = zvz*tnor(ivz) + vlf*bvz(ic)
!-----------------------------------------------

      return
      end
