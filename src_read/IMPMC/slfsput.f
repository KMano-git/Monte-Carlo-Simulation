!***********************************************************************
      subroutine slfsput(ic, ip, rflw, v0slf)
!***********************************************************************
      use cimcns, only : ftom, tomvel
      use cimcom, only : amz, is, vv, yfslf
      use cimntl, only : svx, svy, svz
      use cntpls, only : teme
      use cphcns, only : cev
      implicit none

      integer, intent(in)  :: ic, ip
      real(8), intent(out) :: rflw, v0slf

      integer  ix, kspphy
      real*8   zec, yc
      real*8   zvx, zvy, zvz, zvv
      real*8   xran
! function
      real(8)    random, sputphy
      integer    nsptin

      kspphy = nsptin("C")

      if(is(ip).eq.0)then
        zvx = svx(ip)
        zvy = svy(ip)
        zvz = svz(ip)
        zvv = zvx**2+zvy**2+zvz**2
        zec = 0.5d0*amz*zvv/cev
      else
        zvv = vv(ip)
        zec = 0.5d0*amz*zvv/cev + 3.0d0*teme(ic)*dfloat(is(ip))
      endif
!
!::   sputtering of normal incidence
      yc = sputphy( zec, 0.0d0, kspphy )
!
!::   angular effect of incidence
      yc = yc*yfslf
      rflw = yc

      if(rflw.eq.0.0d0)then
        v0slf = sqrt(zvv)
      else
        xran = random(0)
        ix = int(ftom*xran + 1.0d0)
        v0slf = tomvel(ix)
      endif

      return
      end
