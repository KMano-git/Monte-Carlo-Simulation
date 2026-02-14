!***********************************************************************
      subroutine imscor(ip,dtunt)
!***********************************************************************
!
!        latom = 0
!
!-----------------------------------------------------------------------
      use cimcom, only : denz, ien, il, ir, is, pemt, rr, temz, tt, vv
     >    , wght, zz
      use cimptl, only : wic0, wic1, wrr0, wzz0
      use mod_shexe, only : impmc_model
      implicit none
!
!::argument
      integer, intent(in) :: ip
      real(8), intent(in) :: dtunt
!
!::local variables for fndway
      integer ndmv; parameter (ndmv=201)
      integer nmv, icmv(ndmv)
      real*8  dtmv(ndmv)
!
!::local variables for scoreing
      integer ndam; parameter (ndam=201)
      integer ndsc; parameter (ndsc=201)
!
!::local variables
      real(8)  wrr1, wzz1, ftmv, dtmad, dt0, zwght
      integer iz, ic, i
      integer lon
!
      ic   = ir(ip)
      iz   = is(ip)
!
!::trace
      wrr1 = rr(ip)
      wzz1 = zz(ip)
      if( iz.eq.0 ) then
        wrr1 = wrr0
        wzz1 = wzz0
      endif
      call fndway(wrr0,wzz0,wrr1,wzz1,wic0,wic1,ndmv,nmv,icmv,dtmv)
      if( nmv.eq.0 ) then          ! error
        ien(ip) = 4
        return
      endif
      lon = 0
      if( wic1.le.0 ) lon = wic1
!
      ftmv = 0.0d0
      do i = 1, nmv
        ftmv = ftmv + dtmv(i)
      enddo
!
!::final state
      ic = icmv(nmv)
      dtmad  = dtunt*ftmv
      tt(ip) = tt(ip) + dtmad
      rr(ip) = wrr1
      zz(ip) = wzz1
      ir(ip) = ic
      is(ip) = iz
      if( lon.lt.0 ) then
        ien(ip) = 3
        il(ip)  = lon
      endif
!
!::scoreing (not change charge state)
      if( impmc_model == 0 ) then
        zwght = wght(ip)
      else
        zwght = pemt(ip)*wght(ip)
      endif
      do i = 1, nmv
      dt0   = dtunt*dtmv(i)
      ic = icmv(i)
      denZ(iz,ic) = denZ(iz,ic) + zwght*dt0
      temZ(iz,ic) = temZ(iz,ic) + zwght*dt0*vv(ip)
      enddo
!
      if( is(ip).eq.0 ) ien(ip) = 1  ! C+ ==> C0
      return
      end
