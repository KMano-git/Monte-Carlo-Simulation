!***********************************************************************
      subroutine imp_print(typ)
!***********************************************************************
      use cimcom, only : hgrpw, hregw, i6_pcnt, i6_pcon, i6_wrad, igemt
     >    , igexh, igpmp, igtbl, igten, igtst, itmmc, npmax, sitmz
     >    , scpum, stimz, tmmcz
      use cimden, only : twrd
      implicit none

!::argument
      character*(*), intent(in) :: typ

!::local variables
      real(8) :: tot_regw, wr_tot, wr_reg(10)
      integer :: irg, i6

      select case(typ)

!::pcon
      case("pcon")
      if( i6_pcon > 0 ) then
      tot_regw = 0.0d0
      do irg = 1, 8
        tot_regw = tot_regw + hregw(irg)
      enddo

      i6 = i6_pcon
      write(i6,'(1x,i6,1pe12.4,1p2e10.2,1p6e11.3,1p9e11.3)')
     >  itmmc, tmmcz, dble(npmax), scpum,
     >  hgrpw(igtst), hgrpw(igten), hgrpw(igtbl),
     >  hgrpw(igemt), hgrpw(igexh), hgrpw(igpmp),
     >  tot_regw, hregw(1:8)
      endif

!::pcnT
      case("pcnT")
      if( i6_pcnT > 0 ) then
      tot_regw = 0.0d0
      do irg = 1, 8
        tot_regw = tot_regw + hregw(irg)
      enddo

      i6 = i6_pcnT
      write(i6,'(1x,i6,1pe12.4,1p2e10.2,1p6e11.3,1p9e11.3)')
     >  int(sitmz+1.0d-5), stimz, dble(npmax), scpum,
     >  hgrpw(igtst), hgrpw(igten), hgrpw(igtbl),
     >  hgrpw(igemt), hgrpw(igexh), hgrpw(igpmp),
     >  tot_regw, hregw(1:8)
      endif

!::wrad
      case("wrad")
      if( i6_wrad > 0 ) then
      i6 = i6_wrad
      call intgsc(twrd,wr_tot,wr_reg)
      write(i6,'(1x,i6,1p2e12.4,1p9e12.4)')
     >  int(sitmz+1.0d-5), stimz, dble(npmax), wr_tot, wr_reg(1:8)
      endif

      case default
      end select

      return
      end
