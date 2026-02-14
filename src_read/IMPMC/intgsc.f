!***********************************************************************
      subroutine intgsc(var,tot,tot_rg)
!***********************************************************************
!)
!)      integration in each section
!)       (1:odv/2:sol/3:idv/4:opv/5:ipv/6:edg/7:cor/8:vac)
!)
!)      tot : region 1-6  except for core (7)
!)                            2013/01/18  K. Shimizu
!)
!)      more occurate cal.  irg = mrgnp(ic)
!)       But no data of mrgnp in IMPMC code
!)
!)       totrgn => intgsc   at the some future date
!)
!-----------------------------------------------------------------------
      use cntcom, only : mrgn, ncmax2, volm
      use csize,  only : ndmc
      implicit none

!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/07/04
!ik   real(8) ::  var(ndmc)
!ik   real(8) ::  tot, tot_rg(10)
      real(8), intent(in)  :: var(ndmc)
      real(8), intent(out) :: tot, tot_rg(10)
!
!::local variables
      integer ic, irg

      tot_rg(1:10) = 0.0d0
      do ic = 1, ncmax2  ! incl. vac-region
        irg = mrgn(ic)
        tot_rg(irg) = tot_rg(irg) + var(ic)*volm(ic)
      enddo

      tot = 0.0d0
      do irg = 1, 6
        tot = tot + tot_rg(irg)
      enddo

      return
      end
