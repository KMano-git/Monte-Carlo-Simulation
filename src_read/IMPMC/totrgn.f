!***********************************************************************
      subroutine totrgn(var,tot,tot_rg,tot_it)
!***********************************************************************
!
!      tot : region 1-6  except for core (7)
!                            2013/01/18  K. Shimizu
!
!-----------------------------------------------------------------------
      use csize,  only : ndmc, ndy
      use cplmet, only : icel, itmpe, itmps, itsls, jcel, jtmax, jtmin
      use cntcom, only : mcel, mrgn, ncmax, volm
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8  var(ndmc)
!ik   real*8  tot, tot_rg(10), tot_it(ndy)
      real(8), intent(in)  :: var(ndmc)
      real(8), intent(out) :: tot, tot_rg(10), tot_it(ndy)
!
!::local variables
      integer it, jt, jts, jte, j, i, ic, irg
      real*8  zsum
!
      tot = 0.0d0
      do j = 1, 10
      tot_rg(j) = 0.0d0
      enddo
!
      do ic = 1, ncmax
      irg = mrgn(ic)
        if( irg.ne.7 ) then
      tot = tot + var(ic)*volm(ic)
      endif
      tot_rg(irg) = tot_rg(irg) + var(ic)*volm(ic)
      enddo
!
      do it = itsls, itmpe
      jts = jtmin(it)
      jte = jtmax(it)
      if( it.ge.itmps ) then
      jts = jts + 1
      jte = jte - 1
      endif
      zsum = 0.0d0
      do jt = jts, jte
      j  = jcel(jt,it)
      i  = icel(jt,it)
      ic = mcel(j,i)
      if( ic.eq.0 ) cycle
      zsum = zsum + var(ic)*volm(ic)
      enddo
      tot_it(it) = zsum
      enddo
!
      return
      end
