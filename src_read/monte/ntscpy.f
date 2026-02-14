!**********************************************************************
      subroutine ntscpy
!**********************************************************************
!
!      select source terms    track or collision estimator
!
!-----------------------------------------------------------------------
      use cntcom, only : flxin, lsrmd
      use cntscg, only : sn_mc, sn_mt, sp_mc, sp_mt, we_mc, we_mt, wi_mc
     >    , wi_mt
      use cntwcn, only : wcest, wflx, wssn, wssp, wswe, wswi
      use csize,  only : ndgs, ndmc
      implicit none
!
!::local variables
! modified 2/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer  ig, ic, ia, nsiz, ir, iw, i
!ik   real*8   zdotn, zfc, zn0, ze0, zv0, zmas, zvl
      integer  ig, ic
!
!::track
      if( lsrmd.ne.2 ) then
      wcest = "track"
      do ig = 1, ndgs
      do ic = 0, ndmc
      wssn(ic,ig) = sn_mt(ic,ig)
      wssp(ic,ig) = sp_mt(ic,ig)
      enddo
      enddo
      do ic = 0, ndmc
      wswe(ic) = we_mt(ic)
      wswi(ic) = wi_mt(ic)
      enddo
!
!::collision
      else
      wcest = "colli"
      do ig = 1, ndgs
      do ic = 0, ndmc
      wssn(ic,ig) = sn_mc(ic,ig)
      wssp(ic,ig) = sp_mc(ic,ig)
      enddo
      enddo
      do ic = 0, ndmc
      wswe(ic) = we_mc(ic)
      wswi(ic) = wi_mc(ic)
      enddo
      endif
!
!::flux
      wflx = flxin
!
      return
      end
