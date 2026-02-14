!**********************************************************************
      subroutine ntwcon
!**********************************************************************
!
!::wtion ==> wreg(kr)
!::wssn  ==> swreg(kr)
!
!  kr  soldor : odiv(1), sol(2), idiv(3), oprv(4), iprv(5),
!               main(6), core(7)
!      other  : void(8), error(0)
!
!   normalization factor
!      dotn/wtot/volm(ic) ==> dotn/swnrm/volm(ic)       03/08/01
!
!       sum(wtion) in folow .ne. sum(wssn) in track-length
!          wtion = dlt-weit due to ionization
!          wssn  = weit*dtim*Ne*<sigv>ion      confirm  05/01/28
!
!----------------------------------------------------------------------
      use cntcom, only : mrgn, ipmp, ncmax2, ngas, npmp, npwl2
      use cntwcn, only : swabs, swerr, swion, swnrm, swpmp, swreg, swsum
     >    , wabs, wend, werr, wion, wnrm, wpmp, wreg, wssn, wsum, wtion
     >    , wwal
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer   nsiz, ig, ic, ir, iw, i
      integer   ig, ic, ir, iw, i
!
!----------------------------------------------------------------------
!::monte carlo event (wtion)
!----------------------------------------------------------------------
      do ir = 0, 10
      wreg(ir) = 0.0d0
      enddo
      do ig = 1, ngas
      do ic = 1, ncmax2
      ir   = mrgn(ic)
      wreg(ir) = wreg(ir) + wtion(ic,ig)
      enddo
      enddo
!
      wabs = 0.0d0
      do iw = 1, npwl2
      wabs = wabs + wwal(iw)
      enddo
!
      wpmp = 0.0d0
      do i = 1, npmp
      iw = ipmp(i)
      wpmp = wpmp + wwal(iw)
      enddo
!
      wion = wreg(1)+wreg(2)+wreg(3)+wreg(4)+wreg(5)
     >      +wreg(6)+wreg(7)
      werr = wend(2)+wend(3)+wend(4)+wend(5)
     >      +wreg(0)+wreg(8)+wreg(9)
      wsum = wion + wabs + werr
      wnrm = wion + wabs
!
!----------------------------------------------------------------------
!::source term (wssn)
!----------------------------------------------------------------------
      do ir = 0, 10
      swreg(ir) = 0.0d0
      enddo
      do ig = 1, ngas
      do ic = 1, ncmax2
      ir = mrgn(ic)
      swreg(ir) = swreg(ir) + wssn(ic,ig)
      enddo
      enddo
!
      swabs = wabs
      swpmp = wpmp
!
      swion = swreg(1)+swreg(2)+swreg(3)+swreg(4)+swreg(5)
     >       +swreg(6)+swreg(7)
      swerr = wend(2)+wend(3)+wend(4)+wend(5)
     >       +swreg(0)+swreg(8)+swreg(9)
      swsum = swion + swabs + swerr
      swnrm = swion + swabs
!
!::debug write (%2008/08/10)
!      write(n6,'(2x,"itim = ",i7,"  wion,wabs,werr    =",1p3e14.6,
!     >   "  wnrm  =",1pe14.6)') itim,wion,wabs,werr,wnrm
!      write(n6,'(2x,"       ",7x,"  swion,swabs,swerr =",1p3e14.6,
!     >   "  swnrm =",1pe14.6)') swion,swabs,swerr,swnrm
!      write(n6,'(2x,"swreg =",1p15e14.6)') swreg
!
!::?? (information from SGI)
!      if( wpmp.ne.wabs ) then
!      write(n6,'(/2x,"ntwcon   wpmp.ne.wabs  ",1p2e12.3)') wpmp,wabs
!      call wexit("ntwcon","wpmp.ne.wabs")
!      endif
!
      return
      end
