!***********************************************************************
      subroutine set_ntl
!***********************************************************************
!
!     PLS         NTL            IMP
!      flx     <==(N^0, T^0)     
!      N0, T0  ===========>    imp_pls
!
!----------------------------------------------------------------------
      use csize
      use cntcom
      use cntmnt
      use cntsrc
      use cplcom
      use cplmet
      use cimctl
      use csonic
      use cunit
      implicit none
!
!::local variables
      integer  lenx, i, mj, ic, ig, ir, sv_n6
      integer  ia, isrc, iflx
      character  cmsg*80, cneut*80, dsn*80

      integer :: tbic(4), k, ix, iy, imox, imoy
      integer, save :: lp = 0
   

      lp = lp + 1
      write(n6,'(/2x,"*** set_ntl ***  mygrp,mype =",2i4)') mygrp,mype
!
      sv_n6  = n6
      n6_src = n6
      i6_src = n6
      n6     = n6_src
      n6_src = 0
      i6_src = 0
!
!----------------------------------------------------------------------
!::source
!----------------------------------------------------------------------
!::neutral source      
!   selection of neutral source   see pltcnt/aapst2dp.f
!
!::clear common variables (plntsr_plt.f)
      do ic = 0, ndmc
      do ig = 1, ndgs
      tden0(ic,ig) = 0.0d0
      teng0(ic,ig) = 0.0d0
      tvlp0(ic,ig) = 0.0d0
      enddo
      do ig = 1, 2
      tdeng(ic,ig) = 0.0d0
      tengg(ic,ig) = 0.0d0
      enddo
      do ia = 1, ndsp
      tssn(ic,ia) = 0.0d0
      tssp(ic,ia) = 0.0d0
      enddo
      tswi(ic) = 0.0d0
      tswe(ic) = 0.0d0
      enddo
!
      do ir = 1, 10
      do ia = 1, ndsp
      trgsn(ir,ia) = 0.0d0
      trgsp(ir,ia) = 0.0d0
      enddo
      trgwi(ir) = 0.0d0
      trgwe(ir) = 0.0d0
      enddo     
!
!::see sub. plsorc
      sdotn = 0.0d0
      sdotn2= 0.0d0
      tflex = 0.0d0
      tflpm = 0.0d0
      tnpuf = 0.0d0
!
!::source terms   KSFUJI  0519  use nzmx in plcrad_mc
      call plsorc(2)  ! at  monte cell
!
!::temperature
      do ic = 1, ncmax2  ! include Vac-region  2014/04/04
      do ia = 1, nion
      ig = ia
      if( tden0(ic,ig).gt.0.0d0 ) then
      teng0(ic,ig) = teng0(ic,ig)/tden0(ic,ig)
      tvlp0(ic,ig) = tvlp0(ic,ig)/tden0(ic,ig)
      endif
      enddo
      if( tdeng(ic,1).gt.0.0d0 ) tengg(ic,1) = tengg(ic,1)/tdeng(ic,1)
      if( tdeng(ic,2).gt.0.0d0 ) tengg(ic,2) = tengg(ic,2)/tdeng(ic,2)
      enddo

 100  continue
      ic = 67
      write(n6,'(1x,a,1pe16.8,i4,i5,1p8e12.4, 1p5e12.4)')
     >  "#set_ntl", time, itim, ic, snflx(1:8),
     >   tN0(ic,1), tT0(ic,1), tV0(ic,1), tNg(ic,1), tTg(ic,1)
      write(n6,'(1x,a,1pe16.8,i4,i5,1p8e12.4, 1p5e12.4)')
     >  "#set_ntl(tden0)", time, itim, ic, snflx(1:8),
     >   tden0(ic,1), teng0(ic,1), tvlp0(ic,1), tdeng(ic,1), tengg(ic,1)

      call plprof(itsle-6)
      call plprof(itsle-4)
      call plprof(itsle-2)
      call plprof(itsle)
      call plprof(itmps)
      call plprof(29)
      call plprof(26)
      call plprof(6)
!
!::flux to wall
      call ntwflx
      call plwflx
!
      return
      end
