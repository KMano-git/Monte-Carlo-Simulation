!***********************************************************************
      subroutine sgvcxr_set
!***********************************************************************
!::CXR rate  <sv>(Er)  [m3/s]    Er = 1/2*mH*vrel^2
!-----------------------------------------------------------------------
      use catcom, only : ip_cxr
      use cimcom, only : ismax, ndis
      use cimcxr, only : cxr_edlt, cxr_emax, cxr_emin, cxr_erel
     >    , cxr_nemx, cxr_sigv, cxr_vrel, ndcxre, ndcxrz
      use cphcns, only : cev, cmp
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  ie, iz, imd, izmx
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   real*8   zerl, zne, zn0, zte
      real*8   zerl, zn0, zte
      real*8   acof(ndis)
!
      write(n6,'(/2x,"*** sgvcxr_set ***  ADAS-data ")')
      write(n6,'(2x,"ismax =",i3,"  ndis, ndcxrZ =",2i3)')
     >   ismax, ndis, ndcxrZ
      if( ismax.gt.ndis .or. ismax.gt.ndcxrZ ) then
        call wexit("sgvcxr_set","ismax > ndiz, ndcxrZ")
      endif
!
      cxr_nemx = ndcxrE
      cxr_emin = 0.1d0
      cxr_emax = 10.0d3
      cxr_edlt = dlog(cxr_emax/cxr_emin)/dfloat(cxr_nemx-1)
!
      do ie = 1, cxr_nemx
      zerl = cxr_edlt*dfloat(ie-1)
      zerl = cxr_emin*dexp(zerl)
      cxr_erel(ie) = zerl
      cxr_vrel(ie) = dsqrt(2.0d0*zerl*cev/cmp)
!-----
      zn0 = 1.0d17   ! <== independent on N0
!-----
      zte = zerl
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   if( ip_cxr.eq.0 ) then
      if( ip_cxr(1) == 0 ) then
        cxr_sigv(ie,iz) = 0.0d0
      else
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik     call atm_eval2(ip_cxr,zn0,zte,acof,ismax,ndis)
        call atm_eval2( ip_cxr(1), zn0, zte, acof, ismax, ndis, 1 )
        do iz = 1, ismax
          cxr_sigv(ie,iz) = acof(iz)
        enddo
      endif
      enddo
!
!::log-world
!     cxr_edlt is unit in log-world
      do ie = 1, cxr_nemx
      cxr_erel(ie) = dlog( cxr_erel(ie) )
      cxr_vrel(ie) = dlog( cxr_vrel(ie) )
      do iz = 1, ismax
      cxr_sigv(ie,iz) = dlog( cxr_sigv(ie,iz) )
      enddo
      enddo
!
!::debug write
      write(n6,'(2x,"nemx =",i4,"  emin =",1p2e12.4,"  emax =",
     >  1p2e12.4)') cxr_nemx, cxr_emin, dexp(cxr_erel(1)),
     >   cxr_emax, dexp(cxr_erel(cxr_nemx))
!
      izmx = min0( ismax, 10 )
      imd = cxr_nemx/5
      if(imd.le.0 ) imd = 1
      write(n6,'(4x,"ie",3x,"erel",8x,10("sigv_",i2.2,5x,:))')
     >    (iz,iz=1,izmx)
      do ie = 1, cxr_nemx, imd
      write(n6,'(2x,i4,1p12e12.3)')
     >   ie, dexp(cxr_erel(ie)),
     >  (dexp(cxr_sigv(ie,iz)),iz=1,izmx)
      enddo
!
      return
      end
!
!***********************************************************************
      subroutine sgvcxr(ic,iz,pvel2,sgv)
!***********************************************************************
!
!       ic    (i) :  cell number
!       iz    (i) :  charge state
!       pvel2 (i) :  impurity ion velocity^2
!       sgv   (o) :  N0*<sig*v>_cxr
!
!-----------------------------------------------------------------------
      use catcom, only : ip_cxr
      use cimcxr, only : cxr_edlt, cxr_emax, cxr_emin, cxr_nemx
     >    , cxr_sigv
      use cntcom, only : rmas
      use cphcns, only : cev, cmp, cpi
      use cplcom, only : tn0, tt0
      implicit none
!
!::argument
      integer, intent(in)  :: ic, iz
      real(8), intent(in)  :: pvel2
      real(8), intent(out) :: sgv
!
!::local variables
      integer  ig
      real*8   zt0, zn0, zvr2, zerl
      real*8   dxx, yy, xx1, xx2, xx3
      real*8   deps     ! KSFUJI
      integer  ixa, ixb
!
      if( ip_cxr(1) == 0 ) then
        sgv = 0.0d0
        return
      endif
!
      if( iz.le.0 ) then
        sgv = 0.0d0
        return
      endif
!
!::relative energy in amu
      zn0  = 1.0d0
      zerl = pvel2
!
      if( ic.gt.0 ) then
        ig = 1
        zt0  = tT0(ic,ig)
        zn0  = tN0(ic,ig)
        zvr2 = 8.0d0*cev/(cpi*rmas(ig))*zt0 + pvel2
        zerl = 0.5d0*cmp*zvr2/cev
      endif
!
!::interpolation in log-world
      deps = 1.0d-7                       ! for handling boudaries of index
      xx1 = zerl
      xx1 = dmax1( xx1, cxr_emin+deps )
      xx1 = dmin1( xx1, cxr_emax-deps )
      xx2 = dlog( xx1/cxr_emin )
      xx3 = xx2/cxr_edlt + 1.0d0
      ixa  = int(xx3)
      ixb  = min0( ixa+1, cxr_nemx )  ! KSFUJI
      dxx  = xx3-dfloat(ixa)
      yy   = cxr_sigv(ixa,iz)*(1.0d0-dxx) + cxr_sigv(ixb,iz)*dxx
!
      sgv = dexp(yy)
      sgv = zn0*sgv
      return
      end
