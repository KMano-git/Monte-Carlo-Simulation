!***********************************************************************
      subroutine set_scatt
!***********************************************************************
      use cimcom, only : aimas, azmas, dfm, fxm1, fxm2, fxm3, fxam, fxvu
     >    , fxxs, ifm, slnv, slw0, slw1
      use cntcom, only : fvnth, ncmax, nfai, nthe, nthe, tcfi, tcth
     >    , tsfi, tsth, tvcth, tvsth
      use cntpls, only : dene, teme, temi
      use cphcns, only : cev, cme, cmp, cpi
      use csize,  only : ndmc
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  ic, ms, iu, j, imd
      real*8   fac, du, u, uu, f, g, fx1, fx2
      real*8   derf, slo, zden, ztem, x, y
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   real*8   zdmax, ztmin, zxmax, s0max, s1max
!
      write(n6,'(2x,a)') "*** set_scatt ***"
      write(n6,'(2x,"aimas =",f8.4,"  azmas =",f8.4)') aimas, azmas
!
!::clear
      do ic = 0, ndmc
      do ms = 1, 2
      slw0(ic,ms) = 0.0d0
      slw1(ic,ms) = 0.0d0
      slnv(ic,ms) = 0.0d0
      enddo
      enddo
!
!::scat const
      fac  = 3.0d0*sqrt(cpi)/2.0d0
      dfm  = 100.0
      ifm  = 1000
      du   = 1.0/dfm
      fxam(1) = cme/cmp
      fxam(2) = aimas
!
      do ms = 1, 2
      fxxs(ms) = sqrt(fxam(ms))
      do iu = 1, ifm
      u  = du*iu*fxxs(ms)
      uu = u*u
      f  = fac*derf(u)
      g  =(f-3.*u*exp(-uu))/2.0/uu
      fx1 = g/u
      fx2 = (f-g)/u
      fxm1(iu,ms) = -g*(1.+fxam(ms)/azmas)
      fxm2(iu,ms) = sqrt(fxxs(ms)*fx1)
      fxm3(iu,ms) = sqrt(fxxs(ms)*fx2)
      fxvu(iu,ms) = u
      enddo
      enddo
!
!::debug write
      imd = ifm/5
      if(imd.le.0 ) imd = 1
      write(n6,'(2x,"fxm1, fxm2, fxm3")')
      do iu = 1, ifm, imd
      write(n6,'(2x,i5,1p3e14.5,2x,1p3e14.5)')
     >  iu, (fxm1(iu,ms), fxm2(iu,ms), fxm3(iu,ms), ms= 1,2)
      enddo
!
!::constant
      slo = 1.085d-12/azmas
!
!::slowing down
      do ms = 1, 2
      do ic = 1, ncmax
      zden = dene(ic)
      if( ms.eq.1 ) ztem = teme(ic)
      if( ms.eq.2 ) ztem = temi(ic)
      if(ztem.le.0.0) cycle
      x  = slo*zden/ztem**1.5
      y  = 2.*cev*ztem/cmp
      slw0(ic,ms) = x*sqrt(y)
      slw1(ic,ms) = sqrt(x*y/azmas)
      slnv(ic,ms) = sqrt(cmp/(2.*cev*ztem))
      enddo
      enddo
!
      write(n6,'(2x,"slo0, slo1, ste")')
      imd = ncmax/5
      if( imd.le.0 ) imd = 1
      do ic = 1, ncmax, imd
      write(n6,'(2x,i5,1p3e14.5,2x,1p3e14.5)')
     >  ic, (slw0(ic,ms), slw1(ic,ms), slnv(ic,ms), ms=1,2)
      enddo
!
!::check for tri-function
      write(n6,'(2x)')
      write(n6,'(2x,"tri-function")')
      write(n6,'(4x,"nthe =",i5)') nthe
      write(n6,'(4x,"tcth =",3f8.4,2x,3f8.4)')
     >   (tcth(j),j=1,3),(tcth(j),j=nthe-2,nthe)
      write(n6,'(4x,"tcth =",3f8.4,2x,3f8.4)')
     >   (tsth(j),j=1,3),(tsth(j),j=nthe-2,nthe)
      write(n6,'(4x,"nfai =",i5)') nfai
      write(n6,'(4x,"tcfi =",3f8.4,2x,3f8.4)')
     >   (tcfi(j),j=1,3),(tcfi(j),j=nthe-2,nthe)
      write(n6,'(4x,"tsfi =",3f8.4,2x,3f8.4)')
     >   (tsfi(j),j=1,3),(tsfi(j),j=nthe-2,nthe)
      write(n6,'(4x,"fvnth =",f8.1)') fvnth
      write(n6,'(4x,"tvcth =",3f8.4,2x,3f8.4)')
     >   (tvcth(j),j=1,3),(tvcth(j),j=nthe-2,nthe)
      write(n6,'(4x,"tvsth =",3f8.4,2x,3f8.4)')
     >   (tvsth(j),j=1,3),(tvsth(j),j=nthe-2,nthe)
!
      return
      end
