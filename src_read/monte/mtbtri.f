!**********************************************************************
      subroutine mtbtri
!**********************************************************************
      use cphcns, only : cpi
      use cntcom, only : fnfi, fnth, fvnth, ndfai, ndthe, nfai, nthe
     >    , tcfi, tcth, tsfi, tsth, tvcth, tvsth
      implicit none
!
      real*8  dthe, xran, dfai, th
      integer i, nvnth
!
      nthe = ndthe-1
      fnth = dfloat(nthe)
      dthe = 1.0d0/fnth
      do 15 i = 1, nthe
      xran    = dthe*(dfloat(i)-0.5d0)
      tcth(i) = dsqrt(1.0d0-xran)
      tsth(i) = dsqrt(xran)
   15 continue
      tcth(nthe+1) = tcth(nthe)
      tsth(nthe+1) = tsth(nthe)
!
      nfai = ndfai-1
      fnfi = dfloat(nfai)
      dfai = 2.0d0*cpi/fnfi
      do 20 i = 1, nfai
      xran = dfai*(dfloat(i)-0.5d0)
      tcfi(i) = cos(xran)
      tsfi(i) = sin(xran)
   20 continue
      tcfi(nfai+1) = tcfi(nfai)
      tsfi(nfai+1) = tsfi(nfai)
!
      nvnth = 2000
      fvnth = dfloat(nvnth)
      dthe  = cpi/fvnth
      do 30 i = 1, nvnth
      th = dthe*(dfloat(i)-0.5d0)
      tvcth(i) = cos(th)
      tvsth(i) = sin(th)
   30 continue
      tvcth(nvnth+1)=tvcth(nvnth)
      tvsth(nvnth+1)=tvsth(nvnth)
!
      return
      end
