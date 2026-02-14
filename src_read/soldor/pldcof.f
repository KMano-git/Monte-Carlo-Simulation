!***********************************************************************
      subroutine pldcof
!***********************************************************************
      use cplcom,     only : fcda, fcet, nion, vdda, vdet, vdxe, vdxi
     >    , vlda, vlet, vlxe, vlxi
      use cplmet,     only : icel, icmpe, icmps, itmax, itmpe, itmps
     >    , jcel, jcmax, jtmax, jtmin, kreg
      use cplvpn,     only : vdvp
      use cpmpls,     only : adamp, avpmp, axemp, aximp, jmd1
      use csize,      only : ndsp, ndx, ndy
      use csonic,     only : itim, time
      use topics_mod, only : lgtpc
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jt, jts, jte, jc, ic, irg, ia, iy
      integer  it, jt, jts, jte, jc, ic, irg, ia
      integer  i6
!
!::clear
      vdda(1:ndx,1:ndy,1:ndsp) = 0.0d0
      vdet(1:ndx,1:ndy,1:ndsp) = 0.0d0
      vdxe(1:ndx,1:ndy) = 0.0d0
      vdxi(1:ndx,1:ndy) = 0.0d0
      vdvp(1:ndx,1:ndy) = 0.0d0
!
!::diffusion coef. in cell
      do it = 1, itmax
      jts = jtmin(it)
      jte = jtmax(it)
      if( it.ge.itmps .and. it.le.itmpe ) then
      jts = jts + 1
      jte = jte - 1
      endif
!
      do jt = jts, jte
      jc  = jcel(jt,it)
      ic  = icel(jt,it)
      irg = kreg(jc,ic)
!
      do ia = 1, nion
! added 1 line integrated calculation with TOPICS by kamata 2020/11/16
      if( lgtpc == 0 ) then
        vdda(jc,ic,ia) = vlda(irg)*fcda(ia)
! added 3 lines integrated calculation with TOPICS by kamata 2020/11/16
      else
        vdda(jc,ic,ia) = adamp(icmps,ia)
      endif
      vdet(jc,ic,ia) = vlet(irg)*fcet(ia)
      enddo
! added 1 line integrated calculation with TOPICS by kamata 2020/11/16
      if( lgtpc == 0 ) then
        vdxi(jc,ic) = vlxi(irg)
        vdxe(jc,ic) = vlxe(irg)
! added 4 lines integrated calculation with TOPICS by kamata 2020/11/16
      else
        vdxi(jc,ic) = aximp(icmps)
        vdxe(jc,ic) = axemp(icmps)
      endif
      enddo
      enddo
!
      i6 = 8000
!
!xx      do ic = icmps-2, icmpe
!xx      write(i6,'(i3,1p15e11.3)') ic,(vdxi(jc,ic),jc=1,jcmax,10)
!xx      enddo
!
!::Da, Xi, Xe in main plasma
      do it = itmps, itmpe
      jts = jtmin(it) + 1
      jte = jtmax(it) - 1
      do jt = jts, jte
      jc  = jcel(jt,it)
      ic  = icel(jt,it)
      irg = kreg(jc,ic)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   iy  = ic
!
      do ia = 1, nion
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   vdda(jc,ic,ia) = adamp(iy,ia)
      vdda(jc,ic,ia) = adamp(ic,ia)
! modified 1/1 lines integrated calculation with TOPICS by kamata 2020/11/16
!xx   vdet(jc,ic,ia) = vlet(irg)*fcet(ia)
      vdet(jc,ic,ia) = vlet(irg)*fcet(ia)
      enddo
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   vdxi(jc,ic) = aximp(iy)
!ik   vdxe(jc,ic) = axemp(iy)
      vdxi(jc,ic) = aximp(ic)
      vdxe(jc,ic) = axemp(ic)
! added 1 line addition and reconsider of passed data by kamata 2022/02/23
      vdvp(jc,ic) = avpmp(ic,1)
      enddo
      enddo
!
!::debug write
      write(i6,'(/2x,"pldcof   set Xi ",1pe14.6,i6,
     >  "  icmps,icmpe =",2i4)') time, itim, icmps, icmpe
      do ic = icmps-2, icmpe
      write(i6,'(i3,1p15e11.3)') ic,(vdxi(jc,ic),jc=1,jcmax,10)
      enddo
!
      write(i6,'(2x,"ic, Da, Et, Xi, Xe")')
      ia = 1
      jc = jmd1
      do ic = icmps-2, icmpe
      write(i6,'(i4,1p4e12.3)') ic,vdda(jc,ic,ia), vdet(jc,ic,ia),
     >    vdxi(jc,ic), vdxe(jc,ic)
      enddo
!
      return
      end
