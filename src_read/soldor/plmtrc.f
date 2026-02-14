!**********************************************************************
      subroutine plmtrc
!**********************************************************************
      use cplmet, only : hdsp, icel, itpvs, jcel, jcxp1, jtmax, kce, kcw
      use csize,  only : ndy
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jmax, jw, ii, j, i, jwx, jx, jx0, jx1, jx2, jx3, jx4
      integer  it, jmax, jw, j, i, jwx, jx0, jx1, jx2, jx3, jx4
      real*8   zw0, zw1, zw2, zw3, zw4
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8   zh0, zh1, zh2, zhh, fxmax
      real*8   zh1, zh2, zhh, fxmax
      real*8  zz(ndy)
!
!::input data
      fxmax = 2.0d0
!
      write(n6,'(/2x,"*** plmtrc ***   No effect")')
      return
!
!::X-point
      it = itpvs
      jmax = jtmax(it)
      do jw = 1, jmax
      j = jcel(jw,it)
      if( j.eq.jcxp1 ) then
      jwx = jw
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   jx   = j
      goto 110
      endif
      enddo
 110  continue
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   write(n6,'(2x,"jwx =",i3,"  jx =",i3)') jwx,jx
      write(n6,'(2x,"jwx =",i3,"  j =",i3)') jwx,j
!

      do it = itpvs, itpvs+5
!
      i   = icel(jwx,it)
      jx0 = jcel(jwx,it)
      jx1 = jcel(jwx-2,it)
      jx2 = jcel(jwx-1,it)
      zw0 = hdsp(jx0,i,kce)
      zw1 = hdsp(jx1,i,kce)
      zw2 = hdsp(jx2,i,kce)
      zh1 = dlog(zw1) + (dlog(zw2)-dlog(zw1))*dfloat(jwx-(jwx-2))
      zh1 = dexp(zh1)
      jx3 = jcel(jwx+1,it)
      jx4 = jcel(jwx+2,it)
      zw3 = hdsp(jx3,i,kce)
      zw4 = hdsp(jx4,i,kce)
      zh2 = dlog(zw3) + (dlog(zw4)-dlog(zw3))*dfloat(jwx-(jwx+1))
      zh2 = dexp(zh2)
      zhh = 0.5d0*(zh1+zh2)
!
      if( dabs(zw2/zw0).gt.fxmax .or. dabs(zw3/zw0).gt.fxmax ) then
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   write(n6,'(2x,"hdsp(jx,i,kce/kcw) set new value  fxmax =",f8.3,
      write(n6,'(2x,"hdsp(j,i,kce/kcw) set new value  fxmax =",f8.3,
     >  " it =",i3,"  hdsp =",1p3e12.3,2x,a,1pe12.3)')
     >  fxmax, it, zw2, zw0, zw3, " ==>",zhh
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   hdsp(jx,i,kce) = zhh
!ik   hdsp(jx,i,kcw) = zhh
      hdsp(j,i,kce) = zhh
      hdsp(j,i,kcw) = zhh
      endif
!
      enddo
!
!::debug write
      write(n6,'(/5x,"jw",3x,"j",4x,6(a,i2,6x))')
     >   ( "hdsp",it,it=itpvs,itpvs+5 )
!
      do jw = jwx-10, jwx+10
      do it = itpvs, itpvs+5
      j = jcel(jw,it)
      i = icel(jw,it)
      zz(it) = hdsp(j,i,kce)
      enddo
      write(n6,'(2x,2i5,1p6e12.3)')
     >  jw, j, (zz(it),it=itpvs,itpvs+5)
      enddo
!
      return
      end
