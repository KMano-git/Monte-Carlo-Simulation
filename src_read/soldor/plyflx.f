!***********************************************************************
      subroutine plyflx(crg,it,kdr)
!***********************************************************************
!
!       crg     irg       it     kdr
!       "edg"    6       itmpe    -1
!       "sol"    2       itsle    +1
!       "odv"    1       itsls    +1
!       "idv"    3       itsls    +1
!       "prv"   4/5      itpve    -1
!
!-----------------------------------------------------------------------
      use cgdcom, only : grdx, grdy, mpsp, mqs1
      use cphcns, only : cpi
      use cplcom, only : ama, nion, vne, vte, vti, vva
      use cplmet, only : gare, icel, itmpe, itsle, jcel, jtmax, kreg
      use cplqcn, only : qfy_cd, qfy_df, qfy_vh 
      use cpmpls, only : xmd1, xmd2, ymd1, ymd2
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   character  crg*(*)
!ik   integer    it, kdr
      character, intent(in) :: crg*(*)
      integer,   intent(in) :: it, kdr
!
!::local variables
      character  csgn*1, cdsn*17
      integer  i6, j0, i0, ii, jw, j, i, jmax, isf, irg
      integer  jwst, jwen
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  ia, m1a, m2a, m3, m4
      integer  ia, m1a, m3, m4
      integer  mx1, mx2, mx3, mx4, my1, my2, my3, my4
      real*8   xxp, yxp, xpc, ypc, th0
      real*8   x1, y1, x2, y2, x0, y0, dl, zsv, zsv2, dln
      real*8   fdf3, fdf4, fcd3, fcd4, f1, f3, f4, f34, th
      real*8   smf1, smf3, smf4
!
      character tbrg(7)*3
      data  tbrg/"odv","sol","idv","prv","prv","edg","man"/
!
      csgn = "p"
      if( kdr.lt.0 ) csgn = "m"
      write(cdsn,'(a,i3.3,a,a)') "Yflx_"//crg//"_",it,csgn,".txt"
      i6 = 21
      open(unit=i6,file=cdsn)
!
      write(i6,'(2x,"*** plyflx ***  ",a,2x,"it =",i3,"  kdr =",i2)')
     >  crg, it, kdr
!
!::th0
      j0 = mqs1
      i0 = mpsp
      xxp = grdx(j0,i0)
      yxp = grdy(j0,i0)
      xpc = 0.5d0*(xmd1+xmd2)
      ypc = 0.5d0*(ymd1+ymd2)
      dln = sqrt( (xxp-xpc)**2+(yxp-ypc)**2 )
      th0 = datan2( (yxp-ypc)/dln, (xxp-xpc)/dln )
      th0 = th0/cpi*180.0d0
      write(i6,'(2x,"Xp =",2f9.5,"  Pc =",2f9.5,"  th0 =",f9.3)')
     > xxp, yxp, xpc, ypc, th0
!
      ia = 1
      m1a = 2*ia - 1
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   m2a = 2*ia
      m3  = 2*nion + 1
      m4  = 2*nion + 2
!
!
      smf1 = 0.0d0
      smf3 = 0.0d0
      smf4 = 0.0d0
      ii = 0
!
      write(i6,'(3x,"it",4x,"j",2x,"isf",2x,"x1",7x,"y1",7x,
     >  "are",6x,"gare",5x,"x0",7x,"y0",7x,"th",7x,
     >  "Ne",9x,"Te",9x,"Ti",9x,"Va",9x,
     >  "f1p",8x,"f3p",8x,"f4p",8x,"ftp")')
!
      jmax = jtmax(it)
      jwst = 1
      jwen = jmax
      if( it.ge.itsle .and. it.le.itmpe ) then
      jwst = jwst + 1
      jwen = jwen - 1
      endif
!
      do jw = jwst, jwen
      j = jcel(jw,it)
      i = icel(jw,it)
      isf = i
      if( kdr.lt.0 ) isf = icel(jw,it-1)
      irg = kreg(j,i)
      if( irg.eq.0 ) cycle
      if( crg.ne.tbrg(irg) ) cycle
      ii = ii + 1
!
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
      if( kdr.gt.0 ) then
      x1 = grdx(mx4,my4)
      y1 = grdy(mx4,my4)
      x2 = grdx(mx3,my3)
      y2 = grdy(mx3,my3)
      else
      x1 = grdx(mx1,my1)
      y1 = grdy(mx1,my1)
      x2 = grdx(mx2,my2)
      y2 = grdy(mx2,my2)
      endif
      x0 = 0.5d0*(x1+x2)
      y0 = 0.5d0*(y1+y2)
      dl = sqrt((x2-x1)**2+(y2-y1)**2)
!
      zsv  = 2.0d0*cpi*x0*dl  ! S^psi
      zsv2 = gare(j,isf)      ! S^psi
!
      if( jw.eq.jwst ) then
      if( it.ge.itsle .and. it.le.itmpe ) then
      xxp = x1
      yxp = y1
      dln = sqrt( (xxp-xpc)**2+(yxp-ypc)**2 )
      th0 = datan2( (yxp-ypc)/dln, (xxp-xpc)/dln )
      th0 = th0/cpi*180.0d0
      endif
      endif
!
      fdf3 = qfy_df(j,isf,m3)
      fdf4 = qfy_df(j,isf,m4)
      fcd3 = qfy_cd(j,isf,m3) + qfy_vh(j,isf)
      fcd4 = qfy_cd(j,isf,m4)
      f1   = qfy_df(j,isf,m1a)/ama(ia)
      f3   = fdf3 + fcd3
      f4   = fdf4 + fcd4
      f34  = f3 + f4
      smf1 = smf1 + f1
      smf3 = smf3 + f3
      smf4 = smf4 + f4
!
!::th
      dln = sqrt((x0-xpc)**2+(y0-ypc)**2)
      th  = datan2( (y0-ypc)/dln, (x0-xpc)/dln )
      th  = th/cpi*180.0d0
      th  = th-th0
      if( th.gt.180.0d0 )  th = th - 360.0d0
      if( th.lt.-180.0d0 ) th = th + 360.0d0
!
      if(zsv.ne.0.0d0) then
      write(i6,'(2x,i3,2i5,0p6f9.5,f9.3,1p11e11.3)')
     >  it, j, isf, x1, y1, zsv, zsv2,
     >  x0, y0, th,
     >  vne(j,i), vte(j,i), vti(j,i), vva(j,i,ia),
     >  f1/zsv, f3/zsv, f4/zsv, f34/zsv, th0
      else
      write(i6,'(2x,i3,2i5,0p6f9.5,f9.3,1p11e11.3)')
     >  it, j, isf, x1, y1, zsv, zsv2,
     >  x0, y0, th,
     >  vne(j,i), vte(j,i), vti(j,i), vva(j,i,ia),
     >  0.0d0, 0.0d0, 0.0d0, 0.0d0,th0
      endif
!
      enddo
!
      write(i6,'(2x)')
      write(i6,'(2x,"sum_f1p =",1pe11.3,"  sum_f3p =",1pe11.3,
     >  "  sum_f4p =",1pe11.3)') smf1, smf3, smf4
      close(i6)
!
!::summary
      write(n6,'(2x,a,2x,i3,2x,i2,2x,1p4e12.3)')
     >  "Yflx_"//crg, it, kdr, smf1, smf3+smf4, smf3, smf4
!
      return
      end
