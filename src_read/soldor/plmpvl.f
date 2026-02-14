!**********************************************************************
      subroutine plmpvl
!**********************************************************************
!
!       vlmp  : volume from axis to i+1/2 surface
!       sfmp  : surface area
!       armp  : radius from axis
!       arhmp : radius from axis
!       romp  : normalized effective radius
!       rohmp : normalized effective radius
!       dvmp  : volume of i-tube
!       drmp  : width of i-tube
!
!       wdmp : width of i-tube
!       wfmp : width of i-tube weighted surface
!       wlmp : length along b  between ixp and oxp
!
!----------------------------------------------------------------------
      use cgdcom, only : grdx, grdy, mqs1, npmx
      use cphcns, only : cpi
      use cplmet, only : hgdy, hpit, icaxs, icmpe, icmps, icspx, icwl1
     >    , jcxp1, jcxp2
      use cpmpls, only : arhmp, armp, drmp, dvmp, imd1, imd2, jmd1, jmd2
     >    , r0mp, ramp, rohmp, romp, sfmp, vlmp, wdmp, wfmp, wlmp, xmd1
     >    , xmd2, ymd1, ymd2
      use csize,  only : ndy
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  ldbg
      parameter ( ldbg = 1 )
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  ic, jc, jcs, jce, icp, nsiz, ix, iy, jcm, id, mc
      integer  ic, jc, jcs, jce, icp, ix, iy, jcm, id
      integer  mx1,mx2,mx3,mx4,my1,my2,my3,my4
      integer  mx5,mx6,mx7,mx8,my5,my6,my7,my8
      real*8   r1,r2,r3,r4,r5,r6,r7,r8,z1,z2,z3,z4,z5,z6,z7,z8
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8   ra,rb,rc,rd,ro,re,rw,rn,rs,za,zb,zc,zd,zo,ze,zw,zn,zs
      real*8   ra,rb,rc,rd,ro,re,rw,za,zb,zc,zd,ze,zw
      real*8   sf, vl, ar, dwd, smw, dls, sps, sfp, sfm
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8   smv, sms, sml, smr, smrh, smvp, smvm, smwp, smwm
!ik   real*8   zr1, za1, xaxs, yaxs, eps, dl, pt, bx, by, bz
      real*8   smv, sms, sml, smr, smrh
      real*8   zr1, za1, xaxs, yaxs, eps, dl, pt
      real*8   zsuf(ndy)
! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer   nft
!ik   character cdsn*80
!ik   logical   lex
!
      write(n6,'(/2x,"*** plmpvl ***")')
!
      jcs = jcxp1 + 1
      jce = jcxp2 - 1
!
!----------------------------------------------------------------------
!::separatirx
!----------------------------------------------------------------------
      ic  = icmps
      smv = 0.0d0
      sms = 0.0d0
      do jc = jcs, jce
      call mcpnt(jc,ic, mx1,mx2,mx3,mx4,my1,my2,my3,my4)
      ra = grdx(mx1,my1)
      rb = grdx(mx2,my2)
      za = grdy(mx1,my1)
      zb = grdy(mx2,my2)
      smv = smv + 0.5d0*(rb*rb+ra*ra)*(zb-za)
      sms = sms + 0.5d0*(rb+ra)*sqrt((rb-ra)**2+(zb-za)**2)
      enddo
      smv = smv*cpi
      sms = sms*2.0d0*cpi
      write(n6,'(2x,"spx-Volm =",1pe12.3,"  Are =",1pe12.3)') smv,sms
!
!----------------------------------------------------------------------
!::clear
!----------------------------------------------------------------------
      armp (1:ndy) = 0.0d0
      arhmp(1:ndy) = 0.0d0
      romp (1:ndy) = 0.0d0
      rohmp(1:ndy) = 0.0d0
      vlmp (1:ndy) = 0.0d0
      sfmp (1:ndy) = 0.0d0
      dvmp (1:ndy) = 0.0d0
      drmp (1:ndy) = 0.0d0
      wdmp (1:ndy) = 0.0d0
      wfmp (1:ndy) = 0.0d0
      wlmp (1:ndy) = 0.0d0
!
      zsuf (1:ndy) = 0.0d0
!
!----------------------------------------------------------------------
!::volume, surface, length
!----------------------------------------------------------------------
      do ic = icwl1, icaxs
      smv = 0.0d0
      sms = 0.0d0
      sml = 0.0d0
      do jc = jcs, jce
      call mcpnt(jc,ic, mx1,mx2,mx3,mx4,my1,my2,my3,my4)
      r1 = grdx(mx1,my1);  z1 = grdy(mx1,my1)
      r2 = grdx(mx2,my2);  z2 = grdy(mx2,my2)
      r3 = grdx(mx3,my3);  z3 = grdy(mx3,my3)
      r4 = grdx(mx4,my4);  z4 = grdy(mx4,my4)
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   ro = (r1+r2+r3+r4)/4.0d0;  zo = (z1+z2+z3+z4)/4.0d0
      ro = (r1+r2+r3+r4)/4.0d0
      sf = 0.5d0*((r3-r1)*(z4-z2)-(z3-z1)*(r4-r2))
      id = (mx1-mx3)**2 + (my1-my3)**2 + (mx2-mx4)**2 + (my2-my4)**2
      if( id.eq.0 ) sf = 0.0d0
      vl = 2.0d0*cpi*ro*sf
      ar = 2.0d0*cpi*0.5d0*(r3+r4)*dsqrt((r4-r3)**2+(z4-z3)**2)
      smv = smv + vl
      sms = sms + ar
!
      if( ic.gt.icwl1 .and. ic.le.icmpe ) then
      ra = 0.5d0*(r1+r4);  za = 0.5d0*(z1+z4)
      rb = 0.5d0*(r2+r3);  zb = 0.5d0*(z2+z3)
      dl = sqrt((rb-ra)**2+(zb-za)**2)
      pt = 0.5d0*(hpit(jc,ic)+hpit(jc+1,ic))
      sml = sml + dl/pt
      endif
!
      enddo
!
      dvmp(ic) = smv
      sfmp(ic) = sms
      wlmp(ic) = sml
      enddo
!
      do ic = icwl1+1, icaxs
      drmp(ic) = dvmp(ic)/(0.5d0*(sfmp(ic)+sfmp(ic-1)))
      enddo
!
      ic   = icaxs
      smv  = 0.0d0
      smr  = 0.0d0
      smrh = 0.5d0*drmp(ic)
      do ic = icaxs, icwl1, -1
      vlmp(ic)  = smv
      armp(ic)  = smr
      arhmp(ic) = smrh
      if( ic.eq.icwl1 ) cycle
      smv  = smv + dvmp(ic)
      smr  = smr + drmp(ic)
      smrh = smr + 0.5d0*drmp(ic-1)
      enddo
!
!----------------------------------------------------------------------
!::effective width for flux
!----------------------------------------------------------------------
      do ic = icwl1, icaxs-1
      icp = ic + 1
      sms = 0.0d0
      smw = 0.0d0
!
      do jc = jcs, jce
!
      call mcpnt(jc,ic, mx1,mx2,mx3,mx4,my1,my2,my3,my4)
      r1 = grdx(mx1,my1);  z1 = grdy(mx1,my1)
      r2 = grdx(mx2,my2);  z2 = grdy(mx2,my2)
      r3 = grdx(mx3,my3);  z3 = grdy(mx3,my3)
      r4 = grdx(mx4,my4);  z4 = grdy(mx4,my4)
!
      call mcpnt(jc,icp, mx5,mx6,mx7,mx8,my5,my6,my7,my8)
      r5 = grdx(mx5,my5);  z5 = grdy(mx5,my5)
      r6 = grdx(mx6,my6);  z6 = grdy(mx6,my6)
      r7 = grdx(mx7,my7);  z7 = grdy(mx7,my7)
      r8 = grdx(mx8,my8);  z8 = grdy(mx8,my8)
!
      ra = 0.5d0*(r1+r4);  za = 0.5d0*(z1+z4)
      rb = 0.5d0*(r2+r3);  zb = 0.5d0*(z2+z3)
      rc = 0.5d0*(r6+r7);  zc = 0.5d0*(z6+z7)
      rd = 0.5d0*(r8+r5);  zd = 0.5d0*(z8+z5)
!
      re = r3;  ze = z3
      rw = r4;  zw = z4
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   ro = 0.5d0*(re+rw);  zo = 0.5d0*(ze+zw)
!ik   rn = 0.5d0*(rc+rd);  zn = 0.5d0*(zc+zd)
!ik   rs = 0.5d0*(ra+rb);  zs = 0.5d0*(za+zb)
      ro = 0.5d0*(re+rw)
!
      dls = dsqrt((re-rw)**2+(ze-zw)**2)
      sps = 2.0d0*cpi*ro*dls
      sfp = 0.5d0*((rc-rw)*(zd-ze)-(zc-zw)*(rd-re))
      sfm = 0.5d0*((re-ra)*(zw-zb)-(ze-za)*(rw-rb))
      dwd = (sfp+sfm)/dls
      sms  = sms + sps
      smw  = smw + sps/dwd
      enddo
!
      zsuf(ic) = sms
      wdmp(ic) = 0.5d0*(drmp(ic)+drmp(ic+1))
      wfmp(ic) = sms/smw
      enddo
!
!----------------------------------------------------------------------
!::major radius & minor radius
!----------------------------------------------------------------------
      zr1 = sfmp(icspx)**2/(8.0d0*cpi**2*vlmp(icspx))
      za1 = 2.0d0*vlmp(icspx)/sfmp(icspx)
      write(n6,'(2x,"sum-Volm =",1pe12.3,"  Are =",1pe12.3,"  R0 ="
     > ,0pf8.4,"  a =",0pf8.4)') vlmp(icspx),sfmp(icspx),zr1,za1
!
!::axis
      ix = mqs1
      iy = npmx
      xaxs = grdx(ix,iy)
      yaxs = grdy(ix,iy)
!
!::outer-mid-plane
      ic  = icspx
      eps = 1.0e20
      jcs = jcxp1
      jce = (jcxp1+jcxp2)/2
      do jc = jcs, jce
      if( dabs(hgdy(jc,ic)-yaxs).lt.eps ) then
      eps  = dabs(hgdy(jc,ic)-yaxs)
      jcm = jc
      endif
      enddo
      xmd1 = grdx(jcm,ic); ymd1 = grdy(jcm,ic)
      jmd1 = jcm; imd1 = ic
!
!::inner-mid-plane
      ic  = icspx
      eps = 1.0e20
      jce = jcxp2
      jcs = (jcxp1+jcxp2)/2
      do jc = jcs, jce
      if( dabs(hgdy(jc,ic)-yaxs).lt.eps ) then
      eps  = dabs(hgdy(jc,ic)-yaxs)
      jcm = jc
      endif
      enddo
      xmd2 = grdx(jcm,ic); ymd2 = grdy(jcm,ic)
      jmd2 = jcm; imd2 = ic
!
      r0mp = 0.5d0*(xmd1+xmd2)
      ramp = dsqrt(vlmp(icspx)/(2.0d0*cpi**2*r0mp))
!
      write(n6,'(2x,"axis point =",2f8.4,"  out-mid =",2f8.4,2i5,
     >  "  in-mid =",2f8.4,2i5)') xaxs,yaxs,xmd1,ymd1,jmd1,imd1
     >  ,xmd2,ymd2,jmd2,imd2
      write(n6,'(2x,"R0 =",f8.4,"  a =",f8.4)') r0mp,ramp
      write(n6,'(2x,"icwl1 =",i3,"  icspx =",i3,"  icmps =",i3
     > ,"  icmpe =",i3,"  icaxs =",i3)') icwl1,icspx,icmps,icmpe,icaxs
!
!----------------------------------------------------------------------
!::normalized roh = sqrt(vol(i)/vol(N))
!----------------------------------------------------------------------
      do ic = icwl1, icaxs
      romp(ic) = dsqrt(vlmp(ic)/vlmp(icspx))
      enddo
      do ic = icwl1+1, icaxs
      rohmp(ic) = 0.5d0*(romp(ic)+romp(ic-1))
      enddo
      rohmp(icwl1) = rohmp(icwl1+1)
!
!----------------------------------------------------------------------
!::debug write
!---------------------------------------------------------------------
      write(n6,'(3x,"ic",4x,"ro",10x,"roh",9x,"vol",9x,"dvl",9x,"r",
     >  11x,"rh",10x,"dr",10x,"wd",10x,"wf",10x,"sf",10x,"sf",10x,
     >  "ln")')
      do ic = icwl1, icaxs
      write(n6,'(2x,i3,1p12e12.3)')
     >  ic,romp(ic),rohmp(ic),vlmp(ic)
     >  ,dvmp(ic),armp(ic),arhmp(ic)
     >  ,drmp(ic),wdmp(ic),wfmp(ic),sfmp(ic),zsuf(ic),wlmp(ic)
      enddo
!
      return
      end
