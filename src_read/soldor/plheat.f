!**********************************************************************
      subroutine plheat
!**********************************************************************
!
!      additional heating
!
!---------------------------------------------------------------------
      use cntcom,     only : mcel, volm
      use cplcom,     only : nlp, ssnc, ssnv, swec, swev, swic, swiv
      use cplmet,     only : icaxs, icmpe, icmps, jcxp1, jcxp2
      use csize,      only : ndy
      use csonic,     only : itim, time
      use ctopics,    only : ktop_cr, ktop_ese, ktop_esi, ktop_ps
     >    , ktop_sn, ktop_we, ktop_wi
      use cunit,      only : cdgr, lmspe, lmype, mygrp, mype
      use topics_mod, only : gmion, gvr, khsc, khsn, khsr, kpsc, kpsn
      implicit none
!
!::local variables
      real*8   siad(ndy), wiad(ndy), wead(ndy), tsiad, twiad, twead
      real*8   ztvl, zdvl
      real*8   zsi, zwi, zwe, zsic, zsiv, zwic, zwiv, zwec, zwev
      integer  iy, i, j, ic, ia, i6, jcxst, jcxen
      integer  ldbg
!
      if( mygrp.ne.cdgr(1) ) return
      if( lmype.ne.lmspe )   return
!
      i6 = 2000 + mype
      ldbg = 0
      if( mod(itim,1000).eq.0 .and. nlp.eq.1 ) ldbg = 1
!
      if( ldbg.eq.1 ) then
      write(i6,'(2x,"*** plheat ***  time =",1pe14.6,"  itim =",i6)')
     >   time, itim
      endif
!
!::mesh number
      jcxst = jcxp1 + 1
      jcxen = jcxp2 - 1
!
!::clear
      siad(1:ndy) = 0.0d0
      wiad(1:ndy) = 0.0d0
      wead(1:ndy) = 0.0d0
      tsiad = 0.0d0
      twiad = 0.0d0
      twead = 0.0d0
!
!::source terms
      ztvl = 0.0d0
      iy = 0
      do i = icaxs, icmps, -1
      iy = iy + 1
!
      zdvl = 0.0d0
      do j = jcxst, jcxen
      ic = mcel(j,i)
      zdvl = zdvl + volm(ic)
      enddo
      ztvl = ztvl + zdvl
!
! modified 3/8 lines addition and reconsider of passed data by kamata 2022/02/23
!ik   siad(i) =  gvr(iy,1,kpsb)/zdvl
!ik   wiad(i) = (gvr(iy,1,khsb)+gvr(iy,2,khsb))/zdvl
!ik   wead(i) = (gvr(iy,0,khsb)+gvr(iy,0,khsj))/zdvl
      siad(i) = siad(i) + ( ktop_ps  * gvr(iy,1,kpsc)
     >                    + ktop_sn  * gvr(iy,1,kpsn) ) / zdvl
      wead(i) = wead(i) + ( ktop_we  * gvr(iy,0,khsn)
     >                    + ktop_cr  * gvr(iy,0,khsr)
     >                    + ktop_ese * gvr(iy,0,khsc) ) / zdvl
      wiad(i) = wiad(i) + ( ktop_wi  * sum( gvr(iy,1:gmion,khsn) )
     >                    + ktop_esi * sum( gvr(iy,1:gmion,khsc) ) )
     >                    / zdvl
      do j = jcxst, jcxen
      ic = mcel(j,i)
      tsiad = tsiad + siad(i)*volm(ic)
      twiad = twiad + wiad(i)*volm(ic)
      twead = twead + wead(i)*volm(ic)
      enddo
!
      if( ldbg.eq.1 ) then
      write(i6,'(2x,i3,1p5e12.3)')
     >   iy, zdvl, ztvl, tsiad, twiad, twead
      endif
      enddo
!
!::source terms in soldor
      do i = icmps, icmpe
      zsi  = siad(i)
      zwi  = wiad(i)
      zwe  = wead(i)
      zsic = zsi
      zsiv = 0.0d0
      zwic = zwi
      zwiv = 0.0d0
      zwec = zwe
      zwev = 0.0d0
!
      ia = 1
      do j = jcxst, jcxen
      ssnc(j,i,ia) = ssnc(j,i,ia) + zsic
      ssnv(j,i,ia) = ssnv(j,i,ia) + zsiv
!
      swic(j,i) = swic(j,i) + zwic
      swiv(j,i) = swiv(j,i) + zwiv
!
      swec(j,i) = swec(j,i) + zwec
      swev(j,i) = swev(j,i) + zwev
      enddo
      enddo
!
      return
      end
