!**********************************************************************
      subroutine plntsr(tflux)
!**********************************************************************
!
!       source terms of plasma
!
!----------------------------------------------------------------------
      use cntcom, only : cstyp, iplx, iply, mrgn, ncmax, ncmax2, tflex
     >    , tflpm, volm
      use cntmnt, only : dotn, dotn2, flexh, flpmp, ift => i6_src, sdotn
     >    , sdotn2, sn0, ssn, ssp, sumsn, sumsp, sumwe, sumwi, swe, swi
     >    , totmsn, totmsp, totmwe, totmwi, totsn, totsp, totwe, totwi
      use cntwcn, only : swabs, swerr, swnrm, swpmp, swreg, wabs, wcsrc
     >    , wden, werr, wion, wisrc, wnrm, wpmp, wreg, wssn, wssp, wsum
     >    , wswe, wswi, wtot
      use cplcom, only : ama, nion, nlp,q2a, q3, q4, ssnc, ssnv, sspc
     >    , sspv, swec, swev, swic, swiv, vna, vne,vte
      use cplmet, only : hvol, icel, itpve, itmax, jcel, jtmax, kreg
      use cplqcn, only : mrgnp
      use csize,  only : ndsp
      use csonic, only : itim, time
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8  tflux
      real(8), intent(in) :: tflux
!
!::local variables
      real*8  stalf, stalfp, zbn
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer nsza, nsiz, ia, ig, it, i, jmax, jt, j, ic
      integer nsza, nsiz, ia, it, i, jmax, jt, j, ic
      integer ir, ir0, ir1, ir2
      real*8  zfc, zn0, zsn, zsp, zwe, zwi
      real*8  zsnc, zsnv, zspc, zspv, zwic, zwiv, zwec, zwev
      real*8  zvte, zvnm, zqnm, zab, zbi, zbe
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8  zte, fac, zsmtot, zsmreg(0:10), zfac
!ik   integer kdbg_rec
!ik   integer :: i6, ii, ift, ift2
! modified 1/1 lines replace all include files with module files by kamata 2021/08/18
!ik   integer :: ii, ift, ift2
      integer :: ii, ift2
!
! deleted 1 line replace all include files with module files by kamata 2021/08/18
!ik   ift = i6_src
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   kdbg_rec = 0
!
!x    write(n6,'(/2x,"***  plntsr  ***")')
!
!-----------------------------------------------------------------------
!::linearlized coef.
!-----------------------------------------------------------------------
      stalf  = 1.0d0
      stalfp = 1.0d0
!
!-----------------------------------------------------------------------
!::flux
!-----------------------------------------------------------------------
      dotn = tflux
!
!-----------------------------------------------------------------------
!::total source terms for each source
!-----------------------------------------------------------------------
      nsza = 10*ndsp
      nsiz = 10
      call setd( sumsn, nsza, 0.0d0 )
      call setd( sumsp, nsza, 0.0d0 )
      call setd( sumwi, nsiz, 0.0d0 )
      call setd( sumwe, nsiz, 0.0d0 )
!
!-----------------------------------------------------------------------
!::source
!-----------------------------------------------------------------------
      do ic = 1, ncmax
      ir  = mrgn(ic)
      ir2 = mrgnp(ic)     !  tube (itmpe) is out of system
      j   = iplx(ic)
      i   = iply(ic)
      zfc = dotn/swnrm/volm(ic)
!
!::particle & momentum source
      do ia = 1, nion
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   ig = ia                      !  temporary  2003/02/23
!ik   zn0 = zfc*wden(ic,ig)
      zn0 = zfc*wden(ic,ia)
      zsn = zfc*wssn(ic,ia)
      zsp = zfc*wssp(ic,ia)
!
!::total
      sumsn(ir2,ia) = sumsn(ir2,ia) + zsn*volm(ic)
      sumsp(ir2,ia) = sumsp(ir2,ia) + zsp*volm(ic)
!
!::in plasma (tube it=itmax is out of system)
      if( j.gt.0 .and. i.gt.0 .and. ir2 <= 6) then
!
!::neutral density
      sn0(j,i,ia) = sn0(j,i,ia) + zn0
!
!::linearlized particle source
      zbn = stalf
      if( zsn.lt.0.0d0 ) zbn = -stalf-1.0d0
      zsnc = + (1.0d0+zbn)*zsn
      zsnv = - zbn*zsn/vna(j,i,ia)
      ssnc(j,i,ia) = ssnc(j,i,ia) + zsnc
      ssnv(j,i,ia) = ssnv(j,i,ia) + zsnv
      ssn (j,i,ia) = ssn (j,i,ia) + zsn
!
!::linearlized momentum source  according to B2-code
      stalfp = 0.0d0    ! %KS
      zvte = 4.19d5*dsqrt(vte(j,i))
      zvnm = 0.01d0*zvte
      zqnm = ama(ia)*vne(j,i)*zvnm
      zab  = abs(zsp)
      zspc = + zsp + stalfp*zab*q2a(j,i,ia)/zqnm
      zspv = - stalfp*zab/zqnm
      sspc(j,i,ia) = sspc(j,i,ia) + zspc
      sspv(j,i,ia) = sspv(j,i,ia) + zspv
      ssp (j,i,ia) = ssp (j,i,ia) + zsp
      endif
!
      enddo  ! loop(ia)
!
!::ion & ele energy source
      zwi = zfc*wswi(ic)
      zwe = zfc*wswe(ic)
!
!::total
      sumwi(ir2) = sumwi(ir2) + zwi*volm(ic)
      sumwe(ir2) = sumwe(ir2) + zwe*volm(ic)
!
!::in plasma
!  Sn, Wi (nteraction of neutrals), system excludes the tube it=itmax
!  N0,E0 are defined in whole region (including hot core)
      if( j.gt.0 .and. i.gt.0 .and. ir2<=6 ) then
!
!
!::linearlized ion energy source
      zbi = stalf
      if( zwi.lt.0.0d0 ) zbi = -stalf-1.0d0
      zwic = + (1.0d0+zbi)*zwi
      zwiv = - zbi*zwi/q3(j,i)
      swic(j,i) = swic(j,i) + zwic
      swiv(j,i) = swiv(j,i) + zwiv
      swi (j,i) = swi (j,i) + zwi
!
!::linearlized ele energy source
      zbe = stalf
      if( zwe.lt.0.0d0 ) zbe = -stalf-1.0d0
      zwec = + (1.0d0+zbe)*zwe
      zwev = - zbe*zwe/q4(j,i)
      swec(j,i) = swec(j,i) + zwec
      swev(j,i) = swev(j,i) + zwev
      swe (j,i) = swe (j,i) + zwe
!xx      call check_NaN(ic,j,i,swev(j,i),"plntsr B")
      endif
!
      enddo
!
!::hot core region
      ir = 7
      do ia = 1, nion
      totmsn(ia) = totmsn(ia) + sumsn(ir,ia)
      totmsp(ia) = totmsp(ia) + sumsp(ir,ia)
      enddo
      totmwi = totmwi + sumwi(ir)
      totmwe = totmwe + sumwe(ir)
!
!----------------------------------------------------------------------
!::conv-check
!----------------------------------------------------------------------
      nsza = 10*ndsp
      nsiz = 10
      call setd( totsn, nsza, 0.0d0 )
      call setd( totsp, nsza, 0.0d0 )
      call setd( totwi, nsiz, 0.0d0 )
      call setd( totwe, nsiz, 0.0d0 )
!
!::integrate in system region
!  except tube (itmax) ( ir=1,2,..6, )
      do it = 2, itmax-1
      if( it.eq.itpve ) cycle
      jmax = jtmax(it)
      do jt = 2, jmax-1
      j  = jcel(jt,it)
      i  = icel(jt,it)

!--it=itmax is excluded because the tube is out of system
      ir = kreg(j,i)
      totwi(ir) = totwi(ir) + (swic(j,i)+swiv(j,i)*q3(j,i))*hvol(j,i)
      totwe(ir) = totwe(ir) + (swec(j,i)+swev(j,i)*q4(j,i))*hvol(j,i)
      do ia = 1, nion
      totsn(ir,ia) = totsn(ir,ia) +
     >              (ssnc(j,i,ia)+ssnv(j,i,ia)*vna(j,i,ia))*hvol(j,i)
      totsp(ir,ia) = totsp(ir,ia) +
     >              (sspc(j,i,ia)+sspv(j,i,ia)*q2a(j,i,ia))*hvol(j,i)
      enddo
      enddo
      enddo
!
!::hot core region (plasma center)
      ir = 7
      do ia = 1, nion
      totsn(ir,ia) = totmsn(ia)
      totsp(ir,ia) = totmsp(ia)
      enddo
      totwi(ir) = totmwi
      totwe(ir) = totmwe
!
!::whole region
      do ir = 1, 7
      do ia = 1, nion
      sumsn(10,ia) = sumsn(10,ia) + sumsn(ir,ia)
      sumsp(10,ia) = sumsp(10,ia) + sumsp(ir,ia)
      totsn(10,ia) = totsn(10,ia) + totsn(ir,ia)
      totsp(10,ia) = totsp(10,ia) + totsp(ir,ia)
      enddo
      sumwi(10) = sumwi(10) + sumwi(ir)
      sumwe(10) = sumwe(10) + sumwe(ir)
      totwi(10) = totwi(10) + totwi(ir)
      totwe(10) = totwe(10) + totwe(ir)
      enddo
!
!----------------------------------------------------------------------
!::pumped flux
!----------------------------------------------------------------------
      flexh = dotn*wabs/swnrm
      flpmp = dotn*wpmp/swnrm
      tflex = tflex + flexh
      tflpm = tflpm + flpmp
!
!----------------------------------------------------------------------
!::actual dotn
!----------------------------------------------------------------------
      ia = 1
      dotn2  = sumsn(10,ia) + flexh
      sdotn  = sdotn  + dotn
      sdotn2 = sdotn2 + dotn2
!
!::debug write
!      if( nlp.eq.5 ) then
!      write(n6,'(2x,"dotn =",1pe14.6,"  dotn2 =",1pe14.6,"  flexh =",
!     >  1pe14.6 )')  dotn, dotn2, flexh
!      ia = 1
!      write(n6,'(2x,"sumsn = ",1p10e14.6)') (sumsn(i,ia),i=1,10)
!      endif

!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
!xx      i6 = 0
!xx      if( itim == 20 .and. nlp > nlpmn ) i6 = 21000 + mype
!xx      ift = i6
      if( ift.eq.0 ) return

!::fort.21001  totSn(irg),totWi in plntsr
      ia = 1
      write(ift,'(/2x,"*** plntsr ***  itim =",i7,"  time =",
     >  1pe14.6,"  nlp =",i3,2x,a)') itim, time, nlp, cstyp
      write(ift,'(2x,i2,2x,a,"  dotn =",1pe12.4,
     > "  wtot =",1p2e14.6,"  werr =",1pe14.6,"  wion/wabs =",
     > 1p2e14.6)') wisrc,wcsrc,dotn,wtot,wsum,werr,wion,wabs
      write(ift,'(2x,"typ",2x,"src",5x,"dotn",7x,"sum",11x,"odp",8x,
     > "sol",8x,"idp",8x,"opv",8x,"ipv",8x,"edg",8x,"hot",8x,
     > "abs",8x,"pmp",8x,"err")')
 602  format(2x,a6,2x,a3,1pe11.3,1pe14.6,1x,1p12e11.3)

      write(ift,602) cstyp(1:6),"Wmn",dotn,wnrm,
     >   (wreg(ir)/wnrm,ir=1,7),wabs/wnrm,wpmp/wnrm,werr/wnrm
      write(ift,602) cstyp(1:6),"Wsr",dotn,swnrm,
     >   (swreg(ir)/swnrm,ir=1,7),swabs/swnrm,swpmp/swnrm,swerr/swnrm
      write(ift,602) cstyp(1:6),"Sna",dotn,dotn2,
     >   (sumsn(ir,ia),ir=1,7),flexh,flpmp
      write(ift,602) cstyp(1:6),"Spa",dotn,
     >              sumsp(10,ia),(sumsp(ir,ia),ir=1,7)
      write(ift,602) cstyp(1:6),"Swi",dotn,sumwi(10),(sumwi(ir),ir=1,7)
      write(ift,602) cstyp(1:6),"Swe",dotn,sumwe(10),(sumwe(ir),ir=1,7)

      write(ift,602) "total","Sna",sdotn,sdotn2,
     >   (totsn(ir,ia),ir=1,7), tflex,tflpm
      write(ift,602) "total","Spa",sdotn2,
     >                         totsp(10,ia),(totsp(ir,ia),ir=1,7)
      write(ift,602) "total","SWi",sdotn2,totwi(10),(totwi(ir),ir=1,7)
      write(ift,602) "total","SWe",sdotn2,totwe(10),(totwe(ir),ir=1,7)

!::fort.21051  N0(ic),Sn,Wi in plntsr
      ift2 = 0
      if( ift2 == 0 ) return
!
      ift2 = ift + 50
      write(ift2,'(/2x,"*** plntsr ***  itim =",i7,"  time =",
     >  1pe14.6,"  nlp =",i3,2x,a)') itim, time, nlp, cstyp
      write(ift2,'(4x,a)')
     >  "ic  ir      N0(j,i)    Sn(j,i)   Sp(j,i)   Wi(j,i)   We(j,i)"
      ii = 0
      do ic = 1, ncmax2
        j = iplx(ic)
        i = iply(ic)
        ir0 = kreg(j,i)
        ir1 = mrgn(ic)
        ir2 = mrgnp(ic)

        if( j <= 0 .or. i <= 0 ) cycle
        ii = ii + 1
        if( mod(ii,200) == 0 ) then
          write(ift2,'(i6,3i3,1p5e12.4)') ic, ir0, ir1, ir2,
     >      sn0(j,i,ia), ssn(j,i,ia), ssp(j,i,ia), swi(j,i), swe(j,i)
      endif
      enddo

      return
      end
