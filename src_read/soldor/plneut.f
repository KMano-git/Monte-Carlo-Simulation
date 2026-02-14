!**********************************************************************
      subroutine plneut
!**********************************************************************
!
!       analyical model for neutral transport
!
!     plpflx    nw = 1(sol), 2(idp), 3(prv), 4(odp), 5(man)
!     plneut    neutral source from idp & odp
!
!   iflx       : 1:(sol)/2:(idp)/3:(prv)/4:(odp)/5:(man)
!              : 6:(puf)/7:(vol)      see  sub. ntnflx
!
!----------------------------------------------------------------------
      use cntcom, only : cstyp, eion, eioni, icps, iplx, iply, npep
     >    , npsp
      use cntmnt, only : flexh, flpmp, mfmax, pfl_abs, pfl_err, pfl_ion
     >    , pfl_man, pfl_ntl, pfl_pmp, pfl_src,  sn0, sumsn, sumsp
     >    , sumwe, sumwi
      use cphcns, only : cev
      use cplcom, only : ama, flps, nion, ssnc, sspc, swec, swic, tfps
     >    , vne, vte
      use cplmet, only : hare, hvol, jcel, jtmax, kreg
      use csize,  only : ndsp, ndx
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  ia, iw, ic, j, i, jmax, jmax1, inc, jt, jw
      integer  m, jtm, isrc, iflx, iout
      integer  nsza, nsiz, ir, la
      real*8   alnx(ndx), wn0(ndx), wsi(ndx)
      real*8   tflx, eflx, tsi, twi, twe, esi, zsi, ewi, ewe
      real*8   zar, zf0, zv0, zn0, zfa, tln, trm, dln, zsp
      real*8   zne, zte, zal, zal2, zsgv, rmd, zlmt, zfb
      real*8   famp, zeta, znrm, zwi, zwe, zvl
      real*8   zsum1, zsum2
!::const
!   Note   not send famp & eion data from NET2D
      zlmt = 20.0d0
      famp = 25.0d0   ! flux amplification factor
      zeta = 1.0d0 - 1.0d0/famp
      eion = 5.0d0
      eioni = 40.0d0
!
!::pfl_ion/ntl for multi ion species
      do m = 1, mfmax
      zsum1 = 0.0d0
      do ia = 1, nion
      zsum1 = zsum1 + tfps(m,ia)
      enddo
      pfl_ion(m) = zsum1
      pfl_ntl(m) = 0.0d0
      if( m.eq.2 .or. m.eq.4 ) then
      pfl_ntl(m) = zsum1
      endif
      enddo
!
!::ion species
      do ia = 1, nion
!
!::inner/outer plates
      isrc = 0
      do iflx = 2, 4, 2
      isrc = isrc + 1
      if( iflx.eq.2 ) write(cstyp,'(a,i2.2)') "idp_",ia
      if( iflx.eq.4 ) write(cstyp,'(a,i2.2)') "odp_",ia
!
!::total source terms for each source
      nsza = 10*ndsp
      nsiz = 10
      call setd( sumsn, nsza, 0.0d0 )
      call setd( sumsp, nsza, 0.0d0 )
      call setd( sumwi, nsiz, 0.0d0 )
      call setd( sumwe, nsiz, 0.0d0 )
!
      tflx = tfps(iflx,ia)
      tsi = 0.0d0
      twi = 0.0d0
      twe = 0.0d0
!::wall segment
      do iw = npsp(iflx), npep(iflx)-1
      eflx = flps(iw,ia)    !  > 0.0
      ic = icps(iw)
      j  = iplx(ic)  !  polidal number
      i  = iply(ic)  !  tube number
      jmax  = jtmax(i)
      jmax1 = jmax-1
!
      iout = 0
      call plenb(i,alnx)
!
!::plate information
      if( iflx.eq.4 ) then
        inc = +1
        jt  =  1
        zar = hare(1,i)
      elseif( iflx.eq.2 ) then
        inc = -1
        jt  = jmax
        zar = hare(jmax1,i)
      endif
      zf0 = eflx
      zv0 = dsqrt(10.0d0*cev/ama(ia))
      zn0 = zf0/(zv0*zar)
!
      call setd(wn0,ndx,0.0d0)
      call setd(wsi,ndx,0.0d0)
      zfa = zf0    !  S*N0*V0 at before surface
      tln = 0.0d0  !  poloidal length
      trm = 0.0d0  !  integ(dl/rmd)
      esi = 0.0d0  !  integ(Ne*N0*<sigv>)
!
      do jw = 2, jmax1
      jt  = jt + inc
      jtm = jt - 1
      j   = jcel(jt,i)
      dln = alnx(jt) - alnx(jtm)
      tln = tln + dln
      zne = vne(j,i)
      zte = vte(j,i)
!-----
      zal = zte/10.0d0
      zal2= zal**2
      zsgv= 3.0d-14*zal2/(3.0d0+zal2)
!-----
      rmd = zv0/(zne*zsgv)
      trm = trm + dln/rmd
      if( trm.gt.zlmt ) goto 110
      zfb = zf0*dexp(-trm)
      zvl = hvol(j,i)  ! soldor
      zsi = zfa - zfb
      zsi = zsi/zvl
      zn0 = zsi/(zne*zsgv)
      zfa = zfb
      esi = esi + zsi*zvl
      wn0(j) = zn0
      wsi(j) = zsi
      enddo  ! loop_(jw)
 110  continue
!
!::normalization for each flux tube
      znrm = zeta/esi*(eflx/tflx)
!
!::source terms
!
!::debug write
      if(iout.eq.1) then
      write(n6,'(2x)')
      write(n6,'(2x,1x,"iflx",3x,"iw",4x,"j",4x,"i",3x,"jt",4x,
     >"alnx",8x,"dln",9x,"N0",10x,"Ne",10x,"zsi",9x,"esi",9x,"eflx")')
      endif
!
      esi = 0.0d0
      ewi = 0.0d0
      ewe = 0.0d0
      do jt = 2, jmax1
      j  = jcel(jt,i)
      wn0(j) = wn0(j)*znrm
      wsi(j) = wsi(j)*znrm
!-----
      zn0 = wn0(j)*tflx
      zsi = wsi(j)*tflx
      zsp = 0.0d0
      zwi = -zsi*eioni*cev
      zwe = -zsi*eion *cev
      esi = esi + zsi*hvol(j,i)
      ewi = ewi + zwi*hvol(j,i)
      ewe = ewe + zwe*hvol(j,i)
!-----
      sn0(j,i,ia) = sn0(j,i,ia) + zn0
      ssnc(j,i,ia) = ssnc(j,i,ia) + zsi
      sspc(j,i,ia) = sspc(j,i,ia) + zsp
      swic(j,i) = swic(j,i) + zwi
      swec(j,i) = swec(j,i) + zwe
!-----
      ir = kreg(j,i)
      sumsn(ir,ia) = sumsn(ir,ia) + zsi*hvol(j,i)
      sumsp(ir,ia) = sumsp(ir,ia) + zsp*hvol(j,i)
      sumwi(ir) = sumwi(ir) + zwi*hvol(j,i)
      sumwe(ir) = sumwe(ir) + zwe*hvol(j,i)
!-----
      if( iout.eq.1 ) then
      dln = alnx(jt)-alnx(jt-1)
      write(n6,'(2x,5i5,1p12e12.3)')
     >  iflx, iw, j, i, jt, alnx(jt), dln, zn0, zne, zsi, esi, eflx
     > ,ewi, ewe
      endif
      enddo  ! loop_(jt)
!
      if( iout.eq.1 ) then
      write(n6,'(2x,"it =",i3,"  eflx =",1pe12.3,"  esi =",1pe12.3,
     >   "  eta =",1pe12.3,"  ewi,ewe =",1p2e12.3)')
     >   i, eflx, esi, esi/eflx, ewi, ewe
      endif
!
      tsi = tsi + esi
      twi = twi + ewi
      twe = twe + ewe
      enddo  ! loop_(iw)
!::conv-check
      call plneut_conv(iflx,ia)
!
      zsum1 = 0.0d0
      zsum2 = 0.d0
      do la = 1, nion
      do ir = 1, 6
      zsum1 = zsum1 + sumsn(ir,la)
      enddo
      zsum2 = zsum2 + sumsn(7,la)
      enddo
!
!::save
      pfl_src(iflx) = pfl_src(iflx) + zsum1
      pfl_man(iflx) = pfl_man(iflx) + zsum2
      pfl_abs(iflx) = pfl_abs(iflx) + flexh
      pfl_pmp(iflx) = pfl_pmp(iflx) + flpmp
!
      enddo  ! loop_(m)   [iflx]
      enddo  ! loop_(ia)
!
      do m = 1, mfmax
      pfl_err(m) = pfl_ntl(m)-(pfl_src(m)+pfl_man(m)+pfl_abs(m))
      enddo
!
      end
!
!**********************************************************************
      subroutine plneut_conv(nw,ka)
!**********************************************************************
      use cntcom, only : cstyp,  tflex, tflpm
      use cntmnt, only : dotn, dotn2, flexh, flpmp, ift => i6_src
     >    , sdotn, sdotn2, sumsn, sumsp, sumwe, sumwi, totmsn, totmsp
     >    , totmwe, totmwi, totsn, totsp, totwe, totwi
      use cplcom, only : nion, q2a, q3, q4, ssnc, ssnv, sspc, sspv, swec
     >    , swev, swic, swiv, tfps, vna
      use cplmet, only : hvol, icel, itmax, itpve, jcel, jtmax, kreg
      use csize,  only : ndsp
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nw, ka
      integer, intent(in) :: nw, ka
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nsza, nsiz, it, jmax, jt, j, i, ir, ia, ift
! modified 1/1 lines replace all include files with module files by kamata 2021/08/18
!ik   integer  nsza, nsiz, it, jt, j, i, ir, ia, ift
      integer  nsza, nsiz, it, jt, j, i, ir, ia
      real*8   tflx
!
      tflx = tfps(nw,ka)
! deleted 1 line replace all include files with module files by kamata 2021/08/18
!ik   ift = i6_src
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
!::except for core region ( ir=1,2,..6 )
      do it = 2, itmax-1
      if( it.eq.itpve ) cycle
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   jmax = jtmax(it)
!ik   do jt = 2, jmax-1
      do jt = 2, jtmax(it)-1
      j  = jcel(jt,it)
      i  = icel(jt,it)
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
!::hot core region
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
      dotn  = tflx
      flexh = dotn - sumsn(10,ka)
      flpmp = flexh
      tflex = tflex + flexh
      tflpm = tflpm + flpmp
!
!----------------------------------------------------------------------
!::actual dotn
!----------------------------------------------------------------------
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   ia = ka
!ik   dotn2  = sumsn(10,ia) + flexh
      dotn2  = sumsn(10,ka) + flexh
      sdotn  = sdotn  + dotn
      sdotn2 = sdotn2 + dotn2
!
!----------------------------------------------------------------------
!::debug write
!---------------------------------------------------------------------
      if( ift.eq.0 ) return
!
      do ia = 1, nion
!x      write(ift,'(2x,"=== plntsr Now ===",i2,2x,a,"  dotn =",1pe12.4,
!x     > "  wtot =",1p2e14.6,"  werr =",1pe14.6,"  wion/wabs =",
!x     > 1p2e14.6)') wisrc,wcsrc,dotn,wtot,wsum,werr,wion,wabs
!x
      write(ift,'(2x,"=== plntsr Now ===  nw =",i2,"  ka =",i2," ia =",
     >     i2)') nw, ka, ia
      write(ift,'(2x,"typ",2x,"src",5x,"dotn",7x,"sum",11x,"odp",8x,
     > "sol",8x,"idp",8x,"opv",8x,"ipv",8x,"edg",8x,"hot",8x,
     > "abs",8x,"pmp",8x,"err")')
 602  format(2x,a6,2x,a3,1pe11.3,1pe14.6,1p12e11.3)
!
      write(ift,602) cstyp(1:6),"Sna",dotn,dotn2,
     >   (sumsn(ir,ia),ir=1,7),flexh,flpmp
      write(ift,602) cstyp(1:6),"Swi",dotn,sumwi(10),(sumwi(ir),ir=1,7)
      write(ift,602) cstyp(1:6),"Swe",dotn,sumwe(10),(sumwe(ir),ir=1,7)
!
      write(ift,602) "total","Sna",sdotn,sdotn2,
     >   (totsn(ir,ia),ir=1,7), tflex,tflpm
      write(ift,602) "total","SWi",sdotn2,totwi(10),(totwi(ir),ir=1,7)
      write(ift,602) "total","SWe",sdotn2,totwe(10),(totwe(ir),ir=1,7)
      enddo
!
      return
      end
