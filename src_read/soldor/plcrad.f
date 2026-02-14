!**********************************************************************
      subroutine plcrad
!**********************************************************************
!
!   non-corona model   No output of crad for @outntl_xx @outimp_xx
!
!----------------------------------------------------------------------
      use cntcom, only : iplx, iply, mrgn, ncmax
      use cntmnt, only : dotn, i6_src, totwe
      use cplcom, only : mdl_wrd, nlp, nlpmn, q3, q4, swec, swev, swic
     >    , swiv, trad, trad_max, trad_min, tradrg, wdnz, wrad
      use cplmet, only : hvol, icel, itmax, itpve, jcel, jtmax, kreg
      use cplqcn, only : mrgnp
      use cplwrd, only : wcre, wfac, wime, wimi, wmci, witim, wnlp
     >    , wtime
      use csize,  only : ndmc, ndx, ndy
      use csonic, only : itim, time
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 16/3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nsiz, it, jmax, icx, icy, ir, i, ic, icxd, icxd1, icxd2
!ik   integer  icxt1, icxt2, lprv, is1(ndy), is2(ndy), icm
!ik   integer  jw, j, ii, j1, j2, j3, jws, jwe, jwsv(ndx), l, itm, jt
!ik   integer  ift, ldbg
!ik   real*8   zte, zne, zcn, zned, zsum, fac, fcmd
!ik   real*8   ztem1, ztem2, ztem, zpt1, zpt2, zpt3, hh
!ik   real*8   zft, cfl, bet
!ik   real*8   zwrc, zwrv, zrad
!ik   real*8   zne_min(10), zne_opv(ndy), zne_ipv(ndy)
!ik   real*8   zx, zy, zlz, zne_nz, zne_wr
!ik   real*8   cfl1, znimp1, zwr1, cfl2, znimp2, zwr2, znimp, zwr
!ik   real(8) :: zwe, zwec, zwev, zwi, hq4, hwr1, hwr2
!ik   real(8), dimension(ndmc) :: wcr1, wcr2
!ik   real(8) :: smwrd1, smwrd2, hwsm1, hwsm2
!ik   real(8), dimension(10) :: hwrg1(10), hwrg2(10)
!ik   real(8) :: swrd(ndmc)
      integer  i, ic, ift, ir, ir2, it, j, jmax, jt, ldbg, nsiz
      real(8)  fac, hq4, hwrg1(10), hwsm1, tot_nte, tot_nti
     >       , znimp, zsum, zwe, zwec, zwev
      real(8),dimension(10) :: wrg_nti, wrg_nte  ! ntl
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8   tot_nti, tot_nte
!ik   integer :: ir2
!
      ift = i6_src
!
      nsiz = ndx*ndy
      call setd( wrad, nsiz, 0.0d0 )
      call setd( wdnz, nsiz, 0.0d0 )
      call setd( tradrg, 10, 0.0d0 )
!
!::clear  radiation in MC-code
!xxx             see sub. CHK_plcrad_mc.f
      wcre(1:ndmc) = 0.0d0
      wmci(1:ndmc) = 0.0d0
!KH      wmce(1:ndmc) = 0.0d0
      wimi(1:ndmc) = 0.0d0
      wime(1:ndmc) = 0.0d0
!KH      twrd(1:ndmc) = 0.0d0
      wrg_nti(1:10) = 0.0d0
      tot_nti = 0.0d0
      tot_nte = 0.0d0
!
!KH191120 copy from plcrad_mcN.f
!::energy loss due to neutral particles
      do ic = 1, ncmax    ! KSFUJI
       j  = iplx(ic)
       i  = iply(ic)
       if( j.le.0 .or. i.le.0 ) cycle
       ir2 = mrgnp(ic)
       if( ir2.le.0 ) cycle
       wrg_nti(ir2)=wrg_nti(ir2)+(swic(j,i)+swiv(j,i)*q3(j,i))*hvol(j,i)
       wrg_nte(ir2)=wrg_nte(ir2)+(swec(j,i)+swev(j,i)*q4(j,i))*hvol(j,i)
      enddo
!
!::corona Model
      call plwrdcrG(wcre)               ! general version
!
!::debug write
!xx      if( mod(itim,100) == 0 .and. nlp == nlpmn ) then
!xx        write(n6,'(/2x,"itim =",i6,"  nlp =",i3)') itim, nlp
!xx        do ic = 1, ncmax, 100
!xx        write(n6,'(2x,"CHK W,Ar ",i7,1p4e16.8,1p3e11.3)')
!xx     >     ic, wcr1(ic), wcr_wrd(ic,1), wcr2(ic), wcr_wrd(ic,2)
!xx     >       , wcr1(ic)-wcr_wrd(ic,1), wcr2(ic)-wcr_wrd(ic,2)
!xx     >       , wcr1(ic)+wcr2(ic)-swrd(ic)
!xx        enddo
!xx      endif
!
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   smwrd1 = 0.0d0  ! corona Carbon cimprg
!ik   smwrd2 = 0.0d0  ! corona Argon  cimprg2
!
!::radiation  < 0.0 in soldor
      do ic = 1, ncmax
      j = iplx(ic)
      i = iply(ic)
      ir  = mrgn(ic)
      ir2 = mrgnp(ic)
      if( j.le.0 .or. i.le.0 ) cycle
!-----
!xx   wcre(ic) = wcr1(ic)
!xx   wmce(ic) = wcr2(ic)
!xx   wcre(ic) = wcr1(ic)+wcr2(ic)  ! old varsion
!-----
!xx    wcre(ic) = wcr_wrd(ic,1)  ! W
!xx    wmce(ic) = wcr_wrd(ic,2)  ! Ar
!xx    wcre(ic) = swrd(ic)
!-----
!KH      wmce(ic) = 0.0d0
      znimp = 0.0d0          ! temporary
!
!KH      wrad(j,i) = wcre(ic) + wmce(ic)
      wrad(j,i) = wcre(ic)
      wdnz(j,i) = znimp
      enddo
!
!::total radiation
      call wrintg(wcre,hwsm1,hwrg1)
!KH      call wrintg(wmce,hwsm2,hwrg2)
!KH      zsum = - (hwsm1+hwsm2)
      zsum = -hwsm1
!
!----------------------------------------------------------------------
!::normalization
!----------------------------------------------------------------------
      fac = 1.0d0
      if( mdl_wrd.eq.1 ) then
      if( zsum.le.trad_min ) fac = trad_min/zsum
      if( zsum.ge.trad_max ) fac = trad_max/zsum
      endif
!
!----------------------------------------------------------------------
!::assign Wrad values in common variables
!----------------------------------------------------------------------
      do ic = 1, ncmax
      j = iplx(ic)
      i = iply(ic)
      ir  = mrgn(ic)
      ir2 = mrgnp(ic)
      if( j.le.0 .or. i.le.0 ) cycle
!
!::normalization
      wcre(ic) = wcre(ic)*fac
!KH      wmce(ic) = wmce(ic)*fac
!KH      wime(ic) = wcre(ic) + wmce(ic)
      wime(ic) = wcre(ic)
!KH      twrd(ic) = -wime(ic)  ! twrd > 0  in IMPMC
!
      wrad(j,i) = wrad(j,i)*fac
      tradrg(ir) = tradrg(ir) + wrad(j,i)*hvol(j,i)
!
!::radiation
      if(ir2 < 7)then
      hq4  = q4(j,i)
      zwe  = wime(ic)
      zwec = 0.0d0
      if(hq4 .eq. 0) then
        if(zwe .eq. 0) then
          zwev = 0.0
        else
          call wexit("plcrd","hq4=0")
        endif
      else
        zwev = zwe/hq4
      endif
      swec(j,i)  = swec(j,i) + zwec
      swev(j,i)  = swev(j,i) + zwev
      endif
!
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   zwi = 0.0d0
      enddo
!
      wtime = time
      witim = itim
      wnlp  = nlp
      wfac  = fac
!
!----------------------------------------------------------------------
!::total radiation
!----------------------------------------------------------------------
      call wrintg(wcre,hwsm1,hwrg1)
!KH      call wrintg(wmce,hwsm2,hwrg2)
!KH      trad = hwsm1 + hwsm2
      trad = hwsm1
!
!::conv-check
      nsiz = 10
      call setd( totwe, nsiz, 0.0d0 )
      do it = 2, itmax-1
      if( it.eq.itpve ) cycle
      jmax = jtmax(it)
      do jt = 2, jmax-1
      j  = jcel(jt,it)
      i  = icel(jt,it)
      ir = kreg(j,i)
      totwe(ir) = totwe(ir) + (swec(j,i)+swev(j,i)*q4(j,i))*hvol(j,i)
      enddo
      enddo
!
      totwe(10) = 0.0d0
      do ir = 1, 7
      totwe(10) = totwe(10) + totwe(ir)
      enddo

      do ir2 = 1, 6  !  6:core edge 7:hot core
      tot_nti  = tot_nti  + wrg_nti(ir2)
      tot_nte  = tot_nte  + wrg_nte(ir2)
      enddo
!
!::debug write
      ift = i6_src
      ift = n6
      if( ift.eq.0 ) return
!
      ldbg = 0
      if( mod(itim,50).eq.0 .and. nlp.eq.nlpmn ) ldbg = 1
!
      if( ldbg.eq.1 ) then
      write(ift,602) "NEUT","Swi","--",dotn,1.0d0,tot_nti,
     >                                  (wrg_nti(ir2),ir2=1,7)
      write(ift,602) "NEUT","Swe","--",dotn,1.0d0,tot_nte,
     >                                  (wrg_nte(ir2),ir2=1,7)
      write(ift,602) "cradCR","Swe","--",dotn,fac,hwsm1,
     >                                  (hwrg1(ir),ir=1,7)
!xx   write(ift,602) "radMC","Swe",dotn,fac,hwsm2,(hwrg2(ir),ir=1,7)
 602  format(2x,a6,2x,a3,2x,a2,x,1pe11.3,1pe11.3,1pe14.6,1x,1p12e11.3)
      endif
!
!::debug write  (see plsorc)
!x      if( ift.gt.0 ) then
!xxx  write(ift,'(2x,a,2x,1p10e11.2)') "cimprg", cimprg
!x      write(ift,602) "cradLz","Swe",dotn,trad,(tradrg(ir),ir=1,7)
!x      write(ift,602) "total","SWe",sdotn2
!x     >  ,crad_totwe(10),(crad_totwe(ir),ir=1,7)
!x 602  format(2x,a6,2x,a3,1pe11.3,1pe14.6,1x,1p12e11.3)
!x      endif
!
      return
      end
