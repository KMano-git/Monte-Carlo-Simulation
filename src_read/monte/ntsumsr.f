!**********************************************************************
      subroutine ntsumsr(tflux)
!**********************************************************************
!
!       source terms of plasma in neut2d
!                                                   2006/04/16
!
!----------------------------------------------------------------------
      use cntcom, only : cstyp, den0, deng, eng0, engg, iplx, iply, mrgn
     >    , ncmax, ncmax2, tflex, tflpm, vlp0, volm
      use cntmnt, only : dotn, dotn2, i6_src, flexh, flpmp, sdotn
     >    , sdotn2, sumsn, sumsp, sumwe, sumwi
      use cntsrc, only : tden0, tdeng, teng0, tengg, trgsn, trgsp, trgwe
     >    , trgwi, tssn, tssp, tswe, tswi, tvlp0
      use cntwcn, only : swabs, swerr, swnrm, swpmp, swreg, wabs, wcsrc
     >    , werr, wion, wisrc, wnrm, wpmp, wreg, wssn, wssp, wsum, wswe
     >    , wswi, wtot
      use cplcom, only : nion, nlp
      use cplmet, only : kreg
      use cplqcn, only : mrgnp
      use csize,  only : ndsp
      use csonic, only : itim, time
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik    real*8  tflux
      real(8), intent(in) :: tflux
!
!::local variables
      integer  nsza, nsiz, ic, ia, ir, ir0, ir1, ir2, ig
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik    real*8   zn0, zsn, zsp, zwi, zwe, zfc
!ik    integer :: i6, ii, j, i, ift, ift2
      real*8   zsn, zsp, zwi, zwe, zfc
      integer :: ii, j, i, ift, ift2
!
      ift = i6_src
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
!::neutral density and temperature
!-----------------------------------------------------------------------
      do ic = 1, ncmax2
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   do ia = 1, nion
!ik   ig = ia
      do ig = 1, nion
      tden0(ic,ig) = tden0(ic,ig) + den0(ic,ig)
      teng0(ic,ig) = teng0(ic,ig) + den0(ic,ig)*eng0(ic,ig)
      tvlp0(ic,ig) = tvlp0(ic,ig) + den0(ic,ig)*vlp0(ic,ig)
      enddo
      tdeng(ic,1) = tdeng(ic,1) + deng(ic,1)
      tengg(ic,1) = tengg(ic,1) + deng(ic,1)*engg(ic,1)
      tdeng(ic,2) = tdeng(ic,2) + deng(ic,2)
      tengg(ic,2) = tengg(ic,2) + deng(ic,2)*engg(ic,2)
      enddo
!
!::debug write
!x    write(i6,'(2x,"ic     =",10(i7,5x))') (ic,ic=1000,ncmax2,1000)
!x    write(i6,'(2x,"den0   =",1p10e12.3)')
!x   >       (den0(ic,1),ic=1000,ncmax2,1000)
!x    write(i6,'(2x,"tden0  =",1p10e12.3)')
!x   >       (tden0(ic,1),ic=1000,ncmax2,1000)
!
!-----------------------------------------------------------------------
!::source
!-----------------------------------------------------------------------
      do ic = 1, ncmax
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   ir  = mrgn(ic)
      ir2 = mrgn(ic)
!xx   ir2 = mrgnp(ic)     !  tube (itmpe) is out of system
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   ir2 = ir            !  edit mrgn  see sub. ntmrgn
      zfc = dotn/swnrm/volm(ic)
!
!::particle & momentum source
      do ia = 1, nion
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   ig = ia                      !  temporary  2003/02/23
!ik   zn0 = zfc*wden(ic,ig)
      zsn = zfc*wssn(ic,ia)
      zsp = zfc*wssp(ic,ia)
      tssn(ic,ia) = tssn(ic,ia) + zsn
      tssp(ic,ia) = tssp(ic,ia) + zsp
!
!::total
      sumsn(ir2,ia) = sumsn(ir2,ia) + zsn*volm(ic)
      sumsp(ir2,ia) = sumsp(ir2,ia) + zsp*volm(ic)
      trgsn(ir2,ia) = trgsn(ir2,ia) + zsn*volm(ic)
      trgsp(ir2,ia) = trgsp(ir2,ia) + zsp*volm(ic)
      enddo  ! loop(ia)
!
!KH160315      write(n6,*)"ntsumsr",ic,tssn(ic,1)
!
!::ion & ele energy source
      zwi = zfc*wswi(ic)
      zwe = zfc*wswe(ic)
      tswi(ic) = tswi(ic) + zwi
      tswe(ic) = tswe(ic) + zwe
!
!::total
      sumwi(ir2) = sumwi(ir2) + zwi*volm(ic)
      sumwe(ir2) = sumwe(ir2) + zwe*volm(ic)
      trgwi(ir2) = trgwi(ir2) + zwi*volm(ic)
      trgwe(ir2) = trgwe(ir2) + zwe*volm(ic)
      enddo
!
!::whole region
      do ir = 1, 7
      do ia = 1, nion
      sumsn(10,ia) = sumsn(10,ia) + sumsn(ir,ia)
      sumsp(10,ia) = sumsp(10,ia) + sumsp(ir,ia)
      enddo
      sumwi(10) = sumwi(10) + sumwi(ir)
      sumwe(10) = sumwe(10) + sumwe(ir)
      enddo
!
      do ia = 1, nion
      zsn = 0.0d0
      zsp = 0.0d0
      do ir = 1, 7
      zsn = zsn + trgsn(ir,ia)
      zsp = zsp + trgsp(ir,ia)
      enddo
      trgsn(10,ia) = zsn
      trgsp(10,ia) = zsp
      enddo
      zwi = 0.0d0
      zwe = 0.0d0
      do ir = 1, 7
      zwi = zwi + trgwi(ir)
      zwe = zwe + trgwe(ir)
      enddo
      trgwi(10) = zwi
      trgwe(10) = zwe
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

!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
!xx      i6 = 0
!xx      if( itim == 20 .and. nlp > nlpmn ) i6 = 22000 + mype
!xx      ift = i6
      if( ift.eq.0 ) return

!::fort.22001  totSn(irg),totWi in ntsumsr
      ia = 1
      write(ift,'(/2x,"*** ntsumsr(NTL) ***  itim =",i7,"  time =",
     >  1pe14.6)') itim, time
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
     >   (trgsn(ir,ia),ir=1,7), tflex,tflpm
      write(ift,602) "total","Spa",sdotn2,
     >                         trgsp(10,ia),(trgsp(ir,ia),ir=1,7)
      write(ift,602) "total","SWi",sdotn2,trgwi(10),(trgwi(ir),ir=1,7)
      write(ift,602) "total","SWe",sdotn2,trgwe(10),(trgwe(ir),ir=1,7)

!::fort.21051  N0(ic),Sn,Wi in ntsumsr
      ift2 = 0
      if( ift2 == 0 ) return

      ift2 = ift + 50
      write(ift2,'(/2x,"*** ntsumsr ***  itim =",i7,"  time =",
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
     >     tden0(ic,ia), tssn(ic,ia), tssp(ic,ia), tswi(ic), tswe(ic)
      endif
      enddo

      return
      end
