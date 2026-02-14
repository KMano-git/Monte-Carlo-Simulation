!***********************************************************************
      subroutine plsorc(kgrp)
!***********************************************************************
!
!    source-term linearlization
!
!         ssn = ssnc + ssnv*na
!         ssp = sspc + sspv*q2a
!         swi = swic + swiv*q3*cev
!         swe = swec + swev*q4*cev
!
!    references
!
!          Patankar "Numerical Heat Transfer and Fluid Flow"
!                 p.48 and p.143
!
!          B2-EIRENE code
!              j108b/b2/src/neutrals/neutrals.F     snic,sniv,rf0,t0
!
!    Note ssp   neutrals/neutrals.F ( see 511 )
!
!      ssp = {Sp+bet*|Sp|} - bet*|Sp|/va(l)*va(l+1)
!          =  sspc       + sspv*va(l+1)
!      bet = alf*va(l)/vthe
!      vthe = sqrt(Te/me)/100 = 4.19e5*sqrt(Te_eV)/100
!
!        ss = ssc + ssv*Q              !   2002/11/13 (plsorc)
!        recomb-sorce                  !   2003/01/29 (plrecm)
!        linearlized for eqch sources  !   2003/02/21
!
!   Note
!     neutral source calculated using by pfl_ntl(iflx)
!     n6_src   = 0 (no output)  = 81 (file)  ==> ift
!
!     kgrp = 1  (soldor) plntsr  :  Sn, Sp, Wi, We
!          = 2  (neut2d) ntsumsr :  N0, E0  (Sn, Sp, Wi, We)
!          = 3  (plotp)
!     add the function of N0 and E0 in sub. plntsr
!                        2015/06/25  K.Shimizu
!
!::KSFUJI
!     kgrp = 0  (first cal. of impmc)
!
!::recombination
!      vli         vol
!    H+ + e    =>   H0
!
!   reaction       monte carlo cal.
!   no neutral      N0,  Si, Wi, We
!    N0 = 0         zsnA = dotN/Wemit/Vol(ic) * wssn(ic) = Ne*N0*<sigv>
!                   zwiA = dotN/Wemit/Vol(ic) * wswi(ic)
!
! wssn(ic) = zsnA /dotN*Wemit*Vol(ic)   zsnB = Ne*N0*<sigv>rec*Vol(ic)
!   = zsnA*Vol(ic) /dotN*Wemit                 flvl(ic) = zsnB in plpflx
!   = zsnB /dotN*Wemit
!   = flvl(ic)/dotN*Wemit
!             in plsorc
!
!   Sn = Ne*N0*<sigv> = dotN/Wemit/Vol(ic)*wssn(ic)
!                            in plntsr
!
!   total neutral flux at cal. of source terms (ssnc ssnv) /cmsorc/
!     local tfsrc => snflx(iflx) /cmnflx
!
!----------------------------------------------------------------------
      use cntcom,     only : cstyp, lntmd, tflex, tflpm
      use cntmnt,     only : dotn, i6_src, flexh, flpmp, mcflx, mfmax
     >    , pfl_abs, pfl_err, pfl_ion, pfl_man, pfl_ntl, pfl_pmp
     >    , pfl_src, psm_abs, psm_dt, psm_err, psm_fl, psm_ion, psm_man
     >    , psm_ntl, psm_pmp, psm_si, psm_src, sdotn, sdotn2, sn0, ssn
     >    , ssp, sumsn, swe, swi, tnpuf, totmsn, totmsp, totmwe, totmwi
     >    , vcsrc, vflux, visrc, vnsmp
      use cntsrc,     only : tden0, tdeng, teng0, tengg, trgsn, trgsp
     >    , trgwe, trgwi, tssn, tssp, tswe, tswi, tvlp0
      use cntwcn,     only : wcsrc, wisrc
      use cplcom,     only : cimp, crad_totwe, nion, nlp, rcydt, rcysp
     >    , rcytm, snflx, ssnc, ssnv, sspc, sspv, swec, swev, swic, swiv
     >    , trad, tradrg
      use csize,      only : ndgs, ndmc, ndmfl, ndsp, ndx, ndy
      use csonic,     only : dtim, itim, limp, lntl, time
      use cunit,      only : n6
      use topics_mod, only : lgtpc
      implicit none
!
!::argument
      integer, intent(in) :: kgrp
!
!::local variables
      integer  k1, k2, k3
      integer  ia, i, ift
      integer  m, isrc, iflx, ir
      real*8   tfion, tfntl, tfsrc
      real*8   zsum1, zsum2
!
!::Note
!     please set varibales  pfl_ion & pfl_ntl
!
      ift = i6_src
!
!::output file
      if( ift.gt.0 ) then
      write(ift,'(/2x,"*** plsorc ***  START   time =",1pe14.6,
     >  "  dtim =",1pe11.2,"  itim =",i6,"  nlp =",i3)')
     >   time,dtim,itim,nlp
      endif
!
!
!----------------------------------------------------------------------
!::zero clear
!----------------------------------------------------------------------
      do k3=1,nion
        do k2=1,ndy
          do k1=1,ndx
            sn0 (k1,k2,k3) = 0.0d0
            ssn (k1,k2,k3) = 0.0d0
            ssp (k1,k2,k3) = 0.0d0
            ssnc(k1,k2,k3) = 0.0d0
            ssnv(k1,k2,k3) = 0.0d0
            sspc(k1,k2,k3) = 0.0d0
            sspv(k1,k2,k3) = 0.0d0
          enddo
        enddo
      enddo
!
      do k2=1,ndy
        do k1=1,ndx
          swi (k1,k2) = 0.0d0
          swe (k1,k2) = 0.0d0
          swic(k1,k2) = 0.0d0
          swiv(k1,k2) = 0.0d0
          swec(k1,k2) = 0.0d0
          swev(k1,k2) = 0.0d0
        enddo
      enddo
!
      sdotn = 0.0d0
      sdotn2= 0.0d0
      tflex = 0.0d0
      tflpm = 0.0d0
      tnpuf = 0.0d0
!
!::total source in hot core region
      do k1=1,nion
        totmsn(k1) = 0.0d0
        totmsp(k1) = 0.0d0
      enddo
      totmwi = 0.0d0
      totmwe = 0.0d0
!
!::clear /cntpfl/
      do k1=1,ndmfl
        pfl_src(k1) = 0.0d0
        pfl_man(k1) = 0.0d0
        pfl_abs(k1) = 0.0d0
        pfl_pmp(k1) = 0.0d0
        pfl_err(k1) = 0.0d0
      enddo
      psm_ion = 0.0d0
      psm_ntl = 0.0d0
      psm_src = 0.0d0
      psm_man = 0.0d0
      psm_abs = 0.0d0
      psm_pmp = 0.0d0
      psm_err = 0.0d0
!
!----------------------------------------------------------------------
!::zero clear  /cntsrc/
!----------------------------------------------------------------------
      if( kgrp.eq.2 ) then
        do k1 = 0, ndmc
          do k2 = 1, ndgs
            tden0(k1,k2) = 0.0d0
            teng0(k1,k2) = 0.0d0
            tvlp0(k1,k2) = 0.0d0
          enddo
          tdeng(k1,1) = 0.0d0
          tdeng(k1,2) = 0.0d0
          tengg(k1,1) = 0.0d0
          tengg(k1,2) = 0.0d0
        enddo
!
        do k1 = 0, ndmc
          do k2 = 1, ndsp
            tssn(k1,k2) = 0.0d0
            tssp(k1,k2) = 0.0d0
          enddo
          tswi(k1) = 0.0d0
          tswe(k1) = 0.0d0
        enddo
!
        do k1 = 1, 10
          do k2 = 1, ndsp
            trgsn(k1,k2) = 0.0d0
            trgsp(k1,k2) = 0.0d0
          enddo
          trgwi(k1) = 0.0d0
          trgwe(k1) = 0.0d0
        enddo
      endif
!
!::KSFUJI
      if( kgrp.eq.0 ) then
        write(n6,'(2x,"KSFUJI  call plsorc(0) clear  swic,swiv")')
        return
      endif
!
!----------------------------------------------------------------------
!::Neutral Source  analytical model
!----------------------------------------------------------------------
      if( lntl.eq.0 ) return
      if( lntmd.eq.0 ) then
        call plneut   ! analytical model
        goto 100
      endif
!
!----------------------------------------------------------------------
!::Neutral Source  monte-carlo model
!----------------------------------------------------------------------
      do m = 1, mfmax
        iflx  = m
        isrc  = visrc(iflx)
        cstyp = mcflx(iflx)
        tfion = pfl_ion(iflx)
        tfntl = pfl_ntl(iflx)
!-----
!::suppress the initial high recycling
        rcytm = 1.0d0
        if( time.le.rcydt ) then
          rcytm =  dmin1( 1.0d0, rcysp+(1.0d0-rcysp)*(time/rcydt) )
        endif
        tfntl = tfntl*rcytm
        tfion = tfion*rcytm
!-----
        tfsrc = tfntl
        if( iflx.eq.8 ) tfsrc = tfion
!-----
        if(rcytm.ne.1.0d0 .and. mod(itim,20).eq.0 .and. nlp.eq.6) then
          write(n6,'(2x,"plsorc",2x,i7,1pe12.3,"  iflx =",
     >    i2,"  tfion,tfntl,tfsrc =",1p3e12.3,"  rcytm =",1pe12.3)')
     >    itim,time,iflx,tfion,tfntl,tfsrc,rcytm
        endif
!--------
        if( ift.gt.0 ) then
          write(ift,'(2x)')
          if( isrc.gt.0 ) then
            write(ift,'(2x,"=== plsorc ===   iflx/isrc =",
     >      i2,"/",i2,2x,a,2x,a,"  flux =",1pe14.6," =>",1pe14.6,
     >      "  rcytm =",0pf6.3,"  nsmp =",i7)')  iflx,isrc,
     >      mcflx(iflx),vcsrc(isrc),vflux(isrc),tfsrc,rcytm,vnsmp(isrc)
          else
            write(ift,'(2x,"=== plsorc ===   iflx/isrc =",
     >      i2,"/",i2,2x,a,2x,a,"  flux =",1pe14.6," =>",1pe14.6,
     >      "  rcytm =",0pf6.3,"  nsmp =",i7)')  iflx,isrc,
     >      mcflx(iflx),"---",0.0d0,tfsrc,rcytm,0
          endif
        endif
!
!::total neutral flux
        snflx(iflx) = tfsrc
        if( tfsrc.le.0.0d0 ) cycle

!::monte result
        if( isrc.gt.0 ) then
          wisrc = isrc
          wcsrc = vcsrc(isrc)
          call setvmwork( 1, isrc )
          call ntwcon
        endif

!::vli  recomb.  D+ + e => D0
        if( iflx == 8 ) then
          call plrecm(tfsrc)
          call ntdens(0.0d0)  ! <== No neutral ( set den0 = 0.0 )
          if( kgrp.eq.1 ) then
            call plntsr(tfsrc)
          else
            call ntsumsr(tfsrc)
          endif
!
!::sol/idp/../vol
        else
          call ntdens(tfsrc)
          if( kgrp.eq.1 ) then
            call plntsr(tfsrc)
          else
            call ntsumsr(tfsrc)
          endif
        endif
!
        zsum1 = 0.0d0
        zsum2 = 0.d0
        do ia = 1, nion
          do ir = 1, 6
            zsum1 = zsum1 + sumsn(ir,ia)
          enddo
          zsum2 = zsum2 + sumsn(7,ia)
        enddo
!
!::save flux
        if( iflx.eq.6 ) tnpuf = tfsrc
        pfl_ntl(iflx) = tfntl
        pfl_src(iflx) = zsum1
        pfl_man(iflx) = zsum2
        pfl_abs(iflx) = flexh
        pfl_pmp(iflx) = flpmp
        pfl_err(iflx) = tfsrc - (zsum1+zsum2+flexh)
        if( iflx.eq.8 ) then
          pfl_err(iflx) =-tfsrc - (zsum1+zsum2+flexh)
        endif
!
        if( dabs(pfl_err(iflx)/tfsrc).le.1.0d-14 ) then
          pfl_err(iflx) = 0.0d0
        endif
!
      enddo ! loop(iflx)
!
!----------------------------------------------------------------------
!::radiation  ! No need for neut2d & impmc
!----------------------------------------------------------------------
 100  continue
      if( limp.eq.0 ) then
        if( cimp.gt.0.0d0 ) call plcrad
      else
        call plcrad_mc
      endif
!
!----------------------------------------------------------------------
!::addiational heating : NBI and Joule heating
!----------------------------------------------------------------------
      if( lgtpc == 1 ) call plheat
!
      if( ift.gt.0 ) then
        write(ift,602) "srcrad","Swe",dotn,trad,(tradrg(ir),ir=1,7)
        write(ift,602) "total","SWe",sdotn2
     >  ,crad_totwe(10),(crad_totwe(ir),ir=1,7)
 602    format(2x,a6,2x,a3,1pe11.3,1pe14.6,1x,1p12e11.3)
      endif
!
!----------------------------------------------------------------------
!::partcile conservation
!----------------------------------------------------------------------
      do i = 1, mfmax
        psm_ion = psm_ion + pfl_ion(i)
        psm_ntl = psm_ntl + pfl_ntl(i)
        psm_src = psm_src + pfl_src(i)
        psm_man = psm_man + pfl_man(i)
        psm_abs = psm_abs + pfl_abs(i)
        psm_pmp = psm_pmp + pfl_pmp(i)
        psm_err = psm_err + pfl_err(i)
      enddo
!
      psm_Fl = pfl_ion(1)+pfl_ion(2)+pfl_ion(3)+pfl_ion(4)+pfl_ion(5)
      psm_Si = psm_src
      psm_Dt = -pfl_ion(5)+pfl_ntl(6)-(psm_man+psm_abs)
!
      if( ift.eq.0 ) return
!
      write(ift,'(2x,"=== plsorc ===  partcile conservation ")')
      write(ift,'(4x,"i",1x,"flx",6x,"ion",12x,"ntl",12x,"src",12x,
     >  "hot",12x,"abs",12x,"pmp",12x,"err")')
      do i = 1, mfmax
      write(ift,'(2x,i3,1x,a,1p7e15.6)')
     >  i,mcflx(i),pfl_ion(i),pfl_ntl(i),pfl_src(i),pfl_man(i),
     >  pfl_abs(i),pfl_pmp(i),pfl_err(i)
      enddo
      write(ift,'(2x,i3,1x,a,1p7e15.6)') 0,"total ",
     > psm_ion,psm_ntl,psm_src,psm_man,psm_abs,psm_pmp,psm_err
!
      write(ift,'(2x,1x,a,1pe15.6,2x,a,1pe15.6,2x,a,1p2e15.6)')
     >  "Fl =",psm_Fl, "Si =",psm_Si,"Dt =",-psm_Fl+psm_Si,psm_Dt
!
      write(ift,'(2x,"*** plsorc ***  END   time =",1pe14.6,
     >  "  itim =",i6,"  nlp =",i3)') time,itim,nlp
!
      end