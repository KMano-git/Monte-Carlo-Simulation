!**********************************************************************
      subroutine plinpt(nft)
!**********************************************************************
!
!::dummy variables
!xx      dtrmp, dtrsl, dtrwl, dtrod, dtrid
!xx      dtna, dtva, dtte, dtti
!
!        test  svn  2012/07/03     K. Shimizu     A
!----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : aion, aza, beta, bmpfn, bmpni, bmpte, bmpti
     >    , bwpfn, bwpni, bwpte, bwpti, bwsfn, bwsni, bwste, bwsti, caim
     >    , cdle, cdli, cevpr, cfth, cftp, cftr, cftw, chfcv, chfpr
     >    , cimp, cimprg, cimprg2, clmdnz, clzef, cprm, dtfn, dtlmt
     >    , dtmax, dtmin, dtrateq, dttb, edtb, edtmn, edtmx, elpmx, eltb
     >    , factor_bal, fai, faim, faip, fcda
     >    , fcet, fcna, fcbcvl, fdps, flime, flimi, flimv, gcse
     >    , gcse_pv, gcse_sl, gcsi, gcsi_pv, gcsi_sl, gwpni, gwpte
     >    , gwpti, gwsni, gwste, gwsti, ibcyl, itcn, itfix, ittb
     >    , jen_bal, jst_bal, lbcgd, lbcpw, lbcsw, ldluc, ldps, lfbbc
     >    , lmuscl, lordr, lttb, lttbr, lupwd, mdl_cgen, mdl_bal
     >    , mdl_bale, mdl_eqp, mdl_edt, mdl_fimp, mdl_hcv, mdl_ini
     >    , mdl_srp, mdl_srw, mdl_vis, mdl_wrd, nad0, nas0, nfth, nftp
     >    , nftr, nftw, nimin_aux, nion, nlpmn, nlpmx, nsmp_wrad, ntfn
     >    , nttb, rcydt, rcysp, relmt, rxflp, ted0, temin_aux, temn_dprv
     >    , temn_dpsl, temndf, temndp, tes0, tid0, timin_aux, timin_vis
     >    , timndf, timneta, tis0, tmtb, vlda, vlet, vlxe, vlxi, wfac_lv
     >    , xplm_ic, xplm_jc, xplm_ne, xplm_opt, xplm_te, fdeg2
      use cplhst, only : ndht
      use cplmet, only : icel, icmpe, icmps, icspx, icwl1, itmax, itmps
     >    , itpve, itpvs, itsle, itsls, jcel, jtmax, jtmin, kreg
      use cplqcn, only : itqcn
      use cplvpn, only : vlvp, vpn_man, vpn_sol
      use csize,  only : ndeq, ndsp, ndx, ndy
      use csonic, only : dtim, itend, lpost, mstop, mxdsk, mxhst, mxprf
     >    , ndsk, nhsav, nhst, npcn, nprf, tend, time
      use cunit,  only : n6
!     for soldor time series
      use mod_soldorTimeSeries, only : timeNum, interNum 
      implicit none
!
!::argument
      integer, intent(in) :: nft
!
!::local variables
      integer  lenx, i, ii, it, jt, j, irg, ia, n, ic
      integer  jtst, jten, ndp
      integer  m1a, m2a, m3, m4
      real*8   tfac
      integer, parameter :: rgsize=10, dpsize=2, xplmt_len=80
      real(8), dimension(rgsize) :: dtq1rg, dtq2rg, dtq3rg, dtq4rg
      real(8), dimension(dpsize)  :: dtq1dp, dtq2dp, dtq3dp, dtq4dp
      character xplmt(xplmt_len)
      integer  ierr
      character cdsn*80,clin*120
!
      call trmark("plinpt","start")
!
      dtq1rg(1:10) = 0.0d0
      dtq2rg(1:10) = 0.0d0
      dtq3rg(1:10) = 0.0d0
      dtq4rg(1:10) = 0.0d0
      dtq1dp(1:2)  = 0.0d0
      dtq2dp(1:2)  = 0.0d0
      dtq3dp(1:2)  = 0.0d0
      dtq4dp(1:2)  = 0.0d0
      dtrateq(1:ndx,1:ndy,1:ndeq) = 0.0d0
!
      inquire( unit=nft, name=cdsn )
      write(n6,'(/2x,"*** plinpt ***  cdsn=",a)') cdsn(1:lenx(cdsn))
      call flush(n6)
!
!::input /uplinp/
      call read_uplinp(nft,dtq1rg, dtq2rg, dtq3rg, dtq4rg, rgsize
     > , dtq1dp, dtq2dp, dtq3dp, dtq4dp, dpsize
     > , xplmt, xplmt_len, timeNum, interNum, vlda)
      call flush(n6)
!
!     error check for input parameter
      if(fdeg2.le.0.0) then
        call wexit("plinpt","fdeg2 in inppls must non-negative")
      endif
!
!::test calculation
      if( itend.le.5 ) then
         nfth = 0
         npcn = 1
         dttb(1) = 1.0d-12
         dtmin   = 0.5d-12
      endif
!
!::file
      if( cftr(1:1).eq." " ) nftr = 0
      if( cftw(1:1).eq." " ) nftw = 0
      if( cfth(1:1).eq." " ) nfth = 0
      if( cftp(1:1).eq." " ) nftp = 0
!
!::output
      if( nion.eq.1 )  fcna(1) = 1.0d0
      if( mxdsk.gt.0 ) ndsk = itend/mxdsk
      if( mxhst.gt.0 ) nhst = itend/mxhst
      if( mxprf.gt.0 ) nprf = itend/mxprf
      if( ndsk.le.0 )  ndsk = 100
      if( nhst.le.0 )  nhst = 1
      if( nprf.le.0 )  nprf = 1
      if( nhsav.le.0 ) nhsav = ndht
      if( nhsav.gt.ndht ) nhsav = ndht
      if( npcn.le.0 )  npcn = 100
      if( itfix.eq.99 ) itfix = itend + 10
!
!::time
      time = 0.0d0
!
!::time step control
      if( nftr.gt.0 ) then
      do 110 i = 1, 10
      lttb(i) = lttbr(i)
 110  continue
      endif
!
      if( nttb.le.0 )  nttb = 1
      if( nttb.ge.21 ) nttb = 20
      do 120 i = nttb, 1, -1
      ii = i
      if( dttb(i).gt.0.0d0 .and. eltb(i).gt.0.0d0 .and.
     >    edtb(i).gt.0.0d0 .and. lttb(i).gt.0 ) goto 125
 120  continue
 125  continue
      nttb = ii
      if( dtlmt.le.0.0d0 ) dtlmt = dttb(nttb)*10.0d0
      nttb = nttb + 1
      if( nttb.gt.21 ) goto 910
      tmtb(nttb) = 1.0e10
      dttb(nttb) = dttb(nttb-1)
      eltb(nttb) = eltb(nttb-1)
      edtb(nttb) = edtb(nttb-1)
      lttb(nttb) = itend + 8000
!
      ittb  = 1
      itcn  = lttb(ittb)
      dtim  = dttb(ittb)
      dtmax = dttb(ittb)
      elpmx = eltb(ittb)
      edtmx = edtb(ittb)
      tfac  = 1.0d0
      write(n6,'(2x,"<dtcntl-in>  set  dtim =",1pe11.3,"  tfac ="
     >  ,0pf8.3,"  elpmx,edtmx =",1p2e11.3/)') dtim,tfac,elpmx,edtmx
      if( dtim.le.0.0 .or. elpmx.le.0.0 .or. edtmx.le.0.0 ) goto 920
!
!::fine time step just afer IMPMC cal. (see pldtfn.f)
      ntfn = -1
      do i = 0, 30
        if( dtfn(i).eq.0.0d0 ) exit
        ntfn = i
      enddo
      write(n6,'(2x,"*** pldtfn ***   ntfn =",i3)') ntfn
      if( ntfn.ge.0 ) then
      write(n6,'(2x,"dtfn =",1p10e12.3)') (dtfn(i),i=0,ntfn)
      endif
!
!::ratio of time step   New plinpt 2011/10/13
      do 130 it = 1, itmax
      jtst = jtmin(it)
      jten = jtmax(it)
!
      do 140 jt = jtst, jten
      i = icel(jt,it)
      j = jcel(jt,it)
      irg = kreg(j,i)
      if( irg.lt.1 .or. irg.gt.6 ) goto 925
!-----
!::near dp
      ndp = 0
      if( irg.eq.1 .or. irg.eq.4 ) then
        if( jt.le.jtst+4 ) ndp = 1
      elseif( irg.eq.3 .or. irg.eq.5 ) then
        if( jt.ge.jten-4 ) ndp = 2
      endif
!-----
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      dtrateq(j,i,m1a) = dtq1rg(irg)
      dtrateq(j,i,m2a) = dtq2rg(irg)
      if( ndp.gt.0 ) then
      dtrateq(j,i,m1a) = dtq1dp(ndp)
      dtrateq(j,i,m2a) = dtq2dp(ndp)
      endif
      enddo
      m3 = 2*nion + 1
      m4 = 2*nion + 2
      dtrateq(j,i,m3) = dtq3rg(irg)
      dtrateq(j,i,m4) = dtq4rg(irg)
      if( ndp.gt.0 ) then
      dtrateq(j,i,m3) = dtq3dp(ndp)
      dtrateq(j,i,m4) = dtq4dp(ndp)
      endif
!-----
 140  continue
 130  continue
!
!::input /upldif/
      call read_upldif(nft,ierr)
!
!::anomalous diffusion
      do 220 n = 1, nion
      if( fcda(n).le.0.0 ) fcda(n) = 1.0d0
      if( fcet(n).le.0.0 ) fcet(n) = 1.0d0
 220  continue
      do 225 n = 2, 6
      if( vlda(n).le.0.0 ) vlda(n) = vlda(1)
      if( vlet(n).le.0.0 ) vlet(n) = vlet(1)
      if( vlxi(n).le.0.0 ) vlxi(n) = vlxi(1)
      if( vlxe(n).le.0.0 ) vlxe(n) = vlxe(1)
 225  continue
!
!::vpinch
      do ic = 1, ndy
      vlvp(ic) = 0.0d0
      enddo
      i = 0
      do ic = icspx, icwl1, -1
      if( ic.eq.icwl1 ) cycle
      i = i + 1
      vlvp(ic) = vpn_sol(i)
      enddo
      i = 0
      do ic = icmps, icmpe, +1
      if( ic.eq.icmpe ) cycle
      i = i + 1
      vlvp(ic) = vpn_man(i)
      enddo
!
!::boundary condition
      do 230 n = 1, nion
      if( bwsfn(n).le.0.0 ) bwsfn(n) = fcna(n)
      if( bwpfn(n).le.0.0 ) bwpfn(n) = fcna(n)
      if( bmpfn(n).le.0.0 ) bmpfn(n) = fcna(n)
 230  continue
!
      do 240 it = 1, itmax
      ibcyl(it) = 0
      if( it.ge.itmps ) ibcyl(it) = 1
 240  continue
!
!::third order scheme fai=1/3 , second order fai=-1
      if( lordr.ne.3 ) lordr = 2
      fai  = 1.0d0/3.0d0
      if( lordr.eq.2 ) fai = -1.0d0
      faip = (1.0d0+fai)/4.0d0
      faim = (1.0d0-fai)/4.0d0
      beta = (3.0d0-fai)/(1.0d0-fai)
!
!::sound velocity
      do it = 1, itmax
      gcsi(it) = 1.0d0
      gcse(it) = 1.0d0
      enddo
      do it = itsls, itsle
      if( gcsi_sl.gt.0.0d0 ) gcsi(it) = gcsi_sl
      if( gcse_sl.gt.0.0d0 ) gcse(it) = gcse_sl
      enddo
      do it = itpvs, itpve
      if( gcsi_pv.gt.0.0d0 ) gcsi(it) = gcsi_pv
      if( gcse_pv.gt.0.0d0 ) gcse(it) = gcse_pv
      enddo
!
!::relaxzation factor
!
!::impurity concentration
      call set_wrdcr
!
!::Temin
      do i = 1, ndy
      temndp(i) = 0.0d0
      if( i.ge.itsls .and. i.le.itsle ) temndp(i) = temn_dpsl
      if( i.ge.itpvs .and. i.le.itpve ) temndp(i) = temn_dprv
      enddo
!
!::convection flux  c52 => c32
      chfpr = 1.0d0
      if( mdl_hcv.eq.1 ) chfpr = 0.0d0
      chfcv = 1.5d0 + chfpr
      cevpr = cev*chfpr
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      write(n6,'(2x,"comment")')
      write(n6,'(5x,"caim : ",a)') caim(1:60)
      write(n6,'(5x,"cprm : ",a)') cprm(1:60)
      write(n6,'(/2x,"restart file")')
      write(n6,'(5x,"read  file :",i3,"  cftr =",a)') nftr,cftr(1:50)
      write(n6,'(5x,"write file :",i3,"  cftw =",a)') nftw,cftw(1:50)
      write(n6,'(5x,"hist  file :",i3,"  cfth =",a)') nfth,cfth(1:50)
      write(n6,'(5x,"prof  file :",i3,"  cftp =",a)') nftp,cftp(1:50)
!
      write(n6,'(2x,"output data")')
      write(n6,'(5x,"restart     ndsk =",i4,"  mxdsk =",i4)') ndsk,mxdsk
      write(n6,'(5x,"hist data   nhst =",i4,"  mxhst =",i4,"  nhsav =",
     >  i3)') nhst,mxhst,nhsav
      write(n6,'(5x,"prof data   nprf =",i6,"  mxprf =",i4)') nprf,mxprf
!
      write(n6,'(2x,"ion species  nion =",i3,"  ndsp =",i3)') nion,ndsp
!
      do 310 ia = 1, nion
      write(n6,'(5x,"aion =",f5.2,"  achg =",f5.2,"  fcna =",f5.2)')
     >   aion(ia), aza(ia), fcna(ia)
 310  continue
!
      write(n6,'(2x,"variables related to time")')
      write(n6,'(5x,"time =",1pe11.3,"  tend =",1pe11.3,
     >   "  itend =",i6,"  itfix =",i6)') time, tend, itend, itfix
      write(n6,'(5x,"dtmax =",1pe11.3,"  dtmin =",1pe11.3,"  dtlmt =",
     >  1pe11.3)') dtmax,dtmin,dtlmt
      write(n6,'(2x,"step size & convergence  nttb  =",i4)') nttb
      write(n6,'(5x,"lttb =",10(1x,i6,4x))')  (lttb(n),n=1,nttb)
      write(n6,'(5x,"dttb =",1p11e11.3)')     (dttb(n),n=1,nttb)
      write(n6,'(5x,"eltb =",1p11e11.3)')     (eltb(n),n=1,nttb)
      write(n6,'(5x,"edtb =",1p11e11.3)')     (edtb(n),n=1,nttb)
!
      write(n6,'(2x,"ratio of time step size")')
      write(n6,'(5x,"dtq1rg =",10f8.2)') (dtq1rg(i),i=1,8)
      write(n6,'(5x,"dtq2rg =",10f8.2)') (dtq2rg(i),i=1,8)
      write(n6,'(5x,"dtq3rg =",10f8.2)') (dtq3rg(i),i=1,8)
      write(n6,'(5x,"dtq4rg =",10f8.2)') (dtq4rg(i),i=1,8)
      write(n6,'(5x,"dtq1dp =",10f8.2)') (dtq1dp(i),i=1,2)
      write(n6,'(5x,"dtq2dp =",10f8.2)') (dtq2dp(i),i=1,2)
      write(n6,'(5x,"dtq3dp =",10f8.2)') (dtq3dp(i),i=1,2)
      write(n6,'(5x,"dtq4dp =",10f8.2)') (dtq4dp(i),i=1,2)
      call dbg_dtrateq
      if( mstop > 0 ) return
!
      write(n6,'(5x,"nlpmn =",i3,"  nlpmx =",i3)') nlpmn, nlpmx
      write(n6,'(5x,"rxflp =",20f6.2)') (rxflp(n),n=1,nlpmx)
!
      write(n6,'(2x,"numerical scheme: lmuscl/lordr =",i2,"/",i2,
     > "  fai =",f8.3,"  (1/3:third, -1:second)")') lmuscl,lordr,fai
      write(n6,'(2x,"up-wind scheme  : lupwd =",i2)') lupwd
!
      write(n6,'(2x,"classical diffusion")')
      write(n6,'(5x,"clzef =",f6.3)') clzef
      write(n6,'(5x,"ldluc =",i3)') ldluc
!
      write(n6,'(2x,"anomalous diffusion")')
      write(n6,'(5x,"ldps =",i2,"  0 (const) / 1 (D_psi:flux scaled)"
     >    )') ldps
      write(n6,'(5x,"fdps =",1pe12.3)') fdps
      write(n6,'(5x,"fcda(ion) =",1p6e11.3)') (fcda(n),n=1,nion)
      write(n6,'(5x,"fcet(ion) =",1p6e11.3)') (fcet(n),n=1,nion)
      write(n6,'(5x,"vlda(reg) =",1p6e11.3)') (vlda(n),n=1,6)
      write(n6,'(5x,"vlet(reg) =",1p6e11.3)') (vlet(n),n=1,6)
      write(n6,'(5x,"vlxi(reg) =",1p6e11.3)') (vlxi(n),n=1,6)
      write(n6,'(5x,"vlxe(reg) =",1p6e11.3)') (vlxe(n),n=1,6)
      write(n6,'(5x,"sol icwl1,icspx =",2i5,"  man icmps,icmpe =",2i5)')
     >   icwl1, icspx, icmps, icmpe
      write(n6,'(5x,"vlvp(sol) =",10f9.3)')   (vlvp(i),i=icwl1,icspx)
      write(n6,'(5x,"vlvp(man) =",10f9.3)')   (vlvp(i),i=icmps,icmpe)
      write(n6,'(5x,"temndf =",f9.3,"  timndf =",f9.3)') temndf, timndf
      write(n6,'(5x,"timneta =",f9.3)') timneta
      write(n6,'(5x,"temn_dpsl =",f9.3,"  temn_dprv =",f9.3)')
     >   temn_dpsl, temn_dprv
      write(n6,'(5x,"temndp =",10f9.3)') (temndp(i),i=1,itpve)
!
      write(n6,'(2x,"boundary condition")')
      write(n6,'(5x,"ibcyl=",20i3:/(5x,6x,20i3))')
     >   (ibcyl(it),it=1,itmax)
      write(n6,'(5x,"cdli =",1pe11.3,"  cdle =",1pe11.3)') cdli,cdle
      write(n6,'(5x,"gcsi_sl =",2f6.3,"  gcse_sl =",2f6.3)')
     >     gcsi(itsls),gcsi(itsle),gcse(itsls),gcse(itsle)
      write(n6,'(5x,"gcsi_pv =",2f6.3,"  gcse_pv =",2f6.3)')
     >     gcsi(itpvs),gcsi(itpve),gcse(itpvs),gcse(itpve)
      write(n6,'(5x,"lbcgd =",i3)') lbcgd
      if( lbcsw.eq.0 ) then
      write(n6,'(5x,"S-Wall lbcsw =",i2,"  (gni,gti,gte) =",1p3e11.3)')
     >  lbcsw,gwsni,gwsti,gwste
      else
      write(n6,'(5x,"S-Wall lbcsw =",i2,"  (bni,bti,bte) =",1p3e11.3)')
     >  lbcsw,bwsni,bwsti,bwste
      endif
      if( bwsni.gt.0.0 )
     > write(n6,'(5x,"        bfn(ia) =",15f6.3)') (bwsfn(n),n=1,nion)
      if( lbcpw.eq.0 ) then
      write(n6,'(5x,"P-Wall lbcpw =",i2,"  (gni,gti,gte) =",1p3e11.3)')
     >  lbcpw,gwpni,gwpti,gwpte
      else
      write(n6,'(5x,"P-Wall lbcpw =",i2,"  (bni,bti,bte) =",1p3e11.3)')
     >  lbcpw,bwpni,bwpti,bwpte
      endif
!KH 20121111
      if( lfbbc.ne.0)then
      write(n6,'(5x,"P-Wall feedback boudary condition is used")')
      endif
      if( bwpni.gt.0.0 )
     > write(n6,'(5x,"        bfn(ia) =",15f6.3)') (bwpfn(n),n=1,nion)
      write(n6,'(5x,"Mplas   (bni,bti,bte) =",1p3e11.3)')
     >  bmpni,bmpti,bmpte
      write(n6,'(5x,"         bfn(ia) =",15f6.3)') (bmpfn(n),n=1,nion)
!
      if( bwsti.le.0.0 .or. bwste.le.0.0 ) goto 930
      if( gwsni.le.0.0 .or. gwsti.le.0.0 .or. gwste.le.0.0 ) goto 940
      if( gwpni.le.0.0 .or. gwpti.le.0.0 .or. gwpte.le.0.0 ) goto 940
!
      write(n6,'(2x,"time step control")')
      write(n6,'(5x,"nlpmn =",i3,"  nlpmx =",i3,"  elpmx =",1pe11.3)')
     >   nlpmn, nlpmx, elpmx
      write(n6,'(5x,"edtmx =",1pe11.3,"  edtmn =",1pe11.3,"  relmt =",
     >   1pe11.3)') edtmx,edtmn,relmt
      if( elpmx.le.0.0 .or. edtmx.le.0.0 ) goto 950
!
      write(n6,'(2x,"simple radiation model  cimp =",es12.3)') cimp
      write(n6,'(5x,"cimprg  =",6f8.3)') (cimprg(i),i=1,6)
      write(n6,'(5x,"cimprg2 =",6f8.3)') (cimprg2(i),i=1,6)
!
      write(n6,'(2x,"recycling coef.")')
      write(n6,'(5x,"rcydt =",1pe12.3,"  rcysp =",1pe12.3)') rcydt,rcysp
!
      write(n6,'(2x,"tube number of flux conservation")')
      write(n6,'(5x,"itqcn =",10i5)') (itqcn(i),i=1,10)
!
!::radiation model to avoid Xp-MARFE
      write(n6,'(2x,"radiation model")')
      write(n6,'(5x,"xplmt =",a)') xplmt(1:lenx(xplmt))
      read(xplmt,*) xplm_opt, xplm_jc, xplm_ic, xplm_ne, xplm_te
      write(n6,'(5x,"_opt =",i3," _jc =",i3,"  _ic =",i3,
     >  "  _ne =",1pe12.3,"  _te =",1pe12.3)')
     >  xplm_opt, xplm_jc, xplm_ic, xplm_ne, xplm_te
!
!::MODEL for slimCS
      if( lpost.eq.1 ) mdl_edt = 0
      write(n6,'(2x,"slimcs model for high Ti in outer divertor")')
      write(n6,'(5x,"mdl_wrd =",i3,"  mdl_eqp =",i3,"  mdl_edt =",
     >  i3,"  mdl_srw =",i3,"  mdl_srp =",i3,"  mdl_vis =",i3)')
     >   mdl_wrd, mdl_eqp, mdl_edt, mdl_srw, mdl_srp, mdl_vis
      write(n6,'(5x,"mdl_ini =",i3,"  mdl_cgen =",i3," mdl_fimp =",i3)')
     >   mdl_ini, mdl_cgen, mdl_fimp
      write(n6,'(5x,"temin_aux =",f8.4,"  timin_aux =",f8.4,
     >   "  timin_vis =",f8.4,"  nimin_aux =",f8.4)')
     >   temin_aux, timin_aux, timin_vis, nimin_aux/1.0d19
      write(n6,'(2x,"init sol  nas0 =",1pe12.3,"  tis0 =",1pe12.3,
     >   "  tes0 =",1pe12.3)') nas0, tis0, tes0
      write(n6,'(2x,"init dpl  nad0 =",1pe12.3,"  tid0 =",1pe12.3,
     >   "  ted0 =",1pe12.3)') nad0, tid0, ted0
      write(n6,'(5x,"convection flux: mdl_hcv =",i2,"  chfpr =",f8.3,
     >  "  chfcv =",f8.3,"  cevpr =",1pe12.4)')
     >   mdl_hcv, chfpr, chfcv, cevpr
!
!::flux limiter
      write(n6,'(/2x,"flux limiter")')
      write(n6,'(5x,"flimi =",f8.4,"  flime =",f8.4,"  flimv =",f8.4)')
     >   flimi, flime, flimv
      write(n6,'(5x,"fcbcvl=",f8.4)') fcbcvl
      call flush(n6)
!
      write(n6,'(/2x,"clmdnz = ",f8.3)') clmdnz
!
!
!::YH  Heat diffusivity poloidal distribution setting
      write(n6,'(/2x,"mdl_bal = ",i2)')mdl_bal
      write(n6,'(/2x,"  jst_bal, jen_bal, factor_bal: ",2i4,f9.3)')
     >       jst_bal, jen_bal, factor_bal
!
      write(n6,'(/2x,"mdl_bale = ",i2)')mdl_bale
      if(mdl_bale.eq.1)then
         write(n6,'(/2x,"  Using Balescu model")')
      else
         write(n6,'(/2x,"  Using Braginskii model ")')
      end if
      return
!
!::error
 910  continue
      write(clin,'("nttb =",i3," > 20")') nttb
      mstop = 1
      return

 920  continue
      write(clin,'("invarid time step control  dtim,elpmx,edtmx =",
     >  1p3e12.3)') dtim,elpmx,edtmx
      mstop = 1
      return
!
 925  continue
      write(clin,'("invarid irg =",i2," jt,it =",2i5,"  j,i =",2i5)')
     >   irg, jt, it, j, i
      mstop = 1
      return
!
 930  continue
      write(clin,'("bwsni, bwsti, bwste are not specified.")')
      write(n6,'(2x,a)') trim(clin)
      mstop = 1
      return
!
 940  continue
      write(clin,'("invarid boundary condition  gwsni,gwpni,...")')
      write(n6,'(2x,a)') trim(clin)
      mstop = 1
      return
!
 950  continue
      write(clin,'("input error of time step control")')
      write(n6,'(2x,a)') trim(clin)
      mstop = 1
      return

 960  continue
      write(n6,'(2x,"error at set_atomz in plinpt")')
      mstop = 1
      return
      end
!
!**********************************************************************
      subroutine pldfin
!**********************************************************************
      use cplcom, only : bmpfn, bmpni, bmpte, bmpti, bwpfn, bwsfn
     >    , cimprg, cimprg2, clmdnz, dtfn, dtlmt, dtrid, dtrod, edtmn
     >    , edtmx, elpmx, fcbcvl, flime, flimi, flimv, gcse_pv, gcse_sl
     >    , gcsi_pv, gcsi_sl, factor_bal, heat_flux_para_by_kappa_para
     >    , heat_flux_para_elec, itfix, jen_bal, jst_bal, lbcpw, lbcsw
     >    , lfbbc, lmuscl, lupwd, mdl_bal, mdl_bale, mdl_cgen, mdl_edt
     >    , mdl_eqp, mdl_fimp, mdl_hcv, mdl_ini, mdl_srp, mdl_srw
     >    , mdl_vis, mdl_wrd, nad0, nas0, nimin_aux, nlpmn, nlpmx
     >    , nsmp_wrad, rcydt, rcysp, relmt, rxflp, ted0, temin_aux
     >    , temn_dprv, temn_dpsl, temndf, tes0, tid0, timin_aux
     >    , timin_vis, timndf, timneta, tis0, wfac_lv, fdeg2
      use cplqcn, only : itqcn
      use cplvpn, only : vpn_man, vpn_sol   ! Fuji 2010/03/04
      use csize,  only : ndsp, ndy
!::   for soldor time Series
      use mod_soldorTimeSeries, only : timeNum, interNum 
      implicit none
!
!::local variables
      integer i
!
!::fine time step
      dtfn(0:30) = 0.0d0
!
      dtlmt = 0.0d0
      nlpmn = 2
      nlpmx = 6
      elpmx = 5.0d-4
      edtmx = 4.0d-2
      edtmn = 4.0d-3
      relmt = 0.0d0
      lmuscl= 1
      lupwd = 1
      lbcsw = 0
      lbcpw = 0
      itfix = 99
!
      dtrod = 1.0d0
      dtrid = 1.0d0
!
      do i = 1, 51
      rxflp(i) = 1.0d0
      enddo
!
      rcydt = 1.0d-3
      rcysp = 0.8d0
!
      do i = 1, 10
      cimprg(i) = 0.0d0
      enddo
!
      do i = 1, 10
      itqcn(i) = 0
      enddo
!
!::minimun Te and Ti for diffusion coefficients
      temndf = 0.03d0
      timndf = 0.03d0
!
!  minimun Ti for eta timneta  by hoshino  140514
      timneta = 0.01d0
!
!::minimum Te and Ti for parallel diffusion at divertor plate
      temn_dpsl = 0.00d0    ! off
      temn_dprv = 0.05d0    ! on
!
!--------------Fuji 2010/03/04 ---------pldfin----------
!::pinch velocity
      vpn_sol(1:ndy) = 0.0d0
      vpn_man(1:ndy) = 0.0d0
!
!::boundary condition
      bwsfn(1:ndsp) = 0.0d0
      bwpfn(1:ndsp) = 0.0d0
      bmpfn(1:ndsp) = 0.0d0
!
!::Mach-number
      gcsi_sl = 0.0d0
      gcse_sl = 0.0d0
      gcsi_pv = 0.0d0
      gcse_pv = 0.0d0
!
!::imp concentration
      cimprg (1:10) = 0.0d0
      cimprg2(1:10) = 0.0d0
!
!::boundary condition
      bmpni = 0.0d0
      bmpti = 0.0d0
      bmpte = 0.0d0
      bmpfn(1:ndsp) = 0.0d0
!--------------Fuji 2010/03/04 ---------pldfin----------
!
!::MODEL
      mdl_wrd = 0
      mdl_eqp = 0
      mdl_edt = 0
      mdl_srw = 0   ! No use
      mdl_srp = 0
      mdl_vis = 0
      mdl_ini = 0
      mdl_hcv = 0
      mdl_cgen = 0
      mdl_fimp = 0
!
      temin_aux = 0.01d0
      timin_aux = 0.01d0
      timin_vis = 0.01d0
      nimin_aux = 0.05d19
!
!::initial profile
      nas0 = 0.0d0
      tis0 = 0.0d0
      tes0 = 0.0d0
      nad0 = 0.0d0
      tid0 = 0.0d0
      ted0 = 0.0d0
!
!::flux limiter
      flimi = -0.5d0   ! 2011/05/13  KH111108
      flime = 0.2d0
      flimv = 0.5d0
!
!::Heat flux density vector at cell sides
      heat_flux_para_by_kappa_para = 0.d0
      heat_flux_para_elec = 0.d0
!
!::boundary condition
      fcbcvl = 0.02d0   ! 1/50
!
!::impurity denisty
      clmdnz = 100.0d0
!
!:: feedback boundary condition at p-wall KH20121111
      lfbbc = 0
!
!::  number of test particle in pst_wrad KH150519
      nsmp_wrad = 250000
!
!::  delay of wfac valiation KH150803
      wfac_lv=0.0d0
!
!::  ballooning like diffusion model
      mdl_bal = 0
      jst_bal = -99
      jen_bal = -99
      factor_bal = 0.0d0
!
!    Balescu model for classical transport (From SOLPS5.X)
      mdl_bale = 0

!----- 22/11/2 yamamoto time series of Qdpl_i,o.txt (defoult = 1, no time average)
      timeNum = 1
      interNum = 1

      fdeg2 = 1.0

      return
      end
!
!**********************************************************************
      subroutine dbg_dtrateq
!**********************************************************************
      use cplcom, only : dtrateq, nion
      use cplmet, only : icel, icmpe, icwl1, itmax, itmps, itpvs, itsle
     >    , jcel, jcxp1, jtmax, jtmin
      use csonic, only : mstop
      use cunit,  only : n6
      implicit none
!
      integer :: ia, m1a, m2a, m3, m4
      integer :: it, jt, m, mmx
      integer :: j, i, ist, ien, jtst, jten
      integer :: ner, ier
!
      ia = 1
      m1a = 2*ia - 1
      m2a = 2*ia
      m3  = 2*nion + 1
      m4  = 2*nion + 2
!
      it = itsle
      write(n6,'(5x,"dtrateq  itsle = ",i3)') it
      write(n6,'(5x,"   m1a =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m1a),jt=1,jtmax(it))
      write(n6,'(5x,"   m2a =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m2a),jt=1,jtmax(it))
      write(n6,'(5x,"   m3  =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m3),jt=1,jtmax(it))
      write(n6,'(5x,"   m4  =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m4),jt=1,jtmax(it))
!
      it = itpvs
      write(n6,'(5x,"dtrateq  itpvs = ",i3)') it
      write(n6,'(5x,"   m1a =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m1a),jt=1,jtmax(it))
      write(n6,'(5x,"   m2a =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m2a),jt=1,jtmax(it))
      write(n6,'(5x,"   m3  =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m3),jt=1,jtmax(it))
      write(n6,'(5x,"   m4  =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m4),jt=1,jtmax(it))
!
      it = itmps
      write(n6,'(5x,"dtrateq  itmps = ",i3)') it
      write(n6,'(5x,"   m1a =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m1a),jt=1,jtmax(it))
      write(n6,'(5x,"   m2a =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m2a),jt=1,jtmax(it))
      write(n6,'(5x,"   m3  =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m3),jt=1,jtmax(it))
      write(n6,'(5x,"   m4  =",10f8.2)')
     >     (dtrateq(jcel(jt,it),icel(jt,it),m4),jt=1,jtmax(it))
!
      j = jcxp1 + 10
      ist = icwl1
      ien = icmpe
      write(n6,'(5x,"dtrateq(j,i)  j = ",i3,"  icwl1,icmpe =",2i4)')
     >   j, icwl1, icmpe
      write(n6,'(5x,"   m1a =",10f8.2)') (dtrateq(j,i,m1a),i=ist,ien)
      write(n6,'(5x,"   m2a =",10f8.2)') (dtrateq(j,i,m2a),i=ist,ien)
      write(n6,'(5x,"   m3  =",10f8.2)') (dtrateq(j,i,m3),i=ist,ien)
      write(n6,'(5x,"   m4  =",10f8.2)') (dtrateq(j,i,m4),i=ist,ien)
!
!::undefine dtrateq(j,i,m)
      write(n6,'(5x,"check dtrateq(j,i,m) = 0.0d0")')
      ner = 0
      do it = 1, itmax
      jtst = jtmin(it)
      jten = jtmax(it)
!
      do jt = jtst, jten
      i = icel(jt,it)
      j = jcel(jt,it)
      ier = 0
      mmx = 2*nion + 2
      do m = 1, mmx
      if( dtrateq(j,i,m).le.0.0d0 ) ier = ier + 1
      enddo
!
      if( ier.gt.0 ) then
      write(n6,'(2x,"find cell  jt,it =",2i5,"  j,i =",2i5,
     >    "  dtrateq(j,i,m) =",10f6.2)') jt, it, j, i,
     >    (dtrateq(j,i,m),m=1,max0(mmx,10))
      ner = ner + 1
      endif
      enddo
      enddo
      write(n6,'(5x,"number of cells  ner =",i3)') ner
!
      call flush(n6)
      if( ner.gt.0 ) then
      mstop = 1
      return
      endif
!
      return
      end
