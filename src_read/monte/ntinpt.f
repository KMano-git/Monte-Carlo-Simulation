!***********************************************************************
      subroutine ntinpt(nft)
!***********************************************************************
      use chrt,   only : cdhrt
      use cmolhy, only : lmolhcr
      use cntcom, only : cerel, cftnr, cftnw, cvrel, e0in, e0pfa, e0pfm
     >    , eion, eioni, famp1, famp2, faref, fawid, flxin, ipfnp
     >    , iransd, ldpmd, lelmd, lemd2, lemdp, lempf, lemwl, lnnel
     >    , lntmd, lrcmd, lscrg, lsmax, lsrmd, ltrc, lwlmd, ndrmx, nftnr
     >    , nftnw, ngas, noia, npf, nptldp, nptlvl, nptlwl, pfmag, pfpx1
     >    , pfpx2, pfpy1, pfpy2, rcfc, rmas, scfc, temin_ion, temin_rec
     >    , tlmt_el, v0in, void_ne, void_ni, void_te, void_ti, wmate
     >    , wtmin
     >    , use_gaspufftime, pufftime, pfmag_time, pufftime_max
     >    , pufftime_size
      use cntctl, only : dtntl, hrt_dbg, imxnt, mdl_hrt, mdl_ntsm, nhfl
     >    , nhunt, nini_cal, nntl
      use cntpfctl, only : fac_pfctl, lpfctl, lpfctl_rst, nsep_pfctl
     >    , rst_pfctl
      use cphcns, only : cev, cmp, cpi
      use cplcom, only : ama
      use csonic, only : lntl
      use cunit,  only : n6
      use cfatmhy, only : latomhcr
      use dbg_mod, only : dtmdbgneu, timdbgneu
      implicit none
!
!::argument
      integer, intent(in) :: nft
!
!::local variables
      character cdsn*80
      integer  i, ig, ia,ierr
! function
      integer  lenx
!
      namelist /untinp/
     >   lntmd, nntl, dtntl, imxnt, nhunt, nhfl
     >  ,eioni, eion, famp1, famp2, faref, fawid
     >  ,iransd
     >  ,wmate
     >  ,npf, ipfnp, pfmag, pfpx1, pfpx2, pfpy1, pfpy2
     >  ,nptldp, nptlwl, nptlvl
     >  ,ldpmd, lwlmd, lrcmd, lelmd, lsrmd
     >  ,lemwl, lemdp, lempf, e0pfa, e0pfm
     >  ,lsmax, ltrc
     >  ,nftnr, cftnr, nftnw, cftnw
     >  ,scfc, rcfc
     >  ,temin_ion, temin_rec, tlmt_el
     >  ,nini_cal, mdl_ntsm    ! KH111108
     >  ,latomhcr, lmolhcr !KH160812, SY211112
     >  ,mdl_hrt, cdhrt, hrt_dbg
     >  ,lscrg
     >  ,lemd2 !KH191211
     >  ,lnnel
     >  ,lpfctl, nsep_pfctl, fac_pfctl, rst_pfctl, lpfctl_rst
! added 1 line periodic invocation of IMPMC/figdat by kamata 2023/08/04
     >  , dtmdbgneu, timdbgneu

!::time depending puff rate setting
      namelist /gaspufftime/ use_gaspufftime, pufftime, pfmag_time
!
      inquire( unit=nft, name=cdsn )
      write(n6,'(/2x,"***  ntinpt  ***  cdsn=",a)') cdsn(1:lenx(cdsn))
!
!::default value
      lemwl = 1       !  d0 3 ev outward
      lemdp = 0       !  dummy
      lsmax = 10000   !  max event number in cell
      ltrc  = 0       !  flag of tracing particle
      iransd = 98815  !  initial random number
      temin_ion = -0.7d0  ! ionization
      temin_rec = -0.7d0  ! recombination
      e0pfa = 3.0d0       ! energy of puff d0 (flanck-condon)
      e0pfm = 0.0388d0    ! energy of puff d2
!xx   mdl_elcx = 0        ! reduction of EL & CX cross section
      tlmt_el = 0.9d0     ! limit for el_collission  see ntsctr_el.f
      nini_cal = 1        ! initical cal. < 10 KH111108
      mdl_ntsm = 1        ! smothing source    KH111108
      latomhcr = 0        ! atomic data by 0:DEGAS and 1:CR mode (AMJUEL)
      lmolhcr = 0         ! molecular data by 0:DEGAS and 1:CR mode (AMJUEL)
      mdl_hrt = 0         ! H radiation trapping
      cdhrt = ""          ! location of hrt.dat
      hrt_dbg = 0
      lscrg = 0           ! ntscrg
      lnnel = 0
      lemd2 = 0           !  0:trim, 1:all particles are emitted/reflected as D2
      lpfctl = 0          ! select ipfnp for nsep control by gas puff. 0: off
      fac_pfctl = 1.0d5
      lpfctl_rst = 0
      rst_pfctl = 0.0d0
      void_ne = 0.0d0
      void_ni = 0.0d0
      void_te = 0.0d0
      void_ti = 0.0d0
!
!-------- fuji 2010/03/04 ------ ntinpt ----------------
      rcfc(1:ndrmx) = 1.0d0
      scfc(1:ndrmx) = 1.0d0
!-------- fuji 2010/03/04 ------ ntinpt ----------------
!
!::input /untinp/
      read(nft,untinp)
!
!::check
      if( lntmd.eq.0 ) lrcmd = 0
      if( lntmd.eq.0 ) npf = 0
      if( lntl.ne.1 ) cftnr = " " // cftnr(1:lenx(cftnr))
      if( cftnr(1:1).eq." " ) nftnr = 0
      if( cftnw(1:1).eq." " ) nftnw = 0
      if( ldpmd.eq.0 ) nptldp = 0
      if( lwlmd.eq.0 ) nptlwl = 0
      if( lrcmd.eq.0 ) nptlvl = 0
      if( ltrc.gt.0 )  ltrc = 1
!      if( ltrc.gt.0 )  lsmax = 200

!::input /gaspufftime/
      rewind (nft)
      read(nft,nml=gaspufftime,iostat=ierr)
      rewind (nft)
      if(use_gaspufftime) then
        call check_t_order(pufftime,pufftime_max,pufftime_size)
        write(n6,'(5x,"time depending puff-rate is used")')
      endif
!
!:debug write
      write(n6,'(5x,"lntmd =",i3,2x,"neutral model ",
     >   "(0:analytical model/1:monte carlo)")') lntmd
      write(n6,'(5x,"nntl  =",i3,"  imxnt =",i3,"  dtntl =",1pe12.3)')
     >     nntl, imxnt, dtntl
      write(n6,'(5x,"nhunt =",i3,"  nhfl  =",i3)') nhunt, nhfl
      write(n6,'(5x,"ldpmd =",i3,"  lwlmd =",i3,"  lrcmd =",i3)')
     >     ldpmd,lwlmd,lrcmd
      write(n6,'(5x,"lelmd =",i3,"  lnnel =",i3,"  lsrmd =",i3)')
     >     lelmd,lnnel,lsrmd
      write(n6,'(5x,"lsmax =",i8,"  ltrc =",i2)') lsmax, ltrc
      write(n6,'(5x,"lemwl =",i3,"  lemdp =",i3,"  lempf =",i3)')
     >  lemwl, lemdp, lempf
      write(n6,'(5x,"nftnr =",i3,  "  cftnr =",a)') nftnr,cftnr(1:50)
      write(n6,'(5x,"nftnw =",i3,  "  cftnw =",a)') nftnw,cftnw(1:50)
      write(n6,'(5x,"eion  =",f8.3,"  eioni =",f8.3/5x,"famp1 =",f8.3,
     >  "  famp2 =",f8.3,"  faref =",f8.3,"  fawid =",f8.3)')
     >   eion, eioni, famp1, famp2, faref, fawid
      write(n6,'(5x,"wall material =",a)') wmate
      write(n6,'(5x,"nptldp=",i6,"  nptlwl=",i6,"  nptlvl=",i6)')
     >    nptldp,nptlwl,nptlvl
      write(n6,'(5x,"scfc  =",15f5.2)') (scfc(i),i=1,ndrmx)
      write(n6,'(5x,"rcfc  =",15f5.2)') (rcfc(i),i=1,ndrmx)
      write(n6,'(5x,"temin_ion =",f5.2,"  temin_rec =",f5.2,
     >   "  tlmt_el =",f5.2)') temin_ion, temin_rec, tlmt_el
      write(n6,'(5x,"e0pfa =",f8.4,"  e0pfm =",f8.4)') e0pfa, e0pfm
      write(n6,'(5x,"tlmt_el =",f5.2)')  tlmt_el
      write(n6,'(5x,"nini_cal =",i3,"  mdl_ntsm =",i3)')
     >    nini_cal, mdl_ntsm
      write(n6,'(5x,"mdl_hrt =",i2,"  hrt.dat = ",a)')mdl_hrt, cdhrt
      write(n6,'(5x,"ne/ni/te/ti in void =",4e12.3)')
     >    void_ne, void_ni, void_te,void_ti
      write(n6,'(5x,"lscrg =",i2 )') lscrg
      write(n6,'(5x,"lemd2 =",i2 )') lemd2
      write(n6,'(5x,"lpfctl =",i2,",  nsep_pfctl =",1pe14.7,
     >    ",  rst_pfctl = ",1pe14.7 )')
     >    lpfctl, nsep_pfctl, rst_pfctl
!
      if( famp1.lt.famp2 .or. famp2.le.1.0d0 .or. fawid.le.0.0 ) then
        call wexit("ntinpt","famp1<famp2 famp2<=1.0 fawid.<=0.0")
      endif
!
!::neutral particle generation model
      call ntmset
!
!-----    from ntplas
!
!::number of gasses
      ngas = 1
      ig   = 1
      noia(ig) = 1
!
!::atomic mass of neutral particles
      do 110 ig = 1, ngas
      ia = noia(ig)
      rmas(ig) = ama(ia)
 110  continue
!
!::const in erel of sig_cx
      do ig = 1, ngas
      cvrel(ig) = 8.0d0*cev/(cpi*rmas(ig))
      enddo
      cerel = 0.5d0*cmp/cev
!
!::neutral energy
      e0in = 3.0d0
      v0in = dsqrt(2.0d0*e0in*cev/rmas(1))
!
!::neutral flux
      flxin = 1.0d23
!
!::minimum weight
      wtmin = 1.0d-4
!-----
!
      return
      end