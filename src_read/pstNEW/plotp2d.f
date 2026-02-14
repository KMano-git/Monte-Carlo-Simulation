!***********************************************************************
      subroutine plotp2d(plot_mode)
!***********************************************************************
!
!     monte     Ni   : deni(0:ndmc,ndgs)
!
!     soldor    vni  :  vni(ndx,ndy)
!     DTPRV     vne  :  vne(ndx,ndy)
!               vti  :  vti(ndx,ndy)
!               vte  :  vte(ndx,ndy)
!               vva  :  vva(ndx,ndy,ndsp)
!
!     DTWRD     vwrt :  total radiation
!               vwrm :  radiation calculated with MC
!               vwrc :  non coronal model
!
!     neut2d    N0   :  tden0(0:ndmc,ndgs)   density of D0
!     DTNTL     E0   :  teng0(0:ndmc,ndgs)   energy of  D0
!               P0   :  tden0*teng0*qe       pressure of D0
!               Ng   :  tdeng(0:ndmc,2)      denisty of D2
!               Eg   :  tengg(0:ndmc,2)      energy of  D2
!               Pg   :  tdeng*tengg*qe       pressure of D2
!
!     impmc     Wrd  :  twrd(ndmc)   radiation calculated by IMPMC
!     DTIMP     Nz   :  tdnz(0:ndis2,ndmc)   
!                   select charage state     0, a, i, 1-5   
!               rt0  :  Nz(is=0->ismax)/Ni
!               rt1  :  Nz(is=1->ismax)/Ni 
!               Zef  :  Zeff =  [Ni+sum Z**2*Nz]/[Ni+sum Z*Nz]
!               Zav  :  <Z>  =  sum[Z*Nz]/sum[Nz]
!     IMPSY     Vz   : Averaged velocity parallel to the field line 
!               Ffr  : Friction Force [N]
!               Fth  : Thermal force [N]
!               ion  : ionization source [m-3 s-1]
!               rec  : recombination source [m-3 s-1]
!      bug  setv2d_mc(,,wime,,)  dim(ndmc) wime
!           sub  setv2d_mc(,,vmc,,) dim(0:ndmc) vmc
!           ==> setv2d_mc(,,wime(1),,)  setv2d_mc(,,dene(1),,)
!                        2012/12/13  Kawashima found
!
!----------------------------------------------------------------------
      use csize
      use cimcom
      use cimden
      use cntcom
      use cntpls
      use cntsrc
      use cplcom
      use cplimp
      use cplwrd
      use cplwrd2
      use cunit
      use com_impvar, only : zav, zef
      use cplimp_plimp, only : twrdL
      use catcom, only : catmz
      implicit none
!
      integer plot_mode

!::local variables
      integer :: mji, iz, ic
      character(80) :: cinp, chrg
      character(20) :: cnam
      character(30) :: cvnm
      integer :: ia, ig, nv, kcn, is0, klog, nty
      real(8), dimension(ndmc) :: vtmp
      real(8) :: zsm, sum
      real(8), dimension(ndx,ndy) :: pto !total pressure
      integer :: mj1, mj2
!
      real(4), dimension(ndmc) :: var, varTI
      real(8), dimension(ndmc) :: var2
!
      write(n6,'(/2x,"*** plotp2d ***")')
      ia = 1
      ig = 1
!
 100  continue
      call gdput("  ")
      call gdput("MONTE   : Ni, Ne, Ti, Te, Vf")
      call gdput("SOLDOR  : vni,vne,vti,vte,vva,vwrt,vwrm,vwrc,vwrC")
      call gdput("NEUT2D  : N0, E0, Ng, Eg, P0, Pg, Sn, Sp, Wi, We, Ha")
      call gdput("IMPMC   : Wrd, rt0, rt1, rta, rtn, Zef, Zav")
      call gdput("IMPMC   : Nzi, Nz0, Nza, Nz5-9, Nz10-15")
      call gdput(" log(Wrd), Wrd, rt0=Nz(0-18)/Ni, Vz, Ffr, Fth")
 120  continue
!
      klog = 0
      call gdget(">> enter (Ti/Wrd/?/q) ==> ",cinp,mji)
      if( cinp(1:1).eq."?" ) goto 100
      if( mji.eq.0 .or. cinp(1:1).eq."q" ) goto 150
      cnam = trim(cinp)
      if(cinp(1:3).eq."log" ) then
      mj1 = index(trim(cinp),"(")
      mj2 = index(trim(cinp),")")
      if( mj1.gt.0 .and. mj2.gt.0 ) then
        cnam = cinp(mj1+1:mj2-1)
        klog = 1
      else
        call gdput("wrong input "//trim(cinp))
        goto 120
      endif
      endif
      cvnm = trim(cinp)
!
!----------------------------------------------------------------------
!::MONTE
!----------------------------------------------------------------------
!::Ni
      if( trim(cnam).eq."Ni" ) then
      nv = ncmax
      call setv2d_mc("Ni",1.0d-19,nv,deni(1,ia),var,kcn)
!::Ne
      elseif( trim(cnam).eq."Ne" ) then
      nv = ncmax
      call setv2d_mc("Ne",1.0d-19,nv,dene(1),var,kcn)
!::Ti
      elseif( trim(cnam).eq."Ti" ) then
      nv = ncmax
      call setv2d_mc("Ti",1.0d0,nv,temi(1),var,kcn)
!::Te
      elseif( trim(cnam).eq."Te" ) then
      nv = ncmax
      call setv2d_mc("Te",1.0d0,nv,teme(1),var,kcn)
!::Vf
      elseif( trim(cnam).eq."Vf" ) then
      nv = ncmax
      call setv2d_mc("Vf",1.0d0,nv,vflw(1,ig),var,kcn)
!
!----------------------------------------------------------------------
!::SOLDOR
!----------------------------------------------------------------------
!::vni
      elseif( trim(cnam).eq."vni" ) then
      nv = ncmax
      call setv2d_pl("vni",1.0d-19,nv,vni,var,kcn)
!::vne
      elseif( trim(cnam).eq."vne" ) then
      nv = ncmax
      call setv2d_pl("vne",1.0d-19,nv,vne,var,kcn)
!::vti
      elseif( trim(cnam).eq."vti" ) then
      nv = ncmax
      call setv2d_pl("vti",1.0d0,nv,vti,var,kcn)
!::vte
      elseif( trim(cnam).eq."vte" ) then
      nv = ncmax
      call setv2d_pl("vte",1.0d0,nv,vte,var,kcn)
!::vva
      elseif( trim(cnam).eq."vva" ) then
      nv = ncmax
      call setv2d_pl("vva",1.0d0,nv,vva,var,kcn)
!::ptot
      elseif( trim(cnam).eq."pto" ) then
            nv = ncmax
            pto(:,:) = 1.6d-19*(vni(:,:)*vti(:,:)+vne(:,:)*vte(:,:))!+vni(:,:)*3.35e-27*vva(:,:,1)**2
            call setv2d_pl("pto",1.0d0,nv,pto,var,kcn)
!::vwrt < 0  t: total
      elseif( trim(cnam).eq."vwrt" ) then
      nv = ncmax
      call setv2d_mc("vwrt",-1.0d-6,nv,wime(1),var,kcn)
!::vwrc < 0  c: non corona
      elseif( trim(cnam).eq."vwrc" ) then
      nv = ncmax
      call setv2d_mc("vwrc",-1.0d-6,nv,wcre(1),var,kcn)
!::vwrm < 0  m: monte
      elseif( trim(cnam).eq."vwrm" ) then
      nv = ncmax
      call setv2d_mc("vwrm total",-1.0d-6,nv,wmce(1),var,kcn)
!::vwrm < 0  m: monte 1
      elseif( trim(cnam).eq."vwrm1" ) then
      nv = ncmax
      write(cvnm,'("vwrm Z=",i2)') ismaxL(1)
      call setv2d_mc(cvnm,-1.0d-6,nv,wmc_zwe(1,1),var,kcn)
!::vwrm < 0  m: monte 2
      elseif( trim(cnam).eq."vwrm2" ) then
      nv = ncmax
      write(cvnm,'("vwrm Z=",i2)') ismaxL(2)
      call setv2d_mc(cvnm,-1.0d-6,nv,wmc_zwe(1,2),var,kcn)
!
!::vwrfn < 0  ! fn: function
      elseif( trim(cnam).eq."vwrC" ) then
      call set_corona(vtmp)
      nv = ncmax
      call setv2d_mc("wcrt",-1.0d-6,nv,vtmp(1),var,kcn)
!
!----------------------------------------------------------------------
!::NEUT2D including Vac-region
!----------------------------------------------------------------------
!::N0
      elseif( trim(cnam).eq."N0" ) then
      nv = ncmax2
      call setv2d_mc("N0",1.0d-19,nv,tden0(1,ig),var,kcn)
!::E0
      elseif( trim(cnam).eq."E0" ) then
      nv = ncmax2
      call setv2d_mc("E0",1.0d0,nv,teng0(1,ig),var,kcn)
!::P0
      elseif( trim(cnam).eq."P0" ) then
      nv = ncmax2
      vtmp(1:nv) = tden0(1:nv,ig) * teng0(1:nv,ig) * 1.60210d-19
      call setv2d_mc("P0",1.0d0,nv,vtmp(1),var,kcn)
!::Ng
      elseif( trim(cnam).eq."Ng" ) then
      nv = ncmax2
      call setv2d_mc("Ng",1.0d-19,nv,tdeng(1,ig),var,kcn)
!::Eg
      elseif( trim(cnam).eq."Eg" ) then
      nv = ncmax2
      call setv2d_mc("Eg",1.0d0,nv,tengg(1,ig),var,kcn)
!::Pg
      elseif( trim(cnam).eq."Pg" ) then
      nv = ncmax2
      vtmp(1:nv) = tdeng(1:nv,ig) * tengg(1:nv,ig) * 1.60210d-19
      call setv2d_mc("Pg",1.0d0,nv,vtmp(1),var,kcn)
!::Sn
      elseif( trim(cnam).eq."Sn" ) then
          nv = ncmax
         call setv2d_mc("Sn",1.0d0,nv,tssn(1,ia),var,kcn)
!::Sp
      elseif( trim(cnam).eq."Sp" ) then
          vtmp(1:ncmax)=-1.0*tssp(1:ncmax,ia)
          nv = ncmax
          call setv2d_mc("Sp",1.0d0,nv,vtmp,var,kcn)
!::Wi
      elseif( trim(cnam).eq."Wi" ) then
          vtmp(1:ncmax)=-1.0*tswi(1:ncmax)
          nv = ncmax
          call setv2d_mc("Wi",1.0d0,nv,vtmp,var,kcn)
!::We
      elseif( trim(cnam).eq."We" ) then
          vtmp(1:ncmax)=-1.0*tswe(1:ncmax)
          nv = ncmax
          call setv2d_mc("We",1.0d0,nv,vtmp,var,kcn)
!::Ha
      elseif( trim(cnam).eq."Ha" ) then
          call Ha_dist(ndmc, vtmp, nv)
          nv = ncmax
          call setv2d_mc("Ha",1.0d0,nv,vtmp,var,kcn)
!
!
!----------------------------------------------------------------------
!::IMPMC
!----------------------------------------------------------------------
!::Wrd > 0
      elseif( trim(cnam).eq."Wrd" ) then
      nv = ncmax
      call selnty(wmc_nty,nty)
      vtmp(:)=twrdL(:,nty)
!      cvnm = trim(cnam)//" "//trim(catmzL(nty))
      cvnm = trim(cnam)//" "//trim(catmz(nty))
      call setv2d_mc(trim(cnam),1.0d-6,nv,vtmp(1),var,kcn)
!::Nz0
      elseif( cnam(1:2).eq."Nz") then
      chrg = "  "
      if( len_trim(cnam).gt.2 ) chrg = cnam(3:)
! ###2020-09-15 update NAIS [display vac region. (nv = ncmax -> ncmax2)] start
!      nv = ncmax
      nv = ncmax2
      call selnty(wmc_nty,nty)
      call selchrg(cnam,chrg,nv,vtmp,kcn,nty)
!      cvnm = trim(cnam)//" "//trim(catmzL(nty))//" "//trim(chrg)
      cvnm = trim(cnam)//" "//trim(catmz(nty))//" "//trim(chrg)
!      nv = ncmax
      nv = ncmax2
! ###2020-09-15 update NAIS [display vac region. (nv = ncmax -> ncmax2)]  end 
      call setv2d_mc(cvnm,1.0d-19,nv,vtmp(1),var,kcn)
!::ionization
      elseif( cnam(1:3).eq."ion") then
      chrg = "  "
      if( len_trim(cnam).gt.3 ) chrg = cnam(4:)
      nv = ncmax
      call selnty(wmc_nty,nty)
      call selchrg(cnam,chrg,nv,vtmp,kcn,nty)
!      cvnm = trim(cnam)//" "//trim(catmzL(nty))//" "//trim(chrg)
      cvnm = trim(cnam)//" "//trim(catmz(nty))//" "//trim(chrg)
      nv = ncmax
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!::recombination
      elseif( cnam(1:3).eq."rec") then
      chrg = "  "
      if( len_trim(cnam).gt.3 ) chrg = cnam(4:)
      nv = ncmax
      call selnty(wmc_nty,nty)
      call selchrg(cnam,chrg,nv,vtmp,kcn,nty)
!      cvnm = trim(cnam)//" "//trim(catmzL(nty))//" "//trim(chrg)
      cvnm = trim(cnam)//" "//trim(catmz(nty))//" "//trim(chrg)
      nv = ncmax
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!::Vz
      elseif( cnam(1:2).eq."Vz") then
      chrg = "  "
      if( len_trim(cnam).gt.2 ) chrg = cnam(3:)
      nv = ncmax
      call selnty(wmc_nty,nty)
!      call selchrg_avg(cnam,chrg,nv,vtmp,kcn,nty)
!    debug count particle number
      call selchrg(cnam,chrg,nv,vtmp,kcn,nty)
!      cvnm = trim(cnam)//" "//trim(catmzL(nty))//" "//trim(chrg)
      cvnm = trim(cnam)//" "//trim(catmz(nty))//" "//trim(chrg)
      nv = ncmax
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!::Ffr
      elseif( cnam(1:3).eq."Ffr") then
      chrg = "  "
      if( len_trim(cnam).gt.3 ) chrg = cnam(4:)
      nv = ncmax
      call selnty(wmc_nty,nty)
      call selchrg_avg(cnam,chrg,nv,vtmp,kcn,nty)
!      cvnm = trim(cnam)//" "//trim(catmzL(nty))//" "//trim(chrg)
      cvnm = trim(cnam)//" "//trim(catmz(nty))//" "//trim(chrg)
      nv = ncmax
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!::Ffrt
      elseif( cnam(1:4).eq."Ffrt") then
      chrg = "  "
      if( len_trim(cnam).gt.4 ) chrg = cnam(5:)
      nv = ncmax
      call selnty(wmc_nty,nty)
      call selchrg_avg(cnam,chrg,nv,vtmp,kcn,nty)
!      cvnm = trim(cnam)//" "//trim(catmzL(nty))//" "//trim(chrg)
      cvnm = trim(cnam)//" "//trim(catmz(nty))//" "//trim(chrg)
      nv = ncmax
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!::Fth
      elseif( cnam(1:3).eq."Fth") then
      chrg = "  "
      if( len_trim(cnam).gt.3 ) chrg = cnam(4:)
      nv = ncmax
      call selnty(wmc_nty,nty)
      call selchrg_avg(cnam,chrg,nv,vtmp,kcn,nty)
!      cvnm = trim(cnam)//" "//trim(catmzL(nty))//" "//trim(chrg)
      cvnm = trim(cnam)//" "//trim(catmz(nty))//" "//trim(chrg)
      nv = ncmax
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!::Ftht
      elseif( cnam(1:4).eq."Ftht") then
      chrg = "  "
      if( len_trim(cnam).gt.4 ) chrg = cnam(5:)
      nv = ncmax
      call selnty(wmc_nty,nty)
      call selchrg_avg(cnam,chrg,nv,vtmp,kcn,nty)
!      cvnm = trim(cnam)//" "//trim(catmzL(nty))//" "//trim(chrg)
      cvnm = trim(cnam)//" "//trim(catmz(nty))//" "//trim(chrg)
      nv = ncmax
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!::force balance
      elseif( cnam(1:3).eq."FBl") then
      chrg = "  "
      if( len_trim(cnam).gt.3 ) chrg = cnam(4:)
      nv = ncmax
      call selnty(wmc_nty,nty)
      call selchrg_avg(cnam,chrg,nv,vtmp,kcn,nty)
!      cvnm = trim(cnam)//" "//trim(catmzL(nty))//" "//trim(chrg)
      cvnm = trim(cnam)//" "//trim(catmz(nty))//" "//trim(chrg)
      nv = ncmax
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!::force balance
      elseif( cnam(1:4).eq."Fnet") then
      chrg = "  "
      if( len_trim(cnam).gt.4 ) chrg = cnam(5:)
      nv = ncmax
      call selnty(wmc_nty,nty)
      call selchrg_avg(cnam,chrg,nv,vtmp,kcn,nty)
!      cvnm = trim(cnam)//" "//trim(catmzL(nty))//" "//trim(chrg)
      cvnm = trim(cnam)//" "//trim(catmz(nty))//" "//trim(chrg)
      nv = ncmax
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!::Zef
      elseif( trim(cnam) == "Zef" ) then
      do ic = 1, ncmax
        call impvar(ic,nty)
        vtmp(ic) = zef
      enddo
      cvnm = trim(cnam)
      nv = ncmax
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!::<Z>
      elseif( trim(cnam) == "Zav" ) then
      do ic = 1, ncmax
        call impvar(ic)
        vtmp(ic) = zav
      enddo
      cvnm = trim(cnam)
      nv = ncmax
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!
!::rt0/rt1
      elseif( cnam(1:2).eq."rt" ) then
      call selnty(wmc_nty,nty)
      call set_rt(cnam, nty, vtmp, nv, cvnm)
      call setv2d_mc(cvnm,1.0d0,nv,vtmp(1),var,kcn)
!
!::??
      else
        kcn = 1
        write(n6,'(/2x,"invarid cnam  ",a)') trim(cnam)
        call gdput("  ")
        call gdput("************************************")
        call gdput("Error: invarid value name, try again")
        call gdput("************************************")
        goto 100
      endif
!
      if( kcn.ne.0 ) then
        write(n6,'(/2x,"invalid data  ",a)') trim(cnam)
        call gdput("  ")
        call gdput("************************************************")
        call gdput("Error: invarid plot data (e.g. all value <= 0)")
        call gdput("************************************************")
        goto 100
      endif
!
!::check
!      if( cnam(1:2).eq."rt" ) then
!      do ic = 1, ncmax, 500
!      write(6,'(2x,"CHK rt  ic =",i6,2x,1pe12.3)')
!     >  ic, var(ic)
!      enddo
!      endif
!
!::log(Wrd)
      if( klog.eq.1 ) then
      do ic = 1, nv
        if( var(ic).gt.0.0d0 ) then
          var(ic) = alog10(var(ic))
        else
          var(ic) = 0.0d0
        endif
      enddo
      endif
!
!::check
!      if( cnam(1:2).eq."rt" ) then
!      do ic = 1, ncmax, 500
!      write(6,'(2x,"CHK log10(rt)  ic =",i6,2x,1pe12.3)')
!     >  ic, var(ic)
!      enddo
!      endif
!
      if(plot_mode==1) then
!     ::void  ( var < -0.1d30 )
      if( nv.gt.ncmax) then
        call gdget(">> display vac region (y:<rtn>/n) ==> ",cinp,mji)
        if( mji.eq.0 .or. cinp(1:1).eq."y" ) then
          write(n6,'(2x,"  disp vac region  nv, ncmax, ncmax2 =",3i6)')
     >        nv, ncmax, ncmax2
        else
          nv = ncmax2
          call setv2d_mc("Ti",1.0d0,nv,temi(1),varTi,kcn)
          do ic = 1, nv
            if( varTi(ic).le.0.0d0 ) var(ic) = -1.0d30
          enddo
        endif
      else
        call clear_void(var,ncmax2)
      endif
!
        call gdpl2dN(cvnm,nv,var,"  ")  
!
      elseif(plot_mode==2) then    
        var2 = dble(var)
        call chord_view(cvnm, nv, var2)
      endif
!
      goto 120
!
 150  continue
      return
      end



!***********************************************************************
      subroutine clear_void(var,ncmax2)
!***********************************************************************
      use csize
      use cimcom
      use cimden
      use cntpls
      use cntsrc
      use cplcom
      use cplimp
      use cplwrd
      use cplwrd2
      implicit none

      integer ncmax2
      real(4), dimension(ndmc) :: var, varTI
      integer nv, ic, kcn

      nv = ncmax2
      call setv2d_mc("Ti",1.0d0,nv,temi(1),varTi,kcn)
      do ic = 1, nv
        if( varTi(ic).le.0.0d0 ) var(ic) = -1.0d30
      enddo

      return
      end
