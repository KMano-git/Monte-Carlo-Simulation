!**********************************************************************
      subroutine imp_ini
!**********************************************************************
      use cimcns,    only : iseed, ndseed, nseed
      use cimcom,    only : bkflux, bkfrac, lwrad, ndbkf, timez
      use cimctl,    only : icalz, nprg, tmimp
      use cimden,    only : nzmx, tdnz, tfrz, tionz, tradiz
     >    , tradliz,  tradrz, trecz,  tthz, tvlz, twci, twrd
      use cimi60,    only : i60_imcndv, i60_imhist, i60_tmwgtn
      use csize,     only : ndmc, ndmis
      use csonic,    only : lpost, lrand
      use cunit,     only : cdnm, lmspe, lmype, mygrp, mype, n5, n6
     >    , stcpy, ufl6
      use mod_shexe, only : impmc_model, INP_IMP, use_namelist, MESH_NTL
      use mod_externalgrid, only:use_exdata
      implicit none
!
!::local variables
      integer   nft
      character cdsn*80
      logical   lex
      integer   ic, iz
      integer,parameter::  ndim = 30
!
      call trmark("imp_ini","start")
      write(n6,'(/2x,30("="),"  imp_ini")')
!
      call phcnst
!
      write(n6,*) "set nprg, mygrp, cdnm  =", nprg, mygrp,cdnm(mygrp)
!
!::input data
      nft = 21
      close(nft)
      if(use_namelist) then
         !.. use variable in inppls
         call ncopy(nft,INP_IMP,"text",cdsn,lex)
      else
         !.. use environment variable
         call ncopy(nft,"INP_IMP","text",cdsn,lex)
      endif
      open(unit=n5, file=cdsn)
      call iminpt(n5)
      call imfile(n5)
      if( impmc_model == 1 ) call imp_file("setunt") ! profile/history
      close(n5, status=stcpy)

      ! use external grid data in Gmesh
      if(use_exdata) then
        !::KSFUJI  need data of mcel for plauxv
        if(use_namelist) then
        !..   use variable in inppls
          call nopen(nft,MESH_NTL,"binary",cdsn,lex)
        else ! use environment variable
          call nopen(nft,"MESH_NTL","binary",cdsn,lex)
        endif
        call ntgdsk(nft,"read")
        close (nft)
      endif

!::lwrad <= mdl_wrd
!      lwrad = mdl_wrd  !--> dtypeset
!
      write(n6,'(2x,"call imcnst in sub. imp_ini for 3 Groups")')
      call imcnst
!
!::area & thdeg for output of sputtering
      call ntwdeg
!
!::set variables
      call set_bfld
      call set_tabl
!
      call set_roh
      call set_edrfl
      call set_adas
!
!::add mesh data
      call dstset
      if(.not.use_exdata) then
        call lstgate
      endif
      call set_puff
      call set_bkgt
      call set_limt
!
      if( lrand.eq.1 ) then
      nseed = ndseed
      call rnseed(nseed,iseed)
      endif
!
      write(n6,'(2x,"imp_ini  No call Bcast")')
      write(n6,'(2x,"lwrad = ",i2)') lwrad
!
!::control of imp-cal
      tmimp = -999.0
      icalZ = 0
      timeZ = 0.0d0
!
!::clear
      call imclear("particle")
      call imclear("score")
      call imclear("floss")
!
!::KSFUJI   twrd in plcrad_mc
!::clear /cimden/
      nzmx = 0
      do ic = 1, ndmc
      twrd(ic) = 0.0d0
      twci(ic) = 0.0d0
      do iz = 0, ndmis
      tdnz(iz,ic) = 0.0d0
      tfrz(iz,ic) = 0.0d0
      tthz(iz,ic) = 0.0d0
      tvlz(iz,ic) = 0.0d0
      tionZ(iz,ic) = 0.0d0
      trecZ(iz,ic) = 0.0d0
      tradiZ(iz,ic) = 0.0d0
      tradliZ(iz,ic) = 0.0d0
      tradrZ(iz,ic) = 0.0d0
      enddo
      enddo
!
!::conductance of gate  calculate at sub. imcndv
      bkfrac(0:ndbkf) = 0.0d0
      bkflux(0:ndbkf) = 0.0d0
!
!::diffusion coef
!     This part is moved to imp_pls.
!
!::define i60  edit inc/cimi60
      if( lpost.eq.0 ) then
      i60_imcndv = 360000 + mype  !
      i60_imhist = 360000 + mype  !
      i60_tmwgtn = 370000 + mype  !
      endif
!
      if( trim(ufl6).eq."/dev/null" ) then
      i60_imcndv = 0
      i60_imhist = 0
      i60_tmwgtn = 0
      endif
!
      if( lmype.ne.lmspe ) then
      endif
!
      write(n6,'(2x,"i60_imcndv =",i8,2x,a)')
     >   i60_imcndv, "evaluation of impurity back flow"
      write(n6,'(2x,"i60_imhist =",i8,2x,a)')
     >   i60_imhist, "number of stick particles"
      write(n6,'(2x,"i60_tmwgtn =",i8,2x,a)')
     >   i60_tmwgtn, "position and veloicty at back flow partciles"
!
 100  continue
      call trmark("imp_ini","return")
      return
      end
