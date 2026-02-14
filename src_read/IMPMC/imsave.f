!***********************************************************************
      subroutine imsave(cgty,cact,k_inout)
!***********************************************************************
      use cimcom,      only : denz, friz, hcnpt, hcnwt, hgrpp, hgrpw
     >    , hregp, hregw, hrgpa, hrgpc, idcal, idemt, iderr, idetc
     >    , idexh, idstk, igcal, igemt, igerr, igexh, igmin, ignul
     >    , igpmp, igrun, igsce, igsgt, ionz, ndis, ndmp, phdnz, phfrz
     >    , phionz, phrecz, phthz, phvlz, phwci, phwrd, recz, sdmy1
     >    , sdmy2, sdmy3, sdtmz, sflux, sitmp, sptcl, sptyc, stimp
     >    , stimz, swtot, temz, thfz, vzpz, wsct
     >    , sflux_save, swtot_save, sptyc_save
     >    , denZ_save, temZ_save, wsct_save, hrgpc_save, phwrd_save
     >    , phwci_save, phdnz_save, phvlz_save,hgrpp_save, hgrpw_save
     >    , hregp_save, hregw_save, friZ_save, thfZ_save
     >    , vzpZ_save, ionZ_save,recZ_save
     >    , phfrz_save, phthz_save, phionz_save, phrecz_save
     >    , stimz_save, hcnwt_save, hrgpa_save, sptcl_save
     >    , hcnpt_save, stimp_save, sitmp_save, sdtmz_save
     >    , sdmy1_save, sdmy2_save, sdmy3_save
      use cimctl,      only : nprg
      use csize,       only : ndmc
      use cunit,       only : lmspe, lmype, lnope, mywld, n6
      use mod_dtypdef, only : nddt, typ_dnam, typ_itag
      use mod_shexe,   only : impmc_model
      use mpi!,         only : mpi_bcast
      implicit none
!
!::argument
      character, intent(in) :: cgty*(*), cact*(*)
      integer, intent(in) :: k_inout
!
!::local variables
      integer  nsizp, nsizs, nsizc
      integer  nft
      integer  ierr, itag, ityp
      integer    i
      real(8) :: totp, totw

!-----------------------------------------------------------------------
!::write
!-----------------------------------------------------------------------
      if( cact(1:1).eq."w" ) then
        if( lmype.ne.lmspe ) return
        sflux_save(k_inout) = sflux
        swtot_save(k_inout) = swtot
        stimz_save(k_inout) = stimz
        sptyc_save(k_inout) = sptyc
        denZ_save(:,:,k_inout) = denZ
        temZ_save(:,:,k_inout) = temZ
        wsct_save(:,k_inout) = wsct
        ! following values are for output and/or TOPICS, not used for calculation in SONIC
        sptcl_save(k_inout) = sptcl
        sitmp_save(k_inout) = sitmp
        stimp_save(k_inout) = stimp
        sdtmz_save(k_inout) = sdtmz
        sdmy1_save(k_inout) = sdmy1
        sdmy2_save(k_inout) = sdmy2
        sdmy3_save(k_inout) = sdmy3
        if(impmc_model == 0) then
          hrgpc_save(:,k_inout) = hrgpc
        ! following values are just for log output, not used in calculation
          hrgpa_save(:,k_inout) = hrgpa
          hcnpt_save(:,k_inout) = hcnpt
          hcnwt_save(:,k_inout) = hcnwt
        else
          phwrd_save(:,k_inout) = phwrd
          phwci_save(:,k_inout) = phwci
          phdnz_save(:,:,k_inout) = phdnz
          phvlz_save(:,:,k_inout) = phvlz
          hgrpp_save(:,k_inout) = hgrpp
          hgrpw_save(:,k_inout) = hgrpw
          hregp_save(:,k_inout) = hregp
          hregw_save(:,k_inout) = hregw
        endif
        ! debug write
         write(n6,'(2x,"imsave WRITE   sptyc =",a,
     >     2x,i6,"  time(pls) =",1pe12.4,"  time(imp) =",1pe12.4,
     >     "  sflux =",1pe12.4,"  swtot =",1pe12.4)')
     >     sptyc, int(sitmp+1.0d-8), stimp, stimz,
     >     sflux, swtot
!-----------------------------------------------------------------------
!::write Yamoto
!-----------------------------------------------------------------------
      elseif( cact(1:1).eq."y" ) then
        if( lmype.ne.lmspe ) return
        if(impmc_model == 0) then
          friZ_save(:,:,k_inout) = friZ
          thfZ_save(:,:,k_inout) = thfZ
          vzpZ_save(:,:,k_inout) = vzpZ
          ionZ_save(:,:,k_inout) = ionZ
          recZ_save(:,:,k_inout) = recZ
        else
          phfrz_save(:,:,k_inout) = phfrz
          phthz_save(:,:,k_inout) = phthz
          phionz_save(:,:,k_inout) = phionz
          phrecz_save(:,:,k_inout) = phrecz
        endif
        ! debug write
        write(n6,'(2x,"imsave YAMOTO WRIT  sptyc =",a,
     >     2x,i6,"  time(pls) =",1pe12.4,"  time(imp) =",1pe12.4,
     >     "  sflux =",1pe12.4,"  swtot =",1pe12.4)')
     >     sptyc, int(sitmp+1.0d-8), stimp, stimz,
     >     sflux, swtot
!-----------------------------------------------------------------------
!::read
!-----------------------------------------------------------------------
      elseif( cact(1:1).eq."r" ) then
        if( lmype.eq.lmspe ) then
          sflux = sflux_save(k_inout)
          swtot = swtot_save(k_inout)
          stimz = stimz_save(k_inout)
          sptyc = sptyc_save(k_inout)
          denZ = denZ_save(:,:,k_inout)
          temZ = temZ_save(:,:,k_inout)
          wsct = wsct_save(:,k_inout)
          ! following values are just for log output, not used in calculation
          sptcl = sptcl_save(k_inout)
          sitmp = sitmp_save(k_inout)
          stimp = stimp_save(k_inout)
          sdtmz = sdtmz_save(k_inout)
          sdmy1 = sdmy1_save(k_inout)
          sdmy2 = sdmy2_save(k_inout)
          sdmy3 = sdmy3_save(k_inout)
          if(impmc_model == 0) then
            hrgpc = hrgpc_save(:,k_inout)
          ! following values are just for log output, not used in calculation
            hrgpa = hrgpa_save(:,k_inout)
            hcnpt = hcnpt_save(:,k_inout)
            hcnwt = hcnwt_save(:,k_inout)
          else
            phwrd = phwrd_save(:,k_inout)
            phwci = phwci_save(:,k_inout)
            phdnz = phdnz_save(:,:,k_inout)
            phvlz = phvlz_save(:,:,k_inout)
            hgrpp = hgrpp_save(:,k_inout)
            hgrpw = hgrpw_save(:,k_inout)
            hregp = hregp_save(:,k_inout)
            hregw = hregw_save(:,k_inout)
          endif
          ! debug write
          write(n6,'(2x,"imsave READ    sptyc =",a,
     >     2x,i6,"  time(pls) =",1pe12.4,"  time(imp) =",1pe12.4,
     >     "  sflux =",1pe12.4,"  swtot =",1pe12.4)')
     >     sptyc, int(sitmp+1.0d-8), stimp, stimz,
     >     sflux, swtot
        endif
        if( lnope.gt.1 ) then
          if( impmc_model == 0 ) then
            call tbfind( nddt, typ_dnam, ityp, 'CIMCOM_11' )
          else
            call tbfind( nddt, typ_dnam, ityp, 'CIMCOMT11' )
          endif
          itag = typ_itag(ityp)
          call MPI_Bcast( sflux, 1, itag, lmspe, mywld, ierr )
          call cimcom_lsr( 1 )
        endif
!-----------------------------------------------------------------------
!::read Yamoto
!-----------------------------------------------------------------------
      elseif( cact(1:1).eq."z" ) then
        if( lmype.eq.lmspe ) then
          if(impmc_model == 0) then
            friZ = friZ_save(:,:,k_inout)
            thfZ = thfZ_save(:,:,k_inout)
            vzpZ = vzpZ_save(:,:,k_inout)
            ionZ = ionZ_save(:,:,k_inout)
            recZ = recZ_save(:,:,k_inout)
          else
            phfrz = phfrz_save(:,:,k_inout)
            phthz = phthz_save(:,:,k_inout)
            phionz = phionz_save(:,:,k_inout)
            phrecz = phrecz_save(:,:,k_inout)
          endif
          ! debug write
          write(n6,'(2x,"imsave YAMOTO READ   sptyc =",a,
     >     2x,i6,"  time(pls) =",1pe12.4,"  time(imp) =",1pe12.4,
     >     "  sflux =",1pe12.4,"  swtot =",1pe12.4)')
     >     sptyc, int(sitmp+1.0d-8), stimp, stimz,
     >     sflux, swtot
        endif
        if( lnope.gt.1 ) then
          call cimcom_lsr( 2 )
        endif
      else
        call wexit("imsave","invalid input")
      endif
!
!::debug write
      if( impmc_model == 0 ) then
        write(n6,'(2x,"sptyc =",a,"  sflux =",1pe12.3,"  swtot =",
     >    1pe12.3,"   sptcl =",1pe12.3,"  itim =",i6,"  time =",
     >    1pe12.3,"  timeZ =",1pe12.3)')
     >   sptyc, sflux, swtot, sptcl, int(sitmp+1.0d-8), stimp, stimz
!
        write(n6,'(2x,"emt =",1pe12.4,"  cal =",1pe12.4,
     >    "  exh/stk =",1p2e12.4,"  err/etc =",1p2e12.4)')
     >     hcnpt(idemt), hcnpt(idcal), hcnpt(idexh), hcnpt(idstk),
     >     hcnpt(iderr), hcnpt(idetc)
      else
        totp = 0.0d0
        totw = 0.0d0
        do i = 2, 10   ! exclude emt (i=1) and pmp(i=8) (exh)
          if( i == igemt .or. i == igpmp ) cycle
          totp = totp + hgrpp(i)
          totw = totw + hgrpw(i)
        enddo

        write(n6,'(2x,"imsave",4x,2x,a,1pe12.4,"  em/cl",1p2e12.4
     >    "  gt/ce",1p2e12.4,"  mn/er/ex/pm",1p4e12.4,"  nl/rn",
     >    1p2e12.4)') cact(1:1),
     >    totw, hgrpw(igemt), hgrpw(igcal), hgrpw(igsgt), hgrpw(igsce),
     >   hgrpw(igmin), hgrpw(igerr), hgrpw(igexh), hgrpw(igpmp),
     >   hgrpw(ignul), hgrpw(igrun)
      endif
!
      ! debug Yamoto
      if( impmc_model == 1 ) call dbg_denz("imsavr")
      call dbg_vzpz("imsavr")

      return
      end