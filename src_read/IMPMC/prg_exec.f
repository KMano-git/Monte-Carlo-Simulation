!***********************************************************************
      subroutine prg_exec(kitr,kcnd)
!***********************************************************************
      use cimcom,      only : amz, azmas, denflxztowall, e0bkz, e0emtz
     >    , e0inz, e0tmz, emflx, eneflxztowall, ismax, ltmax, npemt
     >    , npmax, sptyc, v0emtz, sflux
      use cimctl,      only : kimp
      use cimden,      only : csput, eipot, fsput, npsput, nsput
     >    , nsrc_spt, nwspt, nzmx
      use cphcns,      only : cev
      use csonic,      only : itend, itim, lstop
      use cunit,       only : n6
      use dbg_mod,     only : dbgcntl
      use mod_keylist, only : klast, knorm, kstop
      use mod_shexe,   only : impmc_model
      implicit none

!::argument
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
      integer :: k, ii, i
      character :: cmsg*80

      select case(kitr)
!--------------------------------------------------------------------
      case(1)
!--------------------------------------------------------------------
      kimp = 1
      call imp_recd
      call dbg_snflx("prgexec-at_MPprmP13")
      call ntl_pls
      call imp_pls
!  SOLDOR-IMPMC data transfer test  2018/05/16 Y.Homma
!  !!! Erase when test is done, this call could slow down the code.
!      call yh_heat_flux_plot
!
      if( impmc_model == 0 ) then
        call dbg_snflx("prgexec-bf_impcal")
        call imp_cal
        call imp_end
        call dbg_snflx("prgexec-af_impend")
      else
        call dbg_snflx("prgexec-bf_impexec")
        ii = 0
        do k= 1, nsput
          nsrc_spt = k
          nwspt  = csput(k) !hardcoded at the moment
          npemt  = npsput(k)
          write(n6,*) 'Source', nwspt, npemt
          e0emtZ = e0inZ
          if( nwspt.eq."Arpf" ) e0emtZ = e0inZ
          if( nwspt.eq."Arbk" ) e0emtZ = e0bkZ
          if( nwspt.eq."Cphy" ) e0emtZ = e0tmZ
          v0emtZ = sqrt(2.0d0*e0emtZ*cev/amz)

          denFlxZtoWall = 0.0d0
          eneFlxZtoWall = 0.0d0

          call iminit(nwspt)
          if(sflux.le.0.0d0 .and. npmax <= 0) then
            call imsave(nwspt,"w",k)
            call imsave(nwspt,"y",k)
            cycle
          endif
          nzmx = ismax ! moved from immont.f
          ii = ii + 1
          sptyc = csput(ii)
          fsput(ii) = emflx

          write(n6,'(2x,"impurity generation ",i2,2x,a,2x,1pe12.3,
     >      "  azmas =",0pf7.3,"  ltmax =",i6,"  npmax =",i6)')
     >     ii, csput(ii), fsput(ii), azmas, ltmax, npmax
          write(n6,'(2x,"eipot =",10f10.3)') (eipot(i),i=0,ismax)

          call imp_exec
          call imsave(nwspt,"w",k)
          call imsave(nwspt,"y",k)
          call dbg_denz("imsave")
        end do
        call imp_wrd !imp_end replaced
        call dbg_snflx("prgexec-af_impwrd")
      endif
      call debg_imp_cal
      call dbgcntl( 1 )

      if( lstop.eq.1 ) goto 200
      if( mod(itim,100).eq.0 ) call flush(n6)
      if( itim.eq.itend ) goto 200

!::loop
      kcnd = knorm
      return

 200  continue
      if( impmc_model == 1 ) call imp_rstfl(1) !YAMOTO
      kcnd = klast
      if( lstop /= 0 ) kcnd = kstop
      return

!--------------------------------------------------------------------
      case default
!--------------------------------------------------------------------
      write(cmsg,'("kitr =",i2," > 1")') kitr
      call wexit("prg_init(monte)",trim(cmsg))
      end select
      end
