!**********************************************************************
      subroutine prg_prof(kitr,kcnd)
!**********************************************************************
!)
!)      sequential type    2015/07/23
!)
!)   Mnkntl_pls (vna,vti), (flps)  => NEUT2D
!)   Mnkimp_pls (vna,vti), (flps), (q1a), (N0, T0) => IMPMC
!)
!)
!)    order  code     IN-dtyp      => OUT-dtyp
!)   ----------------------------------------------------------
!)      1    soldor                => PLS_2B, PLS_2
!)      2    neut2d   PLS_2B       => NTL_2
!)      3    IMPMC    PLS_2        => IMP_2
!)      4    soldor   NTL_2, IMP_2 => PLS_2B, PLS_2
!)
!::SOLDOR
!)    PLS_2, PLS_2B contain flps calculated in sub. plpflx
!)    plmrgn Not depend on PLS_2(B)
!)    plwtfy depend on q1a
!)    plpflx depend on plwtfy and send to neut2d IMPMC
!)    plmflx tflne tflqi tflqe for output
!)
!----------------------------------------------------------------------
      use cimcom,      only : amz, azmas, denflxztowall, eneflxztowall
     >    , e0bkz, e0emtz, e0inz, e0tmz, emflx, ismax, ltmax, npemt
     >    , npmax, npmaxl, sptyc, tfbz, v0emtz
      use cimctl,      only : kimp, lstdy, icalZ
      use cimden,      only : csput, eipot, fsput, npsput, nsput
     >     , nsrc_spt, nsrc_spt_dummy, nwspt, nzmx
      use cntctl,      only : kntl
      use cphcns,      only : cev
      use cputim,      only : cputot, eltm0, eltm00, eltm5
      use csonic,      only : is_rsta, kpcn
      use cunit,       only : n6
      use mod_dtflw,   only : rst_drr
      use mod_keylist, only : knorm
      use mod_mpicomm, only : m6
      use mod_shexe,   only : impmc_model
      use mpi!,         only : mpi_wtime
      implicit none

!::argument
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
      character(80) :: cmsg
      integer :: k, ii, i

      nsrc_spt = 0
      nsrc_spt_dummy = 0

      if( impmc_model == 0 .or.
     >    ( impmc_model == 1 .and. kitr /= 2 ) ) then
      select case(kitr)
!--------------------------------------------------------------------
      case( 1, 2 )
!--------------------------------------------------------------------
!::option => simple
        is_rsta = .false.
        kntl = 1
        kimp = 1
        kpcn = 1

        write(m6,'(2x,"prg_prof IMPMC STA")')
        ! SYamoto ======
        write(n6,*) 'RSTA-RSTF: check imp_cal will be skipped or not'
        if (rst_drR(1:1) == ".")then
          call rstfile(2) ! SYamoto read rst-file(if exist)
          call debg_imp_cal
          if (is_rsta) then
            call imp_recd
            if(impmc_model == 0) icalZ = 0 ! set for set bkflw distribution
            write(n6,*) 'RSTA-RSTF: skip imp_cal'
          endif
        endif
        if(.not.is_rsta)then
          write(n6,*) 'RSTA-RSTF: no rst file. calling imp_cal'
          call flush(m6)
          call imp_recd   ! No call MPI routine
          call dbg_snflx("prg_mcini-at_MPprmP13")
          call ntl_pls    ! <== add
          call imp_pls
          call dbg_snflx("prgexec-bf_impcal")
          if( impmc_model == 0 ) then
            call imp_cal
            call imp_end
          else
          ii = 0
          npmax = 0 ! SY initial set of npmax.
          npmaxL(:) = 0 ! SY initial set of npmaxL

          do k= 1, nsput
            nsrc_spt = k
            nwspt  = csput(k) !hardcoded at the moment
            npemt  = npsput(k)
            write(n6,*) 'Source', nwspt, npemt
!           npmax = 0 ! SY for initial benchmarking
!           npmax  = npsput(k) SY npmax may not be determined here
            e0emtZ = e0inZ
            if( nwspt.eq."Arpf" ) e0emtZ = e0inZ
            if( nwspt.eq."Arbk" ) e0emtZ = e0bkZ
            if( nwspt.eq."Cphy" ) e0emtZ = e0tmZ
            v0emtZ = sqrt(2.0d0*e0emtZ*cev/amz)

            denFlxZtoWall = 0.0d0
            eneFlxZtoWall = 0.0d0

            call iminit(nwspt)
            if(tfbz.le.0.0d0)then
              call imsave(nwspt,"w",k) !write dummy imsave
              call imsave(nwspt,"y",k) !write dummy imsave
              if(nsrc_spt.eq.nsput) call imp_reset
              cycle
            endif

            nzmx = ismax !moved from immont.f
            ii = ii + 1
            sptyc = csput(ii)
            fsput(ii) = emflx

            write(n6,'(2x,"impurity generation ",i2,2x,a,2x,1pe12.3,
     >        "  azmas =",0pf7.3,"  ltmax =",i6,"  npmax =",i6)')
     >        ii, csput(ii), fsput(ii), azmas, ltmax, npmax
            write(n6,'(2x,"eipot =",10f10.3)') (eipot(i),i=0,ismax)

            call imp_prof
            call imsave(nwspt,"w",k)
            call imsave(nwspt,"y",k)
            call dbg_denz("imsave")
          enddo !k

          call imp_wrd
        endif

        call dbg_snflx("prgexec-af_impend")
        call debg_imp_cal
        endif
        ! ====== SYamoto
        write(m6,'(2x,"prg_prof IMPMC END")')
        call flush(m6)
!
!::stedy cal.
         if( lstdy.ne.0 ) then
           eltm5 = MPI_WTIME()
           cputot = eltm5-eltm00
           write(n6,'(/2x," total cpu = ", f10.3,
     >      " sec.",f6.2," hour")') cputot,  cputot/3600.0
           call wexit("prg_prof","lstdy == 1")
         endif
!--------------------------------------------------------------------
      case default
!--------------------------------------------------------------------
        write(cmsg,'("kitr =",i2," > 1")') kitr
        call wexit("prg_prof",trim(cmsg))
      end select
      endif

      if( kitr == 2 ) then
        eltm0 = MPI_WTIME()
        write(n6,'(2x," initial cpu ",f10.3,"  sec.")') eltm0-eltm00
      endif

      kcnd = knorm
      return
      end
