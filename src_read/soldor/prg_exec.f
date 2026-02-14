!***********************************************************************
      subroutine prg_exec(kitr,kcnd)
!***********************************************************************
!)
!)::OLD typ (wxdr_103)
!)  exec 1 : soldor                 =>  MST_2  PLS_2B  PLS_2
!)  exec 2 : neut2d  MST_2  PLS_2B  =>  NTL_2
!)  exec 2 : IMPMC   MST_2  PLS_2   =>  IMP_2
!)  exec 3 : soldor  NTL_2  IMP_2
!)
!)::CHK_1 type (wxdr_106) completely same as wxdr_103
!)  exec 1 : soldor (only plmdprf)  =>  MST2 PLS_2B PLS_2
!)  exec 2 : neut2d  MST_2  PLS_2B  =>  NTL_2
!)  exec 2 : IMPMC   MST_2  PLS_2   =>  IMP_2
!)  exec 3 : soldor  NTL_2  IMP_2
!)
!)::CHK_3 type (wxdr_107)
!)  exec 1 : soldor (only plmdprf)  =>  MST2 PLS_2B PLS_2
!)  exec 2 : neut2d  MST_2  PLS_2B  =>  NTL_2
!)  exec 2 : IMPMC   MST_2  PLS_2   =>  IMP_2
!)  exec 3 : soldor  NTL_2  IMP_2
!)
!)::NEW type
!)  exec 1 : neut2d  MST_2  PLS_2B  =>  NTL_2
!)  exec 1 : IMPMC   MST_2  PLS_2   =>  IMP_2
!)  exec 2 : soldor  NTL_2  IMP_2
!)
!)
!)  MST_2   time, dtim, kntl, kimp   Mnkall_tim
!)  PLS_2B  vni, vte, vti, flux      Mnkntl_pls
!)  PLS_2   vti, animp, q3, tN0      Mnkimp_pls
!)  NTL_2   vflux, vnsrc, hti ?      Mnkntl_cal
!)  IMP_2   nsput, twrd, tdnz        Mnkimp_cal
!)
!)  NTL_1   bvx, bvy, mcel, cswl     Mnkntl_iniSnd
!)  IMP_1   xdaty,ip_ion,amz,ismax   Mnkimp_iniSnd
!)
!)  define msts (0:continue/1:last exec/2:normal end/-1:stop
!)
!)   SPMD-SONIC  neut2d IMPMC soldor(itim=itim+1)
!)   MPMD-SONIC  soldor_1 neut2d IMPMC soldor_2(itim=itim+1)
!)
!)
!)   prg_exec
!)       *------ pls_cal
!)                  *------ plstep
!)                            *----- time = time + dtim
!)               itim = itim + 1
!)
!----------------------------------------------------------------------
      use cntctl,      only : itmnt, kntl, nclnt
      use cntmnt,      only : vsmty
      use cplcom,      only : dlpq, dtimb, dtmax, dtmin, dtnxt, dttb
     >    , edtb, edtmx, edtq,  elpmx, eltb, itcn, ittb, lttb, nttb
     >    , snflx
      use cputim,      only : cntstp, cpugrp, cpustp, eltm0, eltm1
     >    , eltm2, eltm3, eltm4
      use csonic,      only : dtim, itend, itim, kdsk, khst, kpcn, lcstg
     >    , lstop, mxcpu, ndsk, nhst, npcn, tend, time
      use cunit,       only : n6
      use dbg_mod,     only : dbgcntl
      use mod_keylist, only : klast, knorm, kstop
      use mod_loc_tim, only : msts
      use mpi!,         only : mpi_wtime
      use topics_mod,  only : dtcal, dtduc, dtim_s, kcon, lexdt, lgtpc
      implicit none

!::argument
      integer, intent(in)  :: kitr
      integer, intent(out) :: kcnd

!::local variables
      integer :: ichk, kset, itend2
      logical :: lex
      real(8) :: tfac, cpmm, fclp, fcdt
      character :: cmsg*80
      integer :: l_heat_flux_output = 0
      real(8) :: dtminh = 0.0_8
      integer :: ksre = 0
! dtminh : time 0.0 judgment threshold
! ksre   : end time flag ( = 1 : end time, = 0 : running )

      dtminh = dtmin * 0.5_8

      select case(kitr)
!--------------------------------------------------------------------
      case(1)
!--------------------------------------------------------------------
        if( itim.gt.itend ) goto 200
        if( time.ge.tend  ) goto 200

!::flag
        kpcn = 0
        kdsk = 0
        khst = 0
        if( mod(itim,npcn).eq.0 ) kpcn = 1
        if( mod(itim,ndsk).eq.0 ) kdsk = 1
        if( mod(itim,nhst).eq.0 ) khst = 1
        if( mod(itim+1,npcn) == 0 ) kpcn = 1

!::cancel
        if( mod(itim,50).eq.0 ) then
          inquire( file='.quit', exist=lex )
          if( lex ) then
            write(n6,'(/2x,"*** Job was cancelled. ***  time =",
     >        1pe12.3,"  itim =",i6)') time, itim
            itend = itim+1
          endif
        endif

        if( itim.eq.itend )   kpcn = 1
        if( itim+1.eq.itend ) kpcn = 2

!::profile in main plasma   ichk: control of output
      ichk = 0
      if( kpcn.eq.1 ) ichk = 1
      kset = 0
      if( lgtpc == 0 ) call plmdprf(kset,ichk)

!::time step control
        if( itim.gt.itcn ) then
          ittb = ittb + 1
          if( ittb.gt.nttb ) ittb = nttb
          itcn  = lttb(ittb)
          dtim  = dttb(ittb)
          dtmax = dttb(ittb)
          elpmx = eltb(ittb)
          edtmx = edtb(ittb)
          dtim  = dmin1(dtim,1.2d0*dtimb)
          tfac  = dtim/dtimb
          write(n6,'("<dtcntl-so> ",i6," set","  dtim =",1pe11.3,
     >     "  tfac =",0pf8.3,"  dtimb,elpmx,edtmx =",1p3e11.3)')
     >     itim,dtim,tfac,dtimb,elpmx,edtmx
          if( dtim.le.0.0.or.elpmx.le.0.0.or.edtmx.le.0.0 ) goto 290
          dtim_s = dtim
          if( lgtpc == 1 .and. dtim > dtduc ) then
            if( dtduc > 0.0_8 ) then
              dtim   = dtduc
            else
              dtim   = min( dtim, dtcal )
            endif
          endif
        endif


      call ntl_skip(kntl)

      eltm1 = MPI_WTIME()

      call plnflx(snflx)   ! use snflx evaluated in pwstep
      call pldens("NTcal")
      call dbg_snflx("prgexec-at_MPcalN21")

      eltm2 = MPI_WTIME()
      cpugrp(2) = eltm2-eltm1
      cntstp = eltm2-eltm1

      call dbg_snflx("prgexec-at_MPprmP13")
      call mcplas ! for dene , etc
      eltm3 = MPI_WTIME()
      cpugrp(3) = eltm3-eltm2

!::time step
!::large change in neut2d loss term
      if( kntl.ge.1 ) then
        if( vsmty.eq.2 ) then
          dtim = dmin1( dtim, 2.0d-7 )
          dtim_s = dtim
          if( lgtpc == 1 .and. dtim > dtduc ) then
            if( dtduc > 0.0_8 ) then
              dtim   = dtduc
            else
              dtim   = min( dtim, dtcal )
            endif
          endif
        endif
      endif

!::advance time step
      call pls_cal

      if( lgtpc == 1 .and. lexdt == 1 ) call plmflx

!::profile in hot core
      ichk = 0
      if( kpcn == 1 ) ichk = 1
      kset = 2
      if( lgtpc == 0 ) call plmdprf(kset,ichk)

      itim  = itim + 1    ! <== Note
!
!::soldor(time,itim) => mod_loc_tim(master)
      call set_wftime

!::wf_tmadvan
!     => mod_glb_tim => Bcast glb_tim => mod_loc_tim

!::cpu time
      eltm4 = MPI_WTIME()
      cpugrp(1) = eltm4-eltm3
      cpustp = eltm4-eltm1

!::cpu time (minutes)
        cpmm = (eltm4-eltm0 )/60.0
        if( cpmm.gt.dfloat(mxcpu) ) then
          itend2 = ((itim+1)/20+1)*20
          if(itend2 .lt. itend ) then
            if( lgtpc == 0 ) itend = itend2
            write(n6,'(/2x,"*** Finished cal. due to cpu-over ***",
     >       "  cpu(min) =",2f12.3,"  itim =",i6,"  itend =",i6/)')
     >      cpmm,dfloat(mxcpu),itim,itend
            mxcpu = 999999
          endif
        endif

!::output
      call pls_out

      call dbgcntl( 3 )

! send edge data to TOPICS
      if( lgtpc == 1 .and. lstop == 0 ) then
        dtduc = dtduc - dtim
        ksre = 0
        if( abs( time - tend ) < dtminh ) ksre = 1
! recieve profile data from TOPICS
        if( ksre == 0 ) call topcntl( 0 )
! send edge data to TOPICS
        call topcntl( 2 )
        if( abs( time - tend ) < dtminh ) ksre = 1
! set the data received from TOPICS to variables
        if( ksre == 0 ) then
          call topcntl( 1 )
        else
          goto 200
        endif
      endif
      if( kcon /= 0 ) lstop = 1
      if( lstop == 1 ) goto 200

!::output
      dtimb = dtim
      dtim_s = dtnxt
      if( lgtpc == 1 .and. dtim_s > dtduc ) then
        if( dtduc > 0.0_8 ) then
          dtim   = dtduc
        else
          dtim   = min( dtim_s, dtcal )
        endif
      else
        dtim  = dtim_s  !   dtim of next loop
      endif
      tfac = dtim/dtimb
      fclp = elpmx/dlpq
      fcdt = edtmx/edtq
      write(n6,'(a,i1,2x,i6,2x,i6,2x,i6,"  dtim =",1pe10.3,"  _bf =",
     >   1pe11.3,"  tfac =",0pf6.3,"  fclp,fcdt =",1p2e10.2," cpmm = "
     >   ,0pf12.1,"m")')"<dtcntl-so> level ", lcstg, itim, nclnt, itmnt,
     >   dtim, dtimb, tfac,fclp,fcdt,cpmm

!::time
      if( lstop.eq.1 ) goto 200
      if( mod(itim,100).eq.0 ) call flush(n6)
      if( itim > itend  ) goto 200

!::loop
      kcnd = knorm
      msts = kcnd
      return

!::itim > itend  execute prg_term
 200  continue

!     heat flux from SOLDOR check plot  2018/05/15 Y.Homma
!!!      call yh_sld_q_para_plot
!     Bp/B pitch check plot @ it=24 (Sepx)  2018/05/16 Y.Homma
!!!      call yh_sld_test_hpit(24)
!     Metrics
!!!      call yh_sld_metric_plot
!     Plasma parameters group 1 check plot 2018/05/22 Y.Homma
!!!      call output_gridpoints
!******************************************************************
!     For final smoothing phase, call this at soldor/dtypset.f and
!     comment out from here. To match plasma between soldor and IMPMC.
!      call yh_sld_plasma_parameters_1_plot
!******************************************************************
!     Kappa interpolation test
!      call yh_sld_kappa_interpolation_test
!     Kappa comparison test
!!!      call kappa_para_write

!!!     2018/12/02 FOR POST-TREATMENT use
!!!     Relative change evaluation
!!!      call relative_change_param_4

!******************************************************************
!     Heat flux along a magnetic flux tube  "it_flux_1"
!      call heat_flux_along_flux_tube(24)
!
!     Heat flux and phys. parameters radial prof. data at "jcel_input"
      if(l_heat_flux_output == 1 )then
      call heat_flux_radial_output(1, 0, 1) ! Odiv
      call heat_flux_radial_output(1, 0, 40) ! OX
      call heat_flux_radial_output(1, 0, 50) ! OX-Omid
      call heat_flux_radial_output(1, 0, 65) ! Omid
      call heat_flux_radial_output(1, 0, 84) ! Imid
      call heat_flux_radial_output(1, 0, 100) ! IX-Imid
      call heat_flux_radial_output(1, 0, 109) ! IX
      call heat_flux_radial_output(1, 0, 148) ! Idiv

!     Heat flux and phys. parameters all cells
      call heat_flux_radial_output(0, 1, 1) ! Write_all, all cells.
      endif
!******************************************************************

      kcnd = klast
      if( lstop /= 0 ) kcnd = kstop
      msts = kcnd
      return

!--------------------------------------------------------------------
      case default
!--------------------------------------------------------------------
      write(cmsg,'("kitr =",i2," > 2")') kitr
      call wexit("prg_init(soldor)",trim(cmsg))
      end select
!
!-----------------------------------------------------------------------
!::error
!-----------------------------------------------------------------------
 290  continue
      write(n6,'("invalid time step control  dtim,elpmx,edtmx =",
     >  1p3e12.3)') dtim,elpmx,edtmx
      call wexit( "aasonic", "invalid time step control" )
      end
