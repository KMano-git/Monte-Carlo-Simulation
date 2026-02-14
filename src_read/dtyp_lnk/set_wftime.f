!**********************************************************************
      subroutine set_wftime
!**********************************************************************
!)
!)    MST : master
!)    SOL : soldor   time control 
!)    NTL : neut2d
!)    IMP : IMPMC
!)
!)    clr_wftime  [mod_loc_tim], [mod_glb_tim]
!)
!)    set_wftime
!)       SOL [inc-sonic]  => SOL [mod_loc_tim]
!)
!)    wf_tmcontl
!)       SOL/IMP  kcal         MPIs  ==>  MST  zkcal, zncal
!)       SOL [mod_loc_tim]     SET   ==>  SOL  [mod_glb_tim]
!)       SOL [mod_glb_tim]     MPIs  ==>  MST  [mod_glb_tim]
!)       MST [mod_glb_tim]     Bcast ==>  SOL/IMP [mod_glb_tim]
!)       SOL/IMP [mod_glb_tim] SET   ==>  SOL/IMP [mod_loc_tim]
!)                              
!)      use [mod_loc_tim] in execution
!)
!)    cstime      : time, tend, dtim, itim, iend
!)    mod_loc_tim : timeM, tendM, dtimM, itimM, itend, mstop, lstop
!)
!)     see WFL/set_wftime
!)
!-----------------------------------------------------------------------
      use csonic,      only : dtim, itend, itim, tend, time
      use mod_loc_tim, only : dcal, dtimm, ical, iend, inxt, itimm, kcal
     >    , mcal, ncal, tcal, tendm, timem, tnxt
      use mod_mpicomm, only : m6
      implicit none

!::mod_loc_tim <= SONIC time
      timeM = time
      tendM = tend
      dtimM = dtim
      itimM = itim
      iend  = itend
      kcal  = 1
      ncal  = 0             ! number of calculations  ncal = ncal + 1

      tcal  = time
      dcal  = dtim
      tnxt  = tcal + dcal
      ical  = itim
      mcal  = 1             ! every mcal times
      inxt  = ical          ! Note error inxt = ical + mcal

      write(m6,'("%%07/27 set_wftime timeM,tendM,dtimM =",1p3e12.4,
     >  "  itimM,iend,kcal,ncal =",2i11,2i3)')
     >  timeM, tendM, dtimM, itimM, iend, kcal, ncal
      write(m6,'("                   tcal, dcal, tnxt  =",1p3e12.4,
     >  "  ical, mcal, inxt    =",i5,i11,i5)')
     >  tcal, dcal, tnxt, ical, mcal, inxt

      return
      end

!**********************************************************************
      subroutine clr_wftime
!**********************************************************************
      use mod_loc_tim, only : dcal, dtimm, ical, iend, inxt, itimm, kcal
     >    , mcal, ncal, tcal, tendm, timem, tnxt
      use mod_glb_tim, only : g_dtim, g_iend, g_itim, g_kcal, g_ncal
     >    , g_tend, g_time, ghd_tim, git_tim, gtm_tim, ndgrp
      implicit none

!::mod_loc_tim
      timeM = 0.0d0
      tendM = 0.0d0
      dtimM = 0.0d0
      itimM = 0
      iend = 0 
      kcal = 0
      ncal = 0

      tcal = 0.0d0
      dcal = 0.0d0
      tnxt = 0.0d0
      ical = 0
      mcal = 0
      inxt = 0
      
!::mod_glb_tim
      Ghd_tim = 0
      Git_tim = 0
      Gtm_tim = 0

      G_time = 0.0d0
      G_tend = 0.0d0
      G_dtim = 0.0d0
      G_itim = 0
      G_iend = 0
      G_kcal(0:ndgrp) = 0
      G_ncal(0:ndgrp) = 0

      return
      end
