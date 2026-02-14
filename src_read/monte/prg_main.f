!**********************************************************************
      program  neut2d
!**********************************************************************
      use cunit,       only : cprg
      use mod_dtflw,   only : mlp, mlpmax, mstp
      use mod_keylist, only : knend, kstop
      use mod_loc_tim, only : iend, itimm, kcal, msts, ncal
      use mod_mpicomm, only : m6
      use mod_shexe,   only : set_shexe_variables,npe_log
      use BGProf, only: allocate_ntbgprof
      implicit none

      integer :: kstp, klp, kcnd
      character,save :: clev*20 = "#= neut2d"

      cprg = 'neut2d'

!::work flow
      call set_shexe_variables
      call wf_defcomm( 1, cprg )
      call wf_deffile(npe_log)
      call wf_dtypsys(kcnd)

      call wf_spcdflw(kcnd)
      call wf_infdflw(kcnd)
      call wf_dtypdef(kcnd)
      call spmd_init
      call stjob
      call setmodv( cprg )
      call rstfile(0)
      call allocate_ntbgprof

!::MPMD sonic
      klp = 0
      mainLoop: do
        klp = klp + 1
        write(m6,'(/80("-"))')
        write(m6,'(a,2x,"klp =",i6,"  itimM =",i6)')
     >    trim(clev),klp,itimM

        call wf_maindrv(mlp,kstp,kcnd)

        write(m6,'(a,2x,"mlp/mlpmax =",2i11,"  mstp/msts =",2i2,
     > "  itimM/iend =",2i11,"  kcal/ncal =",2i6)')
     >  trim(clev), mlp, mlpmax, mstp, msts, itimM, iend, kcal, ncal

        if( msts == knend ) exit
        if( msts == kstop ) exit
        if( mlp > mlpmax )  exit
      enddo mainLoop

      write(m6,'(a,2x,a,2x,"klp =",i6)') trim(clev), "END", mlp

!::last message
      write(m6,'(a,2x,a,2x,"mlp/mlpmax =",2i11,"  msts =",i3)')
     >   trim(clev), "END", mlp, mlpmax, msts
      if( msts == kstop ) then
        call wf_mpiterm("Abnormal end prg. ",trim(cprg))
      else
        call wf_mpiterm("Normal end prg. ",trim(cprg))
      endif

      stop
      end program neut2d
