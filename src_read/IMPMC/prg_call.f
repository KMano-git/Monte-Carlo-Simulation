!**********************************************************************
      subroutine prg_call(kstp)
!**********************************************************************
!)
!)    prg_call : program call   (IMPMC)
!)
!)    prg_parm --- iminpt --- read uiminp (nimp, dtimp)
!)             ==> set_calcnd  (mcal,dcal) in mod_loc_time
!)
!)    see  IMPV5/test_jdgcal.f
!)
!)::mod_loc_tim
!)    time      tcal      tnxt     itim ical inxt jc kc nc
!)JDG 2.587_775 2.587_285 2.587_285  49    0  50   0  0  1
!)JDG 2.587_785 2.587_785 2.587_785  50   50 100   2  1  2
!)JDG 2.587_795 2.587_785 2.587_785  51   50 100   0  0  2
!)
!)::/cimcom_9d/
!)   tmsim = time in sonic  25.87, 30.87, 35.87 msec  time
!)   tmmcz = time in IMPMC  0.0, 5.0, 10.0, ... msec  tnow
!)   tmref = time in IMPMC  25.87 msec
!)
!)   tmsim = tmsim + dtm   tmsim = tmref in prg_prof
!)   check tmsim = time in mod_loc_tim
!)
!----------------------------------------------------------------------
      use cimctl,      only : dtimp, nimp
      use mod_keylist, only : kexec, kinit, kparm, kprof, kterm
      use mod_loc_tim, only : dcal, dtadv, dtimm, ical, iend, inxt
     >    , itimm, kcal, mcal, ncal, tcal, tcalb, tendm, timem, timemb
     >    , tnxt
      use mod_mpicomm, only : m6, n6
      use mod_shexe,   only : impmc_model
      implicit none
!
      integer, intent(in) :: kstp
!
      integer :: jcal
      real(8) :: dtgol, plest !timeMB timeM one step before

      select case(kstp)

!::init
      case(kinit)
      timeM = 0.0d0
      timeMB = 0.0d0
      dtimM = 0.0d0
      tendM = 0.0d0
      itimM = 0

      iend = 0
      ncal = 0
      kcal = 1

!::parm
      case(kparm)
      kcal = 1

!::prof
      case(kprof)
      kcal = 1
      inxt = 0
      mcal = nimp
      dcal = dtimp
      write(m6,'(4x,"SET prg_call(IMPMC)  mcal =",i4,"  dcal =",
     >  1pe12.3)') mcal, dcal

!::exec
!::jcal = 0(no cal)
!::     = 1(time) IMPMC called when reaches time limit
!::     = 2(itim) IMPMC called every nimp steps
      case(kexec)
      dtgol = tnxt-timeM
      tcalB = tcal
      plest = timeM - timeMB
      if( impmc_model == 1 )
     >  write(n6,*) 'plest, timeM, timeMB', plest, timeM, timeMB
      call jdgcal(jcal,plest)
      dtadv = tcal - tcalB
      timeMB = timeM

      kcal = 0
      if( jcal > 0 ) kcal = 1
      ! SY n6 write debug
      if( impmc_model == 1 .and. kcal.eq.1.or.mod(itimM,100).eq.0) then
        write(n6,'(4x,"JDG prg_call(IMPMC) ",1p3e18.10,1p2e12.4,6i8)')
     >    timeM,tcalB,tnxt,dtadv,dtgol, itimM,ical,inxt,jcal,kcal,ncal
      endif
      write(m6,'(4x,"JDG prg_call(IMPMC) ",1p3e18.10,1p2e12.4,6i8)')
     >  timeM,tcalB,tnxt,dtadv,dtgol, itimM,ical,inxt,jcal,kcal,ncal

!::term
      case(kterm)
      kcal = 1

!::error
      case default
      call wf_exit("judgcal","bad kstp")

      end select
      return
      end subroutine prg_call

!**********************************************************************
      subroutine jdgcal(jcal,plest)
!**********************************************************************
      use mod_loc_tim, only : dcal, ical, iend, inxt, itimm, kcal, mcal
     >    , ncal, tcal, timem, tnxt
      use mod_shexe, only : impmc_model
      implicit none

!::argument
      integer, intent(out) :: jcal
      real(8), intent(in)  :: plest
!::local vriables
      real(8)    tcalBL
!
      jcal = 0
      tcalBL = 0.0_8 ! set value for avoid undifined error.

      if( tnxt < 0.0d0 ) tnxt = timeM
      if( inxt < 0 )     inxt = itimM
      if( impmc_model == 1 ) then
        tcalBL = tcal
        if( timeM + plest >= tnxt ) jcal = 1
      else
        if( itimM >= inxt ) jcal = 2
      endif
      if( itimM == iend ) jcal = 2

      if( jcal /= 0 ) then
        tcal = timeM
        ical = itimM
        tnxt = tcal + dcal
        inxt = ical + mcal
        if( impmc_model == 1 .and. tcalBL == 0.0_8 ) jcal = 0 ! Never calc. at the beginning
      endif

      kcal = jcal
      if( jcal > 0 ) kcal = 1
      if( kcal == 1 ) ncal = ncal + 1

      return
      end
