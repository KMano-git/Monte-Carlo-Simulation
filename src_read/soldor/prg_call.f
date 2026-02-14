!**********************************************************************
      subroutine prg_call(kstp)
!**********************************************************************
!
!    prg_call : program call   (soldor)
!
!----------------------------------------------------------------------
      use mod_keylist, only : kexec, kinit, kparm, kprof, kterm
      use mod_loc_tim, only : dcal, dtadv, dtimm, ical, iend, inxt
     >    , itimm, kcal, mcal, ncal, tcal, tcalb, tendm, timem, tnxt
      use mod_mpicomm, only : m6
      implicit none
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: kstp
      integer, intent(in) :: kstp

      integer :: jcal
      real(8) :: dtgol

      select case(kstp)

!::init
      case(kinit)
      timeM = 0.0d0
      dtimM = 0.0d0
      tendM = 0.0d0
      itimM = 0
      iend = 0
      kcal = 0
      ncal = 0
      kcal = 1

!::parm
      case(kparm)
      kcal = 1

!::prof
      case(kprof)
      kcal = 1
      inxt = 0
!xxx  call set_calcnd(mcal,dcal) ! after prg_parm
      mcal = 1
      dcal = 1.0d20

!::exec
!::jcal = 0(no cal)  1(time) 2(itim)
      case(kexec)
      dtgol = tnxt - timeM
      tcalB = tcal
      call jdgcal(jcal)
      dtadv = tcal - tcalB

      kcal = 0
      if( jcal > 0 ) kcal = 1
      write(m6,'(4x,"JDG prg_call(soldor) ",1p3e18.10,1p2e12.4,6i8)')
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
      subroutine jdgcal(jcal)
!**********************************************************************
      use mod_loc_tim, only : dcal, ical, iend, inxt, itimm, kcal, mcal
     >    , ncal, tcal, timem, tnxt
      implicit none

!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: jcal
      integer, intent(out) :: jcal

!::local vriables
      real(8) :: eps = 1.0d-10

      jcal = 0

      if( tnxt < 0.0d0 ) tnxt = timeM
      if( inxt < 0 )     inxt = itimM

      if( timeM + eps >= tnxt ) jcal = 1
      if( itimM >= inxt ) jcal = 2
      if( itimM == iend ) jcal = 2

      if( jcal /= 0 ) then
        tcal = timeM
        ical = itimM
        tnxt = tcal + dcal
        inxt = ical + mcal
      endif

      kcal = jcal
      if( jcal > 0 ) kcal = 1
      if( kcal == 1 ) ncal = ncal + 1

      return
      end
