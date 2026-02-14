!**********************************************************************
      subroutine prg_call(kstp)
!**********************************************************************
!
!    prg_call : program call   (topics)
!
!----------------------------------------------------------------------
      use mod_keylist, only : kexec, kinit, kparm, kprof, kterm
      use mod_loc_tim, only : dcal, dtimm, ical, iend, inxt, itimm
     >   , kcal, mcal, ncal, tcal, tendm, timem, tnxt
      use mod_mpicomm, only : m6
      implicit none
      integer, intent(in) :: kstp

      select case(kstp)
!
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
        mcal = 1
        dcal = 1.0d20

!::exec
      case(kexec)
        kcal = 0
        if( itimM == inxt ) then
          kcal = 1
          tcal = timeM
          ical = itimM
          ncal = ncal + 1
          tnxt = tcal + dcal
          inxt = ical + mcal
          write(m6,'(2x,"prg_call master  itim/inxt/mcal/kcal =",
     >     i7,x,i7,x,i5,x,i5)') itimM, inxt, mcal, kcal
        endif

!::term
      case(kterm)
        kcal = 1

      case default
        call wf_exit("judgcal","bad kstp")

      end select

      return
      end subroutine prg_call
