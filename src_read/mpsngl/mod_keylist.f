      module  mod_keylist

!:: step (kstp)
      integer, parameter :: kinit = 1, kparm = 2
      integer, parameter :: kprof = 3, kexec = 4
      integer, parameter :: kterm = 5, kexit = 6
!
      character(4) :: cstep(6)
     >  = (/"init","parm","prof","exec","term","EXIT"/)
!
!:: return code (kcnd)
      integer, parameter :: knorm =  0  ! init, param, exec
      integer, parameter :: klast =  1  ! exec => term
      integer, parameter :: knend =  2  ! normal end
      integer, parameter :: kstop = -1  ! abnormal end

      integer, parameter :: kjump = -1  ! rapid change, no cal
!
      end module mod_keylist
