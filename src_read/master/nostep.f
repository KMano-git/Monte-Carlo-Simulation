!**********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/06/22
!ik   function nostep(lp)
      integer function nostep(lp)
!**********************************************************************
!
!   mlp => kstep
!                  1: init  2: parm  3: prof  4: exec  5: term
!   msts (status)  0: norm  1: last  2: stop
!
!---------------------------------------------------------------------
      use mod_mpicomm, only : m6, n6
      use mod_keylist, only : cstep, kexec, klast, kterm
      use mod_loc_tim, only : msts
      implicit none
! modified 2/1 lines organize local variables and include files by kamata 2021/06/22
!ik   integer :: lp
!ik   integer :: nostep
      integer, intent(in) :: lp
      integer :: istp

! deleted 1 line organize local variables and include files by kamata 2021/06/22
!ik   integer,save   :: nlev = 2
      character,save :: clev*20 = "#== nostep"

      write(m6,'(a,2x,a)') trim(clev), "STA"

      if( lp <= 0 ) call wf_exit("nostep","bad lp")
      if( lp <= 4 ) then
        istp = lp
      else
        istp = kexec
        if( msts == klast ) istp = kterm
      endif

      write(m6,'(a,2x,"  lp =",i6,"  istp =",i2,2x,a,
     >  "  msts =",2i2,"  =0(loop) =1(last) =2(error)")')
     >   trim(clev), lp, istp, cstep(istp), msts
      write(m6,'(a,2x,a)') trim(clev), "END"

      nostep = istp
      return
      end function nostep
