!***********************************************************************
      subroutine stjob
!***********************************************************************
      use cunit,     only : n6, wh_day, wh_dir, wh_job, wh_lod, wh_tim
      use mod_shexe, only : use_namelist, WHAT_DAY, WHAT_TIM, WHAT_LOD
     >    , WHAT_DIR, WHAT_JOB
      implicit none
! function
      integer    lenx
!
      if(use_namelist) then
         !.. use variables in inppls
         wh_day = WHAT_DAY
         wh_tim = WHAT_TIM
         wh_lod = WHAT_LOD
         wh_dir = WHAT_DIR
         wh_job = WHAT_JOB
      else
         !.. use environment variables
         call getenv( "WHAT_DAY", wh_day )
         call getenv( "WHAT_TIM", wh_tim )
         call getenv( "WHAT_LOD", wh_lod )
         call getenv( "WHAT_DIR", wh_dir )
         call getenv( "WHAT_JOB", wh_job )
      end if
!
      write(n6,'(/2x,"*** stjob ***")')
      write(n6,'(4x,"Date = ",a,4x,"Time = ",a)') wh_day,wh_tim
      write(n6,'(4x,"Dir  = ",a)') wh_dir(1:lenx(wh_dir))
      write(n6,'(4x,"Load = ",a)') wh_lod(1:lenx(wh_lod))
      write(n6,'(4x,"Job  = ",a)') wh_job
!
      return
      end
