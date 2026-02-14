!***********************************************************************
      subroutine gdhedr
!***********************************************************************
      use com_gdhedr, only : hd_fno, hd_day, hd_dir
      implicit none
!
!::local variables
      character  wh_day*80, wh_tim*80, wh_lod*80, wh_dir*80, wh_job*80
      character  cno*5
      integer    mj, lenx
!
      write(6,'(/2x,"*** gdhedr ***")')
!
      call getenv( "WHAT_DAY", wh_day )
      call getenv( "WHAT_TIM", wh_tim )
      call getenv( "WHAT_LOD", wh_lod )
      call getenv( "WHAT_DIR", wh_dir )
      call getenv( "WHAT_JOB", wh_job )
!
      write(6,'(2x,"wh_day =",a)') wh_day
      write(6,'(2x,"wh_tim =",a)') wh_tim
      write(6,'(2x,"wh_lod =",a)') wh_lod
      write(6,'(2x,"wh_dir =",a)') wh_dir
      write(6,'(2x,"wh_job =",a)') wh_job
!
!::fig No
!x    call getenv( "WHAT_FNO", cno )
!x    read( cno, * ) hd_fno
      hd_fno = 0
!
!::date
      hd_day = wh_day(1:lenx(wh_day))//"  "//wh_tim(1:lenx(wh_tim))
!
!::directory
      mj = index( wh_dir, "/exe" )
      hd_dir = wh_dir(mj+4:)
!
      write(6,'(2x,"hd_day = ",a)')  hd_day(1:lenx(hd_day))
      write(6,'(2x,"hd_dir = ",a)')  hd_dir(1:lenx(hd_dir))
      write(6,'(2x,"hd_fno = ",i2)') hd_fno
!
      return
      end
