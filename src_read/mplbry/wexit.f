!***********************************************************************
      subroutine wexit(csub,cmsg)
!***********************************************************************
      use cunit,       only : grnm, mygrp, mype, n6, n7
      use mod_mpicomm, only : nwld_all

! For falcon, ifcore is needed for stack trace output
#ifdef FALCON_I
      use ifcore
#endif

      implicit none
!
!:: argument
      character, intent(in) :: csub*(*), cmsg*(*)
!
!:: local variable
      character dsn*80
      integer  ino, ierr
!
!:: function
      integer lenx
!
      call getenv( "ERROR", dsn )
      if( lenx(dsn).le.0 ) dsn = "erstop"
!
      if( n7.eq.0 ) then
        n7   = 99
        mygrp = 1
        grnm(mygrp) = "mpsta"
        if( n7.ne.n6 ) open( unit=n7, file=dsn, form="formatted" )
      endif
!
!::file error    KSFUJI
      write(n7,'(/2x,"*** wexit ***    GRP =",a,"  mype =",i4)')
     >   grnm(mygrp),mype
      write(n7,'(5x,"stop at sub. ",a)') csub(1:lenx(csub))
      write(n7,'(5x,"because of   ",a)') cmsg(1:lenx(cmsg))
!
!::file 6
      write(6,'(/2x,"*** wexit ***    GRP =",a,"  mype =",i4)')
     >   grnm(mygrp),mype
      write(6,'(5x,"stop at sub. ",a)') csub(1:lenx(csub))
      write(6,'(5x,"because of   ",a)') cmsg(1:lenx(cmsg))
!
!::file n6
      write(n6,'(/2x,"*** wexit ***    GRP =",a,"  mype =",i4)')
     >   grnm(mygrp),mype
      write(n6,'(5x,"stop at sub. ",a)') csub(1:lenx(csub))
      write(n6,'(5x,"because of   ",a)') cmsg(1:lenx(cmsg))
!
!::flush
      if( n6.gt.0 ) call flush(n6)
      if( n7.gt.0 ) call flush(n7)
      call flush(6)
!
      call sleep(10)
!
      ino = 10
!
! For falcon, stack trace does not write in MPI_Abort
#ifdef FALCON_I
      call tracebackqq(user_exit_code=-1)
#endif
!
      call MPI_Abort( nwld_all, ino, ierr )
      call mpiend
!
      stop  ! Abnormal end (wexit)
      end subroutine wexit