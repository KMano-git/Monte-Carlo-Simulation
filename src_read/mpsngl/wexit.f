c***********************************************************************
      subroutine wexit(csub,cmsg)
c***********************************************************************
      implicit none
c
      character csub*(*),cmsg*(*)
c
      character dsn*6, wpg*6
      integer n7, mype, lenx
c
      n7 = 7
      dsn = "erstop"
      wpg = "PROG"
      mype = 0
c
      open( unit=n7, file=dsn )
      write(n7,'(/2x,"*** wexit ***    PRG =",a,"  mype =",i4)')
     >   wpg,mype
      write(n7,'(5x,"stop at sub. ",a)') csub(1:lenx(csub))
      write(n7,'(5x,"because of   ",a)') cmsg(1:lenx(cmsg))
c
      write(6,'(/2x,"*** wexit ***    PRG =",a,"  mype =",i4)')
     >   wpg,mype
      write(6,'(5x,"stop at sub. ",a)') csub(1:lenx(csub))
      write(6,'(5x,"because of   ",a)') cmsg(1:lenx(cmsg))
c
      stop  ! Abnormal end (wexit)
      end
