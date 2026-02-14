!**********************************************************************
      subroutine gdtrm(cspc)
!**********************************************************************
!
!      Use X-term
!
!      cspc = "X"  Xterm   "V"  versaterm
!             "L"  landscape  "P"  portrate  
!
!      call trmset(2)    
!      call trmput(xxxx) ! X trmput
!
!      call gdmod(xx) 
!      call trmput(xxxx)
!
!----------------------------------------------------------------------
      implicit none
      include 'clattr'
!
      character  cspc*(*)
!
!::local variables
      character  ctrm*1, cinp*20, cmsg*80
      integer    lenx, mji, kt, kp
!
      ctrm = "X"
!
!::X/V-term
!      ctrm = " "
!      call getenv("GDDEV",cinp)
!      if( lenx(cinp).gt.0 ) then
!      if( index(cinp,"X").gt.0 ) ctrm = "X"
!      if( index(cinp,"V").gt.0 ) ctrm = "V"
!      endif
!      if( lenx(ctrm).le.0 ) then
!      if( index(cspc,"X").gt.0 ) ctrm = "X"
!      if( index(cspc,"V").gt.0 ) ctrm = "V"
!      endif
!      write(6,'(2x,a)') "200"
!
!      if( lenx(ctrm).le.0 ) then
!      cmsg = ">> Enter device (1:VersaTerm, 2:X-lib) ==> "
!      mji = lenx(cmsg)
!      call trmput( 0, cmsg, mji )
!      call trmget( 0, cinp, 80 )
!      ctrm = "X"
!      if( cinp(1:1).eq."1" ) ctrm = "V"
!      endif
!      write(6,'(2x,a)') "300"
!
!      if( lenx(ctrm).le.0 ) then
!      cmsg = "invalid input for X/V-terminal"
!      mji  = lenx(cmsg)
!      call trmput( 0, cmsg, mji )
!      stop
!      endif
!
      kt = 1
      if( ctrm.eq."X" ) kt = 11
!
!::paper
!      cmsg = ">> Enter paper type (1:landscape/2:portrait) ==> "
!      mji = lenx(cmsg)
!      call trmset( 2 )     ! <== when use trmput before gddev
!      call trmput( 0, cmsg, mji )
!      call trmget( 0, cinp, 80 )
!
!      kp = 0
!      if( cinp(1:1).eq."2" ) kp = 1
!------
      kp =0

!
      call gdmod(kt+kp)
!
!::lsymfg  No define
      write(6,'(2x,"after gdmod  ldev, ltdev, lpaper")')
      write(6,'(2x,6i6)')        ldev, ltdev, lpaper
!
      return
      end


