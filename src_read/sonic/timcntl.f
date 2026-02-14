! added periodic invocation of IMPMC/figdat by kamata 2023/08/04
! determine if the current time is the time to fire an event
      subroutine timcntl( timtbl, dtmtbl, ntim, time, timb, judgement )
      implicit none
! arguments
      integer, intent(in)  :: ntim
      integer, intent(out) :: judgement
      real(8), intent(in)  :: dtmtbl(ntim), timb, time, timtbl(ntim)
! dtmtbl    : output time step size ( dtmtbl(i) between timtbl(i) and timtbl(i+1) )
! judgement : judgment result ( = yes, = 0 : no )
! ntim      : number of time table
! timb      : last output time
! time      : current time
! timtbl    : time to change the output time step size

! local variables
      integer    i, it
      real(8) :: dtminh = 5.0d-11, ztim

! get target index number of timtbl
      it = ntim
      do i = 1, ntim
        if( timtbl(i) > time + dtminh ) then
          it = i - 1
          exit
        endif
      enddo
 
      judgement = 0
      if( it > 0 ) then
! set next check time ztim
        ztim = timtbl(it)
        do
          if( ztim > timb + dtminh ) exit
          ztim = ztim + dtmtbl(it)
! modified 1/2 lines bug that occurs when it = ntim by kamata 2023/09/01
!ik       if( it < ntim .and. ztim + dtminh > timtbl(it+1) ) then
          if( it < ntim ) then
            if( ztim + dtminh > timtbl(it+1) ) then
              ztim = timtbl(it+1)
              it = it + 1
            endif
! added 1 line bug that occurs when it = ntim by kamata 2023/09/01
          endif
        enddo

! judgement
        if( time + dtminh >= ztim ) judgement = 1
      endif

      return
      end
