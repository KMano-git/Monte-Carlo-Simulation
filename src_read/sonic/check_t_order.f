!**********************************************************************
!:: Ensure that the array elements in timearray are set in order from smallest to largest value.
      subroutine check_t_order(timearray,timesize,timesize_eff)
!**********************************************************************
      implicit none
!::arguments
      integer, intent(in) :: timesize
      real(8), intent(in) :: timearray(timesize)
      integer, intent(out) :: timesize_eff
!::local variables
      integer i
!::
      do i = 1, timesize
        ! check the effective data size of timearray
        ! We assume (1.zero-clear for timearray is done),  (2.user never set as timearray(i) < 0)
        if(i==1) cycle
        if(timearray(i).eq.0.0d0) then
          if(i < 3) then
            ! data number of timearray must be more then 1.
            call wexit("check_t_order","len(timearray) < 2")
          endif
          timesize_eff = i-1
          exit
        endif

        if(timearray(i-1) .ge. timearray(i)) then
          call wexit("check_t_order","Set the timearray in order
     > from smallest to largest.")
        endif
      enddo

      end subroutine
