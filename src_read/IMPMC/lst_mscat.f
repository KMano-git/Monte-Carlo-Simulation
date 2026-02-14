!*********************************************************************
      subroutine lst_mscat(tnow)
!*********************************************************************
!)
!)    lst_mscat : list for multi scatter plot
!)
!)-------------------------------------------------------------------
      use cimcom, only : i6_mscat
      use cunit,  only : lmspe, lmype, n6
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/07/04
!ik   real(8) :: tnow
      real(8), intent(in) :: tnow

! modified 1/1 lines organize local variables and include files by kamata 2021/07/04
!ik   character :: ctm*4, dsn*20
      character :: dsn*20
      integer :: ksw, i, i6
      integer :: ii = 0

      integer :: ntb = 12
      save :: ii, ntb

      real(8) :: tmout(12)
      tmout(1:12) = (/0.1d-3, 0.5d-3, 1.0d-3, 2.0d-3, 3.0d-3, 4.0d-3,
     >   5.0d-3, 10.0d-3, 20.0d-3, 30.0d-3, 40.0d-3, 50.0d-3/)

      if( lmype /= lmspe ) return
      if( i6_mscat == 0 ) return

      ksw = 0
      do i = 1, ntb
        if( dabs(tnow-tmout(i)) < 1.0d-6 ) then
          ksw = 1
          exit
        endif
      enddo
      write(n6,'(2x,"lst_mscat  tnow =",1pe12.4,"  ksw =",i3)')
      if( ksw == 0 ) return

      ii = ii + 1
      write(dsn,'(a,i2.2)') "lst_mscat_n",ii
      write(n6,'(2x,"*** lst_mscat ***  tnow =",1pe12.4,2x,a)')
     >   tnow, trim(dsn)

      open(unit=i6_mscat,file=dsn)
      i6 = i6_mscat
      call lst_scat(tnow,1,i6)
      close(i6)

      return
      end
