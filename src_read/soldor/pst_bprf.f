!***********************************************************************
      subroutine pst_bprf
!***********************************************************************
      use csonic, only : lpst_bprf
      implicit none
!
!::local variables
      integer  i, it
!
!
!xx   call plprof(itsle)
!xx   call plprof(itpvs)
!
!::KH111108  define tube number by lpst_bprf
      do i = 1, 20
        it = lpst_bprf(i)
        if( it.le.0 ) exit
        call plprof(it)
      enddo
!
      return
      end
