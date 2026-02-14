!
!***********************************************************************
      subroutine delspc(clin)
!***********************************************************************
      implicit none
      character(*) :: clin
!
!::local variables
      character(240) :: ctmp
      integer :: imx, i, ii
!
      ctmp = " "
!
      imx = len(clin)
      ii = 0
      do i = 1, imx
      if( clin(i:i).eq." " ) cycle
      ii = ii + 1
      ctmp(ii:ii) = clin(i:i)
      enddo
!
      clin = trim(ctmp)
      return
      end
