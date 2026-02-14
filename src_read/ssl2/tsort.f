!***********************************************************************
      subroutine tsort(nn,rel,ipnt)
!***********************************************************************
      implicit none
      integer, intent(in)    :: nn
      integer, intent(out)   :: ipnt(nn)
      real(8), intent(inout) :: rel(nn)
!
!::local variables
      integer  j, md, itp, i, l, isav
      real*8  rsav
!
      do j = 1, nn
      ipnt(j) = j
      enddo
!
      md = 1
      do  while ( md .le. nn )
      md = 2 * md
      enddo
      md = ( md - 1 ) / 2
!
      do while ( md .ne. 0 )
      itp = nn - md
      do i = 1 , itp
      j = i
      do while ( j. gt. 0 )
      l = j + md
      if( rel(l) .ge. rel(j) ) go to 40
      rsav    = rel (j)
      rel(j)  = rel (l)
      rel(l)  = rsav
      isav    = ipnt(j)
      ipnt(j) = ipnt(l)
      ipnt(l) = isav
      j = j - md
      enddo
   40 continue
      enddo
      md = ( md - 1 ) / 2
      enddo
!
      end