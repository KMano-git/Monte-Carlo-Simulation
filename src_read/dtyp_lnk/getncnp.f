! added treat 4 or more impurities with IMPMC by kamata 2022/04/21
! get the numbers before and after '_' after the 4th character.
! get code number(nc) and processing number(np) from ctyp('IMP(nc)_(np)')
      subroutine getncnp( ctyp, nc, np )
      implicit none
! arguments
      character(*), intent(in)  :: ctyp
      character(2), intent(out) :: np
      integer,   intent(out) :: nc
! ctyp : data label
! nc   : code number       ( = 0   : couldn't take it out )
! np   : processing number ( = ' ' : couldn't take it out )

! local variables
      integer    ios, k, mj
! ios : status   
! k   : position
! mj  : number of characters

! get processing number
      np = ' '
      mj = len_trim( ctyp )
      k  = index( ctyp, '_' )
      if( k < mj ) np = ctyp(k+1:mj)

! get code number
      nc = 0
      if( k == 4 ) then
        nc = 1
      elseif( k > 4 ) then
        read(ctyp(4:k-1),*,iostat=ios) nc
        if( ios /= 0 ) nc = 0
      endif

      return
      end
