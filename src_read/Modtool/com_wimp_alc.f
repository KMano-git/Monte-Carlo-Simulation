! allocate wimpd (ndmc is not parameter) 23/1/25
      subroutine com_wimp_alc
      use csize, only : ndmc
      use com_wimp, only : wimpd,wimpm,wimps
      implicit none
      
      integer    istat
      character  cmsg*80
  
      if( .not. allocated( wimpd ) ) then
        allocate(wimpd(ndmc),wimpm(ndmc),wimps(ndmc),
     >   stat = istat)
        if( istat /= 0 ) then
          write(cmsg,'(a,i6)')
     >     'wimpd allocate error in com_wimp_alc
     >     , istat = ', istat
          call wexit( 'com_wimp_alc', trim( cmsg ) )
        endif
        wimpd(1:ndmc) = 0.0_8
        wimpm(1:ndmc) = 0.0_8
        wimps(1:ndmc) = 0.0_8
      endif

      return
      end subroutine