! allocate xwlp (ndmc is not parameter) 23/1/25
      subroutine com_pth_line_alc
      use csize, only : ndwp
      use com_pth_line, only : iwlp, xwlp, ywlp
      implicit none

      integer    istat
      character  cmsg*80

      if( .not. allocated( xwlp ) ) then
        allocate(xwlp(ndwp), ywlp(ndwp),iwlp(ndwp),
     >   stat = istat)
        if( istat /= 0 ) then
          write(cmsg,'(a,i6)')
     >     'xwlp allocate error in com_pth_line_alc,
     >      istat = ', istat
          call wexit( 'com_pth_line_alc', trim( cmsg ) )
        endif
        xwlp(1:ndwp) = 0.0_8
        ywlp(1:ndwp) = 0.0_8
        iwlp(1:ndwp) = 0
      endif

      return
      end subroutine