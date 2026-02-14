! added treat 4 or more impurities with IMPMC by kamata 2022/04/21
! allocate variables in czwflx
      subroutine czwflx_alc( kk )
      use csize,  only : ndmis, ndwp, nzsmx
      use czwflx, only : zw_catmz, zw_eflx, zw_ismax, zw_pflx
      implicit none
! arguments
      integer, intent(in) :: kk
! kk : flag = 1: allocate,    = 2: deallocate

! local variables
      integer    istat
      character  cmsg*80

      select case( kk )
      case( 1 )
! allocate
        if( .not. allocated( zw_catmz ) ) then
          allocate( zw_catmz(nzsmx), zw_eflx(0:ndmis,ndwp,nzsmx)
     >        , zw_ismax(nzsmx), zw_pflx(0:ndmis,ndwp,nzsmx)
     >        , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'zw_catmz allocate error in czwflx_alc, istat = ', istat
            call wexit( 'czwflx', trim( cmsg ) )
          endif

! initial set
          zw_catmz(1:nzsmx) = ' '
          zw_eflx(0:ndmis,1:ndwp,1:nzsmx) = 0.0_8
          zw_ismax(1:nzsmx) = 0
          zw_pflx(0:ndmis,1:ndwp,1:nzsmx) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( zw_catmz ) ) then
          deallocate( zw_catmz, zw_eflx, zw_ismax, zw_pflx
     >        , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'zw_catmz deallocate error in czwflx_alc, istat = ', istat
            call wexit( 'czwflx', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
