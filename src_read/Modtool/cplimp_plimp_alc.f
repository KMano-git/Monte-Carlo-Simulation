! added treat 4 or more impurities with IMPMC by kamata 2022/05/08
! allocate variables in cplimp
      subroutine cplimp_plimp_alc( kk )
      use csize,  only : ndmc, ndmis, nzsmx
      use cplimp_plimp, only : tvlzL,tfrzL,tthzL,tionZL,trecZL
     > , twrdL
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
        if( .not. allocated( tvlzL ) ) then
          allocate( tvlzL(0:ndmis,ndmc,nzsmx), tfrzL(0:ndmis,ndmc,nzsmx)
     >      , tthzL(0:ndmis,ndmc,nzsmx), tionZL(0:ndmis,ndmc,nzsmx)
     >      , trecZL(0:ndmis,ndmc,nzsmx), twrdL(ndmc,nzsmx)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'tvlzL allocate error in cplimp_plimp_alc,
     >         istat = ', istat
            call wexit( 'cplimp_plimp_alc', trim( cmsg ) )
          endif

! initial set
          tvlzL(0:ndmis,ndmc,nzsmx) = 0.0_8
          tfrzL(0:ndmis,ndmc,nzsmx) = 0.0_8
          tthzL(0:ndmis,ndmc,nzsmx) = 0.0_8
          tionZL(0:ndmis,ndmc,nzsmx) = 0.0_8
          trecZL(0:ndmis,ndmc,nzsmx) = 0.0_8
          twrdL(1:ndmc,1:nzsmx) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( tvlzL ) ) then
          deallocate( tvlzL,tfrzL,tthzL,tionZL,trecZL,twrdL
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'tvlzL deallocate error in cplimp_plimp_alc
     >         , istat = ', istat
            call wexit( 'cplimp_plimp_alc', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
