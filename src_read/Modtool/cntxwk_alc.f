! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cntxwk
      subroutine cntxwk_alc( kk )
      use cntxwk, only : ndhfl, xmwork
      use csize,  only : ndmsr, nwkmp_nt
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
        if( .not. allocated( xmwork ) ) then
          allocate( xmwork(nwkmp_nt,ndmsr,ndhfl)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'xmwork allocate error in cntxwk_alc, istat = ', istat
            call wexit( 'cntxwk', trim( cmsg ) )
          endif

! initial set
          xmwork(1:nwkmp_nt,1:ndmsr,1:ndhfl) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( xmwork ) ) then
          deallocate( xmwork
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'xmwork deallocate error in cntxwk_alc, istat = ', istat
            call wexit( 'cntxwk', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
