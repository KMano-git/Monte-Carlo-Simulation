! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cplvpn
      subroutine cplvpn_alc( kk )
      use cplvpn, only : qfy_vp, vdvp, vlvp, vpn_man, vpn_sol
      use csize,  only : nqfm => ndeq, nqfx => ndx, nqfy => ndy
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
        if( .not. allocated( qfy_vp ) ) then
          allocate( qfy_vp(nqfx,nqfy,nqfm), vdvp(nqfx,nqfy), vlvp(nqfy)
     >      , vpn_man(nqfy), vpn_sol(nqfy)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'qfy_vp allocate error in cplvpn_alc, istat = ', istat
            call wexit( 'cplvpn', trim( cmsg ) )
          endif

! initial set
          qfy_vp(1:nqfx,1:nqfy,1:nqfm) = 0.0_8
          vdvp(1:nqfx,1:nqfy) = 0.0_8
          vlvp(1:nqfy) = 0.0_8
          vpn_man(1:nqfy) = 0.0_8
          vpn_sol(1:nqfy) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( qfy_vp ) ) then
          deallocate( qfy_vp, vdvp, vlvp, vpn_man, vpn_sol
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'qfy_vp deallocate error in cplvpn_alc, istat = ', istat
            call wexit( 'cplvpn', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
