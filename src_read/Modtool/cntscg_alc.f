! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cntscg
      subroutine cntscg_alc( kk )
      use cntscg, only : sn_mc, sn_mt, sp_mc, sp_mt, we_mc, we_mt, wi_mc
     >    , wi_mt
      use csize,  only : ndmc, ndsp
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
        if( .not. allocated( sn_mc ) ) then
          allocate( sn_mc(0:ndmc,ndsp), sn_mt(0:ndmc,ndsp)
     >      , sp_mc(0:ndmc,ndsp), sp_mt(0:ndmc,ndsp), we_mc(0:ndmc)
     >      , we_mt(0:ndmc), wi_mc(0:ndmc), wi_mt(0:ndmc)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'sn_mc allocate error in cntscg_alc, istat = ', istat
            call wexit( 'cntscg', trim( cmsg ) )
          endif

! initial set
          sn_mc(0:ndmc,1:ndsp) = 0.0_8
          sn_mt(0:ndmc,1:ndsp) = 0.0_8 
          sp_mc(0:ndmc,1:ndsp) = 0.0_8
          sp_mt(0:ndmc,1:ndsp) = 0.0_8
          we_mc(0:ndmc) = 0.0_8
          we_mt(0:ndmc) = 0.0_8
          wi_mc(0:ndmc) = 0.0_8
          wi_mt(0:ndmc) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( sn_mc ) ) then
          deallocate( sn_mc, sn_mt, sp_mc, sp_mt, we_mc, we_mt, wi_mc
     >      , wi_mt
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'sn_mc deallocate error in cntscg_alc, istat = ', istat
            call wexit( 'cntscg', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
