! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cimp_loc
      subroutine cimp_loc_alc( kk )
      use cimp_loc, only : dlbp, dldf, dlpt
      use csize,  only : ndmc
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
        if( .not. allocated( dlbp ) ) then
          allocate( dlbp(ndmc), dldf(ndmc), dlpt(ndmc)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'dlbp allocate error in cimp_loc_alc, istat = ', istat
            call wexit( 'cimp_loc', trim( cmsg ) )
          endif

! initial set
          dlbp(1:ndmc) = 0.0_8 
          dldf(1:ndmc) = 0.0_8 
          dlpt(1:ndmc) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( dlbp ) ) then
          deallocate( dlbp, dldf, dlpt
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'dlbp deallocate error in cimp_loc_alc, istat = ', istat
            call wexit( 'cimp_loc', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
