! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cimpuf
      subroutine cimpuf_alc( kk )
      use cimpuf, only : bk_are, bk_deg, bk_flx, bk_iiw, bk_ity, pf_are
     >    , pf_deg, pf_flx, pf_iiw, pf_ity
      use csize,  only : ndwp
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
        if( .not. allocated( bk_are ) ) then
          allocate( bk_are(ndwp), bk_deg(ndwp), bk_flx(ndwp)
     >      , bk_iiw(ndwp), bk_ity(ndwp), pf_are(ndwp), pf_deg(ndwp)
     >      , pf_flx(ndwp), pf_iiw(ndwp), pf_ity(ndwp)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'bk_are allocate error in cimpuf_alc, istat = ', istat
            call wexit( 'cimpuf', trim( cmsg ) )
          endif

! initial set
          bk_are(1:ndwp) = 0.0_8
          bk_deg(1:ndwp) = 0.0_8
          bk_flx(1:ndwp) = 0.0_8 
          bk_iiw(1:ndwp) = 0.0_8
          bk_ity(1:ndwp) = 0.0_8
          pf_are(1:ndwp) = 0.0_8
          pf_deg(1:ndwp) = 0.0_8
          pf_flx(1:ndwp) = 0.0_8
          pf_iiw(1:ndwp) = 0.0_8
          pf_ity(1:ndwp) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( bk_are ) ) then
          deallocate( bk_are, bk_deg, bk_flx, bk_iiw, bk_ity, pf_are
     >      , pf_deg, pf_flx, pf_iiw, pf_ity
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'bk_are deallocate error in cimpuf_alc, istat = ', istat
            call wexit( 'cimpuf', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
