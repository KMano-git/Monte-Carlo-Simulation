! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cgdcom
      subroutine cgdcom_alc( kk )
      use cgdcom, only : grdx, grdy, hbr, hbt, hbz, psman, psprv, pssol
      use csize,  only : ndq => ndx, ndp => ndy
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
        if( .not. allocated( grdx) ) then
          allocate( grdx(ndq,ndp), grdy(ndq,ndp), hbr(ndq,ndp)
     >      , hbt(ndq,ndp), hbz(ndq,ndp), psman(ndp), psprv(ndp)
     >      , pssol(ndp)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'grdx allocate error in cgdcom_alc, istat = ', istat
            call wexit( 'cgdcom', trim( cmsg ) )
          endif

! initial set
          grdx(1:ndq,1:ndp) = 0.0_8
          grdy(1:ndq,1:ndp) = 0.0_8
          hbr(1:ndq,1:ndp) = 0.0_8
          hbt(1:ndq,1:ndp) = 0.0_8
          hbz(1:ndq,1:ndp) = 0.0_8
          psman(1:ndp) = 0.0_8
          psprv(1:ndp) = 0.0_8
          pssol(1:ndp) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( grdx ) ) then
          deallocate( grdx, grdy, hbr, hbt, hbz, psman, psprv, pssol
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'grdx deallocate error in cgdcom_alc, istat = ', istat
            call wexit( 'cgdcom', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
