! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cplwrd
      subroutine cplwrd_alc( kk )
      use cplwrd, only : wcre, wime, wimi, wmce, wmci
     > , wcr_wrd, ndzty
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
        if( .not. allocated( wcre ) ) then
          allocate( wcre(ndmc), wime(ndmc), wimi(ndmc), wmce(ndmc)
     >      , wmci(ndmc)
     >      , wcr_wrd(ndmc,ndzty)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'wcre allocate error in cplwrd_alc, istat = ', istat
            call wexit( 'cplwrd', trim( cmsg ) )
          endif

! initial set
          wcre(1:ndmc) = 0.0_8
          wime(1:ndmc) = 0.0_8
          wimi(1:ndmc) = 0.0_8
          wmce(1:ndmc) = 0.0_8
          wmci(1:ndmc) = 0.0_8
          wcr_wrd(ndmc,ndzty) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( wcre ) ) then
          deallocate( wcre, wime, wimi, wmce, wmci
     >      , wcr_wrd
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'wcre deallocate error in cplwrd_alc, istat = ', istat
            call wexit( 'cplwrd', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
