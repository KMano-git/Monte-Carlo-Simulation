! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cntpls
      subroutine cntpls_alc( kk )
      use cntpls, only : dene, deni, teme, temi, vflw, zefm
      use csize,  only : ndgs, ndmc
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
        if( .not. allocated( dene ) ) then
          allocate( dene(0:ndmc), deni(0:ndmc,ndgs), teme(0:ndmc)
     >      , temi(0:ndmc), vflw(0:ndmc,ndgs), zefm(0:ndmc)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'dene allocate error in cntpls_alc, istat = ', istat
            call wexit( 'cntpls', trim( cmsg ) )
          endif

! initial set
          dene(0:ndmc) = 0.0_8
          deni(0:ndmc,1:ndgs) = 0.0_8
          teme(0:ndmc) = 0.0_8 
          temi(0:ndmc) = 0.0_8
          vflw(0:ndmc,1:ndgs) = 0.0_8
          zefm(0:ndmc) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( dene ) ) then
          deallocate( dene, deni, teme, temi, vflw, zefm
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'dene deallocate error in cntpls_alc, istat = ', istat
            call wexit( 'cntpls', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
