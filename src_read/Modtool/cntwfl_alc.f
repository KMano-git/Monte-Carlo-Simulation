! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cntwfl
      subroutine cntwfl_alc( kk )
      use cntwfl, only : dwntl, iywl, jxwl
     >    , xare, xdeg, xeema, xeemm, xees, xehta, xehtm, xfhta, xfhti
     >    , xfhtm
      use csize,  only : ndmsr, ndwp, ndxy
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
        if( .not. allocated( dwntl ) ) then
          allocate( dwntl(ndxy,4), iywl(ndwp), jxwl(ndwp)
     >      , xare(ndwp)
     >      , xdeg(ndwp), xeema(ndwp,0:ndmsr), xeemm(ndwp,0:ndmsr)
     >      , xees(ndwp,0:ndmsr), xehta(ndwp,0:ndmsr)
     >      , xehtm(ndwp,0:ndmsr), xfhta(ndwp,0:ndmsr)
     >      , xfhti(ndwp,0:ndmsr), xfhtm(ndwp,0:ndmsr)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'dwntl allocate error in cntwfl_alc, istat = ', istat
            call wexit( 'cntwfl', trim( cmsg ) )
          endif

! initial set
          dwntl(1:ndxy,1:4) = 0.0_8
          iywl(1:ndwp) = 0
          jxwl(1:ndwp) = 0
          xare(1:ndwp) = 0.0_8
          xdeg(1:ndwp) = 0.0_8
          xeema(1:ndwp,0:ndmsr) = 0.0_8
          xeemm(1:ndwp,0:ndmsr) = 0.0_8
          xees(1:ndwp,0:ndmsr) = 0.0_8
          xehta(1:ndwp,0:ndmsr) = 0.0_8
          xehtm(1:ndwp,0:ndmsr) = 0.0_8
          xfhta(1:ndwp,0:ndmsr) = 0.0_8
          xfhti(1:ndwp,0:ndmsr) = 0.0_8
          xfhtm(1:ndwp,0:ndmsr) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( dwntl ) ) then
          deallocate( dwntl, iywl, jxwl
     >      , xare, xdeg, xeema, xeemm, xees, xehta, xehtm, xfhta, xfhti
     >      , xfhtm
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'dwntl deallocate error in cntwfl_alc, istat = ', istat
            call wexit( 'cntwfl', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
