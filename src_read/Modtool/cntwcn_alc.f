! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cntwcn
      subroutine cntwcn_alc( kk )
      use cntwcn, only : wden, weema, weemm, wees, wehta, wehtm, wemt
     >    , weng, wfhta, wfhtm, wgden, wgeng, wnfl0x, wnfl0y, wnfl0z
     >    , wnflgx, wnflgy, wnflgz, wsbr, wssn, wssp, wswe, wswi, wtion
     >    , wvlp, wwal
      use csize,  only : ndgs, ndmc, ndwp
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
        if( .not. allocated( wden ) ) then
          allocate( wden(0:ndmc,ndgs), weema(ndwp), weemm(ndwp)
     >      , wees(ndwp), wehta(ndwp), wehtm(ndwp), wemt(ndwp)
     >      , weng(0:ndmc,ndgs), wfhta(ndwp), wfhtm(ndwp)
     >      , wgden(0:ndmc,2), wgeng(0:ndmc,2), wnfl0x(0:ndmc,ndgs)
     >      , wnfl0y(0:ndmc,ndgs), wnfl0z(0:ndmc,ndgs), wnflgx(0:ndmc,2)
     >      , wnflgy(0:ndmc,2), wnflgz(0:ndmc,2), wsbr(0:ndmc,ndgs)
     >      , wssn(0:ndmc,ndgs), wssp(0:ndmc,ndgs), wswe(0:ndmc)
     >      , wswi(0:ndmc), wtion(0:ndmc,ndgs), wvlp(0:ndmc,ndgs)
     >      , wwal(ndwp)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'wden allocate error in cntwcn_alc, istat = ', istat
            call wexit( 'cntwcn', trim( cmsg ) )
          endif

! initial set
          wden(0:ndmc,1:ndgs) = 0.0_8 
          weema(1:ndwp) = 0.0_8 
          weemm(1:ndwp) = 0.0_8
          wees(1:ndwp) = 0.0_8
          wehta(1:ndwp) = 0.0_8 
          wehtm(1:ndwp) = 0.0_8
          wemt(1:ndwp) = 0.0_8 
          weng(0:ndmc,1:ndgs) = 0.0_8 
          wfhta(1:ndwp) = 0.0_8 
          wfhtm(1:ndwp) = 0.0_8 
          wgden(0:ndmc,1:2) = 0.0_8 
          wgeng(0:ndmc,1:2) = 0.0_8
          wnfl0x(0:ndmc,1:ndgs) = 0.0_8 
          wnfl0y(0:ndmc,1:ndgs) = 0.0_8 
          wnfl0z(0:ndmc,1:ndgs) = 0.0_8
          wnflgx(0:ndmc,1:2) = 0.0_8 
          wnflgy(0:ndmc,1:2) = 0.0_8 
          wnflgz(0:ndmc,1:2) = 0.0_8
          wsbr(0:ndmc,1:ndgs) = 0.0_8
          wssn(0:ndmc,1:ndgs) = 0.0_8 
          wssp(0:ndmc,1:ndgs) = 0.0_8
          wswe(0:ndmc) = 0.0_8 
          wswi(0:ndmc) = 0.0_8 
          wtion(0:ndmc,1:ndgs) = 0.0_8 
          wvlp(0:ndmc,1:ndgs) = 0.0_8
          wwal(1:ndwp) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( wden ) ) then
          deallocate( wden, weema, weemm, wees, wehta, wehtm, wemt, weng
     >      , wfhta, wfhtm, wgden, wgeng, wnfl0x, wnfl0y, wnfl0z, wnflgx
     >      , wnflgy, wnflgz, wsbr, wssn, wssp, wswe, wswi, wtion, wvlp
     >      , wwal
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'wden deallocate error in cntwcn_alc, istat = ', istat
            call wexit( 'cntwcn', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
