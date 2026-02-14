! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cimden
      subroutine cimden_alc( kk )
      use cimden, only : tdnz, tfrz, tionZ, tradiZ, tradliZ, tradrZ
     >    , trecZ, tthz, tvlz, twci, twne, twrd, tengz, tprz
      use csize,  only : ndmc, ndmis
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
        if( .not. allocated( tdnz ) ) then
          allocate( tdnz(0:ndmis,ndmc), tfrz(0:ndmis,ndmc)
     >      , tionZ(0:ndmis,ndmc), tradiZ(0:ndmis,ndmc)
     >      , tradliZ(0:ndmis,ndmc), tradrZ(0:ndmis,ndmc)
     >      , trecZ(0:ndmis,ndmc), tthz(0:ndmis,ndmc)
     >      , tvlz(0:ndmis,ndmc), twci(ndmc), twne(ndmc), twrd(ndmc)
     >      , tengz(0:ndmis,ndmc), tprz(0:ndmis,ndmc)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'tdnz allocate error in cimden_alc, istat = ', istat
            call wexit( 'cimden', trim( cmsg ) )
          endif

! initial set
          tdnz(0:ndmis,1:ndmc) = 0.0_8
          tfrz(0:ndmis,1:ndmc) = 0.0_8
          tionZ(0:ndmis,1:ndmc) = 0.0_8
          tradiZ(0:ndmis,1:ndmc) = 0.0_8
          tradliZ(0:ndmis,1:ndmc) = 0.0_8
          tradrZ(0:ndmis,1:ndmc) = 0.0_8
          trecZ(0:ndmis,1:ndmc) = 0.0_8
          tthz(0:ndmis,1:ndmc) = 0.0_8
          tvlz(0:ndmis,1:ndmc) = 0.0_8
          twci(1:ndmc) = 0.0_8
          twne(1:ndmc) = 0.0_8
          twrd(1:ndmc) = 0.0_8
          tengz(0:ndmis,1:ndmc) = 0.0_8
          tprz(0:ndmis,1:ndmc) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( tdnz ) ) then
          deallocate( tdnz, tfrz, tionZ, tradiZ, tradliZ, tradrZ, trecZ
     >      , tthz, tvlz, twci, twne, twrd, tengz, tprz
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'tdnz deallocate error in cimden_alc, istat = ', istat
            call wexit( 'cimden', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
