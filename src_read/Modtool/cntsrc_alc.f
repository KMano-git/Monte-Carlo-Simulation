! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cntsrc
      subroutine cntsrc_alc( kk )
      use cntsrc, only : tden0, tdeng, teng0, tengg, tssn, tssp, tswe
     >    , tswi, tvlp0
      use csize,  only : ndgs, ndmc, ndsp
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
        if( .not. allocated( tden0 ) ) then
          allocate( tden0(0:ndmc,ndgs), tdeng(0:ndmc,2)
     >      , teng0(0:ndmc,ndgs), tengg(0:ndmc,2), tssn(0:ndmc,ndsp)
     >      , tssp(0:ndmc,ndsp), tswe(0:ndmc), tswi(0:ndmc)
     >      , tvlp0(0:ndmc,ndgs)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'tden0 allocate error in cntsrc_alc, istat = ', istat
            call wexit( 'cntsrc', trim( cmsg ) )
          endif

! initial set
          tden0(0:ndmc,1:ndgs) = 0.0_8
          tdeng(0:ndmc,1:2) = 0.0_8 
          teng0(0:ndmc,1:ndgs) = 0.0_8
          tengg(0:ndmc,1:2) = 0.0_8
          tssn(0:ndmc,1:ndsp) = 0.0_8
          tssp(0:ndmc,1:ndsp) = 0.0_8
          tswe(0:ndmc) = 0.0_8
          tswi(0:ndmc) = 0.0_8
          tvlp0(0:ndmc,1:ndgs) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( tden0 ) ) then
          deallocate( tden0, tdeng, teng0, tengg, tssn, tssp, tswe, tswi
     >      , tvlp0
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'tden0 deallocate error in cntsrc_alc, istat = ', istat
            call wexit( 'cntsrc', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
