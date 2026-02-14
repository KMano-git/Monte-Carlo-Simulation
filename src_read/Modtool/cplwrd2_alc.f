! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cplwrd2
! nzsmx for SOLDOR
      subroutine cplwrd2_alc( kk )
      use cplwrd2, only : sradi, sradli, sradr, wmc_zwe, wmc_zwi
      use csize,   only : ndmc, ndmis, ndx, ndy, nzsmx
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
        if( .not. allocated( sradi ) ) then
          allocate( sradi(0:ndmis,ndmc,nzsmx)
     >      , sradli(0:ndmis,ndmc,nzsmx), sradr(0:ndmis,ndmc,nzsmx)
     >      , wmc_zwe(ndmc,nzsmx), wmc_zwi(ndmc,nzsmx)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'sradi allocate error in cplwrd2_alc, istat = ', istat
            call wexit( 'cplwrd2', trim( cmsg ) )
          endif

! initial set
          sradi(0:ndmis,1:ndmc,1:nzsmx) = 0.0_8
          sradli(0:ndmis,1:ndmc,1:nzsmx) = 0.0_8
          sradr(0:ndmis,1:ndmc,1:nzsmx) = 0.0_8
          wmc_zwi(1:ndmc,1:nzsmx) = 0.0_8
          wmc_zwe(1:ndmc,1:nzsmx) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( sradi ) ) then
          deallocate( sradi, sradli, sradr, wmc_zwi, wmc_zwe
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'sradi deallocate error in cplwrd2_alc, istat = ', istat
            call wexit( 'cplwrd2', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
