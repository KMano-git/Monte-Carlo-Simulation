! added treat 4 or more impurities with IMPMC by kamata 2022/05/08
! allocate variables in cplimp
      subroutine cplimp_alc( kk )
      use csize,  only : ndgt, ndmc, ndmis, ndwp, nzsmx
      use cplimp, only : aimasL, amiL, amzL, azmasL, cdifL, eipotL
     >    , ismaxL, nzmxL, tdnzL, twciL, tengzL, tprzL
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
        if( .not. allocated( aimasL ) ) then
          allocate( aimasL(nzsmx), amiL(nzsmx), amzL(nzsmx)
     >      , azmasL(nzsmx), cdifL(nzsmx), eipotL(0:ndmis,nzsmx)
     >      , ismaxL(nzsmx), nzmxL(nzsmx)
     >      , tdnzL(0:ndmis,ndmc,nzsmx), twciL(ndmc,nzsmx)
     >      , tengzL(0:ndmis,ndmc,nzsmx), tprzL(0:ndmis,ndmc,nzsmx)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'aimasL allocate error in cplimp_alc, istat = ', istat
            call wexit( 'cplimp', trim( cmsg ) )
          endif

! initial set
          aimasL(1:nzsmx)               = 0.0_8
          amiL(1:nzsmx)                 = 0.0_8
          amzL(1:nzsmx)                 = 0.0_8
          azmasL(1:nzsmx)               = 0.0_8
          cdifL(1:nzsmx)                = 0.0_8
          eipotL(0:ndmis,1:nzsmx)       = 0.0_8
          tdnzL(0:ndmis,1:ndmc,1:nzsmx) = 0.0_8
          twciL(1:ndmc,1:nzsmx)         = 0.0_8
          tengzL(0:ndmis,1:ndmc,1:nzsmx)= 0.0_8
          tprzL(0:ndmis,1:ndmc,1:nzsmx) = 0.0_8

          ismaxL(1:nzsmx) = -1
          nzmxL(1:nzsmx)  =  0
        endif
      case( 2 )
! deallocate
        if( allocated( aimasL ) ) then
          deallocate( aimasL, amiL, amzL, azmasL, cdifL, eipotL, ismaxL
     >      , nzmxL, tdnzL, twciL, tengzL, tprzL
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'aimasL deallocate error in cplimp_alc, istat = ', istat
            call wexit( 'cplimp', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
