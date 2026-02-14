! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cntmnt
      subroutine cntmnt_alc( kk )
      use cntmnt, only : sn0, ssn, ssp, swe, swi, vmwork
      use csize,  only : ndmsr, ndsp, ndx, ndy, nwkmp_nt
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
        if( .not. allocated( sn0 ) ) then
          allocate( sn0(ndx,ndy,ndsp), ssn(ndx,ndy,ndsp)
     >      , ssp(ndx,ndy,ndsp), swe(ndx,ndy), swi(ndx,ndy)
     >      , vmwork(nwkmp_nt,ndmsr)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'sn0 allocate error in cntmnt_alc, istat = ', istat
            call wexit( 'cntmnt', trim( cmsg ) )
          endif

! initial set
          sn0(1:ndx,1:ndy,1:ndsp) = 0.0_8
          ssn(1:ndx,1:ndy,1:ndsp) = 0.0_8
          ssp(1:ndx,1:ndy,1:ndsp) = 0.0_8
          swe(1:ndx,1:ndy) = 0.0_8
          swi(1:ndx,1:ndy) = 0.0_8
          vmwork(1:nwkmp_nt,1:ndmsr) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( sn0 ) ) then
          deallocate( sn0, ssn, ssp, swe, swi, vmwork
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'sn0 deallocate error in cntmnt_alc, istat = ', istat
            call wexit( 'cntmnt', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
