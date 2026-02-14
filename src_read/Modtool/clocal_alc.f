! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in clocal
      subroutine clocal_alc( kk )
      use clocal, only : man_na, man_ne, man_ni, man_te, man_ti, tpsi
     >    , troh, tvol, vden0, vdeng, vpre0, vpreg
      use csize,  only : ndgs, ndsp, ndwp, ndy
      use cunit,  only : cprg
      implicit none
! arguments
      integer, intent(in) :: kk
! kk : flag = 1: allocate,    = 2: deallocate

! local variables
      integer    istat
      character  cmsg*80

      if( cprg(1:5) == 'IMPMC' ) then
        select case( kk )
        case( 1 )
! allocate
          if( .not. allocated( tpsi ) ) then
            allocate( tpsi(ndy), troh(ndy), tvol(ndy)
     >        , stat = istat )
            if( istat /= 0 ) then
              write(cmsg,'(a,i6)')
     >          'tpsi allocate error in clocal_alc, istat = ', istat
              call wexit( 'clocal', trim( cmsg ) )
            endif
! initial set
            tpsi(1:ndy) = 0.0_8 
            troh(1:ndy) = 0.0_8 
            tvol(1:ndy) = 0.0_8
          endif
        case( 2 )
! deallocate
          if( allocated( tpsi ) ) then
            deallocate( tpsi, troh, tvol
     >        , stat = istat )
            if( istat /= 0 ) then
              write(cmsg,'(a,i6)')
     >          'tpsi deallocate error in clocal_alc, istat = ', istat
              call wexit( 'clocal', trim( cmsg ) )
            endif
          endif
        end select

      elseif( cprg(1:6) == 'neut2d' ) then
        select case( kk )
        case( 1 )
! allocate
          if( .not. allocated( vden0 ) ) then
            allocate( vden0(ndwp,ndgs), vdeng(ndwp,2), vpre0(ndwp,ndgs)
     >        , vpreg(ndwp,2)
     >        , stat = istat )
            if( istat /= 0 ) then
              write(cmsg,'(a,i6)')
     >          'vden0 allocate error in clocal_alc, istat = ', istat
              call wexit( 'clocal', trim( cmsg ) )
            endif
! initial set
            vden0(1:ndwp,1:ndgs) = 0.0_8
            vdeng(1:ndwp,1:2) = 0.0_8
            vpre0(1:ndwp,1:ndgs) = 0.0_8
            vpreg(1:ndwp,1:2) = 0.0_8
          endif
        case( 2 )
! deallocate
          if( allocated( vden0 ) ) then
            deallocate( vden0, vdeng, vpre0, vpreg
     >        , stat = istat )
            if( istat /= 0 ) then
              write(cmsg,'(a,i6)')
     >          'vden0 deallocate error in clocal_alc, istat = ', istat
              call wexit( 'clocal', trim( cmsg ) )
            endif
          endif
        end select

      elseif( cprg(1:6) == 'soldor' ) then
        select case( kk )
        case( 1 )
! allocate
          if( .not. allocated( man_na ) ) then
            allocate( man_na(ndy,ndsp), man_ne(ndy), man_ni(ndy)
     >        , man_te(ndy), man_ti(ndy)
     >        , stat = istat )
            if( istat /= 0 ) then
              write(cmsg,'(a,i6)')
     >          'man_na allocate error in clocal_alc, istat = ', istat
              call wexit( 'clocal', trim( cmsg ) )
            endif
! initial set
            man_na(1:ndy,1:ndsp) = 0.0_8
            man_ne(1:ndy) = 0.0_8
            man_ni(1:ndy) = 0.0_8
            man_te(1:ndy) = 0.0_8
            man_ti(1:ndy) = 0.0_8
          endif
        case( 2 )
! deallocate
          if( allocated( man_na ) ) then
            deallocate( man_na, man_ne, man_ni, man_te, man_ti
     >        , stat = istat )
            if( istat /= 0 ) then
              write(cmsg,'(a,i6)')
     >          'man_na deallocate error in clocal_alc, istat = ', istat
              call wexit( 'clocal', trim( cmsg ) )
            endif
          endif
        end select
      
        elseif( cprg(1:5) == 'plimp' ) then
          select case( kk )
          case( 1 )
  ! allocate
            if( .not. allocated( tpsi ) ) then
              allocate( tpsi(ndy), troh(ndy), tvol(ndy)
     >        , stat = istat )
              if( istat /= 0 ) then
                write(cmsg,'(a,i6)')
     >          'tpsi allocate error in clocal_alc, istat = ', istat
                call wexit( 'clocal', trim( cmsg ) )
              endif
  ! initial set
              tpsi(1:ndy) = 0.0_8 
              troh(1:ndy) = 0.0_8 
              tvol(1:ndy) = 0.0_8
            endif
          case( 2 )
  ! deallocate
            if( allocated( tpsi ) ) then
              deallocate( tpsi, troh, tvol
     >        , stat = istat )
              if( istat /= 0 ) then
                write(cmsg,'(a,i6)')
     >          'tpsi deallocate error in clocal_alc, istat = ', istat
                call wexit( 'clocal', trim( cmsg ) )
              endif
            endif
          end select
      endif
      return
      end
