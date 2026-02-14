! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cplpst
      subroutine cplpst_alc( kk )
      use cplpst, only : dwrad, farei, fareo, fbdpi, fbdpo, fldpi, fldpo
     >    , flsmi, flsmo, rhf_idp, rhf_imd, rhf_odp, rhf_omd
     >    , rhf_uidp, rhf_uodp
      use csize,  only : ndeq, ndxy, ndy
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
        if( .not. allocated( dwrad) ) then
          allocate( dwrad(ndxy,4), farei(ndy), fareo(ndy)
     >      , fbdpi(ndy,ndeq), fbdpo(ndy,ndeq), fldpi(ndy,ndeq)
     >      , fldpo(ndy,ndeq), flsmi(ndy,10), flsmo(ndy,10)
     >      , rhf_idp(ndy), rhf_imd(ndy), rhf_odp(ndy), rhf_omd(ndy)
     >      , rhf_uodp(ndy),rhf_uidp(ndy)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'dwrad allocate error in cplpst_alc, istat = ', istat
            call wexit( 'cplpst', trim( cmsg ) )
          endif

! initial set
          dwrad(1:ndxy,1:4) = 0.0_8
          farei(1:ndy) = 0.0_8
          fareo(1:ndy) = 0.0_8
          fbdpi(1:ndy,1:ndeq) = 0.0_8
          fbdpo(1:ndy,1:ndeq) = 0.0_8
          fldpi(1:ndy,1:ndeq) = 0.0_8
          fldpo(1:ndy,1:ndeq) = 0.0_8
          flsmi(1:ndy,1:10) = 0.0_8
          flsmo(1:ndy,1:10) = 0.0_8
          rhf_idp(1:ndy) = 0.0_8
          rhf_imd(1:ndy) = 0.0_8
          rhf_odp(1:ndy) = 0.0_8
          rhf_omd(1:ndy) = 0.0_8

          rhf_uodp(1:ndy) = 0.0_8
          rhf_uidp(1:ndy) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( dwrad ) ) then
          deallocate( dwrad, farei, fareo, fbdpi, fbdpo, fldpi, fldpo
     >      , flsmi, flsmo, rhf_idp, rhf_imd, rhf_odp, rhf_omd
     >      , rhf_uodp, rhf_uidp
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'dwrad deallocate error in cplpst_alc, istat = ', istat
            call wexit( 'cplpst', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
