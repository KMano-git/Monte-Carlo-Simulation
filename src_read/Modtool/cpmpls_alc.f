! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cpmpls
      subroutine cpmpls_alc( kk )
      use cpmpls, only : adamp, ahpmp, alsfr, an0mp, arhmp, armp, asnmp
     >    , avpmp, axemp, aximp, drmp, dvmp, mdp_ni, mdp_rh, mdp_ro
     >    , mdp_te, mdp_ti, rohmp, romp, sfmp, vlmp, wdmp, wfmp, wlmp
      use csize,  only : ndsp, ndy
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
        if( .not. allocated( adamp ) ) then
          allocate( adamp(ndy,ndsp), ahpmp(ndy,ndsp), alsfr(ndy)
     >      , an0mp(ndy), arhmp(ndy), armp(ndy), asnmp(ndy)
     >      , avpmp(ndy,ndsp), axemp(ndy), aximp(ndy), drmp(ndy)
     >      , dvmp(ndy), mdp_ni(ndy), mdp_rh(ndy), mdp_ro(ndy)
     >      , mdp_te(ndy), mdp_ti(ndy), rohmp(ndy), romp(ndy), sfmp(ndy)
     >      , vlmp(ndy), wdmp(ndy), wfmp(ndy), wlmp(ndy)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'adamp allocate error in cpmpls_alc, istat = ', istat
            call wexit( 'cpmpls', trim( cmsg ) )
          endif

! initial set
          adamp(1:ndy,1:ndsp) = 0.0_8
          ahpmp(1:ndy,1:ndsp) = 0.0_8
          alsfr(1:ndy) = 0.0_8
          an0mp(1:ndy) = 0.0_8
          arhmp(1:ndy) = 0.0_8
          armp(1:ndy) = 0.0_8
          asnmp(1:ndy) = 0.0_8
          avpmp(1:ndy,1:ndsp) = 0.0_8
          axemp(1:ndy) = 0.0_8
          aximp(1:ndy) = 0.0_8
          drmp(1:ndy) = 0.0_8
          dvmp(1:ndy) = 0.0_8
          mdp_ni(1:ndy) = 0.0_8
          mdp_rh(1:ndy) = 0.0_8
          mdp_ro(1:ndy) = 0.0_8
          mdp_te(1:ndy) = 0.0_8
          mdp_ti(1:ndy) = 0.0_8
          rohmp(1:ndy) = 0.0_8
          romp(1:ndy) = 0.0_8
          sfmp(1:ndy) = 0.0_8
          vlmp(1:ndy) = 0.0_8
          wdmp(1:ndy) = 0.0_8
          wfmp(1:ndy) = 0.0_8
          wlmp(1:ndy) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( adamp ) ) then
          deallocate( adamp, ahpmp, alsfr, an0mp, arhmp, armp, asnmp
     >      , avpmp, axemp, aximp, drmp, dvmp, mdp_ni, mdp_rh, mdp_ro
     >      , mdp_te, mdp_ti, rohmp, romp, sfmp, vlmp, wdmp, wfmp, wlmp
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'adamp deallocate error in cpmpls_alc, istat = ', istat
            call wexit( 'cpmpls', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
