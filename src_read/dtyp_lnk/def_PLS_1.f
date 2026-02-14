!*********************************************************************
      subroutine def_PLS_1( kerr )
!*********************************************************************
!
!      additional data  /cmtrcy/  rcydt, rcysp, rcytm  lnkimp_rcy
!                       /cmmodel/ mdl_wrd              lnkimp_ini
!
!                                             2015/04/22  K.Shimizu
!---------------------------------------------------------------------
      use cgdcom,      only : mpax, mpman, mpprv, mpsol, mpsp, mpw1
     >    , mpw2, mqd1, mqd2, mqh1, mqh2, mqs1, mqs2, mqx1, mqx2, npmx
     >    , nqmx
      use cplcom,      only : aion, ama, aza, cimp, cimprg, cimprg2
     >    , mdl_wrd, nion, rcydt, rcysp, rcytm, trad, trad_fac, trad_max
     >    , trad_min, tradrg, qtim_ini
      use cplmet,      only : icaxs, icmpe, icmps, icspx, icwl1, icwl2
     >    , itmax, itmpe, itmps, itpve, itpvs, itsle, itsls, jcdp1
     >    , jcdp2, jcmax, jcxp1, jcxp2, kce, kcn, kcs, kcw, nompl, noprv
     >    , nosol
      use cpmpls,      only : fflna, flxni, flxqe, flxqi, fprna, imd1
     >    , imd2, jmd1, jmd2, lnedg, prfna, prfni, prfte, prfti, r0mp
     >    , ramp, xmd1, xmd2, ymd1, ymd2, ltedg
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  :: kerr
      integer, intent(out) :: kerr

!::init
      call dtyp_init( "PLS_1", "soldor" )

!::../../soldor/inc/cplcom  /cmispc/
      call dtyp_addv( aion )
      call dtyp_addv( aza )
      call dtyp_addv( ama )
      call dtyp_addv( nion )

!::../../sonic/inc/size_68/cgdcom  /cgmsh/
! deleted 8 lines dynamic allocation of arrays by kamata 2022/05/29
!ik   call dtyp_addv( grdx )
!ik   call dtyp_addv( grdy )
!ik   call dtyp_addv( hbr )
!ik   call dtyp_addv( hbz )
!ik   call dtyp_addv( hbt )
!ik   call dtyp_addv( pssol )
!ik   call dtyp_addv( psprv )
!ik   call dtyp_addv( psman )
      call dtyp_addv( mqd1 )
      call dtyp_addv( mqx1 )
      call dtyp_addv( mqs1 )
      call dtyp_addv( mqs2 )
      call dtyp_addv( mqh1 )
      call dtyp_addv( mqh2 )
      call dtyp_addv( mqx2 )
      call dtyp_addv( mqd2 )
      call dtyp_addv( nqmx )
      call dtyp_addv( mpw1 )
      call dtyp_addv( mpsp )
      call dtyp_addv( mpw2 )
      call dtyp_addv( mpax )
      call dtyp_addv( npmx )
      call dtyp_addv( mpsol )
      call dtyp_addv( mpprv )
      call dtyp_addv( mpman )

!::../../soldor/inc/cplmet  /cmetrc/
! deleted 18 lines dynamic allocation of arrays by kamata 2022/05/29
!ik   call dtyp_addv( vlmn )
!ik   call dtyp_addv( romn )
!ik   call dtyp_addv( hvol )
!ik   call dtyp_addv( hgdx )
!ik   call dtyp_addv( hgdy )
!ik   call dtyp_addv( hdxm )
!ik   call dtyp_addv( hdxp )
!ik   call dtyp_addv( hdsp )
!ik   call dtyp_addv( hdsv )
!ik   call dtyp_addv( hvsb )
!ik   call dtyp_addv( hare )
!ik   call dtyp_addv( hpit )
!ik   call dtyp_addv( hwtm )
!ik   call dtyp_addv( hwtp )
!ik   call dtyp_addv( gdsv )
!ik   call dtyp_addv( gare )
!ik   call dtyp_addv( gwtm )
!ik   call dtyp_addv( gwtp )
      call dtyp_addv( nosol )
      call dtyp_addv( noprv )
      call dtyp_addv( nompl )
      call dtyp_addv( jcdp1 )
      call dtyp_addv( jcdp2 )
      call dtyp_addv( jcxp1 )
      call dtyp_addv( jcxp2 )
      call dtyp_addv( jcmax )
      call dtyp_addv( icwl1 )
      call dtyp_addv( icspx )
      call dtyp_addv( icwl2 )
      call dtyp_addv( icmps )
      call dtyp_addv( icmpe )
      call dtyp_addv( icaxs )
      call dtyp_addv( itsls )
      call dtyp_addv( itsle )
      call dtyp_addv( itpvs )
      call dtyp_addv( itpve )
      call dtyp_addv( itmps )
      call dtyp_addv( itmpe )
      call dtyp_addv( itmax )
! deleted 11 lines dynamic allocation of arrays by kamata 2022/05/29
!ik   call dtyp_addv( icmin )
!ik   call dtyp_addv( icmax )
!ik   call dtyp_addv( jtmin )
!ik   call dtyp_addv( jtmax )
!ik   call dtyp_addv( kgdx )
!ik   call dtyp_addv( kgdy )
!ik   call dtyp_addv( kreg )
!ik   call dtyp_addv( jcel )
!ik   call dtyp_addv( icel )
!ik   call dtyp_addv( jnxp )
!ik   call dtyp_addv( jnxm )
      call dtyp_addv( kce )
      call dtyp_addv( kcw )
      call dtyp_addv( kcn )
      call dtyp_addv( kcs )

!::../../soldor/inc/cpmpls  /cpmain/
      call dtyp_addv( flxni )
      call dtyp_addv( flxqi )
      call dtyp_addv( flxqe )
      call dtyp_addv( fflna )
      call dtyp_addv( prfna )
      call dtyp_addv( prfni )
      call dtyp_addv( prfti )
      call dtyp_addv( prfte )
      call dtyp_addv( fprna )
      call dtyp_addv( r0mp )
      call dtyp_addv( ramp )
! deleted 11 lines dynamic allocation of arrays by kamata 2022/05/29
!ik   call dtyp_addv( armp )
!ik   call dtyp_addv( arhmp )
!ik   call dtyp_addv( romp )
!ik   call dtyp_addv( rohmp )
!ik   call dtyp_addv( vlmp )
!ik   call dtyp_addv( sfmp )
!ik   call dtyp_addv( dvmp )
!ik   call dtyp_addv( drmp )
!ik   call dtyp_addv( wdmp )
!ik   call dtyp_addv( wfmp )
!ik   call dtyp_addv( wlmp )
      call dtyp_addv( xmd1 )
      call dtyp_addv( ymd1 )
      call dtyp_addv( xmd2 )
      call dtyp_addv( ymd2 )
      call dtyp_addv( jmd1 )
      call dtyp_addv( imd1 )
      call dtyp_addv( jmd2 )
      call dtyp_addv( imd2 )
      call dtyp_addv( lnedg )
      call dtyp_addv( ltedg )

!::../../soldor/inc/cplcom  /cmcrad/
! deleted 2 lines dynamic allocation of arrays by kamata 2022/05/29
!ik   call dtyp_addv( wrad )
!ik   call dtyp_addv( wdnz )
      call dtyp_addv( tradrg )
      call dtyp_addv( trad )
      call dtyp_addv( trad_min )
      call dtyp_addv( trad_max )
      call dtyp_addv( trad_fac )
      call dtyp_addv( cimp )
      call dtyp_addv( cimprg )
      call dtyp_addv( cimprg2 )

!::../../soldor/inc/cplcom /cmtrcy/ rcydt,rcysp,rcytm
!                              lnkimp_rcy
      call dtyp_addv( rcydt )
      call dtyp_addv( rcysp )
      call dtyp_addv( rcytm )

!::../../soldor/inc/cplcom /cmmodel/ mdl_wrd
!                               lnkimp_ini
      call dtyp_addv( mdl_wrd )
!
!:: send qtim_init for time dependent profile in continue calculation
      call dtyp_addv( qtim_ini )
!
!::term
      call dtyp_term( kerr )

      return
      end
