! added replace all include files with module files by kamata 2021/08/18
      subroutine strct_comblk( cknd )
      use catmhy,      only : elhio, elhrc, emhio, emhrc, sghio, sghrc
     >    , sne, snedlt, snemax, snemin, ste, stedlt, stemax, stemin
      use cgdcom,      only : mpax, mpman, mpprv, mpsol, mpsp, mpw1
     >    , mpw2, mqd1, mqd2, mqh1, mqh2, mqs1, mqs2, mqx1, mqx2, npmx
     >    , nqmx
      use cimcom,      only : bkflux, bkfrac, scpum, sdmy1, sdmy2, sdmy3
     >    , sdtmz, sflux, sitmp, sitmz, sptcl, stimp, stimz, swtot
      use cmolhy,      only : cpid, cract, els3, elsm, eng3, engm, jt3e
     >    , jt3l, jtmx, nrc, sgvm, stedlt_m => stedlt, stem
     >    , stemax_m => stemax, stemin_m => stemin, te3e, te3l
      use cntcom,      only : chwl, cwal, grax, grxp, gzax, gzxp
     >    , iaxs, ihgw, impl, irfw, ispx, ivb1, ivb2, ivs1
     >    , ivs2, iwl1, iwl2, jdp1, jdp2, jsl1, jsl2, jvb1, jvb2
     >    , jxp1, jxp2, ncmax, ncmax2, negt, ngmax, ngmax2, nhwl, nogt
     >    , nowl, nowl2, npew, npmp, npsw, npwl, npwl2, nsgt, rcywl
     >    , temwl
      use cntmnt,      only : vcsrc, vflux, visrc, vitim, vitnt, vkflx
     >    , vksty, vnsmp, vnsrc, vsmno, vsmty, vtime
      use cntrfl,      only : dlge, mpe, rfelg, rfeng, rfer, rfrn
      use com_eqdat,   only : dr, dr0, dz, dz0, nr, nz, paxs, psep, psi
     >    , raxs, rbt0, rg, rmax, rmin, rsep, zaxs, zg, zmax, zmin, zsep
      use cplcom,      only : aion, ama, aza, nion
      use cplmet,      only : icaxs, icmpe, icmps, icspx, icwl1, icwl2
     >    , itmax, itmpe, itmps, itpve, itpvs, itsle, itsls, jcdp1
     >    , jcdp2, jcmax, jcxp1, jcxp2, kce, kcn, kcs, kcw, nompl, noprv
     >    , nosol
      use cpmpls,      only : fflna, flxni, flxqe, flxqi, fprna, imd1
     >    , imd2, jmd1, jmd2, lnedg, prfna, prfni, prfte, prfti, r0mp
     >    , ramp, xmd1, xmd2, ymd1, ymd2
      use cxdcom,      only : nxs, nyp, xdaty, xdevl, xdnam, xdrnk
     >    , xjend, xjnum, xjsta, xwmax, xwmin, xwmlt, xwnum, xwspc
     >    , xwunt, xwvar
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      implicit none
! arguments
      character, intent(in) :: cknd*(*)
! cknd : type_struct name

! local variables
      integer    ierr

      select case( cknd )
! IMPMC/inc/cimcom  /cimcom_11/
      case( 'CIMCOM_11' ) ! for IMPMC
        call dtyp_init( trim( cknd ), 'IMPMC' )
        call dtyp_addv( sflux )
        call dtyp_addv( swtot )
        call dtyp_addv( sptcl )
        call dtyp_addv( stimp )
        call dtyp_addv( sitmp )
        call dtyp_addv( stimz )
        call dtyp_addv( sdtmz )
        call dtyp_addv( sdmy1 )
        call dtyp_addv( sdmy2 )
        call dtyp_addv( sdmy3 )
        call dtyp_term( ierr )

! IMPMC/inc/cimcom  /cimcom_12d/
      case( 'CIMCOM_12D' ) ! for IMPMC
        call dtyp_init( trim( cknd ), 'IMPMC' )
        call dtyp_addv( bkfrac )
        call dtyp_addv( bkflux )
        call dtyp_term( ierr )
! IMPMC_TD/inc/cimcom  /cimcom_11/
      case( 'CIMCOMT11' ) ! for IMPMC_TD
        call dtyp_init( trim( cknd ), 'IMPMC' )
        call dtyp_addv( sflux )
        call dtyp_addv( swtot )
        call dtyp_addv( sptcl )
        call dtyp_addv( stimp )
        call dtyp_addv( sitmp )
        call dtyp_addv( stimz )
        call dtyp_addv( sitmz )
        call dtyp_addv( sdtmz )
        call dtyp_addv( scpum )
        call dtyp_addv( sdmy1 )
        call dtyp_addv( sdmy2 )
        call dtyp_addv( sdmy3 )
        call dtyp_term( ierr )

! cdflib/inc/cxdcom  /cxdcom/
      case( 'CXDCOM' ) ! for NEUT2D
        call dtyp_init( trim( cknd ), 'cdflib' )
        call dtyp_addv( nxs )
        call dtyp_addv( nyp )
        call dtyp_addv( xdnam )
        call dtyp_addv( xdevl )
        call dtyp_addv( xwvar )
        call dtyp_addv( xwunt )
        call dtyp_addv( xwspc )
        call dtyp_addv( xdrnk )
        call dtyp_addv( xwnum )
        call dtyp_addv( xjsta )
        call dtyp_addv( xjend )
        call dtyp_addv( xjnum )
        call dtyp_addv( xwmin )
        call dtyp_addv( xwmax )
        call dtyp_addv( xwmlt )
        call dtyp_addv( xdaty )
        call dtyp_term( ierr )
! dtdegas/inc/catmhy  /catmhy/
      case( 'CATMHY' )
        call dtyp_init( trim( cknd ), 'dtdegas' )
        call dtyp_addv( stemin )
        call dtyp_addv( stemax )
        call dtyp_addv( stedlt )
        call dtyp_addv( snemin )
        call dtyp_addv( snemax )
        call dtyp_addv( snedlt )
        call dtyp_addv( ste )
        call dtyp_addv( sne )
        call dtyp_addv( sghio )
        call dtyp_addv( elhio )
        call dtyp_addv( emhio )
        call dtyp_addv( sghrc )
        call dtyp_addv( elhrc )
        call dtyp_addv( emhrc )
        call dtyp_term( ierr )
! dtdegas/inc/cmolhy  /cmolhy/
      case( 'CMOLHY' )
        call dtyp_init( trim( cknd ), 'dtdegas' )
        call dtyp_addv( stemin_m )
        call dtyp_addv( stemax_m )
        call dtyp_addv( stedlt_m )
        call dtyp_addv( stem )
        call dtyp_addv( sgvm )
        call dtyp_addv( elsm )
        call dtyp_addv( engm )
        call dtyp_addv( jtmx )
        call dtyp_addv( nrc )
        call dtyp_addv( te3l )
        call dtyp_addv( els3 )
        call dtyp_addv( te3e )
        call dtyp_addv( eng3 )
        call dtyp_addv( jt3l )
        call dtyp_addv( jt3e )
        call dtyp_term( ierr )
! dtdegas/inc/cmolhy  /cmolhy2/
      case( 'CMOLHY2' )
        call dtyp_init( trim( cknd ), 'dtdegas' )
        call dtyp_addv( cract )
        call dtyp_addv( cpid )
        call dtyp_term( ierr )
! gmequ/inc/com_eqdat  /com_eqdt/
      case( 'CEQDT' )
        call dtyp_init( trim( cknd ), 'gmequ' )
        call dtyp_addv( dr )
        call dtyp_addv( dz )
        call dtyp_addv( rbt0 )
        call dtyp_addv( rmin )
        call dtyp_addv( rmax )
        call dtyp_addv( zmin )
        call dtyp_addv( zmax )
        call dtyp_addv( dr0 )
        call dtyp_addv( dz0 )
        call dtyp_addv( raxs )
        call dtyp_addv( zaxs )
        call dtyp_addv( paxs )
        call dtyp_addv( rsep )
        call dtyp_addv( zsep )
        call dtyp_addv( psep )
        call dtyp_addv( rg )
        call dtyp_addv( zg )
        call dtyp_addv( psi )
        call dtyp_addv( nr )
        call dtyp_addv( nz )
        call dtyp_term( ierr )
! monte/inc/cntcom  /cntgat/
      case( 'CNTGAT' )
        call dtyp_init( trim( cknd ), 'neut2d' )
        call dtyp_addv( nsgt )
        call dtyp_addv( negt )
        call dtyp_addv( nogt )
        call dtyp_term( ierr )
! monte/inc/cntcom  /cntgrd/
      case( 'CNTGRD' )
        call dtyp_init( trim( cknd ), 'neut2d' )
        call dtyp_addv( grax )
        call dtyp_addv( gzax )
        call dtyp_addv( grxp )
        call dtyp_addv( gzxp )
        call dtyp_addv( ngmax )
        call dtyp_addv( ngmax2 )
        call dtyp_addv( ncmax )
        call dtyp_addv( ncmax2 )
        call dtyp_addv( jdp1 )
        call dtyp_addv( jxp1 )
        call dtyp_addv( jsl1 )
        call dtyp_addv( jsl2 )
        call dtyp_addv( jxp2 )
        call dtyp_addv( jdp2 )
        call dtyp_addv( iwl1 )
        call dtyp_addv( ispx )
        call dtyp_addv( iwl2 )
        call dtyp_addv( impl )
        call dtyp_addv( iaxs )
        call dtyp_addv( ivs1 )
        call dtyp_addv( ivs2 )
        call dtyp_addv( jvb1 )
        call dtyp_addv( jvb2 )
        call dtyp_addv( ivb1 )
        call dtyp_addv( ivb2 )
        call dtyp_term( ierr )
! monte/inc/cntcom  /cntwal/
      case( 'CNTWAL' )
        call dtyp_init( trim( cknd ), 'neut2d' )
        call dtyp_addv( rcywl )
        call dtyp_addv( temwl )
        call dtyp_addv( cwal )
        call dtyp_addv( chwl )
        call dtyp_addv( npsw )
        call dtyp_addv( npew )
        call dtyp_addv( irfw )
        call dtyp_addv( ihgw )
        call dtyp_addv( nowl )
        call dtyp_addv( npwl )
        call dtyp_addv( npmp )
        call dtyp_addv( nowl2 )
        call dtyp_addv( npwl2 )
        call dtyp_addv( nhwl )
        call dtyp_term( ierr )
! monte/inc/cntcom  /cntgat/
      case( 'CNTMNT' ) ! for NEUT2D
        call dtyp_init( trim( cknd ), 'neut2d' )
        call dtyp_addv( vtime )
        call dtyp_addv( vflux )
        call dtyp_addv( vitim )
        call dtyp_addv( vitnt )
        call dtyp_addv( vnsrc )
        call dtyp_addv( vsmty )
        call dtyp_addv( vsmno )
        call dtyp_addv( visrc )
        call dtyp_addv( vkflx )
        call dtyp_addv( vksty )
        call dtyp_addv( vnsmp )
        call dtyp_addv( vcsrc )
        call dtyp_term( ierr )
! monte/inc/cntrfl  /cmrefl/
      case( 'CMREFL' ) ! for NEUT2D
        call dtyp_init( trim( cknd ), 'neut2d' )
        call dtyp_addv( rfeng )
        call dtyp_addv( rfelg )
        call dtyp_addv( rfrn )
        call dtyp_addv( rfer )
        call dtyp_addv( dlge )
        call dtyp_addv( mpe )
        call dtyp_term( ierr )
! soldor/inc/cplcom  /cmispc/
      case( 'CMISPC' )
        call dtyp_init( trim( cknd ), 'soldor' )
        call dtyp_addv( aion )
        call dtyp_addv( ama )
        call dtyp_addv( aza )
        call dtyp_addv( nion )
        call dtyp_term( ierr )

! soldor/inc/cplmet  /cmetrc/
      case( 'CMETRC' )
        call dtyp_init( trim( cknd ), 'soldor' )
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
        call dtyp_addv( kce )
        call dtyp_addv( kcw )
        call dtyp_addv( kcn )
        call dtyp_addv( kcs )
        call dtyp_term( ierr )
! soldor/inc/cpmpls /cpmain/
      case( 'CPMAIN' ) ! for NEUT2D
        call dtyp_init( trim( cknd ), 'soldor' )
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
        call dtyp_addv( xmd1 )
        call dtyp_addv( ymd1 )
        call dtyp_addv( xmd2 )
        call dtyp_addv( ymd2 )
        call dtyp_addv( jmd1 )
        call dtyp_addv( imd1 )
        call dtyp_addv( jmd2 )
        call dtyp_addv( imd2 )
        call dtyp_addv( lnedg )
        call dtyp_term( ierr )

! sonic/inc/size_XX/cgdcom /cgmsh/
      case( 'CGMSH' )
        call dtyp_init( trim( cknd ), 'sonic' )
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
        call dtyp_term( ierr )
      end select

      return
      end
