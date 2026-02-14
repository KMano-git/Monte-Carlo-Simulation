!*******************************************************************************
      subroutine def_NTL_1( kerr )
!***********************************************************************
      use cmolhy,      only : lmolhcr
      use cntcom,      only : chwl, cpsf, cwal, e0in, e0pfa, e0pfm
     >    , flxin, grax, grxp, gzax, gzxp, iaxs, ihgw, impl, ipfnp
     >    , ipfmx, irfw, ispx, ivb1, ivb2, ivs1, ivs2, iwl1
     >    , iwl2, jdp1, jdp2, jsl1, jsl2, jvb1, jvb2, jxp1, jxp2, ldpmd
     >    , lelmd, lemd2, lemdp, lempf, lemwl, lnnel, lntmd, lrcmd
     >    , lsmax, lsrmd, lwlmd, nbxp, nbyp, ncmax, ncmax2, negt, ngas
     >    , ngmax, ngmax2, nhwl, nogt, noia, nops, nowl, nowl2, npep
     >    , npew, npf, npmp, npps, npsp, npsw, nptldp, nptlvl, nptlwl
     >    , npwl, npwl2, nsgt, pfmag, pfpx1, pfpx2, pfpy1, pfpy2, rcywl
     >    , rmas, temin_ion, temin_rec, temwl, tlmt_el, v0in, wtmin
      use cntctl,      only : dtntl, hrt_dbg, imxnt, itmnt, itmntb, kntl
     >    , mdl_hrt, mdl_ntsm, nclnt, nhfl, nhunt, nini_cal, nntl, tmntl
      use cntmnt,      only : mcflx, mfmax, mnknd, mnsmp, mnvwk
      use cntpfctl,    only : lpfctl
      use mod_sizedef, only : lnnam
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  :: kerr
      integer, intent(out) :: kerr

!::init
      call dtyp_init( "NTL_1", "neut2d" )

!::../../monte/inc/cntcom  /cntgrd/
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

!::../../monte/inc/cntcom  /cntwal/
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

!::../../monte/inc/cntcom  /cntgat/
!     sub. ntgdsk  Bcast for /cntgrd/, /cntwal/, /cntgat/
      call dtyp_addv( nsgt )
      call dtyp_addv( negt )
      call dtyp_addv( nogt )

!::../../monte/inc/cntcom  /cntpsf/
      call dtyp_addv( cpsf )
      call dtyp_addv( npsp )
      call dtyp_addv( npep )
      call dtyp_addv( nbxp )
      call dtyp_addv( nbyp )
      call dtyp_addv( nops )
      call dtyp_addv( npps )

!::../../monte/inc/cntcom  /cntpuf/
      call dtyp_addv( ipfmx )

!::../../monte/inc/cntctl  /cntcmp/
      call dtyp_addv( tmntl )
      call dtyp_addv( dtntl )
      call dtyp_addv( nntl )
      call dtyp_addv( imxnt )
      call dtyp_addv( nhunt )
      call dtyp_addv( nhfl )
      call dtyp_addv( nini_cal )
      call dtyp_addv( mdl_ntsm )
      call dtyp_addv( kntl )
      call dtyp_addv( nclnt )
      call dtyp_addv( itmnt )
      call dtyp_addv( itmntb )
      call dtyp_addv( mdl_hrt )
      call dtyp_addv( hrt_dbg )

!::../../monte/inc/cntmnt  /cntmfl/
      call dtyp_addv( mfmax )
      call dtyp_addv( mnknd )
      call dtyp_addv( mnsmp )
      call dtyp_addv( mnvwk )
      call dtyp_addv( mcflx )

!::../../monte/inc/cntcom  /cntinp/
      call dtyp_addv( rmas )
      call dtyp_addv( e0in )
      call dtyp_addv( v0in )
      call dtyp_addv( flxin )
      call dtyp_addv( wtmin )
      call dtyp_addv( pfmag )
      call dtyp_addv( pfpx1 )
      call dtyp_addv( pfpx2 )
      call dtyp_addv( pfpy1 )
      call dtyp_addv( pfpy2 )
      call dtyp_addv( e0pfa )
      call dtyp_addv( e0pfm )
      call dtyp_addv( temin_ion )
      call dtyp_addv( temin_rec )
      call dtyp_addv( tlmt_el )
      call dtyp_addv( ngas )
      call dtyp_addv( noia )
      call dtyp_addv( nptldp )
      call dtyp_addv( nptlwl )
      call dtyp_addv( nptlvl )
      call dtyp_addv( lsmax )
      call dtyp_addv( lemwl )
      call dtyp_addv( lemdp )
      call dtyp_addv( lempf )
      call dtyp_addv( lemd2 )
      call dtyp_addv( ipfnp )
      call dtyp_addv( npf )
      call dtyp_addv( lmolhcr )

!::../../monte/inc/cntcom  /cntmdl/
      call dtyp_addv( lntmd )
      call dtyp_addv( ldpmd )
      call dtyp_addv( lwlmd )
      call dtyp_addv( lrcmd )
      call dtyp_addv( lelmd )
      call dtyp_addv( lnnel )
      call dtyp_addv( lsrmd )
!:: module cntpfctl
      call dtyp_addv( lpfctl )

!::term
      call dtyp_term( kerr )

      return
      end
