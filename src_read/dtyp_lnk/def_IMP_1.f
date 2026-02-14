!*******************************************************************************
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   subroutine def_lnkimp_iniSnd( kerr )
      subroutine def_IMP_1( ctyp, code, kerr )
!***********************************************************************
      use cimcom,      only : aimas, ami, amz, azmas, cdif, cdifrg
     >    , ismax, lwrad
      use cimctl,      only : cdirz, dtimp, ftav, icalz, kimp, lstdy
     >    , ndirz, nimp, nsvtm, svdsn, svitm, tmimp
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      use mod_sizedef, only : lnnam
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  :: kerr
! added 1 line treat 4 or more impurities with IMPMC by kamata 2022/04/21
      character, intent(in)  :: code*10, ctyp*(lnnam)
      integer,   intent(out) :: kerr
! added 3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
! code : code name
! ctyp : data label
! kerr : completion code

!::init
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   call dtyp_init( "IMP_1", "IMPMC" )
      call dtyp_init( trim( ctyp ), trim( code ) )

!::../../ATadas/inc/catcom  /catcom/
! deleted 27 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   call dtyp_addv( xdaty )
!ik   call dtyp_addv( xdtx1 )
!ik   call dtyp_addv( xdtx2 )
!ik   call dtyp_addv( xdatz )
!ik   call dtyp_addv( xwmin )
!ik   call dtyp_addv( xwmax )
!ik   call dtyp_addv( xwmlt )
!ik   call dtyp_addv( ip_ion )
!ik   call dtyp_addv( ip_rec )
!ik   call dtyp_addv( ip_cxr )
!ik   call dtyp_addv( ip_rad )
!ik   call dtyp_addv( ip_plt )
!ik   call dtyp_addv( ip_prb )
!ik   call dtyp_addv( ip_prc )
!ik   call dtyp_addv( nxs )
!ik   call dtyp_addv( nyp )
!ik   call dtyp_addv( xdrnk )
!ik   call dtyp_addv( xwnum )
!ik   call dtyp_addv( xjsta )
!ik   call dtyp_addv( xjend )
!ik   call dtyp_addv( xjnum )
!ik   call dtyp_addv( xdsn )
!ik   call dtyp_addv( xdnam )
!ik   call dtyp_addv( xdevl )
!ik   call dtyp_addv( xwvar )
!ik   call dtyp_addv( xwunt )
!ik   call dtyp_addv( xwspc )

!::../../ATadas/inc/catcom  /catcom3/
! deleted 3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   call dtyp_addv( catmz )
!ik   call dtyp_addv( natmz )
!ik   call dtyp_addv( aionz )

!::../../IMPMC/inc/cimctl  /cimcmp/
      call dtyp_addv( tmimp )
      call dtyp_addv( dtimp )
      call dtyp_addv( ftav )
      call dtyp_addv( nimp )
      call dtyp_addv( kimp )
      call dtyp_addv( icalZ )
      call dtyp_addv( lstdy )
      call dtyp_addv( svitm )
      call dtyp_addv( nsvtm )
      call dtyp_addv( ndirZ )
      call dtyp_addv( svdsn )
      call dtyp_addv( cdirZ )

!::../../IMPMC/inc/cimcom  /cimcom_1/
      call dtyp_addv( aimas )
      call dtyp_addv( azmas )
      call dtyp_addv( ami )
      call dtyp_addv( amz )
      call dtyp_addv( ismax )

!::../../IMPMC/inc/cimcom  /cimcom_8/
      call dtyp_addv( cdif )
      call dtyp_addv( cdifrg )
! deleted 1 line dynamic allocation of arrays by kamata 2022/05/29
!ik   call dtyp_addv( dfcf )

!::../../IMPMC/inc/cimcom  /cimcom_92/ lwrad
      call dtyp_addv( lwrad )

!::term
      call dtyp_term( kerr )

      return
      end
