! added integrated calculation with TOPICS by kamata 2020/11/16
! remade addition and reconsider of passed data by kamata 2022/02/23
! set type_struct in SONIC
      subroutine strct_sonic( cknd )
      use cplmet,      only : icaxs, icmpe, icmps
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      use topics_mod,  only : dtcal, gmion, gnion, gtim, kcon, lexdt
     >    , lhist, nro, nroblk
      implicit none
! arguments
      character, intent(in) :: cknd*(*)
! cknd : type_struct name

! local variables
      integer    ierr

      select case( cknd )
      case ( 'NRO' )  ! rho mesh size in TOPICS
        call dtyp_init( 'NRO', 'SONIC' )
! trn_mod
        call dtyp_addv( nro )
        call dtyp_addv( nroblk )
        call dtyp_addv( gnion )
        call dtyp_addv( gmion )
        call dtyp_term( ierr )
      case ( 'MESH' ) ! mesh size in SONIC
        call dtyp_init( 'MESH', 'SONIC' )
        call dtyp_addv( icaxs )
        call dtyp_addv( icmpe )
        call dtyp_addv( icmps )
        call dtyp_term( ierr )
      case ( 'PROF' ) ! profile data in TOPICS
        call dtyp_init( 'PROF', 'SONIC' )
! time [s]
        call dtyp_addv( gtim )
! dtcal [s]
        call dtyp_addv( dtcal )
! lexdt
        call dtyp_addv( lexdt )
! lhist
        call dtyp_addv( lhist )
        call dtyp_term( ierr )
      case ( 'EDGE' ) ! edge data in SONIC
        call dtyp_init( 'EDGE', 'SONIC' )
! time [s]
        call dtyp_addv( gtim )
! kcon
        call dtyp_addv( kcon )
        call dtyp_term( ierr )
      end select

      return
      end
