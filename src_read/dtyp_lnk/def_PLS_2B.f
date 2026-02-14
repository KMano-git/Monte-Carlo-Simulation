!*******************************************************************************
      subroutine def_PLS_2B( kerr )
!***********************************************************************
      use cplcom,      only : gfvl, tfps, tfvl
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      implicit none

      integer, intent(out) :: kerr

!::init
      call dtyp_init( "PLS_2B", "soldor" )

!::../../solodr/inc/cplcom  /cmpflx/
      call dtyp_addv( tfps )
      call dtyp_addv( tfvl )
      call dtyp_addv( gfvl )

!::term
      call dtyp_term( kerr )

      return
      end
