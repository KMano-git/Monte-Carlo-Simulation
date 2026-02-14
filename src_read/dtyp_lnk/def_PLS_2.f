!**********************************************************************
      subroutine def_PLS_2( kerr )
!**********************************************************************
      use cntwfl,      only : xcfl, xsno
      use cplcom,      only : gfvl, qtim, snflx, tfps, tfvl
      use cplwrd,      only : wfac
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      implicit none

      integer, intent(out) :: kerr

!::init
      call dtyp_init( "PLS_2", "soldor" )

!::../../monte/inc/cntmnt  /cplpfl/
!::../../soldor/inc/cplcom  /cmpflx/  2015/04/23
      call dtyp_addv( tfps )
      call dtyp_addv( tfvl )
      call dtyp_addv( gfvl )
      call dtyp_addv( qtim )
      call dtyp_addv( snflx )

!::../../monte/inc/cntwfl
      call dtyp_addv( xsno  )
      call dtyp_addv( xcfl  )

!::../../soldor/inc/cplwrd
      call dtyp_addv( wfac )

!::term
      call dtyp_term( kerr )

      return
      end
