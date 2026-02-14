!**********************************************************************
      subroutine def_PLS_3( kerr )
!**********************************************************************
      use cplwrd,      only : wfac
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  :: kerr
      integer, intent(out) :: kerr

!::init
      call dtyp_init( "PLS_3", "soldor" )

!::../../soldor/inc/cplwrd
      call dtyp_addv( wfac )

!::term
      call dtyp_term( kerr )

      return
      end
