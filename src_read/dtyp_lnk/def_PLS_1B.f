!*******************************************************************************
      subroutine def_PLS_1B( kerr )
!***********************************************************************
      use cplcom,      only : mdl_wrd
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  :: kerr
      integer, intent(out) :: kerr

!::init
      call dtyp_init( "PLS_1B", "soldor" )

!::../../soldor/inc/cplcom  /cmmodel/  mdl_wrd
      call dtyp_addv( mdl_wrd )

!::term
      call dtyp_term( kerr )

      return
      end

