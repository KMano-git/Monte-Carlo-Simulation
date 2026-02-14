!***********************************************************************
      subroutine def_MST_2( kerr )
!***********************************************************************
      use cntctl,      only : kntl
      use cimctl,      only : kimp
      use csonic,      only : dtim, itend, itim, kdsk, khst, kpcn, kpls
     >    , lcstg, lstop, mstop, tend, time
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  :: kerr
      integer, intent(out) :: kerr

!::init
      call dtyp_init( "MST_2", "soldor" )

!::../../sonic/inc/csonic  /cscntl/
      call dtyp_addv( time )
      call dtyp_addv( tend )
      call dtyp_addv( dtim )
      call dtyp_addv( itim )
      call dtyp_addv( itend )
      call dtyp_addv( mstop )
      call dtyp_addv( lstop )
      call dtyp_addv( lcstg )

!::../../monte/inc/cntctl  /cntcmp/
      call dtyp_addv( kntl )

!::../../IMPMC/inc/cimctl  /cimcmp/
      call dtyp_addv( kimp )

!::../../IMPMC/inc/csonic
      call dtyp_addv( kpcn )
      call dtyp_addv( kdsk )
      call dtyp_addv( khst )
      call dtyp_addv( kpls )

!::term
      call dtyp_term( kerr )

      return
      end
