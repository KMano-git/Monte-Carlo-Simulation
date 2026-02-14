!*******************************************************************************
      subroutine def_MST_1( kerr )
!***********************************************************************
      use csonic,      only : dtim, itend, itim, kdsk, khst, kpcn, kpls
     >    , lcode, lcstg, lfdbg, lfopt, limp, lmstd, lntl, lpchk, lpcn
     >    , lpls, lpost, lpst_bprf, lqick, lrand, lstep, lstop, mstop
     >    , mxcpu, mxdsk, mxhst, mxprf, ndsk, nfdbg, nfopt, nhsav, nhst
     >    , npcn, npls, nprf, tend, time
      use cunit,       only : cdrout, wh_day, wh_dir, wh_job, wh_lod
     >    , wh_tim
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  :: kerr
      integer, intent(out) :: kerr

!::init
      call dtyp_init( "MST_1", "soldor" )

!::../../sonic/inc/csonic  /csonic/
      call dtyp_addv( lstep )
      call dtyp_addv( lcode )
      call dtyp_addv( ndsk )
      call dtyp_addv( nhst )
      call dtyp_addv( nhsav )
      call dtyp_addv( nprf )
      call dtyp_addv( mxcpu )
      call dtyp_addv( mxdsk )
      call dtyp_addv( mxhst )
      call dtyp_addv( mxprf )
      call dtyp_addv( lfopt )
      call dtyp_addv( lfdbg )
      call dtyp_addv( lpst_bprf )
      call dtyp_addv( nfopt )
      call dtyp_addv( nfdbg )
      call dtyp_addv( lpls )
      call dtyp_addv( npls )
      call dtyp_addv( kpls )
      call dtyp_addv( lpcn )
      call dtyp_addv( npcn )
      call dtyp_addv( kpcn )
      call dtyp_addv( lpchk )
      call dtyp_addv( kdsk )
      call dtyp_addv( khst )
      call dtyp_addv( lntl )
      call dtyp_addv( limp )
      call dtyp_addv( lmstd )
      call dtyp_addv( lqick )
      call dtyp_addv( lpost )
      call dtyp_addv( lrand )

!::../../sonic/inc/csonic  /cscntl/
      call dtyp_addv( time )
      call dtyp_addv( tend )
      call dtyp_addv( dtim )
      call dtyp_addv( itim )
      call dtyp_addv( itend )
      call dtyp_addv( mstop )
      call dtyp_addv( lstop )
      call dtyp_addv( lcstg )

!::../../sonic/inc/cunit  /comstjb/
      call dtyp_addv( cdrout )
      call dtyp_addv( wh_day )
      call dtyp_addv( wh_tim )
      call dtyp_addv( wh_lod )
      call dtyp_addv( wh_dir )
      call dtyp_addv( wh_job )

!::term
      call dtyp_term( kerr )

      return
      end
