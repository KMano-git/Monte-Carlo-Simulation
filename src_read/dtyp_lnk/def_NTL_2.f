!**********************************************************************
      subroutine def_NTL_2( kerr )
!**********************************************************************
      use cntmnt,      only : vcsrc, vflux, visrc, vitim, vitnt, vkflx
     >    , vksty, vnsmp, vnsrc, vsmno, vsmty, vtime
      use cntpfctl,    only : fpfctl
      use mod_dtypdef, only : dtyp_addv, dtyp_init, dtyp_term
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  :: kerr
      integer, intent(out) :: kerr

!::init
      call dtyp_init( "NTL_2", "neut2d" )

!::../../monte/inc/cntmnt  /cntmnt/
      call dtyp_addv( vtime )
      call dtyp_addv( vflux )
! deleted 1 line dynamic allocation of arrays by kamata 2022/05/29
!ik   call dtyp_addv( vmwork )
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

!::../../monte/inc/cntavp  /com_ntavpl/
! deleted 9 lines dynamic allocation of arrays by kamata 2022/05/29
!ik   call dtyp_addv( hno )
!ik   call dtyp_addv( hni )
!ik   call dtyp_addv( hvi )
!ik   call dtyp_addv( hti )
!ik   call dtyp_addv( hte )
!ik   call dtyp_addv( gni )
!ik   call dtyp_addv( gvi )
!ik   call dtyp_addv( gti )
!ik   call dtyp_addv( gte )
!:: module cntpfctl
      call dtyp_addv( fpfctl )

!::term
      call dtyp_term( kerr )

      return
      end
