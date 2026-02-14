!*******************************************************************************
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   subroutine def_lnkimp_cal( kerr )
      subroutine def_IMP_2( ctyp, code, kerr )
!***********************************************************************
      use cimden,      only : bkspt, csput, eipot, fsput, ncmx, npsput
     >    , nsput, nsrc_spt, nsrc_spt_dummy, nwspt, nzmx, wtsput
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
!ik   call dtyp_init( "IMP_2", "IMPMC" )
      call dtyp_init( trim( ctyp ), trim( code ) )

!::../../IMPMC/inc/cimden  /cimden/
      call dtyp_addv( nsput )
      call dtyp_addv( nzmx )
      call dtyp_addv( ncmx )
      call dtyp_addv( npsput )
! added 2 lines replace all include files with module files by kamata 2021/08/18
      call dtyp_addv( nsrc_spt )
      call dtyp_addv( nsrc_spt_dummy )
      call dtyp_addv( fsput )
      call dtyp_addv( wtsput )
      call dtyp_addv( bkspt )
      call dtyp_addv( eipot )
! deleted 4 lines dynamic allocation of arrays by kamata 2022/05/29
!ik   call dtyp_addv( twrd )
!ik   call dtyp_addv( tdnz )
!ik   call dtyp_addv( twci )
!ik   call dtyp_addv( twne )
      call dtyp_addv( csput )
      call dtyp_addv( nwspt )
!::../../IMPMC/inc/cimcom  /cimcom_12/
!      call dtyp_addv( denFlxZtoWall )
!      call dtyp_addv( eneFlxZtoWall )

!::term
      call dtyp_term( kerr )

      return
      end
