! added replace all include files with module files by kamata 2021/08/18
! set module variables ( allocate and derived data type )
      subroutine setmodv( codenm )
      use mod_dtflw,   only : ntyp
      use mod_dtypdef, only : nddt
      use cunit,       only : lmspe, lmype
      implicit none
! arguments
      character, intent(in) :: codenm*(*)
! codenm : code name ( 'master', 'soldor', 'neut2d', 'IMPMC', 'IMPMC2', 'IMPMC3' )

! set the size required for allocate
      call alc_size
      call gvar_alc ( codenm, 1 )

! make the start position common to all codes 
      nddt = ntyp
!
! for master PE only
      if( lmype.eq.lmspe ) then
        call dbgdtyp( codenm, 1 )
      endif
!
!*e* by kamata 2021/10/14
! definition of derived type of common block
! IMPMC/inc/cimcom  /cimcom_11/
      call strct_comblk( 'CIMCOM_11' )
! IMPMC/inc/cimcom  /cimcom_11y/
      call strct_comblk( 'CIMCOM_11Y' )
! IMPMC/inc/cimcom  /cimcom_12d/
      call strct_comblk( 'CIMCOM_12D' )

! IMPMC_TD/inc/cimcom  /cimcom_11/
      call strct_comblk( 'CIMCOMT11' )
! IMPMC_TD/inc/cimcom  /cimcom_11y/
      call strct_comblk( 'CIMCOMT11Y' )

! cdflib/inc/cxdcom  /cxdcom/
      call strct_comblk( 'CXDCOM' )

! dtdegas/inc/catmhy  /catmhy/
      call strct_comblk( 'CATMHY' )
! dtdegas/inc/cmolhy  /cmolhy/
      call strct_comblk( 'CMOLHY' )
! dtdegas/inc/cmolhy  /cmolhy2/
      call strct_comblk( 'CMOLHY2' )

! gmequ/inc/com_eqdat  /com_eqdt/
      call strct_comblk( 'CEQDT' )

! monte/inc/cntcom  /cntgat/
      call strct_comblk( 'CNTGAT' )
! monte/inc/cntcom  /cntgrd/
      call strct_comblk( 'CNTGRD' )
! monte/inc/cntcom  /cntwal/
      call strct_comblk( 'CNTWAL' )
! monte/inc/cntmnt  /cntmnt/
      call strct_comblk( 'CNTMNT' )
! monte/inc/cntrfl  /cmrefl/
      call strct_comblk( 'CMREFL' )

! soldor/inc/cplcom  /cmispc/
      call strct_comblk( 'CMISPC' )
! soldor/inc/cplcom  /cmqnow/
      call strct_comblk( 'CMQNOW' )
! soldor/inc/cplmet  /cmetrc/
      call strct_comblk( 'CMETRC' )
! soldor/inc/cpmpls /cpmain/
      call strct_comblk( 'CPMAIN' )
! soldor/inc/cpmpls /com_plmdprf/
      call strct_comblk( 'CPLMDPRF' )

! sonic/inc/size_XX/cgdcom /cgmsh/
      call strct_comblk( 'CGMSH' )
!
! for master PE only
      if( lmype.eq.lmspe ) then
        call dbgdtyp( codenm, 2 )
      endif
!
      return
      end

!*s* by kamata 2021/10/14
      subroutine dbgdtyp( codenm, kk )
      use mod_dtflw,   only : ntyp
      use mod_dtypdef, only : nddt, typ_dnam, typ_itag
      implicit none
! arguments
      character, intent(in) :: codenm*(*)
      integer,   intent(in) :: kk ! ( = 1 : ntyp, = 2 : nddt )
! local variables
      integer :: i, io = 10140, nmx
      logical :: kfst = .true.

      if( kfst ) then
        select case( codenm )
        case( 'master' )
        case( 'soldor' )
          io = io + 1
        case( 'neut2d' )
          io = io + 2
        case( 'IMPMC' )
          io = io + 3
        case( 'IMPMC2' )
          io = io + 4
        case( 'IMPMC3' )
          io = io + 5
        end select
        kfst = .false.
      endif

      select case( kk )
      case( 1 )
        nmx = ntyp
      case( 2 )
        nmx = nddt
      end select

      write(io,'(a)') 'code = ' // trim( codenm )
      write(io,'(a,i5)') 'ndata = ', nmx

      do i = 1, nmx
        write(io,'(i5,2x,a,i12)') i, typ_dnam(i), typ_itag(i)
      enddo

      return
      end
!*e* by kamata 2021/10/14
