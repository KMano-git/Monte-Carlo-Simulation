! added replace all include files with module files by kamata 2021/08/18
! data copy of vwork4 ( IMPMC/inc/cimcom /cimcom_12c/ )
      subroutine setvwork4( kk, work4 )
      use cimcom, only : ndext, md => nwkmp_im4, wvaci, wvacn
      use csize,  only : ndgt
      implicit none
! arguments
      integer, intent(in)    :: kk
      real(8), intent(inout) :: work4(md)
! kk    : flag to set data, = 1 : work4 to module variables
!                           = 2 : module variables to work4
! work4 : work area

! local variables
      integer    i0, i1, nd

! initial set
      nd  = ( ndgt + 1 ) * ( ndext + 1 )

! set start posiotion of each variable
      i0 = 1       ! wvacn
      i1 = i0 + nd ! wvaci
      
      select case ( kk )
      case( 1 )
! work4 to module variables
        call datcp( work4(i0), wvacn, nd )
        call datcp( work4(i1), wvaci, nd )
      case( 2 )
! module variables to work4
        call datcp( wvacn, work4(i0), nd )
        call datcp( wvaci, work4(i1), nd )
      end select

      return
      end
