! added replace all include files with module files by kamata 2021/08/18
! data copy of vwork2
!  IMPMC/inc/cimcom /cimcom_12/ or IMPMC_TD/inc/cimcom /cimcom_12/ 
      subroutine setvwork2( kk, work2 )
      use cimcom, only : md => nwkmp_im2, wexh, whit, wpemt, wtemt
      use csize,  only : ndis => ndmis
      implicit none
! arguments
      integer, intent(in)    :: kk
      real(8), intent(inout) :: work2(md)
! kk    : flag to set data, = 1 : work2 to module variables
!                           = 2 : module variables to work2
! work2 : work area

! local variables
      integer    i0, i1, i2, i3, nd1, nd2

! initial set
      nd2 = ( ndis + 1 ) * 30
      nd1 = nd2 * 4

! set start posiotion of each variable
      i0 = 1        ! wtemt
      i1 = i0 + 1   ! wpemt
      i2 = i1 + 1   ! whit
      i3 = i2 + nd1 ! wexh
      
      select case ( kk )
      case( 1 )
! work2 to module variables
        wtemt = work2(i0) 
        wpemt = work2(i1)

        call datcp( work2(i2), whit, nd1 )
        call datcp( work2(i3), wexh, nd2 )
      case( 2 )
! module variables to work2
        work2(i0) = wtemt 
        work2(i1) = wpemt

        call datcp( whit, work2(i2), nd1 )
        call datcp( wexh, work2(i3), nd2 )
      end select

      return
      end
