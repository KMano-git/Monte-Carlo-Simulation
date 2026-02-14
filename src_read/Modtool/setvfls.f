! added replace all include files with module files by kamata 2021/08/18
! data copy of vfls ( IMPMC_com/inc/cimfls /cimedg_2/ )
      subroutine setvfls( kk, wfls, md )
      use cimfls, only : fls_dnz, fls_flxi, fls_flxo
      use csize,  only : ndis => ndmis
      implicit none
! arguments
      integer, intent(in)    :: kk, md
      real(8), intent(inout) :: wfls(md)
! kk   : flag to set data, = 1 : wfls to module variables
!                          = 2 : module variables to wfls
! md   : dimension size
! wfls : work area

! local variables
      integer    i0, i1, i2, nd

! initial set
      nd  = ( ndis + 1 ) * 5

! set start posiotion of each variable
      i0 = 1       ! fls_flxi
      i1 = i0 + nd ! fls_flxo
      i2 = i1 + nd ! fls_dnz
      
      select case ( kk )
      case( 1 )
! wfls to module variables
        call datcp( wfls(i0), fls_flxi, nd )
        call datcp( wfls(i1), fls_flxo, nd )
        call datcp( wfls(i2), fls_dnz,  nd )
      case( 2 )
! module variables to wfls
        call datcp( fls_flxi, wfls(i0), nd )
        call datcp( fls_flxo, wfls(i1), nd )
        call datcp( fls_dnz,  wfls(i2), nd )
      end select

      return
      end
