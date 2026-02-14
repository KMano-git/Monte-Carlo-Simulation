! added replace all include files with module files by kamata 2021/08/18
! data copy of vmwksc ( monte/inc/cntscg /cntscr_tot/ )
      subroutine setvmwksc( kk, worksc )
      use cntcom, only : ndrmx
      use cntscg, only : md => nwkmp_sc, tsn_mc, tsn_mt, tsp_mc, tsp_mt
     >    , twe_mc, twe_mt, twi_mc, twi_mt
      implicit none
! arguments
      integer, intent(in)    :: kk
      real(8), intent(inout) :: worksc(md)
! kk     : flag to set data, = 1 : worksc to module variables
!                            = 2 : module variables to worksc
! worksc : work area

! local variables
      integer    i0, i1, i2, i3, i4, i5, i6, i7, nd

! initial set
      nd = ndrmx * 10

! set start posiotion of each variable
      i0 = 1        ! tsn_mt
      i1 = i0 + nd  ! tsp_mt
      i2 = i1 + nd  ! twe_mt
      i3 = i2 + nd  ! twi_mt
      i4 = i3 + nd  ! tsn_mc
      i5 = i4 + nd  ! tsp_mc
      i6 = i5 + nd  ! twe_mc
      i7 = i6 + nd  ! twi_mc

      select case ( kk )
      case( 1 )
! worksc to module variables
        call datcp( worksc(i0), tsn_mt, nd )
        call datcp( worksc(i1), tsp_mt, nd )
        call datcp( worksc(i2), twe_mt, nd )
        call datcp( worksc(i3), twi_mt, nd )
        call datcp( worksc(i4), tsn_mc, nd )
        call datcp( worksc(i5), tsp_mc, nd )
        call datcp( worksc(i6), twe_mc, nd )
        call datcp( worksc(i7), twi_mc, nd )
      case( 2 )
! module variables to worksc
        call datcp( tsn_mt, worksc(i0), nd )
        call datcp( tsp_mt, worksc(i1), nd )
        call datcp( twe_mt, worksc(i2), nd )
        call datcp( twi_mt, worksc(i3), nd )
        call datcp( tsn_mc, worksc(i4), nd )
        call datcp( tsp_mc, worksc(i5), nd )
        call datcp( twe_mc, worksc(i6), nd )
        call datcp( twi_mc, worksc(i7), nd )
      end select

      return
      end
