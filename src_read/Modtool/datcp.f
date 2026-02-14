! added replace all include files with module files by kamata 2021/08/18
! copy of array data ( for array data with different dimensions )
      subroutine datcp( dsour, ddest, ndt )
      implicit none
! arguments
      integer, intent(in)  :: ndt
      real(8), intent(in)  :: dsour(ndt)
      real(8), intent(out) :: ddest(ndt)
! ddest : copy destination array
! dsour : copy source array 
! ndt   : number of data

      ddest(1:ndt) = dsour(1:ndt)

      return
      end
