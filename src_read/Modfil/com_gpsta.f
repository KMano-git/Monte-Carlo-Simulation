      module com_gpsta
      implicit none
      integer, parameter :: ndvr = 100
      integer :: nvar = 0
      character(12), dimension(ndvr) :: ctnm = ' '  ! var name
      character(80), dimension(ndvr) :: ctcv = ' '  ! const values
      integer, dimension(ndvr) :: tblrg(ndvr) = 0 ! on/off void region
      end module com_gpsta