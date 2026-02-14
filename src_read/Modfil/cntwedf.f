! added replace all include files with module files by kamata 2021/08/18
      module cntwedf
      implicit none
!
      integer, parameter :: nwmx = 100 ! max wall cell number
      integer, parameter :: nsmx = 2 ! index for atom or molecular
      real(8), parameter :: deg_iw(5)
     >  = (/0.0d0, 45.0d0, 90.d0, 135.0d0, 180.0d0/)
      integer :: iw_edf(nwmx) = 0
!
! log version
      real(8), parameter :: enemin = 0.1d0   ! minimum energy 
      real(8), parameter :: enemax = 10.0d2  ! maximum energy
      real(8), parameter :: d_ene  = 0.1d0   ! energy spectrum width in log10(ene)
      integer, parameter :: nene_mx = 41 ! energy interval
      double precision :: wedf(nene_mx, nwmx, nsmx) = 0.0_8
      double precision :: xwedf(nene_mx, nwmx, nsmx) = 0.0_8
!
! liner version
      real(8), parameter :: enemin_line =  0.1d0  ! minimum energy(for liner)
      real(8), parameter :: enemax_line = 39.9d0  ! maximum energy(for liner)
      real(8), parameter :: d_ene_line  =  0.2d0  ! energy spectrum width(for liner)
      integer, parameter :: nene_mx_line = 200 ! energy interval
      real(8) :: wedf_line(nene_mx_line, nwmx, nsmx) = 0.0_8
      real(8) :: xwedf_line(nene_mx_line, nwmx, nsmx) = 0.0_8

      end module cntwedf
