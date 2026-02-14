! added replace all include files with module files by kamata 2021/08/18
      module cntrfl
      implicit none

!::[MPI_Bcast in refle]  cntrfl/cmrefl/ (rfeng,mpe)  10/04/21
      integer, parameter :: npe = 48
      real(8) :: rfeng(npe) = 0.0_8, rfelg(npe) = 0.0_8
     >    , rfrn(npe) = 0.0_8, rfer(npe) = 0.0_8, dlge = 0.0_8
      integer :: mpe = 0

      end module cntrfl
