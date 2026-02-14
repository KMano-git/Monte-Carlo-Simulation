! added replace all include files with module files by kamata 2021/08/18
      module cimctl
      implicit none

      real(8) :: tmimp = 0.0_8, dtimp = 0.0_8
      real(8) :: ftav = 0.0_8  !KH130408
      integer :: nimp = 0, kimp = 0, icalZ = 0, lstdy = 0
      integer :: svitm(30) = 0, nsvtm = 0, ndirZ = 0
      character :: svdsn(30)*11 = ' '
      character :: cdirZ*20 = ' '

!::[MPI_Send in lnkimp_iniSnd]   cimctl/cimcmp/ (tmimp,emrk) 10/04/21
!::[MPI_Recv in lnkimp_iniSnd]   cimctl/cimcmp/ (tmimp,emrk) 10/04/21

!-----
!xx   integer    svitm(30), nsvtm
!xx   character  svdsn(30)*11
!-----

      integer :: nprg = 0

      end module cimctl
