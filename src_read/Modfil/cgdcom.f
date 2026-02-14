! added replace all include files with module files by kamata 2021/08/18
      module cgdcom
      implicit none

      integer :: mqd1 = 0, mqd2 = 0, mqh1 = 0, mqh2 = 0, mqs1 = 0
     >    , mqs2 = 0, mqx1 = 0, mqx2 = 0, nqmx = 0  
      integer :: mpax = 0, mpsp = 0, mpw1 = 0, mpw2 = 0, npmx = 0  
      integer :: mpman = 0, mpprv = 0, mpsol = 0

      real(8), allocatable :: grdx(:,:), grdy(:,:), hbr(:,:), hbt(:,:)
     >    , hbz(:,:), psman(:), psprv(:), pssol(:)

!::[MPI_Bcast in mtrdsk]   in cgdcom (grdx,emrk)  10/04/21
!::[MPI_Bcast in ntl_ini]  cgdcom/cgmsh/ (grdx,emrk)  10/04/21
!::[MPI_Send  in lnkntl_ini]  cgdcom/cgmsh/  (grdx,emrk)  10/04/21
!::[MPI_Recv  in lnkntl_ini]  cgdcom/cgmsh/  (grdx,emrk)  10/04/21

      end module cgdcom
