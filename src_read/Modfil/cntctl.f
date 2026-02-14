! added replace all include files with module files by kamata 2021/08/18
      module cntctl
      implicit none

!----------------------------------------------------------------------
!::New computation
!----------------------------------------------------------------------
!::[MPI_Send  in lnkntl_iniSnd]  cntctl/cntcmp/ (tmntl,emrk)  10/04/21
!::[MPI_Recv  in lnkntl_iniSnd]  cntctl/cntcmp/ (tmntl,emrk)  10/04/21

      real(8) :: tmntl = 0.0_8, dtntl = 0.0_8
      integer :: nntl = 0, imxnt = 0, nhunt = 0, nhfl = 0, nini_cal = 0
     >    , mdl_ntsm = 0
      integer :: kntl = 0, nclnt = 0, itmnt = 0, itmntb = 0
      integer :: mdl_hrt = 0, hrt_dbg = 0

      end module cntctl
