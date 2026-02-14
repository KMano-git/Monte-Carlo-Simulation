! added replace all include files with module files by kamata 2021/08/18
      module celcom
      implicit none
!-----------------------------------------------------------------------
      integer, parameter :: ndel = 2

      character :: el_dnam(ndel)*20 = ' '

!::[MPI_Bcast in elinit]   celcom (el_dnam,emrk)  10/04/21

      integer :: el_atm = 0, el_mol = 0, el_cros(ndel) = 0
     >    , el_rate(ndel) = 0, el_I_10(ndel) = 0, el_I_11(ndel) = 0
     >    , el_I_12(ndel) = 0, el_angl(ndel) = 0, el_sgmx(ndel) = 0
     >    , el_agmi(ndel) = 0

      real(8) :: sgmx(ndel) = 0.0_8, agmi(ndel) = 0.0_8

      integer :: dcl_ipl = 0, dtr_ipl = 0
      real(8) :: pcl_sgv = 0.0_8, ptr_sgv = 0.0_8
      real(8) :: pcl_mvx = 0.0_8, pcl_mvy = 0.0_8, pcl_mvz = 0.0_8
     >    , pcl_mvp = 0.0_8, pcl_eng = 0.0_8
      real(8) :: dcl_mvx = 0.0_8, dcl_mvy = 0.0_8, dcl_mvz = 0.0_8
     >    , dcl_mvp = 0.0_8, dcl_eng = 0.0_8
      real(8) :: ptr_mvx = 0.0_8, ptr_mvy = 0.0_8, ptr_mvz = 0.0_8
     >    , ptr_mvp = 0.0_8, ptr_eng = 0.0_8
      real(8) :: dtr_mvx = 0.0_8, dtr_mvy = 0.0_8, dtr_mvz = 0.0_8
     >    , dtr_mvp = 0.0_8, dtr_eng = 0.0_8

      end module celcom
