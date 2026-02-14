      module com_gmsiz
      implicit none
      
      integer :: mfig = 0   !  <== lfigxx
      real*4 :: rmi = 0.0_4, rmx = 0.0_4, zmi = 0.0_4, zmx = 0.0_4
     > , xln = 0.0_4, yln = 0.0_4
      integer :: ngsz = 0
      real(4) :: pslgx = 0.0_4, pslgy = 0.0_4
      character(80),dimension(10) :: cgsz = ' '
      character(5) ,dimension(10) :: cgnm = ' '
      character :: cplsb*80  = ' '

      end module com_gmsiz