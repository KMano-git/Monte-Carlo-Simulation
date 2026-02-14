! added replace all include files with module files by kamata 2021/08/18
      module cplmet
      implicit none

      real(8), allocatable :: gare(:,:), gdsv(:,:,:), gwtm(:,:)
     >    , gwtp(:,:), hare(:,:), hdsp(:,:,:), hdsv(:,:,:), hdxm(:,:)
     >    , hdxp(:,:), hgdx(:,:), hgdy(:,:), hpit(:,:), hvol(:,:)
     >    , hvsb(:,:), hwtm(:,:), hwtp(:,:), romn(:), vlmn(:)
      integer :: set_hdsv = 1 ! set by inppls, (0:hdsv=0,1:calucate hdsv)
      integer :: nosol = 0, noprv = 0, nompl = 0
      integer :: jcdp1 = 0, jcdp2 = 0, jcxp1 = 0, jcxp2 = 0, jcmax = 0
      integer :: icwl1 = 0, icspx = 0, icwl2 = 0, icmps = 0, icmpe = 0
     >    , icaxs = 0
      integer :: itsls = 0, itsle = 0, itpvs = 0, itpve = 0, itmps = 0
     >    , itmpe = 0, itmax = 0

      integer, allocatable :: icel(:,:), icmax(:), icmin(:), jcel(:,:)
     >    , jnxm(:,:), jnxp(:,:), jtmax(:), jtmin(:), kgdx(:,:,:)
     >    , kgdy(:,:,:), kreg(:,:)
      integer :: kce = 0, kcw = 0, kcn = 0, kcs = 0

!::[MPI_Bcast in mtrdsk]   cplmet (vlmn,emrk)    10/04/21
!::[MPI_Bcast in ntl_ini]   cplmet(cmetrc) (vlmn,emrk)  10/04/21
!::[MPI_Send  in lnkntl_ini]  cplmet/cmetrc/ (vlmn,emrk)  10/04/21
!::[MPI_Recv  in lnkntl_ini]  cplmet/cmetrc/ (vlmn,emrk)  10/04/21

      end module cplmet
