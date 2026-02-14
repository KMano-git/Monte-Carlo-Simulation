! added replace all include files with module files by kamata 2021/08/18
      module com_eqdat
      implicit none

!::eqdisk
!::[MPI_Bcast in eqdisk]     com_eqdat (dr,emrk)    10/04/21
      integer, parameter :: ndr = 400, ndz = 513
      integer :: nr = 0, nz = 0
      real(8) :: rg(ndr) = 0.0_8, zg(ndz) = 0.0_8, psi(ndr*ndz) = 0.0_8
      real(8) :: rbt0 = 0.0_8, raxs = 0.0_8, zaxs = 0.0_8, paxs = 0.0_8
     >    , rsep = 0.0_8, zsep = 0.0_8, psep = 0.0_8
      real(8) :: dr = 0.0_8, dz = 0.0_8
      real(8) :: rmin = 0.0_8, rmax = 0.0_8, zmin = 0.0_8, zmax = 0.0_8
     >    , dr0 = 0.0_8, dz0 = 0.0_8

!::unit b
      real(8) :: ubx(ndr*ndz) = 0.0_8, uby(ndr*ndz) = 0.0_8
     >    , ubz(ndr*ndz) = 0.0_8

!::eqinit
      real(8) :: rax = 0.0_8, zax = 0.0_8, pax = 0.0_8, rxp(2) = 0.0_8
     >    , zxp(2) = 0.0_8, pxp(2) = 0.0_8
      real(8) :: psepi = 0.0_8, psepo = 0.0_8
      integer :: iax = 0, jax = 0, nax = 0, ixp(2) = 0, jxp(2) = 0
     >    , lxp(2) = 0, nxp = 0

      real(8) :: pmin = 0.0_8, pdlt = 0.0_8
      integer :: neqc = 0

      end module com_eqdat
