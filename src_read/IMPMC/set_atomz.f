!***********************************************************************
      subroutine test_set_atomz
!***********************************************************************
      implicit none
      character(80) :: drnm
!
      drnm = "/home/g9/a094069/sonicV2/adas/Ar_Y89"
      call set_atomz(drnm)
      return
      end
!
!***********************************************************************
      subroutine set_atomz(drnm)
!***********************************************************************
!
!       drnm : full path of dir which contains adas-data
!            : specified in input data(uimimp)
!
!----------------------------------------------------------------------
      use catcom, only : aionz, catmz, dratm, natmz
      use cimcom, only : azmas, ismax
      use cimden, only : eipot
      use csize,  only : ndmis
      use csonic, only : mstop
      use cunit,  only : lmspe, lnope, mywld, n6, lmype
      use mpi!,    only : mpi_bcast, mpi_character
      implicit none

      character(*), intent(in) :: drnm
!
!::local variables
      integer :: ierr, m1, m2, ityp, iz
!
!::atomic data for C, Ne, Ar, Kr, Xe
      character(2),dimension(9) :: tbcatmz
      integer,     dimension(9) :: tbnatmz
      real(8),     dimension(9) :: tbaionz
      real(8), dimension(0:6)   :: tbeipot_C
      real(8), dimension(0:10)  :: tbeipot_Ne
      real(8), dimension(0:18)  :: tbeipot_Ar
      real(8), dimension(0:36)  :: tbeipot_Kr
      real(8), dimension(0:54)  :: tbeipot_Xe
      real(8), dimension(0:2)   :: tbeipot_He
      real(8), dimension(0:4)   :: tbeipot_Be
      real(8), dimension(0:74)  :: tbeipot_W
      real(8), dimension(0:7)   :: tbeipot_N
!
      write(n6,'(/2x,"*** set_atomz *** ",a)')drnm
!
!::dratm
      dratm(1) = trim(drnm)
      if( lnope > 1 ) then
      call MPI_Bcast( dratm(1), 80, MPI_CHARACTER, lmspe, mywld, ierr )
      write(n6,'(2x,"MPI_Bcast  lnope =",i5,"  lmype =",i5,
     >  "  lmspe =",i3)') lnope, lmype, lmspe
      else
      write(n6,'(2x,"No pass MPI_Bcast")')
      endif

!
      m1 = index( dratm(1), "/", .true. )
      m2 = index( dratm(1), "_", .true. )
      catmz(1) = dratm(1)(m1+1:m2-1)
!
      write(n6,'(2x,"catmz = ",a,"  dratm = ",a)')
     >     catmz(1), trim(dratm(1))
      call flush(n6)
!
!
!::atomic data
      tbcatmz = (/ "C ",  "Ne",    "Ar",    "Kr",   "Xe",
     &            "He",  "Be",     "W ",     "N " /)
      tbnatmz = (/   6,    10,      18,      36,     54,
     &               2,     4,      74,       7 /)
      tbaionz = (/ 12.011d0,  20.179d0,  39.95d0,  83.80d0,  131.3d0,
     >             4.003d0,   9.012d0,  183.8d0, 14.01d0 /)
!
!::ionization potential
      tbeipot_C = (/
     >    11.0d0,   24.0d0,   48.0d0,   64.0d0,   392.0d0,
     >   490.0d0,    0.0d0 /)
!
      tbeipot_Ne = (/
     >    21.57d0,  40.98d0,  63.48d0,  97.17d0,  126.26d0,
     >   158.00d0, 207.37d0, 239.21d0, 1196.4d0,  1362.8d0,  0.0d0 /)
!
      tbeipot_Ar = (/
     >    15.760,   27.630,  40.741,  59.810,  75.023,
     >    91.010,   124.32,  143.46,  422.45,  478.69,
     >    538.96,   618.26,  686.11,  755.75,  854.78,
     >    918.03,   4120.9,  4426.2,  0.0000 /)
!
      tbeipot_Kr = (/
     >  1.400d+01,  2.436d+01,  3.695d+01,  5.249d+01,  6.469d+01,
     >  7.849d+01,  1.110d+02,  1.258d+02,  2.310d+02,  2.682d+02,
     >  3.082d+02,  3.501d+02,  3.909d+02,  4.466d+02,  4.918d+02,
     >  5.407d+02,  5.915d+02,  6.410d+02,  7.859d+02,  8.330d+02,
     >  8.838d+02,  9.366d+02,  9.976d+02,  1.051d+03,  1.151d+03,
     >  1.205d+03,  2.928d+03,  3.070d+03,  3.227d+03,  3.381d+03,
     >  3.594d+03,  3.760d+03,  3.966d+03,  4.108d+03,  1.730d+04,
     >  1.794d+04,  0.0d0 /)
!
      tbeipot_Xe = (/
     >  1.213d+01,  2.098d+01,  3.105d+01,  4.090d+01,  5.414d+01,
     >  6.670d+01,  9.160d+01,  1.060d+02,  1.798d+02,  2.017d+02,
     >  2.290d+02,  2.635d+02,  2.944d+02,  3.253d+02,  3.583d+02,
     >  3.896d+02,  4.209d+02,  4.522d+02,  5.725d+02,  6.077d+02,
     >  6.429d+02,  6.781d+02,  7.260d+02,  7.624d+02,  8.527d+02,
     >  8.570d+02,  1.490d+03,  1.491d+03,  1.587d+03,  1.684d+03,
     >  1.781d+03,  1.877d+03,  1.987d+03,  2.085d+03,  2.183d+03,
     >  2.281d+03,  2.548d+03,  2.637d+03,  2.726d+03,  2.814d+03,
     >  3.001d+03,  3.093d+03,  3.296d+03,  3.334d+03,  7.224d+03,
     >  7.491d+03,  7.758d+03,  8.024d+03,  8.617d+03,  8.899d+03,
     >  9.607d+03,  9.813d+03,  4.027d+04,  4.130d+04,  0.0d0 /)

      tbeipot_He = (/ 2.459d1,    5.442d1,    0.0d0 /)

      tbeipot_Be = (/9.323d0,     1.821d1,    1.539d2,    2.177d2,
     >                 0.0d0 /)

      tbeipot_W = (/
     >  7.864D+00,  1.637D+01,  2.600D+01,  3.820D+01,  5.160D+01,
     >  6.477D+01,  1.220D+02,  1.412D+02,  1.602D+02,  1.790D+02,
     >  2.089D+02,  2.316D+02,  2.582D+02,  2.907D+02,  3.253D+02,
     >  3.619D+02,  3.879D+02,  4.207D+02,  4.621D+02,  5.026D+02,
     >  5.434D+02,  5.945D+02,  6.407D+02,  6.856D+02,  7.341D+02,
     >  7.844D+02,  8.334D+02,  8.814D+02,  1.132D+03,  1.180D+03,
     >  1.230D+03,  1.283D+03,  1.335D+03,  1.387D+03,  1.460D+03,
     >  1.512D+03,  1.569D+03,  1.622D+03,  1.830D+03,  1.883D+03,
     >  1.941D+03,  1.995D+03,  2.149D+03,  2.210D+03,  2.355D+03,
     >  2.414D+03,  4.057D+03,  4.180D+03,  4.309D+03,  4.446D+03,
     >  4.578D+03,  4.709D+03,  4.927D+03,  5.063D+03,  5.209D+03,
     >  5.348D+03,  5.719D+03,  5.840D+03,  5.970D+03,  6.093D+03,
     >  6.596D+03,  6.735D+03,  7.000D+03,  7.130D+03,  1.557D+04,
     >  1.590D+04,  1.625D+04,  1.659D+04,  1.848D+04,  1.887D+04,
     >  1.936D+04,  1.969D+04,  7.918D+04,  8.076D+04,    0.0D+00 /)

      tbeipot_N = (/  1.453d1,     2.96d1,    4.745d1,    7.747d1,
     >                9.789d1,    5.521d2,    6.671d2,      0.0d0 /)

!
!::atom to be treated in this group PE
      ityp = 0
      if( trim(catmz(1)).eq."C"  ) ityp = 1
      if( trim(catmz(1)).eq."Ne" ) ityp = 2
      if( trim(catmz(1)).eq."Ar" ) ityp = 3
      if( trim(catmz(1)).eq."Kr" ) ityp = 4
      if( trim(catmz(1)).eq."Xe" ) ityp = 5
      if( trim(catmz(1)).eq."He" ) ityp = 6
      if( trim(catmz(1)).eq."Be" ) ityp = 7
      if( trim(catmz(1)).eq."W"  ) ityp = 8
      if( trim(catmz(1)).eq."N"  ) ityp = 9
      if( ityp.eq.0 ) then
      write(n6,'(2x,"No found catmz in tbcatmz  ",a)') trim(catmz(1))
      mstop = 1
      return
      endif
!
!::common variabels in catcom
      catmz(1) = tbcatmz(ityp)
      natmz(1) = min(tbnatmz(ityp), ndmis)
      aionz(1) = tbaionz(ityp)
!
!::common variables in cimcom
      ismax = natmz(1)
      azmas = aionz(1)
!
      eipot(0:ndmis) = 0.0d0
      if( ityp.eq.1 ) eipot(0:natmz(1)) = tbeipot_C(0:natmz(1))
      if( ityp.eq.2 ) eipot(0:natmz(1)) = tbeipot_Ne(0:natmz(1))
      if( ityp.eq.3 ) eipot(0:natmz(1)) = tbeipot_Ar(0:natmz(1))
      if( ityp.eq.4 ) eipot(0:natmz(1)) = tbeipot_Kr(0:natmz(1))
      if( ityp.eq.5 ) eipot(0:natmz(1)) = tbeipot_Xe(0:natmz(1))
      if( ityp.eq.6 ) eipot(0:natmz(1)) = tbeipot_He(0:natmz(1))
      if( ityp.eq.7 ) eipot(0:natmz(1)) = tbeipot_Be(0:natmz(1))
      if( ityp.eq.8 ) eipot(0:natmz(1)) = tbeipot_W(0:natmz(1))
      if( ityp.eq.9 ) eipot(0:natmz(1)) = tbeipot_N(0:natmz(1))
!
      write(n6,'(2x,"catmz = ",a,"  natmz =",i3,"  aionz =",f8.3)')
     >   trim(catmz(1)), natmz(1), aionz(1)
      write(n6,'(2x,"eipot = ",10f10.3)') (eipot(iz),iz=0,natmz(1))
      call flush(n6)

      if( eipot(1).eq.0.0d0 ) then
      write(n6,'(2x,"error occured at set_atomz ",2x,"eipot = 0.0")')
      mstop = 1
      call flush(n6)
      return
      endif
!
      return
      end
