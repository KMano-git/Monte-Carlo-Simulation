! added replace all include files with module files by kamata 2021/08/18
      module catmhy
      implicit none

!----------------------------------------------------------------------
!::atomic data of degas2-code  
!       sghio    H0 + e -> H+ + e + e   ionization
!       sghrc    H+ + e + e -> H0 + e   recombination
!----------------------------------------------------------------------
      integer, parameter :: ndtm = 60, nddn = 15

      real(8) :: stemin = 0.0_8, stemax = 0.0_8, stedlt = 0.0_8
     >    , snemin = 0.0_8, snemax = 0.0_8, snedlt = 0.0_8
     >    , ste(ndtm) = 0.0_8, sne(nddn) = 0.0_8
     >    , sghio(ndtm,nddn) = 0.0_8, elhio(ndtm,nddn) = 0.0_8
     >    , emhio(ndtm,nddn) = 0.0_8, sghrc(ndtm,nddn) = 0.0_8
     >    , elhrc(ndtm,nddn) = 0.0_8, emhrc(ndtm,nddn) = 0.0_8

!::[MPI_Bcast in atomhy]    catmhy  (stemin-emhrc)  10/04/21

      end module catmhy
