! added replace all include files with module files by kamata 2021/08/18
      module cplpst
      implicit none

!----------------------------------------------------------------------
!::radial position for post list
!----------------------------------------------------------------------
      real(8), allocatable :: rhf_odp(:), rhf_idp(:), rhf_omd(:)
     >    , rhf_imd(:)
      real(8), allocatable :: rhf_uodp(:),rhf_uidp(:)

!----------------------------------------------------------------------
!::heat flux to divertor plate
!----------------------------------------------------------------------
      real(8) :: tflsmi(10) = 0.0_8, tflsmo(10) = 0.0_8
      real(8), allocatable :: fbdpi(:,:), fbdpo(:,:), fldpi(:,:)
     >    , fldpo(:,:), flsmi(:,:), flsmo(:,:), farei(:), fareo(:)

!----------------------------------------------------------------------
!::heat load of radiation onto the wall
!----------------------------------------------------------------------
      integer :: nphtn = 0
      real(8), allocatable :: dwrad(:,:)

      end module cplpst
