!**********************************************************************
      subroutine dtypset(dtyp,mtyp)
!**********************************************************************
      use mod_sizedef, only : lnnam
      implicit none

      character(*), intent(in) :: dtyp
      integer, intent(in) :: mtyp ! 1:pre-Bcst,2:post-bcast,3:Send,4:Recv
      character(lnnam) :: ctyp 

      ctyp = trim(dtyp)

!      if(ctyp == "IMP_1")

!      if(mtyp == 1)then     ! pre-Bcast
!      if(mtyp == 2)then     ! post-Bcast
!      if(mtyp == 3)then     ! Send
!      if(mtyp == 4)then     ! Recv
!      endif

      return
      end
