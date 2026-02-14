!**********************************************************************
      subroutine dtypset(dtyp,mtyp)
!**********************************************************************
      implicit none

      character(*), intent(in) :: dtyp ! dummy
      integer,      intent(in) :: mtyp ! 1:pre-Bcst,2:post-bcast,3:Send,4:Recv ! dummy
! deleted 1 line organize local variables and include files by kamata 2021/06/22
!ik   character(lnnam) :: ctyp

! deleted 1 line organize local variables and include files by kamata 2021/06/22
!ik   ctyp = trim(dtyp)

!      if(ctyp == "IMP_1")

!      if(mtyp == 1)then     ! pre-Bcast
!      if(mtyp == 2)then     ! post-Bcast
!      if(mtyp == 3)then     ! Send
!      if(mtyp == 4)then     ! Recv
!      endif

      return
      end
