!**********************************************************************
      subroutine plnflx(anflx)
!**********************************************************************
      use cntmnt, only : mfmax, visrc, vflux
      use csize,  only : ndmfl
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real(8) :: anflx(ndmfl)
      real(8), intent(out) :: anflx(ndmfl)

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: m, iflx, isrc
      integer :: iflx, isrc

      anflx(1:ndmfl) = 0.0d0

! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   do m = 1, mfmax
!ik     iflx = m
      do iflx = 1, mfmax
        isrc = visrc(iflx)
        if( isrc > 0 )  anflx(iflx) = vflux(isrc)
      enddo

!::vli = vol   ! see ntmset.f
       anflx(8) = anflx(7)

!xx      write(n6,'(2x,"itim =",i7,"  plnflx(snflx) =",1p8e14.6)')
!xx     >  itim, anflx(1:mfmax)

      return
      end
