!**********************************************************************
      subroutine ntl_skip(kcal)
!**********************************************************************
      use cntctl, only : itmnt, nntl
      use csonic, only : itend, itim, lntl
      use cunit,  only : n6
      implicit none

!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: kcal
      integer, intent(out) :: kcal

!::local variables
      integer, save :: itcal = -1

!::lntl modelling of neutal code
      if( lntl.eq.0 ) goto 100
      if( lntl.ne.3 ) goto 100

!::flag
      kcal = 0
      if( itcal == -1 )   itcal = itim
      if( itim == itcal ) kcal = 1
      if( itim == itend ) kcal = 1
      if( itim == itend .and. mod(itend,100) /= 0) kcal = 0

      itmnt = itcal
      if( kcal == 0 ) goto 100
      itmnt = itim
      itcal = itcal + nntl

 100  continue
      if( kcal > 0 ) then
      write(n6,'(2x,"ntl_skip  itim/itmnt =",2i6,"  nntl/kcal =",
     >  2i4)') itim, itmnt, nntl, kcal
      endif

      return
      end
