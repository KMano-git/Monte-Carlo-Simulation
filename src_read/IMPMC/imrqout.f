!***********************************************************************
      subroutine test_imrqout
!***********************************************************************
      use cimcom, only : ltdst
      use cunit,  only : n6
      implicit none
!
      integer :: lp
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer :: iprf, ihst, mhst
      integer :: iprf, ihst
!
      ltdst(1:50) = 0
      ltdst(1:7) =(/1, 5, 25, 100, 500, 505, 0/)
!
      n6 = 6
!
      do lp = 1, 300
      call imrqout(lp,iprf,ihst)
      enddo
!
      do lp = 1, 50000
      call imrqout(lp,iprf,ihst)
      enddo
!
      stop
      end
!
!***********************************************************************
      subroutine imrqout(lptm,iprf,ihst)
!***********************************************************************
      use cimcom, only : ltdst
      use cunit,  only : n6
      implicit none
!
!::local common
      integer :: kskp, ksk2, mhst
      save    :: kskp, mhst
!
!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/28
!ik   integer :: lptm, iprf, ihst
      integer, intent(in)  :: lptm
      integer, intent(out) :: iprf, ihst
!
!::local variables
      integer :: i
!
!xx   write(n6,'(2x,"ltdst =",10i5)') ltdst(1:10)
!
      if( lptm.eq.1 ) then
        kskp = 1
        mhst = 0
      endif
!
!::profile data
      iprf = 0
      do i = 1, 50
        if( ltdst(i).le.0 ) exit
        if( lptm.eq.ltdst(i) ) then
           iprf = 1
           exit
        endif
      enddo
!
!::hitory data
      ihst = 0
      if( dfloat(lptm)/dfloat(kskp).gt.10.0d0 ) kskp = kskp*10
      ksk2 = kskp/2
      if( ksk2.le.1 ) ksk2 = 1
      if( mod(lptm,ksk2).eq.0 ) ihst = 1
!
      if( iprf.gt.0 .or. ihst.gt.0 ) then
        if( ihst.gt.0 ) mhst = mhst + 1
        write(n6,'(2x,"=== imrqout","  lt =",i8,"  iprf/ihst =",2i4,
     >     "  mhst =",i5)')
     >      lptm, iprf, ihst, mhst
      endif
!
      return
      end
