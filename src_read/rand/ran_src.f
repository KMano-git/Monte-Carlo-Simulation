!%%mkspt +
!
!    How to set random-seed
!
!::for particle
!
!      integer iseed(1000), jseed
!
!      nptl = 1000
!      call rnseed(nptl,iseed)
!
!      do n = 1, nptl
!      jseed = iseed(n)
!      call srandom(jseed)
!
!      xran = random(0)
!
!
!***********************************************************************
      real*8 function random( ix )
!***********************************************************************
!
!    Note   0.0 <= random < 1.0    ( ran = 0.0 with prob = 1/10E9 )
!           when optimize, make an error of calcualtion [im+1]
!
!-----------------------------------------------------------------------
      use crand, only : ia, ic, iur, iur0, rm, ur
      implicit none
!
      integer, intent(in) :: ix ! dummy
!
      iur0 = iur
      iur = iur * ia + ic
      iur = ibclr( iur, 31 )
!
      ur = iur * rm
      if( iur.eq.0 ) ur = 1.0d-15
!
      random = ur
!
      return
      end
!
!**********************************************************************
      subroutine srandom(iseed)
!**********************************************************************
      use crand, only : iur
      implicit none
!
      integer, intent(in) :: iseed
!
      iur = iseed
!
      return
      end
!
!**********************************************************************
      subroutine rnseed(n,iseed)
!**********************************************************************
      use crand, only : im, rm
      use cunit, only : n6
      implicit none
!
      integer, intent(in)  :: n
      integer, intent(out) :: iseed(n)
!
!::local variables
      integer  i, iran, jseed
      real*8   xran, frac, rnlast
! function
      real(8)  random
!
      write(n6,'(/2x,"*** rnseed ***  n =",i5)') n
!
      do i = 1, n
      frac = (dfloat(i-1)+random(0))/dfloat(n)
      iran = int(dfloat(im)*frac)
      xran = dfloat(iran)*rm
      iseed(i) = iran
      if( iran.le.0 .or. iran.ge.im ) then
      call wexit("rnseed","incorrect iran")
      endif
      enddo
!
!::debug write
      do i = 1, n
      jseed = iseed(i)
      call srandom(jseed)
      xran = random(0)
      if( mod(i,100).eq.0 .or. i.le.5 .or. i.ge.n-5 ) then
      write(n6,'(2x,i6,i14,2f12.9)') i,iseed(i),xran,rnlast(0)
      endif
      enddo
      return
      end
!
!***********************************************************************
      real*8 function rnlast( ix )
!***********************************************************************
      use crand, only : ur
      implicit none
!
      integer, intent(in) :: ix ! dummy
!
      rnlast = ur
!
      return
      end
!
!***********************************************************************
      subroutine rand_last( isd0, rand )
!***********************************************************************
      use crand, only : iur0, ur
      implicit none
!
      integer, intent(out) :: isd0
      real(8), intent(out) :: rand
!
      isd0 = iur0
      rand = ur
!
      return
      end
