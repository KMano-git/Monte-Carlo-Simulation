!***********************************************************************
      subroutine set_tabl
!***********************************************************************
      use catcom, only : catmz
      use cimcns, only : fgar, fgau, gaus0, gaus1, ndgar, ndghf, ngar
     >    , ngau
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer imax, j, i
      real*8  dx, x, y, z
!
      write(n6,'(/2x,"*** set_tabl ***")')
      write(n6,'(2x,"make tables of normal random number")')
!
!::gauss
      imax = ndghf
      ngau = ndghf*2
      fgau = ngau
      j    = 2*ndghf+1
      dx = 0.49977d0/dfloat(imax)
!
      do i = 1, imax
      j = j-1
      x = dx*dfloat(i)
      y = sqrt(2.0d0*dlog(1.0d0/(0.5d0-x)))
      z = y-(2.30753d0+0.27061d0*y)/(1.0d0+y*(0.99229d0+0.04481d0*y))
      gaus0(i) =  z
      gaus0(j) = -z
      enddo
      gaus0(ngau+1) = gaus0(ngau)
!
!::log
      ngar = ndgar-1
      fgar = ngar
      do i = 1, ngar
      x = -dlog(1.0d0-dfloat(i)/dfloat(ngar+1))
      gaus1(i) = sqrt(x)
      enddo
      gaus1(ngar+1) = gaus1(ngar)
!
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   if(catmz.eq."C") then
      if(catmz(1) == "C") then
!::     tompson model for sputtering energy
! deleted 1 line for bug ( Es not defined ) by kamata 2021/08/18
!ik     call set_tomp
!::     physical sputtering
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik     call set_sputphys(catmz)
        call set_sputphys( catmz(1) )
        call test_sputphy
! added 1 line for bug ( Es not defined ) by kamata 2021/08/18
        call set_tomp
      endif
!
      return
      end
