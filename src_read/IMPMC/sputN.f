!***********************************************************************
      real(8) function sputphy(e0,alf,katm)
!***********************************************************************
!  31 July 2017 K.Hoshino
!
!  the physical sputtering yield for  normal incidence
!
!   Ref. W.Eckstein, et alyy., Tech. Rep. IPP 9/82, IPP Garching (1993)
!
!   e0   : incident energy
!   alf  : incident angle
!   katm : incident particle
!           1:H 2:D 3:T 4:He 5:Be 6:C 7:W 8:Ne 9:Ar 10:Kr 11:Xe
!
!-----------------------------------------------------------------------
      use clocal, only : etf, eth, qy
      implicit none
!
!::argument
      real(8), intent(in) :: e0
      real(8), intent(in) :: alf  ! dummy
      integer, intent(in) :: katm
!
!::local variables
      real*8 epp, eps, sn, gep
!
      sputphy = 0.0
!
      if( eth(katm) .le. 0.0d0) return
!
      epp = eth(katm)/e0
      if( epp.ge.1.0d0 ) return
!
      eps = e0/etf(katm)
!
      sn = 0.5d0 * dlog(1.0d0 + 1.2288d0*eps)
     &      / ( eps + 0.1728d0*sqrt(eps) + 0.008d0*eps**0.1504d0 )
!
      gep = (1.0d0-epp**(2.0d0/3.0d0))*(1.0d0-epp)**2.0d0
!
      sputphy  = qy(katm)*sn*gep
!
! Angular dependence under development
!      amax = cpi*0.5d0 - al*n**(1.0d0/3.0d0) * ( 2.0d0 *eps*sqrt(Es/(gk*E )**(-0.5d0)
!      amax =    a1     -         a2          * ( 2      eps*( a3 * sqrt(1/e) )**(-0.5)

!      yd = 1.0d0 / cos(alf)**fy * exp(fy * (1.0d0-1/cos(alf)) cos(amax))
      return
      end
!
!***********************************************************************
      subroutine set_sputphys(catmz)
!***********************************************************************
! Ref. W.Eckstein, et alyy., Tech. Rep. IPP 9/82, IPP Garching (1993)
! Etf is calculated by Eq.(7)
!-----------------------------------------------------------------------
      use clocal, only : csptin, es, etf, eth, fy, qy
      implicit none
!
      character(*), intent(in) :: catmz
!
      integer i
      integer Z(11), Z1, Z2
      real*8  M(11), M1, M2

      csptin = (/  "H ",     "D ",     "T ",    "He",
     &             "Be",     "C ",     "W ",
     &             "Ne",    "Ar",    "Kr",    "Xe" /)

      Z = (/         1,       1,       1,       2,
     &                4,       6,      74,
     &               10,      18,      36,      54  /)

      M = (/   1.008d0, 2.014d0, 3.016d0, 4.003d0,
     &          9.012d0, 12.01d0, 183.8d0,
     &          20.18d0, 39.95d0, 83.80d0, 131.3d0  /)

!-----------------------------------------------------------------------
!::target  carbon
!-----------------------------------------------------------------------
      if( trim(catmz).eq."C"  ) then
      Z2 = Z(6)
      M2 = M(6)
      Es = 7.42d0
!                    H cal        D cal        T cal       He cal
!                   Be cal        C cal        W cal
!                   Ne exp       Ar exp       Kr exp       Xe exp
!
      eth = (/      27.3d0,      24.3d0,      22.4d0,      25.4d0,
     &               0.0d0,      44.0d0,      0.0d0 ,
     &              69.6d0,      71.1d0,      71.3d0,     230.0d0/)

      qy  = (/    0.0271d0,    0.0601d0,    0.0684d0,     0.169d0,
     &               0.0d0,     0.620d0,       0.0d0,
     &              1.60d0,      3.06d0,      4.04d0,      5.90d0/)
!
!-----------------------------------------------------------------------
!:: target Tungsten
!-----------------------------------------------------------------------
      elseif( trim(catmz).eq."W"  ) then
      Z2 = Z(7)
      M2 = M(7)
      Es = 8.68d0
!                    H exp        D exp        T cal       He exp
!                   Be            C exp        W exp
!                   Ne exp       Ar exp       Kr exp       Xe exp
!
      eth = (/     429.0d0,     178.0d0,     129.0d0,      107.0d0,
     &               0.0d0,      27.6d0,      59.0d0,
     &              26.0d0,      36.5d0,      36.0d0,       42.8d0/)

      qy  = (/     0.007d0,    0.0179d0,    0.0654d0,      0.110d0,
     &               0.0d0,     0.780d0,      30.0d0,
     &              2.67d0,      7.60d0,      16.2d0,       22.1d0/)
!
!-----------------------------------------------------------------------
!::target  b4c
!-----------------------------------------------------------------------
      elseif( trim(catmz).eq."Be"  ) then
      Z2 = Z(5)
      M2 = M(5)
      Es = 3.38d0
!                    H exp        D exp        T cal       He exp
!                   Be            C exp        W exp
!                   Ne exp       Ar exp       Kr exp       Xe exp
!
      eth = (/     429.0d0,     178.0d0,      129.0d0,      107.0d0,
     &               0.0d0,      27.6d0,       59.0d0,
     &              26.0d0,      36.5d0,       36.0d0,       42.8d0/)

      qy  = (/     0.007d0,    0.0179d0,     0.0654d0,      0.110d0,
     &               0.0d0,     0.780d0,       30.0d0,
     &              2.67d0,      7.60d0,       16.2d0,       22.1d0/)
!
!-----------------------------------------------------------------------
!:: target undefined
!-----------------------------------------------------------------------
      else
       write(*,*) "no target data"
       call wexit( 'sputN', 'no target data' )
       stop
      endif


      do i = 1, 11
        Z1 = Z(i)
        M1 = M(i)
        etf(i) = 30.74d0 * (M1+M2)/M2 * Z1*Z2
     &               * (Z1**(2.0d0/3.0d0)+Z2**(2.0d0/3.0d0))**0.5d0
        fy(i)  = sqrt(Es) * ( 0.94d0 - 0.00133*M2/M1 )
      enddo
!
!
      return
      end


!***********************************************************************
      integer function nsptin(catmz)
!***********************************************************************
      use clocal, only : csptin, msptin
      implicit none
!
      character(*), intent(in) :: catmz
!
      integer i

!      catmz = trim(catmz)
      do i = 1, msptin
        if(catmz.eq.csptin(i)) then
          nsptin=i
!          write(*,*) "nsptin: catmz = ",catmz," nsptin = ",nsptin
          exit
        endif
      enddo

      return
      end
!
!!***********************************************************************
      subroutine test_sputphy
!***********************************************************************
      use cunit, only : n6
      implicit none
!
!::local variables
      real*8  emin, emax, de, e0, y1, y2, y3, y4, y5
      real*8  alf
      integer i, n

! function
      real(8)    sputphy
!
      write(n6,'(2x,"=== test_sputphy ===")')
      write(n6,'(6x,"i",3x,"E0",10x,"Y_D",8x,"Y_He",9x,"Y_C",9x,"Y_W",
     >  8x,"Y_Ar")')
!
      call set_sputphys("C")
!
      alf = 0.0d0
!
      n = 41
      emin = 10.0d0
      emax = 100.0d3
      de   = dlog(emax/emin)/dfloat(n-1)
      do i = 1, n
      e0  = emin*dexp(de*dfloat(i-1))
      y1  = sputphy(e0,alf,3)
      y2  = sputphy(e0,alf,4)
      y3  = sputphy(e0,alf,6)
      y4  = sputphy(e0,alf,7)
      y5  = sputphy(e0,alf,9)
      write(n6,'(2x,i5,1p6e12.3)')
     >  i, e0, y1, y2, y3, y4, y5
      enddo
!
      return
      end
!
!***********************************************************************
      subroutine cal_etf
!***********************************************************************
      implicit none

      real*8  M(11), Es(3)
      integer Z(11)
      integer i,j,Z1,Z2
      real*8  etf, fy,M1, M2

!               1H,     2D,   3T,   4He,
!              5Be,     6C,   7W,
!              8Ne,    9Ar, 10Kr,   11Xe

      Z = (/     1,     1,     1,     2,
     &           4,     6,    74,
     &          10,    18,    36,     54 /)

      M = (/ 1.008d0, 2.014d0, 3.016d0, 4.003d0,
     &        9.012d0, 12.01d0, 183.8d0,
     &        20.18d0, 39.95d0, 83.80d0, 131.3d0 /)
!
!binding energy  Be     C     W
      Es= (/ 3.38,  7.42, 8.68 /)


      do j = 1, 3
        Z2 = Z(j+4)
        M2 = M(j+4)
        do i = 1, 11
          Z1 = Z(i)
          M1 = M(i)
          etf = 30.74d0 * (M1+M2)/M2 * Z1*Z2
     &         * (Z1**(2.0d0/3.0d0)+Z2**(2.0d0/3.0d0))**0.5d0
          fy = sqrt(Es(j)) * ( 0.94d0 - 0.00133*M2/M1 )

          write(*,'(2i4,2f10.3,f20.5,f20.5)') Z1, Z2, M1, M2,etf, fy
        enddo
      enddo

      return
      end
