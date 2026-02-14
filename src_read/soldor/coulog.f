!xx      call test_coulog
!xx      stop
!xx      end
!
!*********************************************************************
      subroutine test_coulog
!*********************************************************************
      implicit none
!
      real(8), dimension(6)  :: tbne
     >  = ( /1.0d19, 5.0d19, 1.0d20, 5.0d20, 1.0d21, 5.0d21/ )
      real(8), dimension(6)  :: tbclg
      real(8), dimension(50) :: tbte
      integer :: nmax, i, j
      real(8) :: temin, temax, dte
!
      real(8) :: Ai, Zi, zeff, ne, te, ti, coulog
!
      nmax = 31
      temin = 0.5d0
      temax = 10.0d3
      dte = dlog10(temax/temin)/dfloat(nmax-1)
      do i = 1, nmax
      tbte(i) = temin*10.0d0**(dte*dfloat(i-1))
      enddo
!
      write(6,'(2x,5x,"Te",9x,1p10e12.3)') (tbne(j),j=1,6)
!
      zeff = 1.5d0
      Ai = 2.0d0
      Zi = 1.0d0
!
      do i = 1, nmax
      te = tbte(i)
      ti = te
      do j = 1, 6
      ne = tbne(j)
      tbclg(j) = coulog( zeff, ne, te, ti, Ai, Zi )
      enddo
!
      write(6,'(2x,2x,1pe12.3,2x,1p10e12.3)')
     >   tbte(i), (tbclg(j),j=1,6)
      enddo
!
      return
      end
!
!*****************************************************************************
!   Coulomb logarithm Formulae given by Mitsuru HONDA (2011/06/13)
!     References: M. Honda, in preparation
!                 D.V. Sivukhin, Rev. Plasma Phys. Vol. 1 (1966)
!
! Note: Beam temperature defined by Tb = (2/3)Eb, where Eb: injection energy
!*****************************************************************************
!
!     Note: Coulog assumes that the ions share a common temperature, ti.
!           Set "2" to beam ions if "1" is the ion species.
!
!     electron - ion Coulomg logarithm
!
!*********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   function coulog( zeff, ne, te, ti, Ai, Zi )
      real(8) function coulog( zeff, ne, te, ti, Ai, Zi )
!*********************************************************************
      implicit none
!
!::arguments
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik    real(8) :: coulog
       real(8), intent(in) :: zeff, ne, te, ti, Ai, Zi
!
!      zeff   : effective charge
!      ne     : electron density in m^{-3}
!      te, ti : electron and ion temperatures in eV
!      A1, A2 : mass number
!      Z1, Z2 : ABSOLUTE charge number, e.g. Ze = 1
!
!
!::local variables
      real(8), parameter :: Ae = 5.446170219d-4 ! from CODATA 2010
      real(8) :: A1, Z1, A2, Z2
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real(8) :: lambda, CD, t1, t2
      real(8) :: lambda, CD
      real(8) :: coef, coef1, coef2, thres_l, thres_r

      CD = ne * ( 1.d0 / te + zeff / ti )
!
!::electron-ion Coulomb logarithm
      A1 = Ae
      Z1 = 1.0d0
      A2 = Ai
      Z2 = Zi
!
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   t1 = te
!ik   t2 = ti
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   coef1 = ( A1 + A2 ) / ( A2 * t1 + A1 * t2 )
      coef1 = ( A1 + A2 ) / ( A2 * te + A1 * ti )
      coef2 = A1 * A2 / ( A1 + A2 )
      coef  = log( coef1 * sqrt(CD) )

      thres_l = 1.d0 / coef1
      thres_r = 2.45d4 * ( Z1 * Z2 )**2 * coef2
!
      if ( thres_l <= thres_r ) then
!   Classical formula
        lambda = 30.4d0 - log( Z1 * Z2 ) - coef
      else
!   Quantum-mechanical formula
        lambda = 35.4d0 - coef + 0.5d0 * log( coef1 * coef2 )
      end if

      coulog = lambda
      return
      end
