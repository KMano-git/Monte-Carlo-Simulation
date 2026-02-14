!**********************************************************************
      subroutine test_funclz
!**********************************************************************
      implicit none
!
!::local variables
      real*8  temin, temax, dlgte, te, funclz
      integer nmax, i
!
      temin = 0.1d0
      temax = 1.0d3
      nmax  = 41
      dlgte = dlog(temax/temin)/dfloat(nmax-1)
      do i = 1, nmax
      te = temin*exp(dlgte*dfloat(i-1))
      write(6,'(2x,i5,1p2e12.3)') i, te, funclz(te)
      enddo
!
      stop
      end
!
!**********************************************************************
      function funclz(ate)
!**********************************************************************
      implicit none
!
!::argument
      real*8  funclz, ate
!
!::local variables
      real*8  zx, zy, zlz
!
      zx = ate
      if( zx.gt.1.0d3 ) zx = 1.0d3
      zx = dlog10(zx)
      zy = -36.664d0 + 11.499d0*zx - 9.9541d0*zx**2
     >               + 7.2144d0*zx**3 - 4.2755d0*zx**4
     >               + 1.3616d0*zx**5 - 0.16327d0*zx**6
      zlz = 10.0d0**zy
!
!-----------------------------------------------------------------------
!::Lz(Te)
!-----------------------------------------------------------------------
      funclz = zlz
      return
      end
