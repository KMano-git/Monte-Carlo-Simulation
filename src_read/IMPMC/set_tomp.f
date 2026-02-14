!***********************************************************************
      subroutine set_tomp
!***********************************************************************
      use cimcns, only : ftom, ndtom, ntom, tomeng, tomran, tomvel
      use cimcom, only : amz
      use clocal, only : es
      use cphcns, only : cev
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  n, i, n1
      real*8   u0, emin, emax
      real*8   etbl(ndtom), ptbl(ndtom), rtbl(ndtom)
      real*8   sum, e, p, de, eave, dx, x
! finction
      real(8)    dintp
!
      write(n6,'(2x,"=== set_tomps ===")')
!
!::const (Thompson model)
!      u0   = 7.4d0
      u0   = Es
      n    = ndtom-1
      emin = 0.4d0
      emax = 100.0d0
!
!::probability
      de = (emax-emin)/dfloat(n)
      sum = 0.0d0
      do i = 1, n
      e = emin + de*(dfloat(i)-0.5d0)
!      p = (e/(e+u0/3.0d0)**4)**(2.0d0/3.0d0)
      p = (e/(e+u0)**(8.0d0/3.0d0))
      etbl(i) = e
      ptbl(i) = p
      rtbl(i) = sum
      sum = sum + p*de
      enddo
      n1 = n + 1
      i  = n1
      etbl(i) = emax
      ptbl(i) = ptbl(n)
      rtbl(i) = sum
!
!::normalization
      do i = 1, n1
      ptbl(i) = ptbl(i)/sum
      rtbl(i) = rtbl(i)/sum
      enddo
!
!::mean value
      eave = 0.0d0
      do i = 1, n1
      eave = eave + etbl(i)*ptbl(i)*de
      enddo
      write(n6,'(2x,"n =",i5,"  emin =",1pe10.3,"  emax =",1pe10.3,
     >   "  <E> =",1pe10.3)') n, emin, emax, eave
!
!::uniform random number
      ntom = n
      dx = 1.0d0/dfloat(ntom)
      do i = 1, ntom
      x = dx*(dfloat(i)-0.5d0)
      e = dintp(x, n1,rtbl,etbl)
      tomran(i) = x
      tomeng(i) = e
      tomvel(i) = sqrt(2.0d0*cev*e/amz)
      enddo
      i = ntom + 1
      tomran(i) = tomran(ntom)
      tomeng(i) = tomeng(ntom)
      tomvel(i) = tomvel(ntom)
      ftom = ntom
!
!::debug
      call test_set_tomp
!
      return
      end
!
!***********************************************************************
      subroutine test_set_tomp
!***********************************************************************
      use cimcns, only : ftom, ntom, tomeng
      use cunit,  only : n6
      implicit none
!
!::local variables
      real*8   de, emax, x(201), y(201), eave, xran, e, sum
      integer  nmax, nsmp, i, ix, ih
! function
      real(8)    random
!
      nmax = 100
      emax = tomeng(ntom)
      de   = emax/nmax
      do i = 1, nmax+1
      x(i) = 0.0d0
      y(i) = 0.0d0
      enddo
!
      nsmp = 10000
      eave = 0.0d0
      do i = 1, nsmp
      xran  = random(0)
      ix = int(ftom*xran + 1.0d0)
      e  = tomeng(ix)
      ih = int(e/de + 1.0d0)
      y(ih) = y(ih) + 1.0d0
      eave = eave + e
      enddo
      eave = eave/dfloat(nsmp)
      write(n6,'(2x,"<E> =",1pe10.3," for nsmp =",i8)') eave, nsmp
!
      write(n6,'(6x,"i",3x,"E(eV)",7x,"prb",9x,"Iprb")')
      sum = 0.0d0
      do i = 1, nmax+1
      x(i) = de*(dfloat(i)-0.5d0)
      y(i) = y(i)/dfloat(nsmp)/de
      sum  = sum + y(i)*de
      write(n6,'(2x,i5,1p3e12.3)')
     >  i, x(i), y(i), sum
      enddo
!
      return
      end
