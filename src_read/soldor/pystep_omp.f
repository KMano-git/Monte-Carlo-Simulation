!%%mkspt +
!**********************************************************************
      subroutine pystep
!**********************************************************************
!                                   -
!   dq/dt*vl + viscos + pressure = dq*vl
!
!----------------------------------------------------------------------
      use cplcom, only : c00, c11, cc, dd, dq1a, dq2a, dq3, dq4, dtrateq
     >    , ee, ff, lbcgd, nion, ss
      use cplmet, only : hvol, icmax, jcmax
      use csize,  only : ndeq, ndxy
      use csonic, only : dtim, lfopt
      implicit none
!
!::local variables
      integer  j, imax, imax1, nequ, n, m, i, ia
      integer  m1a, m2a, m3, m4
      integer  iwb
      real*8   dtvl, dtyy
      real*8 :: zero = 1.0e-50

!
!----------------------------------------------------------------------
!::plasma parameter at the main plasma
!----------------------------------------------------------------------
      dtyy = dtim
      if( lfopt(2).eq.0 ) dtyy = 0.0d0
!
      call plbman
!
!----------------------------------------------------------------------
!::colum : vertical direction
!----------------------------------------------------------------------
!==May/04    Exclude  do-variables jc,n,m,ia,iw   Error nexy,dtyy
!KH130604    Include n, m, ia, iw
!
!ik!$omp parallel private( nj, j, imax, imax1, nequ, i, iwb,
!ik!$omp& m1a, m2a, m3, m4, nd1, nd2, neq, nmx, dtvl,
!ik!$omp& n, m, iw, ia )
!$omp parallel private( imax, imax1, nequ, i, iwb,
!$omp& m1a, m2a, m3, m4, dtvl, n, m, ia )
!
!$omp do
      do 100 j = 2, jcmax-1
      imax  = icmax(j)
      imax1 = imax-1
!
!----------------------------------------------------------------------
!::yacobian matrix
!----------------------------------------------------------------------
!
!::clear
      call pycler(j)
!
!::viscosity term
      call pyvisc(j)
      call pyvpin(j)
!
!-----------------------------------------------------------------------
!::number of equations
!-----------------------------------------------------------------------
      nequ = 2*nion + 2
      m3   = 2*nion + 1
      m4   = 2*nion + 2
!
!-----------------------------------------------------------------------
!::factor dtim/hvol(j,i) & term I
!-----------------------------------------------------------------------
      do 110 n = 1, nequ
      do 115 m = 1, nequ
      do 120 i = 2, imax1
      dtvl = dtyy/hvol(j,i)*dtrateq(j,i,n)
      dd(n,m,i) = dd(n,m,i) + ss(n,m,i)
      cc(n,m,i) = cc(n,m,i)*dtvl
      dd(n,m,i) = dd(n,m,i)*dtvl
      ee(n,m,i) = ee(n,m,i)*dtvl
      if(abs(dd(n,m,i))<zero ) then
        dd(n,m,i) = 0.0d0
      endif
      if(abs(ee(n,m,i))<zero ) then
        ee(n,m,i) = 0.0d0
      endif
 120  continue
 115  continue
 110  continue
!
      do 130 n  = 1, nequ
      do 135 i = 2, imax1
      dd(n,n,i) = dd(n,n,i) + c11
 135  continue
 130  continue
!
!-----------------------------------------------------------------------
!::dq^bar in RHS
!-----------------------------------------------------------------------
      do 150 ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      do 160 i = 2, imax1
      ff(m1a,i) = dq1a(j,i,ia)
      ff(m2a,i) = dq2a(j,i,ia)
 160  continue
 150  continue
      do 170 i = 2, imax1
      ff(m3,i) = dq3(j,i)
      ff(m4,i) = dq4(j,i)
 170  continue
!
!-----------------------------------------------------------------------
!::boundary condition : iwb = 1
!-----------------------------------------------------------------------
      iwb = 1
      do 210 n = 1, nequ
      do 215 m = 1, nequ
      cc(n,m,iwb) = c00
      dd(n,m,iwb) = c00
      ee(n,m,iwb) = c00
  215 continue
      dd(n,n,iwb) = c11
      ff(n,iwb)   = c00
  210 continue
!
!-----------------------------------------------------------------------
!::boundary condition iwb = imax
!-----------------------------------------------------------------------
      iwb = imax
      do 230 n = 1, nequ
      do 235 m = 1, nequ
      cc(n,m,iwb) = c00
      dd(n,m,iwb) = c00
      ee(n,m,iwb) = c00
  235 continue
      dd(n,n,iwb) = c11
      ff(n,iwb)   = c00
  230 continue
!
!-----------------------------------------------------------------------
!::boundary condition iw = 1 & imax
!-----------------------------------------------------------------------
      if( lbcgd.eq.1 ) then
        call pybcon(j)
      else
        call pybcon2(j)
      endif
!
!-----------------------------------------------------------------------
!::solve
!-----------------------------------------------------------------------
      call tdmslv( ndeq, ndxy, nequ, imax, cc, dd, ee, ff )
!
!-----------------------------------------------------------------------
!::store the solution
!-----------------------------------------------------------------------
      do 320 i = 1, imax
      do 330 ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      dq1a(j,i,ia) = ff(m1a,i)
      dq2a(j,i,ia) = ff(m2a,i)
 330  continue
      dq3(j,i) = ff(m3,i)
      dq4(j,i) = ff(m4,i)
 320  continue
!
 100  continue
!$omp end parallel
!
      return
      end
