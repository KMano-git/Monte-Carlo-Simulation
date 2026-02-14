!***********************************************************************
      subroutine pxbcnv(nt)
!***********************************************************************
!
!     pxbcnv   n : density  v : velocity at d-plate
!
!   boundary condition at divertor plate
!      j=1  (outer plate) / j=jmax (inner plate)
!
!  (1) na   Vj*dQ1a + S^sita*b^sita*Q2a = S^sita*b^sita*Q2a @BC
!  (2) va   va = Csd   sound velocity
!
!  Note : Not include Vj*dQ/dt term for BC condition of va
!         see sub. pxstep_omp.f
!
!-----------------------------------------------------------------------
      use cplcom, only : c12, c13, c32, cc, dd, ee, ff, nion, q1a, q2a
     >    , ss, vcs, vva
      use cplmet, only : icel, jcel, jtmax
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nt
      integer, intent(in) :: nt
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, i, jmax, jwb1, jwb2, mlx, m, n, jw, j, jp, jm, jsf
      integer  i, jmax, jwb1, jwb2, m, n, jw, j
      integer  mlx1, mlx2
      real*8   sigv, zcs, zroi, dcs
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   real*8   flb1a(ndsp)
      integer  ia, m1a, m2a, m3, m4, ib, m1b, m2b
!
! modified 3/2 lines organize local variables and include files by kamata 2021/05/31
!ik   it   = nt
!ik   i    = icel(1,it)
!ik   jmax = jtmax(it)
      i    = icel(1,nt)
      jmax = jtmax(nt)
      m3 = 2*nion + 1
      m4 = 2*nion + 2
!
!-----------------------------------------------------------------------
!::clear  ss,cc,dd,ee,ff for q1a and q2a
!-----------------------------------------------------------------------
!   Do not change for q1a, q3 and q4
      mlx1 = 2*nion + 1
      mlx2 = 2*nion
!
      jwb1 = 1
      jwb2 = jmax
      do m = mlx1, mlx2
      ff(m,jwb1) = 0.0d0
      ff(m,jwb2) = 0.0d0
      do n = 1, m4
      ss(m,n,jwb1) = 0.0d0
      cc(m,n,jwb1) = 0.0d0
      dd(m,n,jwb1) = 0.0d0
      ee(m,n,jwb1) = 0.0d0
      ss(m,n,jwb2) = 0.0d0
      cc(m,n,jwb2) = 0.0d0
      dd(m,n,jwb2) = 0.0d0
      ee(m,n,jwb2) = 0.0d0
      enddo
      enddo
!
!-----------------------------------------------------------------------
!::boundary condition j = 1
!-----------------------------------------------------------------------
      jw   = 1
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j    = jcel(jw,it)
!ik   jp   = jcel(jw+1,it)
!ik   jsf  = j
      j    = jcel(jw,nt)
      sigv = -1.0d0
      zcs  = vcs(j,i)   !  for all species
!
      zroi = 0.0d0
      do ia = 1, nion
      zroi = zroi + q1a(j,i,ia)
      enddo
!
      dcs = sigv*c13/(zcs*zroi)
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::Va    va-sigv*Cs = 0
!
      do ib = 1, nion
      m1b = 2*ib - 1
      m2b = 2*ib
      dd(m2a,m1b,jw) = -dcs*(c12*vva(j,i,ib)**2-c32*zcs**2)
      dd(m2a,m2b,jw) = +dcs*vva(j,i,ib)
      enddo  ! loop(ib)
!
      dd(m2a,m1a,jw) = dd(m2a,m1a,jw) -q2a(j,i,ia)/q1a(j,i,ia)**2
      dd(m2a,m2a,jw) = dd(m2a,m2a,jw) +1.0d0/q1a(j,i,ia)
      dd(m2a,m3,jw)  = -dcs
      dd(m2a,m4,jw)  = -dcs
      ff(m2a,jw)     = -vva(j,i,ia)+sigv*zcs
      enddo  ! loop(ia)
!
!-----------------------------------------------------------------------
!::boundary condition j = jmax
!-----------------------------------------------------------------------
      jw   = jmax
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j    = jcel(jw,it)
!ik   jm   = jcel(jw-1,it)
!ik   jsf  = jcel(jw-1,it)
      j    = jcel(jw,nt)
      sigv = +1.0d0
      zcs  = vcs(j,i)   !  for all species
!
      zroi = 0.0d0
      do ia = 1, nion
      zroi = zroi + q1a(j,i,ia)
      enddo
!
      dcs = sigv*c13/(zcs*zroi)
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::Va    va-sigv*Cs = 0
!
      do ib = 1, nion
      m1b = 2*ib - 1
      m2b = 2*ib
      dd(m2a,m1b,jw) = -dcs*(c12*vva(j,i,ib)**2-c32*zcs**2)
      dd(m2a,m2b,jw) = +dcs*vva(j,i,ib)
      enddo  ! loop(ib)
!
      dd(m2a,m1a,jw) = dd(m2a,m1a,jw) -q2a(j,i,ia)/q1a(j,i,ia)**2
      dd(m2a,m2a,jw) = dd(m2a,m2a,jw) +1.0d0/q1a(j,i,ia)
      dd(m2a,m3,jw)  = -dcs
      dd(m2a,m4,jw)  = -dcs
      ff(m2a,jw)     = -vva(j,i,ia)+sigv*zcs
      enddo  ! loop(ia)
!
      return
      end
!
!***********************************************************************
      subroutine XX_pxbcnv(nt)
!***********************************************************************
!
!     pxbcnv   n : density  v : velocity at d-plate
!
!   boundary condition at divertor plate
!      j=1  (outer plate) / j=jmax (inner plate)
!
!  (1) na   Vj*dQ1a + S^sita*b^sita*Q2a = S^sita*b^sita*Q2a @BC
!  (2) va   va = Csd   sound velocity
!
!  Note : Not include Vj*dQ/dt term for BC condition of va
!         see sub. pxstep_omp.f
!
!-----------------------------------------------------------------------
      use cplcom, only : cc, dd, ee, ff, nion, q1a, q2a, ss, vcs, vva
      use cplmet, only : icel, jcel, jtmax
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nt
      integer, intent(in) :: nt
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, i, jmax, jwb1, jwb2, mlx, m, n, jw, j, jp, jm, jsf
      integer  i, jmax, jwb1, jwb2, m, n, jw, j
      integer  mlx1, mlx2
      real*8   sigv, zcs
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   real*8   flb1a(ndsp)
      integer  ia, m1a, m2a, m4
!
! modified 3/2 lines organize local variables and include files by kamata 2021/05/31
!ik   it   = nt
!ik   i    = icel(1,it)
!ik   jmax = jtmax(it)
      i    = icel(1,nt)
      jmax = jtmax(nt)
      m4 = 2*nion + 2
!
!-----------------------------------------------------------------------
!::clear  ss,cc,dd,ee,ff for q1a and q2a
!-----------------------------------------------------------------------
!   Do not change for q1a, q3 and q4
      mlx1 = 2*nion + 1
      mlx2 = 2*nion
!
      jwb1 = 1
      jwb2 = jmax
      do m = mlx1, mlx2
      ff(m,jwb1) = 0.0d0
      ff(m,jwb2) = 0.0d0
      do n = 1, m4
      ss(m,n,jwb1) = 0.0d0
      cc(m,n,jwb1) = 0.0d0
      dd(m,n,jwb1) = 0.0d0
      ee(m,n,jwb1) = 0.0d0
      ss(m,n,jwb2) = 0.0d0
      cc(m,n,jwb2) = 0.0d0
      dd(m,n,jwb2) = 0.0d0
      ee(m,n,jwb2) = 0.0d0
      enddo
      enddo
!
!-----------------------------------------------------------------------
!::boundary condition j = 1
!-----------------------------------------------------------------------
      jw   = 1
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j    = jcel(jw,it)
!ik   jp   = jcel(jw+1,it)
!ik   jsf  = j
      j    = jcel(jw,nt)
      sigv = -1.0d0
      zcs  = sigv*vcs(j,i)   !  for all species
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::Va
      dd(m2a,m1a,jw) = -q2a(j,i,ia)/q1a(j,i,ia)**2
      dd(m2a,m2a,jw) = 1.0d0/q1a(j,i,ia)
      ff(m2a,jw)     = -vva(j,i,ia)+zcs
      enddo  ! loop(ia)
!
!-----------------------------------------------------------------------
!::boundary condition j = jmax
!-----------------------------------------------------------------------
      jw   = jmax
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j    = jcel(jw,it)
!ik   jm   = jcel(jw-1,it)
!ik   jsf  = jcel(jw-1,it)
      j    = jcel(jw,nt)
      sigv = +1.0d0
      zcs  = sigv*vcs(j,i)   !  for all species
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::Va
      dd(m2a,m1a,jw) = -q2a(j,i,ia)/q1a(j,i,ia)**2
      dd(m2a,m2a,jw) = 1.0d0/q1a(j,i,ia)
      ff(m2a,jw)     = -vva(j,i,ia)+zcs
      enddo  ! loop(ia)
!
      return
      end
