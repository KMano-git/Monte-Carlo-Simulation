!***********************************************************************
      subroutine pxbcon(nt)
!***********************************************************************
!
!   Yacobian at boundary
!
!   boundary condition at divertor plate
!      j=1  (outer plate) / j=jmax (inner plate)
!
!  (1) na   Vj*dQ1a + S^sita*b^sita*Q2a = S^sita*b^sita*Q2a @BC
!  (2) va   va = Csd   sound velocity
!  (3) Ti   Vj*dQ3/dt + qfx_cv + qfx_cd + qfx_vh =
!             S^sita*b^sita*{
!              sum_a{qfx_cva/Ma*1/2*Ma*Va**2} + cdli*flb1ni*Ti*eV}
!  (4) Te   Vj*dQ4/dt + qfx_cv + qfx_cd =
!             S^sita*b^sita*{cdle*flb1ne*Te*eV}
!
!   flux at dplate
!
!   flb1a = S^sita*b^sita * ma*na*va
!   flb2a = S*b * (ma*na*va^2+Pa + ]~[a )
!   flb3  = S*b *{ Sum_a[(1/2*ma*na*va^2+3/2*na*Ti+Pa)*va] + q//i
!                               + Sum_a[(va*]~[)]a }
!   flb4  = S*b *{ (3/2*ne*Te+Pe)*va + q//e }
!
!   flb1ro = sum_a[flb1a]
!   flb1ni = sum_a[1/ma*flb1a]
!   flb1ne = sum_a[za/ma*flb1a]
!
!     Vdq/dt + sigv*Hbc + ...
!     flb1a = sigv*qfx_cv
!
!                                           ! 2002/11/11
!::error trap   Na < 0 due to flb1a         ! 2009/07/16
!::set Na at dplate  when Na < bcna1        ! 2011/05/19
!::Na  include V*dQ1a/dt                    ! 2011/06/16
!::Ti  1/2*Ma*Cs^2 (error) ==> 1/2*Ma*Va^2  ! 2011/06/19
!
!-----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : ama, aza, c12, c13, c23, cc, cdle, cdli, dd, ee
     >    , fcbcvl, ff, nion, nlp, q1a, q2a, ss, vcs, vna, vne, vni, vte
     >    , vti, vva
      use cplmet, only : hvol, hvsb, icel, jcel, jtmax
      use cplqcn, only : qfx_cv
      use csize,  only : ndsp
      use csonic, only : itim
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nt
      integer, intent(in) :: nt
!
!::local variables
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, i, jmax, jwb1, jwb2, jm, jw, j, jp, jsf
!ik   integer  m, n, mlx1, mlx2, kbc
      integer  i, jmax, jwb1, jwb2, jm, jw, j, jp, jsf
      integer  m, n, mlx1, mlx2
      integer  ia, m1a, m2a, m3, m4
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8   sigv, vbc, zsmro, zsmni, zsmne, zroi
      real*8   sigv, vbc, zsmro, zsmni, zsmne
      real*8   d3va, d3ti, d4te
! modified 5/4 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8   flb1a(ndsp), flb2a(ndsp), glb1a(ndsp)
!ik   real*8   flb1ro, flb1ni, flb1ne, flb3, flb4
!ik   real*8   tibc, tebc, fac, fcvol
!ik   integer  i6, lout, l, ierr
!ik   character clin*150
      real*8   flb1a(ndsp), glb1a(ndsp)
      real*8   flb1ni, flb1ne
      real*8   fac, fcvol
      integer  ierr
!
!::dq/dt*Vol  Vol(1) = fcvol*Vol(2)
      fcvol = fcbcvl
      if( itim.le.50 ) fcvol = 0.02d0
!
! modified 3/2 lines organize local variables and include files by kamata 2021/05/31
!ik   it   = nt
!ik   i    = icel(1,it)
!ik   jmax = jtmax(it)
      i    = icel(1,nt)
      jmax = jtmax(nt)
!
      m3   = 2*nion + 1
      m4   = 2*nion + 2
!
!::debug write
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   i6 = 0   !  231
!
!-----------------------------------------------------------------------
!::clear  ss,cc,dd,ee,ff for q1a and q2a
!-----------------------------------------------------------------------
!::Not clear  q1a
      mlx1  = nion + 1
      mlx2  = 2*nion
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
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   j    = jcel(jw,it)
!ik   jp   = jcel(jw+1,it)
      j    = jcel(jw,nt)
      jp   = jcel(jw+1,nt)
      jsf  = j
      sigv = -1.0d0
      vbc  = sigv*vcs(j,i)   !  for all species including sign
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   kbc  = 1
!
!::hvol of dummy cell
      hvol(j,i) = hvol(jp,i)*fcvol
!
!::flux and variables at boundary
      zsmro = 0.0d0
      zsmni = 0.0d0
      zsmne = 0.0d0
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   zroi  = 0.0d0
      ierr = 0
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      flb1a(ia) = qfx_cv(jsf,i,m1a)*sigv
      glb1a(ia) = hvsb(jsf,i)*ama(ia)*vna(j,i,ia)*vbc*sigv
!-----
      fac = flb1a(ia)/glb1a(ia)
      if( fac.lt.0.0 .or. dabs(fac).lt.0.3d0 ) then
        flb1a(ia) = glb1a(ia)
        ierr = ierr + 1
        write(n6,'(2x,"pxbcon error qfx_cv(m1a)",i7,i3,2i5,1p2e11.3)')
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik  >  itim,nlp,it,jw, qfx_cv(jsf,i,m1a)/ama(ia), glb1a(ia)/ama(ia)
     >  itim,nlp,nt,jw, qfx_cv(jsf,i,m1a)/ama(ia), glb1a(ia)/ama(ia)
      endif
!-----
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   flb2a(ia) = qfx_cv(jsf,i,m2a) + qfx_cd(jsf,i,m2a)
      zsmro = zsmro + flb1a(ia)
      zsmni = zsmni + flb1a(ia)/ama(ia)
      zsmne = zsmne + flb1a(ia)/ama(ia)*aza(ia)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   zroi  = zroi + q1a(j,i,ia)
      enddo
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   flb1ro = zsmro
      flb1ni = zsmni
      flb1ne = zsmne
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   flb3 = qfx_cv(jsf,i,m3)+qfx_cd(jsf,i,m3)+qfx_vh(jsf,i)
!ik   flb4 = qfx_cv(jsf,i,m4)+qfx_cd(jsf,i,m4)
!
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   tibc = (flb3-c12*vbc**2*flb1ro)/(cdli*flb1ni*cev)
!ik   tebc = flb4/(cdle*flb1ne*cev)
!
!::Not use Cs
      d3ti = cdli*flb1ni
      d4te = cdle*flb1ne
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::Na
!xx      dd(m1a,m1a,jw) = 1.0d0
!xxxx ff(m1a,jw)     = -q1a(j,i,ia)+ama(ia)*bna(ia)
!xxxx ff(m1a,jw)     = -q1a(j,i,ia)+qfx_cv(jsf,i,m1a)/(hvsb(jsf,i)*vbc)
!xx   ff(m1a,jw)     = -q1a(j,i,ia)+flb1a(ia)/(hvsb(jsf,i)*vbc)
!
!::q1a
      dd(m1a,m2a,jw) = dd(m1a,m2a,jw) + sigv*hvsb(jsf,i)
      ff(m1a,jw) = ff(m1a,jw) - sigv*hvsb(jsf,i)*q2a(j,i,ia)
!
!::Va  see pxbcnv
      dd(m2a,m1a,jw) = -q2a(j,i,ia)/q1a(j,i,ia)**2
      dd(m2a,m2a,jw) = 1.0d0/q1a(j,i,ia)
      ff(m2a,jw)     = -vva(j,i,ia)+vbc
!
!::Hbc_I = sig(va)*(sum_a{flb1a/Ma*1/2*Ma*Va**2}
!                     + cdli*flb1ni*Ti*eV)
      d3va = flb1a(ia)*vva(j,i,ia)
      dd(m3,m1a,jw) = dd(m3,m1a,jw)
     >       - d3va*q2a(j,i,ia)/q1a(j,i,ia)**2
     >       + d3ti*(c13*vva(j,i,ia)**2-vti(j,i)*cev/ama(ia))/vni(j,i)
      dd(m3,m2a,jw) = dd(m3,m2a,jw)
     >       + d3va/q1a(j,i,ia)
     >       - d3ti*c23*vva(j,i,ia)/vni(j,i)
      ff(m3,jw) = ff(m3,jw) - flb1a(ia)*c12*vva(j,i,ia)**2
!
!::Hbc_e = sig(va)*cdle*flb1ne*Te*eV
      dd(m4,m1a,jw) = dd(m4,m1a,jw)
     >      - d4te*aza(ia)/ama(ia)*vte(j,i)*cev/vne(j,i)
      enddo  ! loop(ia)
!
      dd(m3,m3, jw) = dd(m3,m3, jw) + d3ti*c23/vni(j,i)
      ff(m3,jw) =ff(m3,jw) - d3ti*vti(j,i)*cev
!
      dd(m4,m4, jw) = dd(m4,m4, jw) + d4te*c23/vne(j,i)
      ff(m4,jw)     = ff(m4,jw) - d4te*vte(j,i)*cev
!
!-----------------------------------------------------------------------
!::boundary condition j = jmax
!-----------------------------------------------------------------------
      jw   = jmax
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   j    = jcel(jw,it)
!ik   jm   = jcel(jw-1,it)
!ik   jsf  = jcel(jw-1,it)
      j    = jcel(jw,nt)
      jm   = jcel(jw-1,nt)
      jsf  = jcel(jw-1,nt)
      sigv = +1.0d0
      vbc  = sigv*vcs(j,i)   !  for all species
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   kbc  = 2
!
!::hvol of dummy cell
      hvol(j,i) = hvol(jm,i)*fcvol
!
!::flux and variables at boundary
      zsmro = 0.0d0
      zsmni = 0.0d0
      zsmne = 0.0d0
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   zroi  = 0.0d0
      ierr  = 0
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      flb1a(ia) = qfx_cv(jsf,i,m1a)*sigv
      glb1a(ia) = hvsb(jsf,i)*ama(ia)*vna(j,i,ia)*vbc*sigv
!-----
      fac = flb1a(ia)/glb1a(ia)
      if( fac.lt.0.0 .or. dabs(fac).lt.0.3d0 ) then
        flb1a(ia) = glb1a(ia)
        ierr = ierr + 1
        write(n6,'(2x,"pxbcon error qfx_cv(m1a)",i7,i3,2i5,1p2e11.3)')
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik  >  itim,nlp,it,jw, qfx_cv(jsf,i,m1a)/ama(ia), glb1a(ia)/ama(ia)
     >  itim,nlp,nt,jw, qfx_cv(jsf,i,m1a)/ama(ia), glb1a(ia)/ama(ia)
      endif
!-----
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   flb2a(ia) = qfx_cv(jsf,i,m2a) + qfx_cd(jsf,i,m2a)
      zsmro = zsmro + flb1a(ia)
      zsmni = zsmni + flb1a(ia)/ama(ia)
      zsmne = zsmne + flb1a(ia)/ama(ia)*aza(ia)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   zroi  = zroi + q1a(j,i,ia)
      enddo
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   flb1ro = zsmro
      flb1ni = zsmni
      flb1ne = zsmne
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   flb3 = qfx_cv(jsf,i,m3)+qfx_cd(jsf,i,m3)+qfx_vh(jsf,i)
!ik   flb4 = qfx_cv(jsf,i,m4)+qfx_cd(jsf,i,m4)
!
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   tibc = (flb3-c12*vbc**2*flb1ro)/(cdli*flb1ni*cev)
!ik   tebc =  flb4/(cdle*flb1ne*cev)
!
!::Not use Cs
      d3ti = cdli*flb1ni
      d4te = cdle*flb1ne
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::Na
!xx   dd(m1a,m1a,jw) = 1.0d0
!xxxx ff(m1a,jw)     = -q1a(j,i,ia)+ama(ia)*bna(ia)
!xxxx ff(m1a,jw)     = -q1a(j,i,ia)+qfx_cv(jsf,i,m1a)/(hvsb(jsf,i)*vbc)
!xx   ff(m1a,jw)     = -q1a(j,i,ia)+flb1a(ia)/(hvsb(jsf,i)*vbc)
!
!::q1a
      dd(m1a,m2a,jw) = dd(m1a,m2a,jw) + sigv*hvsb(jsf,i)
      ff(m1a,jw) = ff(m1a,jw) - sigv*hvsb(jsf,i)*q2a(j,i,ia)
!
!::Va  see pxbcnv
      dd(m2a,m1a,jw) = -q2a(j,i,ia)/q1a(j,i,ia)**2
      dd(m2a,m2a,jw) = 1.0d0/q1a(j,i,ia)
      ff(m2a,jw)     = -vva(j,i,ia)+vbc
!
!::Hbc_I = sig(va)*(sum_a{flb1a/Ma*1/2*Ma*Va**2}
!                     + cdli*flb1ni*Ti*eV)
      d3va = flb1a(ia)*vva(j,i,ia)
      dd(m3,m1a,jw) = dd(m3,m1a,jw)
     >       - d3va*q2a(j,i,ia)/q1a(j,i,ia)**2
     >       + d3ti*(c13*vva(j,i,ia)**2-vti(j,i)*cev/ama(ia))/vni(j,i)
      dd(m3,m2a,jw) = dd(m3,m2a,jw)
     >       + d3va/q1a(j,i,ia)
     >       - d3ti*c23*vva(j,i,ia)/vni(j,i)
      ff(m3,jw) = ff(m3,jw) - flb1a(ia)*c12*vva(j,i,ia)**2
!
!::Hbc_e = sig(va)*cdle*flb1ne*Te*eV
      dd(m4,m1a,jw) = dd(m4,m1a,jw)
     >      - d4te*aza(ia)/ama(ia)*vte(j,i)*cev/vne(j,i)
      enddo  ! loop(ia)
!
      dd(m3,m3, jw) = dd(m3,m3, jw) + d3ti*c23/vni(j,i)
      ff(m3,jw) =ff(m3,jw) - d3ti*vti(j,i)*cev
!
      dd(m4,m4, jw) = dd(m4,m4, jw) + d4te*c23/vne(j,i)
      ff(m4,jw)     = ff(m4,jw) - d4te*vte(j,i)*cev
!
      return
      end
