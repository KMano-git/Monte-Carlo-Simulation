!%%mkspt +
!***********************************************************************
      subroutine pxstep
!***********************************************************************
!
!   Vij/dt*(q(N+1)-q(N)) + F^sita(q(N+1)) + F^psi(q(N+1)) = 0
!
!   ( I + dt/V*L^sita ) * ( I + dt/V*L^psi ) * dq = R(q(N+1,l))
!
!  # first-step : pxstep
!                            -
!     ( I + dt/V*L^sita ) * dq = dt/V*R
!        L^sita * dq = ss*dq(j) + cc*dq(j-1) + dd*dq(j) + ee*dq(j+1)
!
!         R(q(N+1,l)) = -(q(N+1,l)-q(N)) + dt/V*f^sita
!                                        + dt/V*f^psi
!
!  # second-step : pystep
!                                -
!     ( I + dt/V*L^psi ) * dq = dq
!
!      dq = q(N+1,l+1) - q(N+1,l)    l : iterlation loop
!      diffusion coef. Xicl = X(q(N+1,l))
!
!       Note.   f^psi should be included in R.
!       Note.   include dq/dt term at boundary cell  2011/06/08
!
!-----------------------------------------------------------------------
      use cplcom, only : ibcyl, c11, cc, dd, dq1a, dq2a, dq3, dq4
     >    , dtrateq, ee, ff, gg, nion, nlp, nlpmn, q1a, q2a, q3, q4, ss
     >    , wq1a, wq2a, wq3, wq4
      use cplmet, only : hvol, icel, itmax, itpve, jcel, jtmax
      use cplqcn, only : itqcn, qvl_al, qvl_dt
      use csize,  only : ndeq, ndxy
      use csonic, only : dtim, kpcn, lfopt, lstop
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer mlpmx, nequ, m1a, m2a, m3, m4, lp, nd1, nd2, neq, nmx
      integer nequ, m1a, m2a, m3, m4, nmx
      integer it, nt, i, jmax, jmax1, ipbc, n, m, jw, j, ia, jwst, jwen
      integer nexy
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer kdmy, jw1, i6
      integer kdmy, jw1
      real*8, dimension(ndeq,ndxy) :: xx
!
      real*8  dtvl, dtxx, dtvl1, dtvl2, dtvl3, dtvl4
!
      real*8 cc0(ndeq,ndeq,ndxy), dd0(ndeq,ndeq,ndxy)
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik  >      ,ee0(ndeq,ndeq,ndxy), ff0(ndeq,ndxy), qb0(ndeq,2)
     >      ,ee0(ndeq,ndeq,ndxy), ff0(ndeq,ndxy)
!xx   real*8  dtvbc(ndeq,ndxy)   ! KS  error 2011/09/22
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   real*8  dtvbc(ndxy,ndeq)   ! KS  No use in this sub.
!
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   i6 = 80000 + 30
!
!xx      write(n6,'(/2x,"time =",1pe12.3)') time
!
!-----------------------------------------------------------------------
!::number of equations
!-----------------------------------------------------------------------
!
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   mlpmx = 5
      nequ = 2*nion + 2
      m3   = 2*nion + 1
      m4   = 2*nion + 2
!
      nexy = ndeq*ndxy
      dtxx = dtim
      if( lfopt(1).eq.0 ) dtxx = 0.0d0
!
!-----------------------------------------------------------------------
!::flux tube
!-----------------------------------------------------------------------
!==May/04    Exclude do-variables  it,n,m,jw,ia,lp
!KH130604    Include n,m,jw,ia,lp,kdmy,jw1,xx,dtvbc
!KH130604    Exclude mlpmx, nequ, m3, m4, dtxx, nexy
!
! modified 5/5 lines organize local variables and include files by kamata 2021/05/31
!ik!$omp parallel private( m1a, m2a, nd1, nd2,
!ik!$omp& neq, nmx,  nt, i, jmax, jmax1, ipbc, j, jwst,
!ik!$omp& jwen, dtvl, dtvl1, dtvl2, dtvl3, dtvl4,
!ik!$omp& cc0, dd0, ee0, ff0, qb0,
!ik!$omp& n, m, jw, ia, lp, kdmy, jw1, xx, dtvbc )
!$omp parallel private( m1a, m2a,
!$omp& nmx, i, jmax, jmax1, ipbc, j, jwst,
!$omp& jwen, dtvl, dtvl1, dtvl2, dtvl3, dtvl4,
!$omp& cc0, dd0, ee0, ff0,
!$omp& n, m, jw, ia, kdmy, jw1, xx )
!
!$omp do
! modified 5/4 lines organize local variables and include files by kamata 2021/05/31
!ik   do 100 it = 2, itmax-1
!ik   nt = it
!ik   if( it.eq.itpve ) goto 100
!ik   i     = icel(1,it)
!ik   jmax  = jtmax(it)
      do 100 nt = 2, itmax-1
      if( nt.eq.itpve ) goto 100
      i     = icel(1,nt)
      jmax  = jtmax(nt)
      jmax1 = jmax-1
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   ipbc  = ibcyl(it)
      ipbc  = ibcyl(nt)
!
      jwst = 1
      jwen = jmax
      if( ipbc.eq.1 ) then
      jwst = 2
      jwen = jmax1
      endif
!
!-----------------------------------------------------------------------
!::yacobian matrix
!-----------------------------------------------------------------------
!
!::clear
      call pxcler(nt)
!
!::convection term
      call pxconv(nt)
      if( lstop .eq. 1 ) goto 100   ! Open_Mp  May/4
!
!::viscosity term
      call pxdiff(nt)
      call pxvisc(nt)
!
!::pressuer term
      call pxpres(nt)
!xx      if( nt.le.itpve ) call pxprdf(nt)
!
!::collision term
      call pxcoll(nt)
!
!::thermal force term
!xx   call pxthrm(nt)
!
!::source term
      call pxsorc(nt)
!
!::residual of F^psi -- gg
      call pxrdyv(nt)
      call pxrdvp(nt)
!xx   call pxrdyp(nt)
!
!::bc-conditon
      if( ipbc.eq.0 ) then
      call pxbcon(nt)
      endif
!
!::set dtvbc
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik   do jw = jwst, jwen
!ik   j = jcel(jw,it)
!ik   do ia = 1, nion
!ik   m1a = 2*ia - 1
!ik   m2a = 2*ia
!!      dtvbc(jw,m1a) = dtxx/hvol(j,i)*dtrat(j,i)*dteq(m1a)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   dtvbc(jw,m1a) = dtxx/hvol(j,i)*dtrateq(j,i,m1a)
!!      dtvbc(jw,m2a) = dtxx/hvol(j,i)*dtrat(j,i)*dteq(m2a)
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   dtvbc(jw,m2a) = dtxx/hvol(j,i)*dtrateq(j,i,m2a)
!ik   enddo
!!      dtvbc(jw,m3) = dtxx/hvol(j,i)*dtrat(j,i)*dteq(m3)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   dtvbc(jw,m3) = dtxx/hvol(j,i)*dtrateq(j,i,m3)
!!      dtvbc(jw,m4) = dtxx/hvol(j,i)*dtrat(j,i)*dteq(m4)
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   dtvbc(jw,m4) = dtxx/hvol(j,i)*dtrateq(j,i,m4)
!ik   enddo
!
! deleted 11 lines organize local variables and include files by kamata 2021/05/31
!ik   if( ipbc.eq.0 ) then
!ik   do jw = jwst, jwen
!ik   j = jcel(jw,it)
!ik   do ia = 1, nion
!ik   m1a = 2*ia - 1
!ik   m2a = 2*ia
!ik   dtvbc(jw,m1a) = 0.0d0
!ik   dtvbc(jw,m2a) = 0.0d0
!ik   enddo
!ik   enddo
!ik   endif
!
!::qcon   ! %%% 2002/11/02
!::V/dt*dltQ term
      do jw = jwst, jwen
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j = jcel(jw,it)
      j = jcel(jw,nt)
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!!      dtvl1 = dtxx/hvol(j,i)*dtrat(j,i)*dteq(m1a)
      dtvl1 = dtxx/hvol(j,i)*dtrateq(j,i,m1a)
!!      dtvl2 = dtxx/hvol(j,i)*dtrat(j,i)*dteq(m2a)
      dtvl2 = dtxx/hvol(j,i)*dtrateq(j,i,m2a)
      qvl_al(j,i,m1a) = ff(m1a,jw) + gg(m1a,jw)
      qvl_al(j,i,m2a) = ff(m2a,jw) + gg(m2a,jw)
      qvl_dt(j,i,m1a) = (q1a(j,i,ia)-wq1a(j,i,ia))/dtvl1
      qvl_dt(j,i,m2a) = (q2a(j,i,ia)-wq2a(j,i,ia))/dtvl2
      enddo
!!      dtvl3 = dtxx/hvol(j,i)*dtrat(j,i)*dteq(m3)
      dtvl3 = dtxx/hvol(j,i)*dtrateq(j,i,m3)
!!      dtvl4 = dtxx/hvol(j,i)*dtrat(j,i)*dteq(m4)
      dtvl4 = dtxx/hvol(j,i)*dtrateq(j,i,m4)
      qvl_al(j,i,m3)  = ff(m3,jw) + gg(m3,jw)
      qvl_al(j,i,m4)  = ff(m4,jw) + gg(m4,jw)
      qvl_dt(j,i,m3)  = (q3(j,i)-wq3(j,i))/dtvl3
      qvl_dt(j,i,m4)  = (q4(j,i)-wq4(j,i))/dtvl4
      enddo
!
!-----------------------------------------------------------------------
!::factor dtim/hvol(j,i) & term I
!-----------------------------------------------------------------------
      do 110 n = 1, nequ
      do 115 m = 1, nequ
      do 120 jw = jwst, jwen
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j    = jcel(jw,it)
      j    = jcel(jw,nt)
!!      dtvl = dtxx/hvol(j,i)*dtrat(j,i)*dteq(n)
      dtvl = dtxx/hvol(j,i)*dtrateq(j,i,n)
      dd(n,m,jw) = dd(n,m,jw) + ss(n,m,jw)
      cc(n,m,jw) = cc(n,m,jw)*dtvl
      dd(n,m,jw) = dd(n,m,jw)*dtvl
      ee(n,m,jw) = ee(n,m,jw)*dtvl
 120  continue
 115  continue
 110  continue
!
!::solv-xy
      if( lfopt(1).eq.0 ) then
      call setd( ff, nexy, 0.0d0 )
      endif
      if( lfopt(2).eq.0 ) then
      call setd( gg, nexy, 0.0d0 )
      endif
!
      do 130 n = 1, nequ
      do 135 jw = jwst, jwen
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j = jcel(jw,it)
      j = jcel(jw,nt)
!!      dtvl = dtim/hvol(j,i)*dtrat(j,i)*dteq(n)
      dtvl = dtim/hvol(j,i)*dtrateq(j,i,n)
      dd(n,n,jw) = dd(n,n,jw) + c11
      ff(n,jw)   = (ff(n,jw)+gg(n,jw))*dtvl
 135  continue
 130  continue
!
!-----------------------------------------------------------------------
!::dt/V*R = -(q(N+1,l)-q(N)) + ...
!-----------------------------------------------------------------------
      do 150 ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      do 155 jw = jwst, jwen
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j = jcel(jw,it)
      j = jcel(jw,nt)
      ff(m1a,jw) = ff(m1a,jw)-(q1a(j,i,ia)-wq1a(j,i,ia))
      ff(m2a,jw) = ff(m2a,jw)-(q2a(j,i,ia)-wq2a(j,i,ia))
 155  continue
 150  continue
!
      do 160 jw = jwst, jwen
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j = jcel(jw,it)
      j = jcel(jw,nt)
      ff(m3,jw) = ff(m3,jw)-(q3(j,i)-wq3(j,i))
      ff(m4,jw) = ff(m4,jw)-(q4(j,i)-wq4(j,i))
 160  continue
!
!-----------------------------------------------------------------------
!::boundary condition j = 1, jmax
!-----------------------------------------------------------------------
      if( ipbc.eq.0 ) then
        call pxbcnv(nt)
      else
        call pxbset(nt)
      endif
!
!-----------------------------------------------------------------------
!::solve with periodic boundary  exclude dummy cells
!-----------------------------------------------------------------------
      if( ipbc.eq.1 ) then
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik   nd1 = ndeq
!ik   nd2 = ndxy
!ik   neq = nequ
!ik   nmx = jmax
!ik   kdmy = 1
!
!::exclude dummy cells
      kdmy = 0
      do jw = 2, jmax1
      jw1 = jw - 1
!-----
!xx   cc0(1:ndeq,1:ndeq,jw1) = cc(1:ndeq,1:ndeq,jw)
!xx   dd0(1:ndeq,1:ndeq,jw1) = dd(1:ndeq,1:ndeq,jw)
!xx   ee0(1:ndeq,1:ndeq,jw1) = ee(1:ndeq,1:ndeq,jw)
!xx   ff0(1:ndeq,jw1) = ff(1:ndeq,jw)
!-----  KS 2011/09/22
! modified 4/4 lines organize local variables and include files by kamata 2021/05/31
!ik   cc0(1:neq,1:neq,jw1) = cc(1:neq,1:neq,jw)
!ik   dd0(1:neq,1:neq,jw1) = dd(1:neq,1:neq,jw)
!ik   ee0(1:neq,1:neq,jw1) = ee(1:neq,1:neq,jw)
!ik   ff0(1:neq,jw1) = ff(1:neq,jw)
      cc0(1:nequ,1:nequ,jw1) = cc(1:nequ,1:nequ,jw)
      dd0(1:nequ,1:nequ,jw1) = dd(1:nequ,1:nequ,jw)
      ee0(1:nequ,1:nequ,jw1) = ee(1:nequ,1:nequ,jw)
      ff0(1:nequ,jw1) = ff(1:nequ,jw)
!-----
      enddo
      nmx = jmax - 2
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   call tdmeqp( neq, nmx, kdmy, cc0, dd0, ee0, xx, ff0, ndeq, ndxy )
      call tdmeqp( nequ, nmx, kdmy, cc0, dd0, ee0, xx, ff0, ndeq, ndxy )
!
!-----
!xx      if( itim.eq.1000 .and. nlp.eq.5 ) then
!xx      if( nt.ge.itmps .and. mod(nt,5).eq.0 ) then
!xx      write(i6,'(/2x,"itim =",i7,"  nlp =",i3,"   nt,i =",2i4)')
!xx     >    itim, nlp, nt, i
!xx      call tdmchk( neq, nmx, kdmy, cc0, dd0, ee0, xx, ff0, nd1, nd2 )
!xx      endif
!xx      endif
!-----
      do jw1 = 1, nmx
      jw = jw1 + 1
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   ff(1:neq,jw) = xx(1:neq,jw1)  ! ndeq => neq  KS 2011/09/22
      ff(1:nequ,jw) = xx(1:nequ,jw1)  ! ndeq => nequ  KS 2011/09/22
      enddo
      else
!
!-----------------------------------------------------------------------
!::solve
!-----------------------------------------------------------------------
! modified 5/1 lines organize local variables and include files by kamata 2021/05/31
!ik   nd1  = ndeq
!ik   nd2  = ndxy
!ik   neq  = nequ
!ik   nmx  = jmax
!ik   call tdmslv( nd1, nd2, neq, nmx, cc, dd, ee, ff )
      call tdmslv( ndeq, ndxy, nequ, jmax, cc, dd, ee, ff )
      endif
!
!-----------------------------------------------------------------------
!::store the solution
!-----------------------------------------------------------------------
      jwst = 1
      jwen = jmax
      if( ipbc.eq.1 ) then
      jwst = 2
      jwen = jmax1
      endif
!
      do 520 jw = jwst, jwen
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j = jcel(jw,it)
      j = jcel(jw,nt)
      do 530 ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      dq1a(j,i,ia) = ff(m1a,jw)
      dq2a(j,i,ia) = ff(m2a,jw)
 530  continue
      dq3(j,i) = ff(m3,jw)
      dq4(j,i) = ff(m4,jw)
 520  continue
!
 100  continue
!$omp end parallel
!
 900  continue
!
!::qcon ! %%% 2002/11/3
      if( kpcn.eq.1 .and. nlp.eq.nlpmn+1 ) then
      do i = 1, 10
      it = itqcn(i)
      if( it.le.0 ) exit
      call plqcon_flx(it)
      call plqcon_pcn(it)
      call plqcon_chk(it)
!x    call plqcon_chk2(it)
      enddo
      endif
!
      return
      end
