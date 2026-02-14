!***********************************************************************
      subroutine test_tdmeqp
!***********************************************************************
!
!  test program of tdmeqp (cyclic Tridiagonal System)
!       "p" means periodic boundary condition
!
!  /home/shimizu/sonicV1old/src/ssl2
!   Oct 30  2010  after IAEA conference
!   Sun Jun 2011  restart
!
!-----------------------------------------------------------------------
      implicit none
!
      integer  mdsz, ndsz
      integer  mvar, nmsh
      parameter (mdsz=10,ndsz=121)
      real*8, dimension(mdsz,mdsz,ndsz) :: cc, dd, ee
      real*8, dimension(mdsz,ndsz)      :: xx, ff
!
      real*8, dimension(mdsz,mdsz)      :: c, d, e
      real*8, dimension(mdsz)           :: x
!
      data  c(1:4,1) /-1.0d0,  0.5d0, -0.2d0, 0.8d0/
      data  c(1:4,2) /0.3d0,   1.5d0,  0.3d0, 0.2d0/
      data  c(1:4,3) /-0.2d0,  1.1d0,  1.2d0, 2.3d0/
      data  c(1:4,4) /0.7d0,   0.8d0,  0.5d0, 0.7d0/
!
      data  d(1:4,1) /6.0d0,  0.5d0, -0.2d0, 0.8d0/
      data  d(1:4,2) /0.3d0,   6.5d0,  0.3d0, 0.2d0/
      data  d(1:4,3) /-0.2d0,  1.1d0,  5.2d0, 2.3d0/
      data  d(1:4,4) /0.7d0,   0.8d0,  0.5d0, 4.0d0/
!
      data  e(1:4,1) /-1.0d0,  0.5d0, -0.2d0, 0.8d0/
      data  e(1:4,2) /0.3d0,   1.5d0,  0.3d0, 0.2d0/
      data  e(1:4,3) /-0.2d0,  1.1d0,  1.2d0, 2.3d0/
      data  e(1:4,4) /0.7d0,   0.8d0,  0.5d0, 0.7d0/
!
      data  x(1:4) /13.0d0, 14.0d0,  15.0d0, 16.0d0/
!
      integer  i, n, kdmy, nb1, nb2
      real*8   tx, err
      real*8, dimension(mdsz,ndsz) :: xx0
      real*8, dimension(mdsz)      :: w1, w2, w3
!
      mvar = 4
      nmsh = 10
!
      cc(1:mvar,1:mvar,1:nmsh) = 0.0d0
      dd(1:mvar,1:mvar,1:nmsh) = 0.0d0
      ee(1:mvar,1:mvar,1:nmsh) = 0.0d0
      xx(1:mvar,1:nmsh) = 0.0d0
      ff(1:mvar,1:nmsh) = 0.0d0
!
!::tdma+periodic
!xx   write(6,'(/2x,"birth sol x =",1p4e12.4)') x(1:4)
!
      do n = 1, nmsh
      tx = -n
      cc(1:mvar,1:mvar,n) = c(1:mvar,1:mvar)
      dd(1:mvar,1:mvar,n) = d(1:mvar,1:mvar)
      ee(1:mvar,1:mvar,n) = e(1:mvar,1:mvar)
      xx(1:mvar,n) = x(1:mvar) + tx
      enddo
!
      xx0(1:mvar,1:nmsh) = xx(1:mvar,1:nmsh)
!
!::dummy cells
      kdmy = 1
      nb1 = 1
      nb2 = nmsh
      if( kdmy.eq.1 ) then
      nb1 = 2
      nb2 = nmsh-1
      xx0(1:mvar,1) = xx0(1:mvar,nb2)
      xx0(1:mvar,nmsh) = xx0(1:mvar,nb1)
      endif
!
!::debug
!xx      write(6,'(/2x,"solution")')
!xx      do n = 1, nmsh
!xx      write(6,'(2x,"n =",i3,"  xx0 =",1p4e12.4)')
!xx     >   n, (xx0(i,n),i=1,mvar)
!xx      enddo
!
!::set F = C*X(n-1)+D*X(n)+E*X(n+1)
!--------------------------------------------------
      do n = nb1, nb2
        w1(1:mvar) = 0.0d0
        w2(1:mvar) = 0.0d0
        w3(1:mvar) = 0.0d0
        if( n.ne.nb1 ) then
          call mvmlt(cc(1,1,n),xx(1,n-1),w1,mvar,mdsz)
        else
          call mvmlt(cc(1,1,n),xx(1,nb2),w1,mvar,mdsz) ! p-BC (n=1)
        endif
          call mvmlt(dd(1,1,n),xx(1,n)  ,w2,mvar,mdsz)
        if( n.ne.nb2 ) then
          call mvmlt(ee(1,1,n),xx(1,n+1),w3,mvar,mdsz)
        else
          call mvmlt(ee(1,1,n),xx(1,nb1),w3,mvar,mdsz) ! p-BC (n=nmsh)
        endif
        ff(1:mvar,n) = w1(1:mvar)+w2(1:mvar)+w3(1:mvar)
      enddo
!--------------------------------------------------
!
      xx(1:mvar,1:ndsz) = 1.234567d0  ! for check
!
      call tdmeqp( mvar, nmsh, kdmy, cc, dd, ee, xx, ff, mdsz, ndsz )
      call tdmchk( mvar, nmsh, kdmy, cc, dd, ee, xx, ff, mdsz, ndsz )
!
!::debug write
      write(6,'(/2x,"test_tdmeqp  2011/06/08")')
      write(6,'(2x,"mvar/mdsz =",2i4,"  nmsh/ndsz =",2i4)')
     >   mvar, mdsz, nmsh, ndsz
      do n = 1, nmsh
      write(6,'(/2x,"n =",i3)') n
      do i = 1, mvar
      err = dabs(xx(i,n)-xx0(i,n))
      write(6,'(2x,"x0 =",1pe12.4, "  xx =",1pe12.4, "  err =",
     >   1pe12.4)') xx0(i,n), xx(i,n), err
      enddo
      enddo
!
      return
      end
!
!***********************************************************************
      subroutine tdmeqp( mvar,nmsh,kdmy, c0, d0, e0, x0, f0, mdsz,ndsz)
!***********************************************************************
!
!   periodic Boudary Condition   with No dummy cells   kdmy = 0
!
!        D(1)*X(1) + E(1)*X(2) + C(1)*X(J) = F(1)
!        C(2)*X(1) + D(2)*X(2) + E(2)*X(3) = F(2)
!           ............
!        C(n)*X(n-1) + D(n)*X(n) + E(n)*X(n+1) = F(n)  n = 2, N-1
!           ...........
!        E(N)*X(1) + C(N)*X(N-1) + D(N)*X(N) = F(N)
!
!
!   periodic Boudary Condition   with dummy cells   kdmy = 1
!
!        X(1) = X(N-1)
!        C(2)*X(1) + D(2)*X(2) + E(2)*X(3) = F(2)
!           ............
!        C(n)*X(n-1) + D(n)*X(n) + E(n)*X(n+1) = F(n)  n = 2, N-1
!           ...........
!        X(N) = X(2)
!
!
!       C,D,E(equ,var,mesh)
!       X(var,mesh), F(var,mesh)
!
!       X : vector X(1:mvar,1,1:nmsh)  kvec = 1
!       X : multiple vector  X(1:mvar,1:kvec,1:nmesh)
!
!        mvar :  the number ot equations
!        mvar :  the number of variables
!        nmsh :  mesh number
!        kdmy :  = 0 (No dummy cell) = 1 (dummy cell)
!
!        kvec :  the number of vector
!
!::note
!     LU solve   A*x = b     A:destory
!     tdmeqp   C*X(n-1)+D*X(n)+E*X(n+1) = F   D:destory
!
!-----------------------------------------------------------------------
      implicit none
!::argument
! modified 3/4 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  mvar, nmsh, kdmy, mdsz, ndsz
!ik   real*8, dimension(mdsz,mdsz,ndsz) :: c0, d0, e0
!ik   real*8, dimension(mdsz,ndsz)      :: x0, f0
      integer, intent(in) :: mvar, nmsh, kdmy, mdsz, ndsz
      real(8), dimension(mdsz,mdsz,ndsz), intent(in) :: c0, d0, e0
      real(8), dimension(mdsz,ndsz), intent(in)      :: f0
      real(8), dimension(mdsz,ndsz), intent(out)     :: x0
!
!::local variables
      real*8, dimension(mvar,mvar,ndsz) :: cc, dd, ee
      real*8, dimension(mvar,ndsz)      :: xx, ff
!
!::local variables (LU routine)
      real*8   x(mvar), b(mvar)
      real*8   det
      integer  nmax, mm(mvar)
      real*8, dimension(mvar,mvar) :: wi, wh, wa
!
!::local variables (TDMA routine)
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8, dimension(mvar,mvar,ndsz) :: wx, wy, wz, wf, wu
      real*8, dimension(mvar,mvar,ndsz) :: wy, wz, wf, wu
      real*8, dimension(mvar) :: why
      integer  kvec
!
!::local variables
      integer  i, k, m, n, nb1, nb2
      real*8   sum
! deleted 1 line deleting unused common variables in SONIC by kamata 2021/05/31
!ik   real*8, dimension(mvar,mvar) :: wa0
!
!::dimension size
      if( nmsh.gt.ndsz ) then
      write(6,'(/2x,"stop at sub. tdmequ  nmsh.gt.ndsz  ",
     >   2i5,2x,2i5)')  nmsh, ndsz
      call wexit("tdmeqp","nmsh.gt.ndsz")
      endif
!
!::dummy cell
      nb1 = 1
      nb2 = nmsh
      if( kdmy.eq.1 ) then
      nb1 = 2
      nb2 = nmsh-1
      endif
!
!::coefficients   ndsz => nmsh   KS 2011/09/22
      cc(1:mvar,1:mvar,1:nmsh) = c0(1:mvar,1:mvar,1:nmsh)
      dd(1:mvar,1:mvar,1:nmsh) = d0(1:mvar,1:mvar,1:nmsh)
      ee(1:mvar,1:mvar,1:nmsh) = e0(1:mvar,1:mvar,1:nmsh)
      ff(1:mvar,1:nmsh) = f0(1:mvar,1:nmsh)
!
!::D(1)=D(1)-C(1)
!::D(N)=D(N)-E(N)
      n = nb1
      dd(1:mvar,1:mvar,n) = dd(1:mvar,1:mvar,n) - cc(1:mvar,1:mvar,n)
      n = nb2
      dd(1:mvar,1:mvar,n) = dd(1:mvar,1:mvar,n) - ee(1:mvar,1:mvar,n)
!
!::U
      wu(1:mvar,1:mvar,1:nmsh) = 0.0d0
      n = nb1
      wu(1:mvar,1:mvar,n) = cc(1:mvar,1:mvar,n)
      n = nb2
      wu(1:mvar,1:mvar,n) = ee(1:mvar,1:mvar,n)
!
!::C(1)=0
!::E(N)=0
      n = nb1
      cc(1:mvar,1:mvar,n) = 0.0d0
      n = nb2
      ee(1:mvar,1:mvar,n) = 0.0d0
!
!::A*Y = F   Y(i,n), F(i,n)
      kvec = 1
      wf(1:mvar,1,nb1:nb2) = ff(1:mvar,nb1:nb2)
!
!xx   write(6,'(/2x,"A*Y-F=0  kvec=",i3)') kvec
      call tdmsol( mvar,kvec,nb1,nb2, cc, dd, ee, wy, wf, ndsz )
!
!::A*Z = U
      kvec = mvar
!xx   write(6,'(/2x,"A*Z-U=0  kvec=",i3)') kvec
      call tdmsol( mvar,kvec,nb1,nb2, cc, dd, ee, wz, wu, ndsz )
!
!::(I+Z1+ZN)*H = I
!xx   write(6,'(2x)')
      wi(1:mvar,1:mvar) = 0.0d0
      do i = 1, mvar
      wi(i,i) = 1.0d0
      enddo
      wa(1:mvar,1:mvar) = wi(1:mvar,1:mvar)
     >    + wz(1:mvar,1:mvar,nb1) + wz(1:mvar,1:mvar,nb2)
! deleted 1 line deleting unused common variables in SONIC by kamata 2021/05/31
!ik   wa0(1:mvar,1:mvar) = wa(1:mvar,1:mvar)
!
      nmax = mvar
      call lu ( wa, det, nmax, mm )   ! Note
!
      do k = 1, mvar
      b(1:mvar) = 0.0d0
      b(k) = 1.0d0
      call solve( wa, b, x, nmax, mm )
      wh(1:mvar,k) = x(1:mvar)
      enddo
!
!::check  Wa*Wh-I=0
!xx      call mtmlt(wa0, wh, wb0, mvar, mvar)
!xx      do i = 1, mvar
!xx      do j = 1, mvar
!xx      if( j.ne.i ) cycle
!xx      wb0(i,i) = wb0(i,j) - 1.0d0
!xx      enddo
!xx      enddo
!xx      call mtout("Wa*Wh-I = 0", wb0, mvar)
!
!:X(n) = Y(n) - Z(n)*H*(Y(1)+Y(n))
      do i = 1, mvar
        sum = 0.0d0
        do m = 1, mvar
          sum = sum + wh(i,m)*(wy(m,1,nb1)+wy(m,1,nb2))
        enddo
      why(i) = sum
      enddo
!
      do n = nb1, nb2
      do i = 1, mvar
        sum = 0.0d0
        do m = 1, mvar
          sum = sum + wz(i,m,n)*why(m)
        enddo
      xx(i,n) = wy(i,1,n) - sum
      enddo
      enddo
!
!::solution
      do n = nb1, nb2
      x0(1:mvar,n) = xx(1:mvar,n)
      enddo
!
      if( kdmy.eq.1 ) then
      x0(1:mvar,1) = x0(1:mvar,nb2)
      x0(1:mvar,nmsh) = x0(1:mvar,nb1)
      endif
!
      return
      end
!
!***********************************************************************
      subroutine tdmchk( mvar,nmsh,kdmy, c0, d0, e0, x0, f0, mdsz,ndsz)
!***********************************************************************
      implicit none
!
!::argument
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  mvar, nmsh, kdmy, mdsz, ndsz
!ik   real*8, dimension(mdsz,mdsz,ndsz) :: c0, d0, e0
!ik   real*8, dimension(mdsz,ndsz)      :: x0, f0
      integer, intent(in) :: mvar, nmsh, kdmy, mdsz, ndsz
      real(8), dimension(mdsz,mdsz,ndsz), intent(in) :: c0, d0, e0
      real(8), dimension(mdsz,ndsz), intent(in)      :: x0, f0
!
!::local variables
      integer  i, n, nb1, nb2
      real*8, dimension(mdsz) :: w1, w2, w3, w4
      real*8  :: emx
      integer :: i6 = 80000 + 30
!
      nb1 = 1
      nb2 = nmsh
      if( kdmy.eq.1 ) then
      nb1 = 2
      nb2 = nmsh-1
      endif
!
!::C*X(n-1)+D*X(n)+E*X(n+1)-F=0
!
      write(i6,'(2x,"res  C*X(n-1)+D*X(n)+E*X(n+1)-F = 0  kdmy =",
     >   i2)') kdmy
!
      do n = nb1, nb2
        w1(1:mdsz) = 0.0d0
        w2(1:mdsz) = 0.0d0
        w3(1:mdsz) = 0.0d0
        if( n.ne.nb1 ) then
           call mvmlt(c0(1,1,n),x0(1,n-1),w1,mvar,mdsz)
        else
           call mvmlt(c0(1,1,n),x0(1,nb2),w1,mvar,mdsz)
        endif
           call mvmlt(d0(1,1,n),x0(1,n)  ,w2,mvar,mdsz)
        if(n.ne.nb2) then
           call mvmlt(e0(1,1,n),x0(1,n+1),w3,mvar,mdsz)
        else
           call mvmlt(e0(1,1,n),x0(1,nb1),w3,mvar,mdsz)
        endif
!
        emx = -1.0d20
        do i = 1, mvar
          w4(i) = w1(i)+w2(i)+w3(i)-f0(i,n)
          emx = dmax1( emx, dabs(w4(i)) )
        enddo
!
      write(i6,'(2x,"tdmchk n =",i3,"  x =",1p4e12.4,"  res =",1p4e12.4,
     >  "  emx =",1pe12.4)') n, x0(1:4,n), w4(1:4), emx
      enddo
!
      if( kdmy.eq.1 ) then
      n = 1
      w4(1:mvar) = x0(1:mvar,n) - x0(1:mvar,nb2)
      write(i6,'(2x,"tdmchk n =",i3,"  x =",1p4e12.4,"  res =",1p4e12.4,
     >  "  emx =",1pe12.4)') n, x0(1:4,n), w4(1:4)
      n = nmsh
      w4(1:mvar) = x0(1:mvar,n) - x0(1:mvar,nb1)
      write(i6,'(2x,"tdmchk n =",i3,"  x =",1p4e12.4,"  res =",1p4e12.4,
     >  "  emx =",1pe12.4)') n, x0(1:4,n), w4(1:4)
      endif
!
      return
      end
!
!
!***********************************************************************
      subroutine tdmsol( mvar,kvec,nb1,nb2, cc,dd,ee,xx,ff,ndsz )
!***********************************************************************
!
!     solve   C(n)*X(n-1) + D(n)*X(n) + E(n)*X(n+1) = F(n)
!                                 mesh      n = 1,2,3 ... nmsh
!
!       C,D,E(equ,var,mesh)
!       X(var,vec,mesh), F(var,vec,mesh)
!
!       X : vector X(1:mvar,1,1:nmsh)  kvec = 1
!       X : multiple vector  X(1:mvar,1:kvec,1:nmesh)
!
!        mvar :  the number ot equations
!        mvar :  the number of variables
!        kvec :  the number of vector
!        nmax :  mesh number
!
!     Note.  size = mvar
!
!-----------------------------------------------------------------------
      implicit none
!
!::argument
! modified 3/4 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  mvar, kvec, nb1, nb2, ndsz
!ik   real*8, dimension(mvar,mvar,ndsz) :: cc, dd, ee
!ik   real*8, dimension(mvar,mvar,ndsz) :: xx, ff
      integer, intent(in) :: mvar, kvec, nb1, nb2, ndsz
      real*8, dimension(mvar,mvar,ndsz), intent(in)  :: cc, dd, ee
      real*8, dimension(mvar,mvar,ndsz), intent(in)  :: ff
      real*8, dimension(mvar,mvar,ndsz), intent(out) :: xx
!
!::local variables (LU subroutine)
      real*8   a(mvar,mvar), x(mvar), b(mvar)
      real*8   det
      integer  nmax, mm(mvar)
!
!::TDMA
      real*8, dimension(mvar,mvar,ndsz) :: p
      real*8, dimension(mvar,mvar,ndsz) :: q
!
!::local variables
      integer  i, j, k, m, n, nmsh
      real*8   sum
!
!::dimension size
      nmsh = nb2
      if( nb1.eq.2 ) nmsh = nb2 + 2
      if( nmsh.gt.ndsz ) then
      write(6,'(/2x,"stop at sub. tdmsol  nmsh.gt.ndsz  ",
     >   2i5,2x,2i5)')  nmsh, ndsz
      call wexit("tdmsol","nmsh.gt.ndsz")
      endif
      if( kvec.ne.1 .and. kvec.ne.mvar ) then
      write(6,'(/2x,"stop at sub. tdmsol  kvec,mvar =",2i4)')
     >     kvec, mvar
      call wexit("tdmsol","kvec.ne.mvar")
      endif
!
!::solve D1*P1 =-E1
!        D1*Q1 = F1
      n = nb1
      a(1:mvar,1:mvar) = dd(1:mvar,1:mvar,n)
      nmax = mvar
      call LU ( A, DET, NMAX, MM )
!
      do k = 1, mvar
      b(1:mvar) = -ee(1:mvar,k,n)
      call SOLVE( A, B, X, NMAX, MM )
      p(1:mvar,k,n) = x(1:mvar)
      enddo
!
      do k = 1, kvec
        b(1:mvar) = ff(1:mvar,k,n)
        call SOLVE( A, B, X, NMAX, MM )
        q(1:mvar,k,n) = x(1:mvar)
      enddo
!
!::loop  (n)
      do n = nb1+1, nb2
!xx   write(6,'(/2x,"--- n = ",i3," ----")') n
!
!::solve (C(n)*P(n-1)+D(n))*P(n) = -E(n)
!        (C(n)*P(n-1)+D(n))*Q(n) = F(n)-C(n)*Q(n-1)
!                          (i,k,n) = (i,j,n)*(j,k,n) n:msh
!
      do i = 1, mvar
      do j = 1, mvar
        sum = 0.0d0
        do m = 1, mvar
          sum = sum + cc(i,m,n)*p(m,j,n-1)
        enddo
        a(i,j) = sum + dd(i,j,n)
      enddo
      enddo
!
      nmax = mvar
      call LU ( A, DET, NMAX, MM )
!
      do k = 1, mvar
        b(1:mvar) = -ee(1:mvar,k,n)
        call SOLVE( A, B, X, NMAX, MM )
        p(1:mvar,k,n) = x(1:mvar)
      enddo
!
!--------------------------------------------------
!::check (C(n)*P(n-1)+D(n))*P(n) + E(n) = 0
!xx   call mtmlt(cc(1,1,n),p(1,1,n-1),w1,mvar)
!xx   call mtadd(w1,dd(1,1,n),w2,mvar)
!xx   call mtmlt(w2,p(1,1,n), w1,mvar)
!xx   call mtadd(w1,ee(1,1,n),w2,mvar)
!xx   call mtout("(C(n)*P(n-1)+D(n))*P(n)+E(n)=0",
!xx  >                        w2,mvar)
!--------------------------------------------------
!
      do k = 1, kvec
      do i = 1, mvar
        sum = 0.0d0
        do m = 1, mvar
          sum = sum + cc(i,m,n)*q(m,k,n-1)
        enddo
        b(i) = ff(i,k,n) - sum
      enddo
        call SOLVE( A, B, X, NMAX, MM )
        q(1:mvar,k,n) = x(1:mvar)
      enddo
!
!--------------------------------------------------
!::check (C(n)*P(n-1)+D(n))*Q(n) - F(n)+C(n)*Q(n-1) = 0
!xx   call mtmlt(cc(1,1,n),p(1,1,n-1),w1,mvar)
!xx   call mtadd(w1,dd(1,1,n),w2,mvar)
!xx   call mtmlt(w2,q(1,1,n), w1,mvar)
!xx   call mtsub(w1,ff(1,1,n),w2,mvar)
!xx   call mtmlt(cc(1,1,n),q(1,1,n-1),w3,mvar)
!xx   call mtadd(w2,w3,w1,mvar)
!xx   call mtout("(C*P(n-1)+D)*Q(n)-F+C*Q(n-1)=0",
!xx  >                        w1,mvar)
!--------------------------------------------------
!
      enddo   ! loop(n)
!xx   write(6,'(2x,50("-"))')
!
!::X(nmsh)
      n = nb2
      xx(1:mvar,1:kvec,n) = q(1:mvar,1:kvec,n)
!
!::X(n-1) = P(n-1)*X(n) + Q(n-1)
      do n = nb2, nb1+1, -1
        do k = 1, kvec
        do i = 1, mvar
          sum = 0.0d0
          do m = 1, mvar
            sum = sum + p(i,m,n-1)*xx(m,k,n)
          enddo
          xx(i,k,n-1) = sum + q(i,k,n-1)
        enddo
        enddo
      enddo
!
      return
      end
!
!******************************************************************
      subroutine LU ( A, DET, N, MM )
!******************************************************************
!
!     LU factorization method of solution of a system of
!     linear equations (A) * X = B
!
!    usage: call LU ( A, DET, N, MM )
!           call SOLVE ( A, B, X, N, MM )
!
!    arguments  A ...... (N*N) matrix
!               B ...... constant vector
!               X ...... solution vector
!               DET .... determinant
!               N ...... number of unknowns
!               MM ...... order of pivotting
!               NDIM .... size of A in calling programme
!
!-----------------------------------------------------------------
!
      implicit none

!::argument
      integer, intent(in) :: N
      real*8, intent(inout) :: A(N,N)
      real*8, intent(out) :: DET
      integer, intent(inout) :: MM(N)

!::local variables
      integer i,k,L,j,IW
      real*8  B(N),X(N),PIVOT,W
      real*8 :: zero = 1.0d-50
!
! LU factorization
      do i = 1, N
         MM(i) = i
      end do
      DET = 1.0d0
!
      do k = 1, N-1
!
!:: partial pivotting
      L =k
      PIVOT = A(k,k)
      do j = k+1, N
         if ( abs(PIVOT)<abs(A(j,k)) ) then
            L = j
            PIVOT = A(j,k)
         end if
      end do
!
      if(L /= k) then
        DET = - DET
        do j = 1, N
           W = A(k,j)
           A(k,j) = A(L,j)
           A(L,j) = W
        end do
        IW = MM(k)
        MM(k) = MM(L)
        MM(L) = IW
      end if
!::end of pivotting
!
      do j = k+1, N
         A(k,j) = A(k,j)/PIVOT
      end do
      do i = k+1, N
         do j=k+1, N
            A(i,j) = A(i,j) - A(i,k) * A(k,j)
         end do
      end do
      end do
!
!::determinant
      do k = 1, N
         DET = DET * A(k,k)
      end do
!
      return
!------------------------------------------------------------------
!
! solution by forward and backward substitution
!
      entry SOLVE( A, B, X, N, MM )
!
!::forward substitution
      do k = 1, N
! for avoid underflow error
         if( abs(B(MM(k))) < zero) then
           X(k) = 0.0d0
         else
           X(k) = B(MM(k))
         endif
         do j = 1, k-1
            X(k) = X(k) - A(k,j) * X(j)
         end do
         X(k) = X(k) / A(k,k)
      end do
!
!::backward substitution
      do k = N, 1, -1
      do j = k+1, N
         X(k) = X(k) - A(k,j) * X(j)
      end do
      end do
!
      return
      end
!
!******************************************************************
      subroutine mtmlt(A,B,X,msz,ksz)
!******************************************************************
      implicit none
!
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  msz, ksz
!ik   real*8   A(msz,msz)
!ik   real*8   B(msz,ksz),X(msz,ksz)
      integer, intent(in)  :: msz, ksz
      real(8), intent(in)  :: A(msz,msz), B(msz,ksz)
      real(8), intent(out) :: X(msz,ksz)
!
!::local variables
      integer  i, k, m
      real*8   sum
!
      X(1:msz,1:ksz) = 0.0d0
!
      do i = 1, msz
      do k = 1, ksz
        sum = 0.0d0
        do m = 1, msz
          sum = sum + A(i,m)*B(m,k)
        enddo
        X(i,k) = sum
      enddo
      enddo
!
      return
      end
!
!******************************************************************
      subroutine mvmlt(A,B,X,n,msz)
!******************************************************************
      implicit none
!
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  n, msz
!ik   real*8   A(msz,msz)
!ik   real*8   B(msz),X(msz)
      integer, intent(in)  :: n, msz
      real(8), intent(in)  :: A(msz,msz), B(msz)
      real(8), intent(out) :: X(msz)
!
!::local variables
      integer  i, m
      real*8   sum
!
      X(1:msz) = 0.0d0
!
      do i = 1, n
        sum = 0.0d0
        do m = 1, n
          sum = sum + A(i,m)*B(m)
        enddo
        X(i) = sum
      enddo
!
      return
      end
!
!******************************************************************
      subroutine mtadd(A,B,X,msz)
!******************************************************************
      implicit none
!
! modified 2/3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  msz
!ik   real*8  A(msz,msz),B(msz,msz),X(msz,msz)
      integer, intent(in)  :: msz
      real(8), intent(in)  :: A(msz,msz), B(msz,msz)
      real(8), intent(out) :: X(msz,msz)
!
!::local variables
      integer  i, j
!
      X(1:msz,1:msz) = 0.0d0
!
      do i = 1, msz
      do j = 1, msz
        X(i,j) = A(i,j) + B(i,j)
      enddo
      enddo
!
      return
      end
!
!******************************************************************
      subroutine mtsub(A,B,X,msz)
!******************************************************************
      implicit none
!
! modified 2/3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  msz
!ik   real*8  A(msz,msz),B(msz,msz),X(msz,msz)
      integer, intent(in)  :: msz
      real(8), intent(in)  :: A(msz,msz), B(msz,msz)
      real(8), intent(out) :: X(msz,msz)
!
!::local variables
      integer  i, j
!
      X(1:msz,1:msz) = 0.0d0
!
      do i = 1, msz
      do j = 1, msz
        X(i,j) = A(i,j) - B(i,j)
      enddo
      enddo
!
      return
      end
!
!******************************************************************
      subroutine mtout(cmsg,X,msz)
!******************************************************************
      implicit none
!
!::argument
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   character  cmsg*(*)
!ik   integer  msz
!ik   real*8   X(msz,msz)
      character, intent(in) :: cmsg*(*)
      integer,   intent(in) :: msz
      real(8),   intent(in) :: X(msz,msz)
!
!::local variables
      integer  i, j
      real*8   xmax
!
      xmax = -1.0d20
      do i = 1, msz
      do j = 1, msz
        xmax = dmax1(xmax, dabs(x(i,j)))
      enddo
      enddo
!
      write(6,'(2x,a,"  max =",1pe12.4)') trim(cmsg), xmax
!
      do i = 1, msz
      write(6,'(3x,1p10e12.4)') (x(i,j),j=1,msz)
      enddo
!
      return
      end
