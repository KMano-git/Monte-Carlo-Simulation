!***********************************************************************
      subroutine corona_3(rlos,fn0,ns,alf,bet,gam,anz)
!***********************************************************************
!
!      including recycling and CXR effects
!
!      a(i)*N(i-1) - b(i)*N(i) + c(i)*N(i+1) = d(i)
!                   i = 0, 1, 2, 3, .... N
!
!        a(0) = 0.0d0
!        c(N) = 0.0d0
!        d(i) = 0.0d0    d(0) = 1.0d0   source
!
!        alf : ionization
!        bet : recombination
!        gam : charge exchange recombination with H0
!        rlos :  1/(Ne*tauI)
!        fn0  :  NH0/Ne
!
!         corona_2         corona_3    2009/01/30
!         N(0) = 1.0  ==> S(0) = 1.0
!
!-----------------------------------------------------------------------
      implicit none
!
!::argument
      real*8, intent(in) :: rlos, fn0
      integer, intent(in) :: ns
      real*8, intent(in) :: alf(0:ns)
      real*8, intent(inout) :: bet(0:ns), gam(0:ns),anz(0:ns)
!
!::tdmslv
      integer, parameter:: ndim=80
      real*8 aa(0:ndim), bb(0:ndim), cc(0:ndim), dd(0:ndim)
!
!::local variables
      integer is, i, n
      real*8 sum

!::dimension check
      if( ns > ndim ) then
        write(6,'(/2x,"*** dimension error at sub. corona_3 ***")')
        write(6,'(2x,"ns > ndim ",2i6)') ns, ndim
        stop
      endif
!
!----------------------------------------------------------------------
!::rlos = 0
!----------------------------------------------------------------------
      if ( rlos<=0.0D0 ) then
        anz(0) = 1.0D0
        do is = 1, ns
          anz(is) = alf(is-1)/bet(is)*anz(is-1)
        enddo
        goto 100
      endif
!
!----------------------------------------------------------------------
!::rlos ne 0.0
!----------------------------------------------------------------------
      do i = 0, ndim
        aa(i) = 0.0D0
        bb(i) = 0.0D0
        cc(i) = 0.0D0
        dd(i) = 0.0D0
      enddo
!
      bet(0) = 0.0D0
      gam(0) = 0.0D0
!
      do i = 0, ns
        if ( i==0 ) then
          aa(i) = 0.0D0
        else
          aa(i) = alf(i-1)
        endif
        bb(i) = alf(i) + bet(i) + fn0*gam(i) + rlos
        if ( i==ns ) then
          cc(i) = 0.0D0
        else
          cc(i) = bet(i+1) + fn0*gam(i+1)
        endif
      enddo
      dd(0) = 1.0D0
!
      n = ns
      call tdmslv_3(n,aa,bb,cc,dd)
      do i = 0, ns
        anz(i) = dd(i)
      enddo
!
!----------------------------------------------------------------------
!::normalization
!----------------------------------------------------------------------
 100  continue
      sum = 0.0D0
      do is = 0, ns
        sum = sum + anz(is)
      enddo
      do is = 0, ns
        anz(is) = anz(is)/sum
      enddo
!
      continue
      end subroutine corona_3
!
!***********************************************************************
      subroutine tdmslv_3(n,a,b,c,d)
!***********************************************************************
!
!     solve   a(j)*u(j-1) - b(j)*u(j) + c(j)*u(j+1) = d(j)
!                    j = 0,1,2,3 ... n
!
!       a,b,c(mesh)  d(mesh)
!
!        n :  mesh number
!
!-----------------------------------------------------------------------
      implicit none
!
!::argument
      integer, intent(in) :: n
      real*8, intent(in) :: a(0:n), b(0:n), c(0:n)
      real*8, intent(inout) :: d(0:n)
!
!::local variables
      integer j
      real*8 zaa, e(0:110), f(0:110)
!
      e(n-1) = a(n)/b(n)
      f(n-1) = -d(n)/b(n)
      do j = n - 1, 1, -1
        zaa = b(j) - c(j)*e(j)
        e(j-1) = a(j)/zaa
        f(j-1) = (c(j)*f(j)-d(j))/zaa
      enddo
!
      d(0) = (c(0)*f(0)-d(0))/(b(0)-c(0)*e(0))
      do j = 1, n
        d(j) = e(j-1)*d(j-1) + f(j-1)
      enddo
!
      continue
      end subroutine tdmslv_3