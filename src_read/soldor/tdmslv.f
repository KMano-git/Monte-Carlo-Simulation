!***********************************************************************
!::test program  tdmslv
!***********************************************************************
      subroutine tst_tdmslv
      implicit none
      integer  j, jmax, md, nd, m, n, i6
      real*8   sum1, sum2
!
      real*8  a(5,5,10),b(5,5,10),c(5,5,10),d(5,10)
      real*8  x(5,10)
!
      data (x(1,j),x(2,j),j=1,5)
     >  /  1.0, 12.0, 0.8, 15.0, 0.3, 18.0, 0.9, 12.0, 0.2, 19.0/
      data (a(1,1,j),a(1,2,j),a(2,1,j),a(2,2,j),j=1,5)
     >  /  0.0, 0.0,  0.0, 0.0,
     >     0.8, 0.2, -0.1, 0.3,
     >     0.3, 0.1,  0.5, 0.2,
     >     0.2,-0.2,  0.4, 0.1,
     >     0.1, 0.6,  0.2, 0.3/
      data (b(1,1,j),b(1,2,j),b(2,1,j),b(2,2,j),j=1,5)
     >  /  2.0, 0.2,  0.3, 2.0,
     >     2.3, 0.1,  0.2, 2.5,
     >     2.4, 0.3,  0.1, 2.8,
     >     2.3, 0.1,  0.6, 2.5,
     >     2.7, 0.3,  0.2, 2.7/
      data (c(1,1,j),c(1,2,j),c(2,1,j),c(2,2,j),j=1,5)
     >  /  0.6, 0.8, 0.2, 0.3,
     >     0.4, 0.1, 0.4, 0.5,
     >     0.3, 0.6, 0.2, 0.4,
     >     0.4, 0.1, 0.1, 0.3,
     >     0.0, 0.0, 0.0, 0.0/
!
      jmax=5
!
      j=1
      sum1=
     >      +b(1,1,j)*x(1,j)  +b(1,2,j)*x(2,j)
     >      +c(1,1,j)*x(1,j+1)+c(1,2,j)*x(2,j+1)
      sum2=
     >      +b(2,1,j)*x(1,j)  +b(2,2,j)*x(2,j)
     >      +c(2,1,j)*x(1,j+1)+c(2,2,j)*x(2,j+1)
      d(1,1)=sum1
      d(2,1)=sum2
      do 110 j=2,4
      sum1=  a(1,1,j)*x(1,j-1)+a(1,2,j)*x(2,j-1)
     >      +b(1,1,j)*x(1,j)  +b(1,2,j)*x(2,j)
     >      +c(1,1,j)*x(1,j+1)+c(1,2,j)*x(2,j+1)
      sum2=  a(2,1,j)*x(1,j-1)+a(2,2,j)*x(2,j-1)
     >      +b(2,1,j)*x(1,j)  +b(2,2,j)*x(2,j)
     >      +c(2,1,j)*x(1,j+1)+c(2,2,j)*x(2,j+1)
      d(1,j)=sum1
      d(2,j)=sum2
  110 continue
      j=jmax
      sum1=  a(1,1,j)*x(1,j-1)+a(1,2,j)*x(2,j-1)
     >      +b(1,1,j)*x(1,j)  +b(1,2,j)*x(2,j)
      sum2=  a(2,1,j)*x(1,j-1)+a(2,2,j)*x(2,j-1)
     >      +b(2,1,j)*x(1,j)  +b(2,2,j)*x(2,j)
      d(1,j)=sum1
      d(2,j)=sum2
!
      md = 5
      nd = 10
      m  = 2
      n  = jmax
      call tdmslv( md, nd, m, n, a, b, c, d)
!
      i6 = 6
      do 210 j=1,jmax
      write(i6,602) j,x(1,j),x(2,j),d(1,j),d(2,j)
  602 format(2x,i2,2x,1p2e18.10,2x,1p2e18.10)
  210 continue
!
      stop   ! prog. test_tdmslv
      end
!
!***********************************************************************
      subroutine tdmslv( md, nd, m, n, a, b, c, d)
!***********************************************************************
!
!     solve   a(j)*u(j-1) + b(j)*u(j) + c(j)*u(j+1) = d(j)
!                    j = 1,2,3 ... n
!
!       a,b,c(equ,var,mesh)  d(equ,mesh)
!
!        m :  the number ot equations
!        m :  the number of variables
!        n :  mesh number
!
!-----------------------------------------------------------------------
      implicit none
!
! modified 2/3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer md, nd, m, n
!ik   real*8  a(md,md,nd),b(md,md,nd),c(md,md,nd),d(md,nd)
      integer, intent(in) :: md, nd, m, n
      real(8), intent(inout) :: a(md,md,nd), b(md,md,nd),c(md,md,nd)
     >       , d(md,nd)
!
!::local variables
      integer i, j, k, l
      real*8  x
      real*8 :: zero = 1.0e-30

!
!-----------------------------------------------------------------------
!::lu decomposition
!-----------------------------------------------------------------------
!
      do 510 i=1,n
      do 500 j=1,m
!
!::a(i,j,k)
!
      if(i.gt.1)then
      do 110 k=2,m
      x=0.
      do 100 l=1,k-1
 100  x=x+a(j,l,i)*b(l,k,i-1)
 110  a(j,k,i)=a(j,k,i)-x
      endif
!
!::b(j,k,i) : k <= j
!
      do k=1,j
        x=0.
        if(i.gt.1)then
          do l=1,m
            x=x+a(j,l,i)*c(l,k,i-1)
          enddo !l
        endif
        if(k.gt.1)then
          do l=1,k-1
            x=x+b(j,l,i)*b(l,k,i)
          enddo !l
        endif
        b(j,k,i)=b(j,k,i)-x ! b
        if(abs(b(j,k,i))<zero) then
          b(j,k,i) = 0.0d0
        endif
      enddo !k
      

!
!::b(j,k,i) : k > j
!
      if(j.lt.m)then
        do k=j+1,m
          x=0.
          if(i.gt.1)then
            do l=1,m
              x=x+a(j,l,i)*c(l,k,i-1)
            enddo !l
          endif
          do l=1,j-1
            x=x+b(j,l,i)*b(l,k,i)
          enddo !l
          b(j,k,i)=(b(j,k,i)-x)/b(j,j,i)
          if(abs(b(j,k,i))<zero) then
            b(j,k,i)=0.0d0
          endif
        enddo !k
      endif
!
!::c(i,j,k)
!
      if(i.lt.n)then
        do k=1,m
          x=0.0d0
          if(j.gt.1)then
            do l=1,j-1
              x=x+b(j,l,i)*c(l,k,i)
            enddo
          endif
          c(j,k,i)=(c(j,k,i)-x)/b(j,j,i)
        enddo
      endif
!
 500  continue
 510  continue
!
!::d(j,i)
!
      do 630 i=1,n
      do 620 j=1,m
      x=0.
      if(i.gt.1)then
      do 600 k=1,m
 600  x=x+a(j,k,i)*d(k,i-1)
      endif
!
      if(j.gt.1)then
      do 610 k=1,j-1
 610  x=x+b(j,k,i)*d(k,i)
      endif
      d(j,i)=(d(j,i)-x)/b(j,j,i)
 620  continue
 630  continue
!
!-----------------------------------------------------------------------
!::back substitution
!-----------------------------------------------------------------------
!
      do 730 i=n,1,-1
      do 720 j=m,1,-1
      x=0.
      if(i.lt.n)then
      do 700 k=1,m
 700  x=x+c(j,k,i)*d(k,i+1)
      endif
      if(j.lt.m)then
      do 710 k=j+1,m
 710  x=x+b(j,k,i)*d(k,i)
      endif
      d(j,i)=d(j,i)-x
 720  continue
 730  continue
!
      return
      end
