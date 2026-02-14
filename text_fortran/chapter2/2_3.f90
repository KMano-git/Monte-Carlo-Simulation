program integral
   implicit none
   real(8) :: a, b, h, sum
   integer(8) :: i, n
   real(8), parameter :: pi = 4.0d0*datan(1.0d0)

   print *, 'Enter: n = '
   read *, n
   !sin(x)の積分 (a~b)
   a = 0.0d0
   b = pi
   h = (b-a)/dble(n)
   sum = 0.0d0
   do i = 1, n
      sum = sum + dsin(a + (i-1)*h)
   end do
   sum = sum + (dsin(a)+dsin(b)) / 2.0d0
   sum = sum * h
   print *, 'sin(x)の積分 = ', sum, '理論値 = ', 2.0

   !(x-2)^3の積分(a~b)
   a = 1.0d0
   b = 5.0d0
   h = (b-a)/dble(n)
   sum = 0.0d0
   do i = 1, n
      sum = sum + (a+(i-1)*h-2.0d0)**3
   end do
   sum = sum + ((a-2.0d0)**3 + (b-2.0d0)**3) / 2.0d0
   sum = sum * h
   print *, '(x-2)^3の積分 = ', sum, '理論値 = ', 20.0d0
   ! ちなみに、桁溢れなどを防ぐためには、全部足してから割るのではなく、割りながら足していくのが良い
   ! 計算速度は遅くなるけど。
end program integral
