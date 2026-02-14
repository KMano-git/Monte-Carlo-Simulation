program newton_method
   implicit none
   real(8) :: x, x0, epsilon
   integer :: i, n

   print *, 'ニュートン法の計算 x = '
   read *, x0

   epsilon = 1.0d-7
   ! f(x) = x - cos(x)のニュートン法
   ! f'(x) = 1 + sin(x)
   x = x0
   n = 100000
   do i = 1, n
      x = x - (x - cos(x))/(1+sin(x))
      if (abs((x - cos(x))/(1+sin(x))) < epsilon) exit
   end do
   if (i == n) print *, '収束しません'

   print "(A, F15.8)", 'ニュートン法の解：', x
   print "(A, I15)", 'ニュートン法の反復回数：', i

end program newton_method
