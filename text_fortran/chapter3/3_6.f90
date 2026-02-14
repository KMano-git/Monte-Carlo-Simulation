program main
   implicit none
   real(8) :: x, f, f_prime, epsilon
   integer :: i, n

   ! ニュートン法の計算、2_8.f90の改良版
   ! 使いたい関数に応じて、subroutine fun_diffを書き換える
   print *, 'ニュートン法の計算 x = '
   read *, x

   epsilon = 1.0d-7
   n = 100000
   do i = 1, n
      call fun_diff(x, f, f_prime)
      x = x - f / f_prime
      if (abs(f / f_prime) < epsilon) exit
   end do
   if (i == n) print *, '収束しません'

   print "(A, F15.8)", 'ニュートン法の解：', x
   print "(A, I15)", 'ニュートン法の反復回数：', i

end program main

subroutine fun_diff(x, f, f_prime)
   implicit none
   real(8) :: x, f, f_prime

   f = x - cos(x)
   f_prime = 1 + sin(x)
end subroutine fun_diff
