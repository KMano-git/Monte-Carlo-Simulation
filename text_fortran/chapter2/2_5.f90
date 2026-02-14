program taylor
   implicit none
   real(8) :: x, e, cos_val, sin_val, term
   integer :: i, n

   print *, 'Enter: x, n'
   read *, x, n

   e = 1.0d0
   cos_val = 1.0d0
   sin_val = x
   term = 1.0d0
   !e^x,cosx,sinxをTaylor展開して計算する
   do i = 1, n
      term = term * x / dble(i)
      e = e + term
   end do

   term = 1.0d0
   do i = 1, n
      term = term * (-1.0d0) * x * x / (dble(2*i) * dble(2*i-1))
      cos_val = cos_val + term
   end do

   term = x
   do i = 1, n
      term = term * (-1.0d0) * x * x / (dble(2*i+1) * dble(2*i))
      sin_val = sin_val + term
   end do
   !理論値と計算値を出力する
   print *, 'e^x = ', e, ' (理論値: ', exp(x), ')'
   print *, 'cos(x) = ', cos_val, ' (理論値: ', cos(x), ')'
   print *, 'sin(x) = ', sin_val, ' (理論値: ', sin(x), ')'
end program taylor
