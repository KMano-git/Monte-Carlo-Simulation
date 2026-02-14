program zenkasiki
   implicit none
   integer :: i, n
   integer, allocatable :: a(:)
   real(8), allocatable :: b(:)
   real(8) :: mean_a, mean_b, sigma_a, sigma_b

   print *, 'Enter: n = '
   read *, n
   allocate(a(n))
   allocate(b(n))

   a(1)= -1
   b(1)= 3.0d0
   ! a_(n+1)=3a_n+2を解く　理論解：a_n=-1
   do i = 1, n-1
      a(i+1) = 3 * a(i) + 2
   end do
   ! b_(n+1)=2/3 b_n-4を解く　理論解：b_n=(2/3)^(n-1) * 15 -12
   do i = 1, n-1
      b(i+1) = (2.0d0/3.0d0) * b(i) - 4.0d0
   end do
   ! 平均、標準偏差を求める
   do i = 1, n
      mean_a = mean_a + a(i)/dble(n)
      mean_b = mean_b + b(i)/dble(n)
   end do
   do i = 1, n
      sigma_a = sigma_a + ((a(i) - mean_a)*(a(i) - mean_a))
      sigma_b = sigma_b + ((b(i) - mean_b)*(b(i) - mean_b))
   end do
   sigma_a = sqrt(sigma_a / real(n))
   sigma_b = sqrt(sigma_b / real(n))
   ! 理論解との誤差を出力
   print *, '理論解との誤差：'
   print *, 'a_n = ', a(n)
   print *, 'b_n = ', b(n)
   print *, 'mean_a = ', mean_a
   print *, 'mean_b = ', mean_b
   print *, 'sigma_a = ', sigma_a
   print *, 'sigma_b = ', sigma_b
end program zenkasiki
