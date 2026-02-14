program poisson
   implicit none
   real(8) :: pi, h, V_0, V_N
   integer :: i, N, is_prob
   real(8), allocatable :: x(:), rho(:), V(:), G(:), Hvec(:)
   ! 1次元ポワソン方程式を連立一次方程式として解く
   print *, '空間の分割数 N = '
   read *, N
   h = 1.0d0 / dble(N)
   allocate(x(0:N), rho(0:N), V(0:N), G(0:N), Hvec(0:N))
   do i = 0, N
      x(i) = i * h
   end do

   pi = 4.0d0 * atan(1.0d0)
   ! 問題を選択
   print *, '解く問題を(1)(2)から選ぶ'
   read *, is_prob
   if (is_prob == 1) then
      do i = 0, N
         rho(i) = 0.0d0
      end do
      V_0 = 0.0d0
      V_N = 1.0d0
   else if (is_prob == 2) then
      do i = 0, N
         rho(i) = sin(2 * pi * x(i) / (dble(N) * h))
      end do
      V_0 = 0.0d0
      V_N = 0.0d0
   end if

   V(0) = V_0
   V(N) = V_N

   G(0) = 0.0d0
   Hvec(0) = V_0

   do i = 1, N-1
      G(i) = - 1 / (-2 + G(i-1))
      Hvec(i) = (rho(i) * h * h - Hvec(i-1)) / (-2.0d0 + G(i-1))
   end do

   do i = N-1, 1, -1
      V(i) = G(i) * V(i+1) + Hvec(i)
   end do

   ! データをCSV形式で出力
   open(20, file='poisson_result.dat', status='replace')
   do i = 0, N
      write(20, '(F12.6, A, F12.6)') x(i), ',', V(i)
   end do
   close(20)
   print *, 'データを poisson_result.dat に出力しました'

end program poisson
