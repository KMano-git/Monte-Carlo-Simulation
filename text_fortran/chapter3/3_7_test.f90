program main
   implicit none
   integer :: n
   real(8), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:)
   real(8) :: t1, t2, t3, t4

   n = 2000  ! 大きな行列で比較
   allocate(a(n,n), b(n,n), c(n,n), d(n,n))

   ! ランダムな値で初期化
   call random_number(a)
   call random_number(b)

   ! テスト1: j外側、i内側（Fortran最適）
   call cpu_time(t1)
   call cal_matrix_ji(a, b, c, d, n)
   call cpu_time(t2)
   print '(A, F10.4, A)', 'j外側, i内側 (最適): ', t2 - t1, ' 秒'

   ! テスト2: i外側、j内側（非最適）
   call cpu_time(t3)
   call cal_matrix_ij(a, b, c, d, n)
   call cpu_time(t4)
   print '(A, F10.4, A)', 'i外側, j内側 (非最適): ', t4 - t3, ' 秒'

   print '(A, F10.2, A)', '速度比: ', (t4 - t3) / (t2 - t1), ' 倍'

   deallocate(a, b, c, d)
end program main

! j外側、i内側（Fortran最適）
subroutine cal_matrix_ji(a, b, c, d, n)
   implicit none
   integer, intent(in) :: n
   real(8), intent(in) :: a(n,n), b(n,n)
   real(8), intent(out) :: c(n,n), d(n,n)
   integer :: i, j, k

   do j = 1, n
      do i = 1, n
         c(i,j) = 0.0d0
         d(i,j) = 0.0d0
      end do
   end do
   do j = 1, n
      do i = 1, n
         c(i,j) = a(i,j) + b(i,j)
         do k = 1, n
            d(i,j) = d(i,j) + a(i,k) * b(k,j)
         end do
      end do
   end do
end subroutine cal_matrix_ji

! i外側、j内側（非最適）
subroutine cal_matrix_ij(a, b, c, d, n)
   implicit none
   integer, intent(in) :: n
   real(8), intent(in) :: a(n,n), b(n,n)
   real(8), intent(out) :: c(n,n), d(n,n)
   integer :: i, j, k

   do i = 1, n
      do j = 1, n
         c(i,j) = 0.0d0
         d(i,j) = 0.0d0
      end do
   end do
   do i = 1, n
      do j = 1, n
         c(i,j) = a(i,j) + b(i,j)
         do k = 1, n
            d(i,j) = d(i,j) + a(i,k) * b(k,j)
         end do
      end do
   end do
end subroutine cal_matrix_ij
