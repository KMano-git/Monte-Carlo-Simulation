program main
   implicit none
   integer :: n
   real(8), allocatable :: a(:,:),b(:,:),c(:,:),d(:,:)

   print *, 'Enter the size of the matrix: '
   read *, n
   allocate(a(n,n), b(n,n), c(n,n), d(n,n))
   print *, 'Enter the matrix A: '
   read *, a(1:n,1:n)
   print *, 'Enter the matrix B: '
   read *, b(1:n,1:n)
   call cal_matrix(a, b, c, d, n)
   print *, 'The matrix C is: '
   print *, c(1:n,1:n)
   print *, 'The matrix D is: '
   print *, d(1:n,1:n)
   deallocate(a, b, c, d)
end program main

subroutine cal_matrix(a, b, c, d, n)
   implicit none
   integer, intent(in) :: n
   real(8) :: a(n,n), b(n,n), c(n,n), d(n,n)
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
end subroutine cal_matrix
