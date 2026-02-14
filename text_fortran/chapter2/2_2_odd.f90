program stat_calc_odd
   implicit none
   integer :: n, i
   real :: mean, sigma
   integer, allocatable :: a(:)

   print *, 'Enter: n = '
   read *, n
   allocate(a(n))

   do i = 1, n
      a(i) = 2*i-1
   end do

   do i = 1, n
      mean = mean + real(a(i))
   end do
   mean = mean / real(n)

   do i = 1, n
      sigma = sigma + ((real(a(i)) - mean)*(real(a(i)) - mean))
   end do
   sigma = sqrt(sigma / real(n))

   print *, 'mean = ', mean
   print *, 'sigma = ', sigma
end program stat_calc_odd
