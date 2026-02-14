program successive_assignment
   implicit none
   integer :: i, n
   real(8) :: x

   n = 100000
   ! コサインを対象に実施、
   x = 0.0d0
   do i = 1, n
      x = cos(x)
      if (abs(cos(x)-x) < 1.0d-7) exit
   end do
   if (i == n) print *, '収束しません'

   print *, 'cos(x) = x の解：', x
end program successive_assignment
