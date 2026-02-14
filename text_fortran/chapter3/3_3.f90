program main
   implicit none
   real(8), allocatable :: R(:)
   integer :: n
   real(8) :: R_s, R_p

   print *, 'Enter the number of resistances: '
   read *, n
   allocate(R(n))
   print *, 'Enter the resistances: '
   read *, R(1:n)
   call cal_resist(R, n, R_s, R_p)

   print *, 'The series resistance is: ', R_s
   print *, 'The parallel resistance is: ', R_p

end program main

subroutine cal_resist(R, n, R_s, R_p)
   implicit none
   integer :: n
   real(8) :: R(n)
   real(8) :: R_s, R_p
   integer :: i

   R_s = R(1)
   R_p = R(1)

   do i = 2, n
      R_s = R_s + R(i)
      R_p = 1.0d0 / (1.0d0 / R_p + 1.0d0 / R(i))
   end do

end subroutine cal_resist
