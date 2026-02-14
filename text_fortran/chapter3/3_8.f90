program main
   implicit none
   real(8) :: u(3), v(3), w(3), determinant, cal_determinant

   print *, 'Enter the vector u: '
   read *, u
   print *, 'Enter the vector v: '
   read *, v
   print *, 'Enter the vector w: '
   read *, w

   determinant = cal_determinant(u, v, w)
   print *, 'The determinant of the matrix is: ', determinant
end program main

function cal_determinant(u, v, w) result(determinant)
   implicit none
   real(8) :: u(3), v(3), w(3)
   real(8) :: determinant

   determinant = u(1) * (v(2) * w(3) - v(3) * w(2)) - u(2) * (v(1) * w(3) - v(3) * w(1)) + u(3) * (v(1) * w(2) - v(2) * w(1))
end function cal_determinant

