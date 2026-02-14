program area_equilateral_triangle
   implicit none
   real :: a, area

   print *, 'Enter the side of the equilateral triangle: '
   read *, a

   area = (sqrt(3.0) / 4.0) * a**2

   print *, 'The area of the equilateral triangle is: ', area
end program area_equilateral_triangle
