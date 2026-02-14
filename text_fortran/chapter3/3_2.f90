program main
   implicit none
   real(8) :: a,b,c,area,heron

   print *, 'Enter side lengths a, b, and c: '
   read *, a, b, c

   area = heron(a,b,c)
   print *, 'Area = ', area
end program main

function heron(a,b,c)
   implicit none
   real(8) :: heron, a,b,c,s

   s = (a + b + c) / 2
   heron = sqrt(s * (s - a) * (s - b) * (s - c))
end function heron
