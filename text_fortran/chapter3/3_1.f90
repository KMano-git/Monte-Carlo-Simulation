subroutine quandratic_equation(a,b,c, x1,x2,discriminant)
   implicit none
   real(8) :: a, b, c, x1, x2, discriminant

   discriminant = b*b - 4*a*c

   if (discriminant > 0) then
      x1 = (-b + sqrt(discriminant)) / (2*a)
      x2 = (-b - sqrt(discriminant)) / (2*a)
   elseif (discriminant == 0) then
      x1 = -b / (2*a)
      x2 = 0.0
   else
      x1 = -b / (2*a) ! 実部
      x2 = sqrt(-discriminant) / (2*a) ! 虚部
   end if

end subroutine quandratic_equation

program main
   implicit none
   real(8) :: a, b, c, x1, x2, discriminant

   print *, 'Enter coefficients a, b, and c: '
   read *, a, b, c

   call quandratic_equation(a, b, c, x1, x2, discriminant)

   if (discriminant > 0) print *, '2実数解 x1, x2 = ', x1, x2
   if (discriminant == 0) print *, '重解 x = ', x1
   if (discriminant < 0) print *, '複素数解 x1 +- x2i = ', dcmplx(x1, x2)
end program
