program quandratic_equation
   implicit none
   real :: a, b, c
   real :: discriminant, root1, root2
   complex :: root3, root4

   print *, 'Enter coefficients a, b, and c: '
   read *, a, b, c

   discriminant = b**2 - 4*a*c

   if (discriminant > 0) then
      root1 = (-b + sqrt(discriminant)) / (2*a)
      root2 = (-b - sqrt(discriminant)) / (2*a)
      print *, 'Root 1 = ', root1
      print *, 'Root 2 = ', root2
   elseif (discriminant == 0) then
      root1 = -b / (2*a)
      print *, 'Root = ', root1
   else
      root3 = cmplx(-b / (2*a), sqrt(-discriminant) / (2*a))
      root4 = cmplx(-b / (2*a), -sqrt(-discriminant) / (2*a))
      print *, 'Root 1 = ', root3
      print *, 'Root 2 = ', root4
   end if
   ! かなり先走って、2－6の内容までやってたわ。if文の文法は参考にすること
end program quandratic_equation
