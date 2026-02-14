program do_practice
   implicit none
   integer :: i, n

   print *, 'Enter: n = '
   read *, n
   do i = 1, n
      print *, i
      print *, i*i
      print *, 1.0/real(i)
      print *, sqrt(real(i))
   end do

end program do_practice
