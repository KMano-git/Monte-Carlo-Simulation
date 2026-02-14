program vector_calc
   implicit none
   real(8) :: E(3), B(3), v(3), v_dot_E, v_cross_E(3)
   integer :: i

   E = [1.0d0, 2.0d0, 3.0d0]
   B = [4.0d0, 5.0d0, 6.0d0]
   v = [7.0d0, 8.0d0, 9.0d0]

   !内積
   v_dot_E = 0.0d0
   do i = 1, 3
      v_dot_E = v_dot_E + v(i) * E(i)
   end do
   !外積
   v_cross_E = [0.0d0, 0.0d0, 0.0d0]
   v_cross_E(1) = v(2)*E(3) - v(3)*E(2)
   v_cross_E(2) = v(3)*E(1) - v(1)*E(3)
   v_cross_E(3) = v(1)*E(2) - v(2)*E(1)

   print *, '内積: ', v_dot_E
   print *, '外積: ', v_cross_E

end program vector_calc
