program main
   implicit none
   real(8) :: a(3,3), b(3), x(3)
   integer :: i, j

   ! 行列Aの入力
   print *, 'Enter the 3x3 matrix A (row by row):'
   do i = 1, 3
      print '(A, I1, A)', 'Row ', i, ': '
      read *, (a(i,j), j=1,3)
   end do

   ! ベクトルbの入力
   print *, 'Enter the vector b (3 elements):'
   read *, b

   ! 連立方程式を解く
   call solve_linear(a, b, x)

   ! 結果の出力
   print *, ''
   print *, 'Solution x:'
   do i = 1, 3
      print '(A, I1, A, F12.6)', 'x(', i, ') = ', x(i)
   end do

end program main

! クラメルの公式を使って連立一次方程式 Ax = b を解く
subroutine solve_linear(a, b, x)
   implicit none
   real(8), intent(in) :: a(3,3), b(3)
   real(8), intent(out) :: x(3)
   real(8) :: det_a, det_ai
   real(8) :: u(3), v(3), w(3)
   real(8) :: cal_determinant
   integer :: i

   ! 行列Aの列ベクトルを取り出す
   u = a(:,1)
   v = a(:,2)
   w = a(:,3)

   ! det(A)を計算
   det_a = cal_determinant(u, v, w)

   if (abs(det_a) < 1.0d-10) then
      print *, 'Error: Matrix is singular (det = 0)'
      x = 0.0d0
      return
   end if

   ! x(1): 第1列をbで置き換え
   det_ai = cal_determinant(b, v, w)
   x(1) = det_ai / det_a

   ! x(2): 第2列をbで置き換え
   det_ai = cal_determinant(u, b, w)
   x(2) = det_ai / det_a

   ! x(3): 第3列をbで置き換え
   det_ai = cal_determinant(u, v, b)
   x(3) = det_ai / det_a

end subroutine solve_linear

! 3つの列ベクトルから行列式を計算
function cal_determinant(u, v, w) result(determinant)
   implicit none
   real(8), intent(in) :: u(3), v(3), w(3)
   real(8) :: determinant

   determinant = u(1) * (v(2) * w(3) - v(3) * w(2)) &
      - u(2) * (v(1) * w(3) - v(3) * w(1)) &
      + u(3) * (v(1) * w(2) - v(2) * w(1))
end function cal_determinant
