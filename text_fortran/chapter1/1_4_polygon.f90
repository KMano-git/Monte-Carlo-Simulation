program biot_savart_polygon
   implicit none
   real(8) :: mu0, i_current, a, b_total, b_one_side, pi, r, alpha, beta
   integer :: n, j

   print *, '=== 正n角形導線の中心での磁場計算 ==='
   print *, '（1-4.f90 のビオ・サバールの法則を利用）'
   print *, ''
   print *, 'Enter current I (A), side length a (m), and n (polygon sides): '
   read *, i_current, a, n

   pi = 4.0d0 * datan(1.0d0)
   mu0 = 4.0d0 * pi * 1.0d-7

   ! 正n角形の中心から各辺までの距離
   ! r = a / (2 * tan(pi/n))
   r = a / (2.0d0 * dtan(pi / dble(n)))

   ! 各辺の端点が中心から見える角度（対称なので alpha = beta）
   alpha = pi / dble(n)
   beta = alpha

   ! 1辺からの磁場（1-4.f90 と同じ計算式）
   b_one_side = (mu0 * i_current) / (4.0d0 * pi * r) * (dcos(alpha) + dcos(beta))

   ! n辺すべてからの寄与を足し合わせる
   b_total = 0.0d0
   do j = 1, n
      b_total = b_total + b_one_side
   end do

   print *, ''
   print *, '--- 計算条件 ---'
   print *, '多角形の辺数 n = ', n
   print *, '電流 I = ', i_current, ' A'
   print *, '1辺の長さ a = ', a, ' m'
   print *, '中心までの距離 r = ', r, ' m'
   print *, 'alpha = beta = ', alpha * 180.0d0 / pi, ' 度'
   print *, ''
   print *, '--- 結果 ---'
   print *, '1辺からの磁場 B_1 = ', b_one_side, ' T'
   print *, '中心での全磁場 B = ', b_total, ' T'

end program biot_savart_polygon
