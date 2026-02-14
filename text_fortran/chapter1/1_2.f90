program heron_formula
   implicit none
   real :: a, b, c
   real :: s, area

   !ヘロンの公式(1-2_1の正三角形のプログラムの結果と比較）
   print *, 'Enter the sides of the triangle: '
   read *, a, b, c

   s = (a + b + c) / 2
   area = sqrt(s * (s - a) * (s - b) * (s - c))

   print *, 'The area of the triangle is: ', area
end program heron_formula
