program area_volume_polygon
   implicit none
   real(8) :: a, pi, area_circle, volume_sphere, area_polygon
   integer :: n

   print *, 'Enter a real number and integer: '
   read *, a, n

   pi = 4.0d0 * datan(1.0d0) !円周率
   area_circle = pi * a * a !半径aの円の面積
   volume_sphere = (4.0d0/3.0d0) * pi * a * a * a !半径aの円の球の体積
   area_polygon = (n * a * a) / (4.0d0 * dtan(pi / n)) !1辺aの正n角形の面積

   print *, 'Area of circle: ', area_circle
   print *, 'Volume of sphere: ', volume_sphere
   print *, 'Area of polygon: ', area_polygon
end program area_volume_polygon
