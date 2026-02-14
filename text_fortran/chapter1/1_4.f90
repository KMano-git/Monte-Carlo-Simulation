program biot_savart
   implicit none
   real(8) :: mu0, i, r, alpha, beta, b, pi

   print *, 'Enter current i, distance r, angles (alpha, beta) in degrees: '
   read *, i, r, alpha, beta

   pi = 4.0d0*datan(1.0d0)
   mu0 = 4.0d-7 * pi

   alpha = alpha * pi / 180.0d0
   beta = beta * pi / 180.0d0

   !ビオ・サバールの法則
   b = (mu0 * i) / (4.0d0 * pi * r) * (dcos(alpha) + dcos(beta))

   print *, 'Magnetic field: ', b

end program biot_savart
