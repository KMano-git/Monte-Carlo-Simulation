program energy
   implicit none
   real(8) :: c, me, mi, e, Ee, Ei, v_9, v_1, gamma_9, gamma_1, K_e_9, K_e_1, K_i_9, K_i_1

   c = 2.99792458d8
   me = 9.10938356d-31
   mi = 1.672621898d-27
   e = 1.6021766208d-19

   ! 電子、陽子の静止エネルギーを計算
   Ee = me * c * c / e
   Ei = mi * c * c / e

   ! 速度vで運動しているときの運動エネルギーを計算
   v_9 = 9d-1 * c
   v_1 = 1d-1 * c

   gamma_9 = 1.0d0 / dsqrt(1.0d0 - v_9 * v_9 / (c * c))
   gamma_1 = 1.0d0 / dsqrt(1.0d0 - v_1 * v_1 / (c * c))
   K_e_9 = gamma_9 * me * c * c / e
   K_e_1 = gamma_1 * me * c * c / e
   K_i_9 = gamma_9 * mi * c * c / e
   K_i_1 = gamma_1 * mi * c * c / e

   print "(a,e12.5,a)", 'Electron energy: ', Ee, ' J'
   print "(a,e12.5,a)", 'Ion energy: ', Ei, ' J'
   print "(a,e12.5,a)", 'Electron kinetic energy (0.9c): ', K_e_9, ' J'
   print "(a,e12.5,a)", 'Electron kinetic energy (0.1c): ', K_e_1, ' J'
   print "(a,e12.5,a)", 'Ion kinetic energy (0.9c): ', K_i_9, ' J'
   print "(a,e12.5,a)", 'Ion kinetic energy (0.1c): ', K_i_1, ' J'

end program energy
