program symplectic
   implicit none
   real(8) :: pi, omega, dt, u, v, u_old, N
   integer :: i

   ! シンプレクティック法を用いると、次の形で式変形できる
   ! u = u + omega * v * dt
   ! v = v - omega * u * dt

   print *, '1周期の分割数 N = '
   read *, N
   pi = 4.0d0 * atan(1.0d0)
   omega = 1.0d0
   dt = 2 * pi / (omega * N)
   u = 1.0d0
   v = 0.0d0

   do i = 1, int(N)
      u_old = u
      u = u + omega * v * dt
      v = v - omega * u * dt
   end do

   print *, 'u0 = 1.0, uN =', u
   print *, 'v0 = 0.0, vN =', v

   do i = 1, 9 * int(N)
      u_old = u
      u = u + omega * v * dt
      v = v - omega * u * dt
   end do

   print *, 'u0 = 1.0, u10N =', u
   print *, 'v0 = 0.0, v10N =', v

end program symplectic
