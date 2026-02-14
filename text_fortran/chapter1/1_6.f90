program impedance
   implicit none
   real(8) :: L, R_L, C, R, f, E, I_current_series, I_current_parallel
   real(8) :: theta_series, theta_parallel, P_series, P_parallel, pi, omega
   complex(8) :: Z_R, Z_L, Z_C, Z_series, Z_parallel
   namelist /option/ L, C, R, f, E, R_L

   ! ファイルから変数を読み込む
   open(15, file='variable_1_6')
   read(15, nml=option)


   pi = 4.0d0 * datan(1.0d0)
   omega = 2.0d0 * pi * f  ! 角周波数

   ! 各素子のインピーダンス
   Z_R = dcmplx(R, 0.0d0)                    ! 抵抗: R
   Z_L = dcmplx(0.0d0, omega * L)            ! インダクタ: jωL
   Z_C = dcmplx(0.0d0, -1.0d0 / (omega * C)) ! キャパシタ: -j/(ωC)

   ! 直列に繋いだときの合成インピーダンス
   Z_series = Z_R + Z_L + Z_C

   ! 並列に繋いだときの合成インピーダンス
   Z_parallel = 1.0d0 / (1.0d0/Z_R + 1.0d0/Z_L + 1.0d0/Z_C)

   ! 電流
   I_current_series = E / abs(Z_series)
   I_current_parallel = E / abs(Z_parallel)

   ! 力率と有効電力　（数式の喧騒はしてないので悪しからず）
   theta_series = datan2(dimag(Z_series), dreal(Z_series))
   P_series = E * I_current_series * dcos(theta_series)
   theta_parallel = datan2(dimag(Z_parallel), dreal(Z_parallel))
   P_parallel = E * I_current_parallel * dcos(theta_parallel)

   ! 結果出力
   print *, '--- 電力計算（並列回路）---'
   print '(A,E12.5,A)', ' 電流 I = ', I_current_parallel, ' A'
   print '(A,F8.3,A)', ' 力率 cos(θ) = ', dcos(theta_parallel), ''
   print '(A,E12.5,A)', ' 有効電力 P = ', P_parallel, ' W'
   print *, '--- 電力計算（直列回路）---'
   print '(A,E12.5,A)', ' 電流 I = ', I_current_series, ' A'
   print '(A,F8.3,A)', ' 力率 cos(θ) = ', dcos(theta_series), ''
   print '(A,E12.5,A)', ' 有効電力 P = ', P_series, ' W'
   close(15)
end program impedance
