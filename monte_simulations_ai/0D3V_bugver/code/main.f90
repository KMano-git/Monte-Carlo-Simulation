!===============================================================================
! Program: monte_carlo_0d
! 0D-3V Monte Carlo シミュレーション
! 衝突による速度分布の緩和過程を追跡
!===============================================================================
program monte_carlo_0d
   use constants
   use particle_data
   use random_utils
   use collisions
   use statistics
   use cdf_reader
   implicit none

   ! パラメータ型
   type(sim_params) :: sim
   type(plasma_params) :: plasma
   type(init_params) :: init_p
   type(diag_params) :: diag

   ! 粒子配列
   type(particle_0d), allocatable :: particles(:)

   ! 統計量
   real(dp) :: E_mean, E_std, T_eff
   real(dp) :: time
   real(dp) :: delta_E_total, delta_E_avg
   real(dp) :: delta_E_cx_total, delta_E_el_total
   real(dp) :: delta_E_cx_avg, delta_E_el_avg

   ! ヒストグラム
   integer, allocatable :: hist(:)
   real(dp), allocatable :: E_centers(:)

   ! 断面積（簡略化: 定数と仮定）
   real(dp) :: sigma_cx, sigma_el, sigma_total
   real(dp) :: nu_cx, nu_el, nu_total, P_coll
   real(dp) :: v_mag, E_particle

   ! ループ変数
   integer :: i, istep
   integer :: unit_relax, unit_ntscrg
   real(dp) :: r, r_type
   real(dp) :: v_th_bg, v_rel_max, nu_max    ! Null Collision用
   real(dp) :: vx_bg, vy_bg, vz_bg           ! 背景速度
   real(dp) :: v_rel, P_accept               ! 相対速度と採用確率
   real(dp) :: E_rel                         ! 衝突の相対エネルギー
   real(dp) :: delta_E
   integer :: ierr
   real(dp), parameter :: SIGMA_MAX = 2.0d-18 ! 最大想定断面積 [m^2]

   !-----------------------------------------------------------------------------
   ! 初期化
   !-----------------------------------------------------------------------------
   write(*,'(A)') '============================================'
   write(*,'(A)') '  0D-3V Monte Carlo Simulation'
   write(*,'(A)') '============================================'
   write(*,*)

   ! パラメータ読み込み
   call read_input_file(sim, plasma, init_p, diag)

   ! CDFファイルの読み込み
   if (.not. sim%use_isotropic) then
      call load_elastic_cdf(trim(sim%cdf_file), ierr)
      if (ierr /= 0) then
         write(*,*) 'Error: Failed to load CDF file: ', trim(sim%cdf_file)
         stop
      end if
   end if

   ! 乱数シード初期化
   call init_random_seed(sim%seed)

   ! 粒子配列の割り当て
   allocate(particles(sim%n_particles))

   ! 粒子初期化
   call initialize_particles(particles, sim%n_particles, init_p)

   ! ヒストグラム配列の割り当て
   allocate(hist(diag%n_hist_bins))
   allocate(E_centers(diag%n_hist_bins))

   ! 出力ファイルのオープン
   open(newunit=unit_relax, file='relaxation.csv', status='replace')
   write(unit_relax, '(A)') '# time[s], E_mean[eV], T_eff[eV]'

   open(newunit=unit_ntscrg, file='ntscrg.csv', status='replace')
   write(unit_ntscrg, '(A)') '# time[s], delta_E_avg[eV], delta_E_cx_avg[eV], delta_E_el_avg[eV]'

   ! 初期統計
   call compute_energy_stats(particles, sim%n_particles, E_mean, E_std, T_eff)
   write(*,'(A,ES12.4,A,F10.4,A,F10.4,A)') &
      '  t = ', 0.0d0, ' s:  E_mean = ', E_mean, ' eV,  T_eff = ', T_eff, ' eV'

   call write_relaxation_data(unit_relax, 0.0d0, E_mean, T_eff)
   call write_relaxation_data(unit_relax, 0.0d0, E_mean, T_eff)
   write(unit_ntscrg, '(ES15.6E3,",",ES15.6E3,",",ES15.6E3,",",ES15.6E3)') &
      0.0d0, 0.0d0, 0.0d0, 0.0d0

   !-----------------------------------------------------------------------------
   ! 時間ループ
   !-----------------------------------------------------------------------------
   do istep = 1, sim%n_steps
      time = istep * sim%dt
      delta_E_total = 0.0d0
      delta_E_cx_total = 0.0d0
      delta_E_el_total = 0.0d0

      ! 各粒子について衝突判定
      do i = 1, sim%n_particles
         if (.not. particles(i)%alive) cycle

         ! 粒子速度とエネルギー
         v_mag = sqrt(particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2)
         if (v_mag < 1.0d-10) cycle

         E_particle = 0.5d0 * M_D * v_mag**2 * J_TO_EV

         !-----------------------------------------------------------------------
         ! Null Collision Method (Rejection Sampling)
         !-----------------------------------------------------------------------

         ! 1. 最大衝突頻度の見積もり
         ! 相対速度の最大値を推定 (v_particle + 5*v_th_bg)
         ! 安全率として少し大きめに取る
         v_th_bg = sqrt(plasma%T_bg * EV_TO_J / M_D)
         v_rel_max = v_mag + 5.0d0 * v_th_bg

         ! 最大衝突頻度
         ! 定数SIGMA_MAXを利用して上限を見積もる
         nu_max = plasma%n_bg * SIGMA_MAX * v_rel_max

         ! 仮の衝突確率
         P_coll = 1.0d0 - exp(-nu_max * sim%dt)

         ! 2. 仮の衝突判定
         call random_number(r)
         if (r > P_coll) cycle  ! 衝突なし

         ! 3. ターゲット速度のサンプリングとRejection判定
         call sample_maxwell_velocity(plasma%T_bg, vx_bg, vy_bg, vz_bg)

         ! 相対速度
         v_rel = sqrt((particles(i)%vx - vx_bg)**2 + &
            (particles(i)%vy - vy_bg)**2 + &
            (particles(i)%vz - vz_bg)**2)

         ! 重心系での相対エネルギー [eV]
         ! 換算質量 = M_D / 2.0 (同種粒子)
         E_rel = 0.5d0 * (0.5d0 * M_D) * v_rel**2 * J_TO_EV

         ! 実際の断面積を相対エネルギーで評価
         sigma_cx = get_sigma_cx(E_rel)
         sigma_el = get_sigma_el(E_rel)
         sigma_total = sigma_cx + sigma_el

         ! Rejection判定:
         P_accept = (sigma_total * v_rel) / (SIGMA_MAX * v_rel_max)

         call random_number(r)
         if (r > P_accept) cycle ! 虚の衝突 (Null Collision)

         ! 4. 実際の衝突処理（タイプ選択）
         call random_number(r_type)

         if (r_type < sigma_cx / sigma_total .and. sim%enable_cx) then
            ! 荷電交換
            ! デバッグ: 最初の数回の衝突のみ詳細出力
            if (istep < 5 .and. i == 1) then
               call collision_cx(particles(i), vx_bg, vy_bg, vz_bg, delta_E)
               print *, "DEBUG: CX Collision (step=", istep, " particle=", i, ")"
               print *, "  Before E=", 0.5d0*M_D*(particles(i)%vx**2+particles(i)%vy**2+particles(i)%vz**2)*J_TO_EV - delta_E
               print *, "  After  E=", 0.5d0*M_D*(particles(i)%vx**2+particles(i)%vy**2+particles(i)%vz**2)*J_TO_EV
               print *, "  Delta  E=", delta_E
            else
               call collision_cx(particles(i), vx_bg, vy_bg, vz_bg, delta_E)
            end if
            delta_E_cx_total = delta_E_cx_total + delta_E
            delta_E_total = delta_E_total + delta_E

         else if (sim%enable_elastic) then
            ! 弾性散乱
            call collision_elastic(particles(i), vx_bg, vy_bg, vz_bg, sim%use_isotropic, &
               E_rel, delta_E)
            delta_E_el_total = delta_E_el_total + delta_E
            delta_E_total = delta_E_total + delta_E
         end if

      end do

      ! 平均エネルギー変化
      delta_E_avg = delta_E_total / sim%n_particles
      delta_E_cx_avg = delta_E_cx_total / sim%n_particles
      delta_E_el_avg = delta_E_el_total / sim%n_particles

      ! 診断出力
      if (mod(istep, diag%output_interval) == 0) then
         call compute_energy_stats(particles, sim%n_particles, E_mean, E_std, T_eff)
         call write_relaxation_data(unit_relax, time, E_mean, T_eff)

         write(*,'(A,ES12.4,A,F10.4,A,F10.4,A)') &
            '  t = ', time, ' s:  E_mean = ', E_mean, ' eV,  T_eff = ', T_eff, ' eV'
      end if

      ! ntscrg.csv への書き込み
      write(unit_ntscrg, '(ES15.6E3,",",ES15.6E3,",",ES15.6E3,",",ES15.6E3)') &
         time, delta_E_avg, delta_E_cx_avg, delta_E_el_avg

   end do

   !-----------------------------------------------------------------------------
   ! 最終出力
   !-----------------------------------------------------------------------------
   write(*,*)
   write(*,'(A)') '  Simulation completed!'

   ! 最終エネルギー分布
   call compute_energy_histogram(particles, sim%n_particles, diag%n_hist_bins, &
      diag%E_hist_min, diag%E_hist_max, hist, E_centers)
   call write_histogram_data('energy_dist.csv', hist, E_centers, diag%n_hist_bins)

   ! 最終統計
   call compute_energy_stats(particles, sim%n_particles, E_mean, E_std, T_eff)
   write(*,'(A,F10.4,A)') '  Final E_mean = ', E_mean, ' eV'
   write(*,'(A,F10.4,A)') '  Final T_eff  = ', T_eff, ' eV'
   write(*,'(A,F10.4,A)') '  Target T_bg  = ', plasma%T_bg, ' eV'

   ! クリーンアップ
   close(unit_relax)
   close(unit_ntscrg)
   deallocate(particles, hist, E_centers)

   write(*,*)
   write(*,'(A)') '============================================'

contains

   !-----------------------------------------------------------------------------
   ! 入力ファイル読み込み
   !-----------------------------------------------------------------------------
   subroutine read_input_file(sim, plasma, init_p, diag)
      type(sim_params), intent(out) :: sim
      type(plasma_params), intent(out) :: plasma
      type(init_params), intent(out) :: init_p
      type(diag_params), intent(out) :: diag

      integer :: unit
      logical :: file_exists

      namelist /simulation/ n_particles, n_steps, dt, seed, &
         enable_cx, enable_elastic, use_isotropic, cdf_file
      namelist /plasma_nml/ n_bg, T_bg
      namelist /particle_init/ E_init, T_init, init_mode
      namelist /diagnostics/ output_interval, n_hist_bins, E_hist_min, E_hist_max

      ! デフォルト値
      integer :: n_particles = 10000
      integer :: n_steps = 10000
      real(dp) :: dt = 1.0d-8
      integer :: seed = 12345
      logical :: enable_cx = .true.
      logical :: enable_elastic = .true.
      logical :: use_isotropic = .true.
      character(len=256) :: cdf_file = 'dd_00_elastic.cdf'

      real(dp) :: n_bg = 1.0d21
      real(dp) :: T_bg = 10.0d0

      real(dp) :: E_init = 100.0d0
      real(dp) :: T_init = 2.0d0
      integer :: init_mode = 0

      integer :: output_interval = 100
      integer :: n_hist_bins = 100
      real(dp) :: E_hist_min = 0.0d0
      real(dp) :: E_hist_max = 150.0d0

      ! ファイル確認
      inquire(file='input.nml', exist=file_exists)
      if (.not. file_exists) then
         write(*,'(A)') '  Warning: input.nml not found, using default values'
      else
         open(newunit=unit, file='input.nml', status='old')
         read(unit, nml=simulation, end=100)
100      rewind(unit)
         read(unit, nml=plasma_nml, end=200)
200      rewind(unit)
         read(unit, nml=particle_init, end=300)
300      rewind(unit)
         read(unit, nml=diagnostics, end=400)
400      close(unit)
      end if

      ! パラメータのセット
      sim%n_particles = n_particles
      sim%n_steps = n_steps
      sim%dt = dt
      sim%seed = seed
      sim%enable_cx = enable_cx
      sim%enable_elastic = enable_elastic
      sim%use_isotropic = use_isotropic
      sim%cdf_file = cdf_file

      plasma%n_bg = n_bg
      plasma%T_bg = T_bg

      init_p%E_init = E_init
      init_p%T_init = T_init
      init_p%init_mode = init_mode

      diag%output_interval = output_interval
      diag%n_hist_bins = n_hist_bins
      diag%E_hist_min = E_hist_min
      diag%E_hist_max = E_hist_max

      ! パラメータ表示
      write(*,'(A)') '  Parameters:'
      write(*,'(A,I0)') '    n_particles   = ', sim%n_particles
      write(*,'(A,I0)') '    n_steps       = ', sim%n_steps
      write(*,'(A,ES10.2)') '    dt            = ', sim%dt
      write(*,'(A,L1)') '    enable_cx     = ', sim%enable_cx
      write(*,'(A,L1)') '    enable_elastic= ', sim%enable_elastic
      write(*,'(A,ES10.2)') '    n_bg          = ', plasma%n_bg
      write(*,'(A,F8.2)') '    T_bg          = ', plasma%T_bg
      write(*,*)

   end subroutine read_input_file

   !-----------------------------------------------------------------------------
   ! 粒子初期化
   !-----------------------------------------------------------------------------
   subroutine initialize_particles(particles, n_particles, init_p)
      type(particle_0d), intent(out) :: particles(:)
      integer, intent(in) :: n_particles
      type(init_params), intent(in) :: init_p

      integer :: i

      do i = 1, n_particles
         particles(i)%alive = .true.

         if (init_p%init_mode == 0) then
            ! ビーム分��（+x方向）
            call set_beam_velocity(init_p%E_init, particles(i)%vx, &
               particles(i)%vy, particles(i)%vz)
         else
            ! Maxwell分布
            call sample_maxwell_velocity(init_p%T_init, particles(i)%vx, &
               particles(i)%vy, particles(i)%vz)
         end if
      end do

   end subroutine initialize_particles

   !-----------------------------------------------------------------------------
   ! 荷電交換断面積（Janev/Rabinovitch近似）
   ! D⁰ + D⁺ → D⁺ + D⁰
   !-----------------------------------------------------------------------------
   function get_sigma_cx(E_eV) result(sigma)
      real(dp), intent(in) :: E_eV
      real(dp) :: sigma

      real(dp) :: E_clip

      ! 低エネルギーでの発散を防ぐ
      E_clip = max(E_eV, 0.01d0)

      ! Janev近似式 (D + D+ → D+ + D)
      ! σ_cx ≈ 3.0×10^-19 × (1 - 0.05×log10(E))^2 [m²]
      sigma = 3.0d-19 * (1.0d0 - 0.05d0 * log10(E_clip))**2

      ! 最小値を設定（高エネルギーでの過小評価防止）
      sigma = max(sigma, 1.0d-22)

   end function get_sigma_cx

   function get_sigma_el(E_eV) result(sigma)
      real(dp), intent(in) :: E_eV
      real(dp) :: sigma

      if (cdf_loaded) then
         sigma = get_sigma_elastic(E_eV)
      else
         ! 簡略化: 典型的な値
         sigma = 2.0d-19  ! [m^2]
      end if

   end function get_sigma_el

end program monte_carlo_0d
