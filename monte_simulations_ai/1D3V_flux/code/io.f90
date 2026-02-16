!===============================================================================
! Module: io
! 入出力関連（1D3V フラックス入射モデル）
!===============================================================================
module io
   use constants, only: dp, M_D_kg, J_TO_EV, EV_TO_J
   use particle_data, only: particle_t, sim_params, plasma_params, init_params, diag_params
   use random_utils, only: set_beam_velocity, sample_maxwell_velocity, sample_half_maxwell_velocity
   implicit none

contains

   !-----------------------------------------------------------------------------
   ! 入力パラメータ読み込み
   !-----------------------------------------------------------------------------
   subroutine read_input_file(sim, plasma, init_p, diag)
      type(sim_params), intent(out)  :: sim
      type(plasma_params), intent(out) :: plasma
      type(init_params), intent(out) :: init_p
      type(diag_params), intent(out) :: diag

      integer, parameter :: unit_input = 20
      logical :: file_exist

      namelist /simulation/ n_particles, n_steps, dt, n_grid, x_min, x_max, &
         n_inject_per_step, seed, &
         enable_cx, enable_el, enable_ei, use_isotropic, &
         cdf_file, output_ntscrg, output_hist, output_statx
      namelist /plasma_nml/ n_i, T_i, n_e, T_e
      namelist /particle_init/ E_init, T_init, gamma_in, init_mode
      namelist /diagnostics/ output_interval, n_hist_bins, E_hist_min, E_hist_max, hist_timing

      ! デフォルト値
      integer :: n_particles = 100000
      integer :: n_steps = 5000
      real(dp) :: dt = 1.0d-8
      integer :: n_grid = 100
      real(dp) :: x_min = 0.0d0
      real(dp) :: x_max = 0.1d0
      integer :: n_inject_per_step = 100
      integer :: seed = 12345
      logical :: enable_cx = .true.
      logical :: enable_el = .true.
      logical :: enable_ei = .false.
      logical :: use_isotropic = .false.
      character(len=256) :: cdf_file = 'dd_00_elastic.cdf'
      character(len=256) :: output_ntscrg = 'ntscrg.csv'
      character(len=256) :: output_hist = 'energy_hist.csv'
      character(len=256) :: output_statx = 'statx.csv'

      real(dp) :: n_i = 1.0d21
      real(dp) :: T_i = 10.0d0
      real(dp) :: n_e = 1.0d21
      real(dp) :: T_e = 10.0d0

      real(dp) :: E_init = 3.0d0
      real(dp) :: T_init = 2.0d0
      real(dp) :: gamma_in = 1.0d23
      integer :: init_mode = 0

      integer :: output_interval = 100
      integer :: n_hist_bins = 100
      real(dp) :: E_hist_min = 0.0d0
      real(dp) :: E_hist_max = 150.0d0
      integer :: hist_timing(5) = (/ 100, 500, 1000, 2000, 5000 /)

      ! ファイル確認
      inquire(file = 'input.nml', exist=file_exist)
      if (.not. file_exist) then
         write(*,*) 'Error: input.nml not found'
         stop
      end if

      ! ファイルオープン
      open(unit=unit_input, file='input.nml', status='old')

      ! Namelistの読み込み
      read(unit_input, nml=simulation)
      read(unit_input, nml=plasma_nml)
      read(unit_input, nml=particle_init)
      read(unit_input, nml=diagnostics)

      close(unit_input)

      ! パラメータのセット
      sim%n_particles = n_particles
      sim%n_steps = n_steps
      sim%dt = dt
      sim%n_grid = n_grid
      sim%x_min = x_min
      sim%x_max = x_max
      sim%n_inject_per_step = n_inject_per_step
      sim%seed = seed
      sim%enable_cx = enable_cx
      sim%enable_el = enable_el
      sim%enable_ei = enable_ei
      sim%use_isotropic = use_isotropic
      sim%cdf_file = cdf_file
      sim%output_ntscrg = output_ntscrg
      sim%output_hist = output_hist
      sim%output_statx = output_statx

      plasma%n_i = n_i
      plasma%T_i = T_i
      plasma%n_e = n_e
      plasma%T_e = T_e

      init_p%E_init = E_init
      init_p%T_init = T_init
      init_p%gamma_in = gamma_in
      init_p%init_mode = init_mode

      diag%output_interval = output_interval
      diag%n_hist_bins = n_hist_bins
      diag%E_hist_min = E_hist_min
      diag%E_hist_max = E_hist_max
      diag%hist_timing = hist_timing

   end subroutine read_input_file

   !-----------------------------------------------------------------------------
   ! 粒子配列の初期化（全粒子を dead 状態にする）
   !-----------------------------------------------------------------------------
   subroutine initialize_particles(particles, n_particles)
      type(particle_t), intent(out) :: particles(:)
      integer, intent(in) :: n_particles

      integer :: i

      do i = 1, n_particles
         particles(i)%alive = .false.
         particles(i)%x = 0.0d0
         particles(i)%vx = 0.0d0
         particles(i)%vy = 0.0d0
         particles(i)%vz = 0.0d0
         particles(i)%w = 0.0d0
      end do

   end subroutine initialize_particles

   !-----------------------------------------------------------------------------
   ! 粒子注入: dead スロットに新しい入射粒子を配置
   ! 戻り値: 実際に注入した粒子数
   !-----------------------------------------------------------------------------
   subroutine inject_particles(particles, n_particles, sim, init_p, weight, n_injected)
      type(particle_t), intent(inout) :: particles(:)
      integer, intent(in) :: n_particles
      type(sim_params), intent(in) :: sim
      type(init_params), intent(in) :: init_p
      real(dp), intent(in) :: weight
      integer, intent(out) :: n_injected

      integer :: i, count

      count = 0

      do i = 1, n_particles
         if (count >= sim%n_inject_per_step) exit

         if (.not. particles(i)%alive) then
            particles(i)%alive = .true.
            particles(i)%x = sim%x_min
            particles(i)%w = weight

            if (init_p%init_mode == 0) then
               ! ビーム状（+x方向に E_init のエネルギー）
               call set_beam_velocity(init_p%E_init, particles(i)%vx, &
                  particles(i)%vy, particles(i)%vz)
            else
               ! 半Maxwell分布（vx > 0）
               call sample_half_maxwell_velocity(init_p%T_init, particles(i)%vx, &
                  particles(i)%vy, particles(i)%vz)
            end if

            count = count + 1
         end if
      end do

      n_injected = count

   end subroutine inject_particles

   !-----------------------------------------------------------------------------
   ! 中性粒子のエネルギーヒストグラム出力
   !-----------------------------------------------------------------------------
   subroutine output_energy_histogram(particles, sim, diag, istep)
      type(particle_t), intent(in) :: particles(:)
      type(sim_params), intent(in)  :: sim
      type(diag_params), intent(in) :: diag
      integer, intent(in) :: istep

      integer :: i
      integer, parameter :: unit_hist = 20
      integer :: bin_index, n_alive
      real(dp) :: E_ev
      real(dp) :: bin_width
      real(dp) :: E_min, E_max
      integer, allocatable :: hist(:)

      ! ヒストグラムの範囲とビン幅
      E_min = diag%E_hist_min
      E_max = diag%E_hist_max
      bin_width = (E_max - E_min) / dble(diag%n_hist_bins)

      ! ヒストグラム配列の確保
      allocate(hist(diag%n_hist_bins))
      hist = 0
      n_alive = 0

      ! エネルギーのカウント
      do i = 1, sim%n_particles
         if (.not. particles(i)%alive) cycle
         n_alive = n_alive + 1

         ! 粒子のエネルギーをeVに変換
         E_ev = 0.5d0 * M_D_kg * (particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2) * J_TO_EV

         ! ビン番号を計算
         if (E_ev >= E_min .and. E_ev < E_max) then
            bin_index = int((E_ev - E_min) / bin_width) + 1
            hist(bin_index) = hist(bin_index) + 1
         end if
      end do

      ! ファイルオープン
      open(unit=unit_hist, file=trim(sim%output_hist), status='unknown', position='append')

      ! ヘッダー書き込み（初回のみ）
      if (istep == 0) then
         write(unit_hist, '(A)', advance='no') 'step,E_min,E_max,bin_width,n_alive'
         do i = 1, diag%n_hist_bins
            write(unit_hist, '(A,I0)', advance='no') ',bin_', i
         end do
         write(unit_hist, *)
      end if

      ! データ書き込み
      write(unit_hist, '(I0,",",F0.4,",",F0.4,",",F0.4,",",I0)', advance='no') &
         istep, E_min, E_max, bin_width, n_alive
      do i = 1, diag%n_hist_bins
         write(unit_hist, '(",",I0)', advance='no') hist(i)
      end do
      write(unit_hist, *)

      close(unit_hist)
      deallocate(hist)

   end subroutine output_energy_histogram

   !-----------------------------------------------------------------------------
   ! 粒子の位置ヒストグラム出力
   !-----------------------------------------------------------------------------
   subroutine output_position_histogram(particles, sim, istep)
      type(particle_t), intent(in) :: particles(:)
      type(sim_params), intent(in)  :: sim
      integer, intent(in) :: istep

      integer :: i
      integer, parameter :: unit_pos = 21
      integer :: bin_index, n_alive
      real(dp) :: dx
      real(dp) :: x_min, x_max
      integer, allocatable :: hist(:)

      ! ヒストグラムの範囲とビン幅
      x_min = sim%x_min
      x_max = sim%x_max
      if (sim%n_grid > 0) then
         dx = (x_max - x_min) / dble(sim%n_grid)
      else
         dx = 1.0d0
      endif

      ! ヒストグラム配列の確保
      allocate(hist(sim%n_grid))
      hist = 0
      n_alive = 0

      ! 位置のカウント
      do i = 1, sim%n_particles
         if (.not. particles(i)%alive) cycle
         n_alive = n_alive + 1

         ! ビン番号を計算
         if (particles(i)%x >= x_min .and. particles(i)%x < x_max) then
            bin_index = int((particles(i)%x - x_min) / dx) + 1
            if (bin_index >= 1 .and. bin_index <= sim%n_grid) then
               hist(bin_index) = hist(bin_index) + 1
            end if
         end if
      end do

      ! ファイルオープン
      open(unit=unit_pos, file=trim(sim%output_statx), status='unknown', position='append')

      ! ヘッダー書き込み（初回のみ）
      if (istep == 0) then
         write(unit_pos, '(A)', advance='no') 'step,x_min,x_max,dx,n_alive'
         do i = 1, sim%n_grid
            write(unit_pos, '(A,I0)', advance='no') ',bin_', i
         end do
         write(unit_pos, *)
      end if

      ! データ書き込み
      write(unit_pos, '(I0,",",F0.4,",",F0.4,",",ES12.4,",",I0)', advance='no') &
         istep, x_min, x_max, dx, n_alive
      do i = 1, sim%n_grid
         write(unit_pos, '(",",I0)', advance='no') hist(i)
      end do
      write(unit_pos, *)

      close(unit_pos)
      deallocate(hist)

   end subroutine output_position_histogram

   !-----------------------------------------------------------------------------
   ! 最終統計出力
   !-----------------------------------------------------------------------------
   subroutine output_statistics(particles, n_particles, plasma, time, weight)
      type(particle_t), intent(in) :: particles(:)
      integer, intent(in) :: n_particles
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in) :: time
      real(dp), intent(in) :: weight

      real(dp) :: E_mean, E_std, T_eff
      real(dp) :: E_sum, E2_sum, E_i
      integer  :: i, n_alive

      E_sum = 0.0d0
      E2_sum = 0.0d0
      n_alive = 0

      do i = 1, n_particles
         if (particles(i)%alive) then
            n_alive = n_alive + 1
            E_i = 0.5d0 * M_D_kg * (particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2) * J_TO_EV
            E_sum = E_sum + E_i
            E2_sum = E2_sum + E_i**2
         end if
      end do

      if (n_alive > 0) then
         E_mean = E_sum / n_alive
         E_std = sqrt(max(E2_sum / n_alive - E_mean*E_mean, 0.0d0))
         T_eff = E_mean * 2.0d0 / 3.0d0
      else
         E_mean = 0.0d0
         E_std = 0.0d0
         T_eff = 0.0d0
      end if

      write(*,'(A)')       '  --- Final Statistics ---'
      write(*,'(A,I8)')    '  Alive particles  = ', n_alive
      write(*,'(A,ES12.4)') '  Weight per particle = ', weight
      write(*,'(A,ES12.4)') '  Real particles (n_alive * w) = ', n_alive * weight
      write(*,'(A,F10.4,A)') '  Final E_mean = ', E_mean, ' eV'
      write(*,'(A,F10.4,A)') '  Final E_std  = ', E_std, ' eV'
      write(*,'(A,F10.4,A)') '  Final T_eff  = ', T_eff, ' eV'
      write(*,'(A,F10.4,A)') '  Target T_i   = ', plasma%T_i, ' eV'
      write(*,'(A,ES12.4,A)') '  Total time   = ', time, ' s'

   end subroutine output_statistics

end module io
