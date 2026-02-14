!===============================================================================
! Module: io
! 入出力関連
!===============================================================================
module io
   use constants, only: dp, M_D_kg, J_TO_EV, EV_TO_J
   use particle_data, only: particle_t, sim_params, plasma_params, init_params, diag_params
   use random_utils, only: set_beam_velocity, sample_maxwell_velocity
   implicit none

contains

   subroutine read_input_file(sim, plasma, init_p, diag)
      type(sim_params), intent(out)  :: sim
      type(plasma_params), intent(out) :: plasma
      type(init_params), intent(out) :: init_p
      type(diag_params), intent(out) :: diag

      integer, parameter :: unit_input = 20
      logical :: file_exist

      namelist /simulation/ n_particles, n_steps, dt, n_grid, x_min, x_max, seed, &
         enable_cx, enable_el, enable_ei, use_isotropic, &
         cdf_file, output_ntscrg, output_hist, output_statx
      namelist /plasma_nml/ n_i, T_i, n_e, T_e
      namelist /particle_init/ E_init, T_init, gamma_in, init_mode
      namelist /diagnostics/ output_interval, n_hist_bins, E_hist_min, E_hist_max, hist_timing

      ! デフォルト値
      integer :: n_particles = 10000
      integer :: n_steps = 10000
      real(dp) :: dt = 1.0d-8
      integer :: n_grid = 100
      real(dp) :: x_min = 0.0d0
      real(dp) :: x_max = 0.1d0
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

      real(dp) :: E_init = 100.0d0
      real(dp) :: T_init = 2.0d0
      real(dp) :: n_init = 1.0d19
      integer :: init_mode = 0

      integer :: output_interval = 100
      integer :: n_hist_bins = 100
      real(dp) :: E_hist_min = 0.0d0
      real(dp) :: E_hist_max = 150.0d0
      integer :: hist_timing(5) = (/ 1, 10, 100, 1000, 10000 /) ! デフォルト値

      ! ファイル確認
      inquire(file = 'input.nml', exist=file_exist)
      if (.not. file_exist) then
         write(*,*) 'Error: input.nml not found'
         stop
      end if

      ! ファイルオープン
      open(unit=unit_input, file='input.nml', status='old')

      ! Namelistの読み込み
      ! 順序次第ではrewind関数が必要
      read(unit_input, nml=simulation)
      read(unit_input, nml=plasma_nml)
      read(unit_input, nml=particle_init)
      read(unit_input, nml=diagnostics)

      ! ファイルクローズ
      close(unit_input)

      ! パラメータのセット
      sim%n_particles = n_particles
      sim%n_steps = n_steps
      sim%dt = dt
      sim%n_grid = n_grid
      sim%x_min = x_min
      sim%x_max = x_max
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
      init_p%n_init = n_init
      init_p%init_mode = init_mode

      diag%output_interval = output_interval
      diag%n_hist_bins = n_hist_bins
      diag%E_hist_min = E_hist_min
      diag%E_hist_max = E_hist_max
      diag%hist_timing = hist_timing

   end subroutine read_input_file

   !-----------------------------------------------------------------------------
   ! 粒子初期化
   !-----------------------------------------------------------------------------
   subroutine initialize_particles(particles, n_particles, init_p)
      type(particle_t), intent(out) :: particles(:)
      integer, intent(in) :: n_particles
      type(init_params), intent(in) :: init_p

      integer :: i

      do i = 1, n_particles
         particles(i)%alive = .true.

         if (init_p%init_mode == 0) then
            ! ビー�������（+x方向）
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
   ! 中性粒子のエネルギーヒストグラム出力
   !-----------------------------------------------------------------------------
   subroutine output_energy_histogram(particles, n_particles, sim, diag, init_p, istep)
      type(particle_t), intent(in) :: particles(:)
      integer, intent(in) :: n_particles
      type(sim_params), intent(in)  :: sim
      type(diag_params), intent(in) :: diag
      type(init_params), intent(in) :: init_p
      integer, intent(in) :: istep

      integer :: i
      integer, parameter :: unit_hist = 20
      integer :: bin_index
      real(dp) :: E_ev
      real(dp) :: bin_width
      real(dp) :: E_min, E_max
      integer, allocatable :: hist(:)

      ! 単位体積あたりの粒子数を計算
      real(dp) :: n_per_volume
      n_per_volume = dble(n_particles) * init_p%n_init

      ! ヒストグラムの範囲とビン幅
      E_min = diag%E_hist_min
      E_max = diag%E_hist_max
      bin_width = (E_max - E_min) / dble(diag%n_hist_bins)

      ! ヒストグラム配列の確保
      allocate(hist(diag%n_hist_bins))
      hist = 0

      ! エネルギーのカウント
      do i = 1, n_particles
         if (.not. particles(i)%alive) cycle

         ! 粒子のエネルギーをeVに変換
         E_ev = 0.5d0 * M_D_kg * (particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2) * J_TO_EV

         ! ビン番号を計算
         if (E_ev >= E_min .and. E_ev < E_max) then
            bin_index = int((E_ev - E_min) / bin_width) + 1
            hist(bin_index) = hist(bin_index) + 1
         end if
      end do

      ! ファイルオープン
      open(unit=unit_hist, file=trim(sim%output_hist), status='unknown',position='append')

      ! ヘッダー書き込み（初回のみ）
      if (istep == 0) then
         write(unit_hist, '(A)', advance='no') 'step,E_min,E_max,bin_width,n_per_volume'
         do i = 1, diag%n_hist_bins
            write(unit_hist, '(A,I0)', advance='no') ',bin_', i
         end do
         write(unit_hist, *)
      end if

      ! データ書き込み
      write(unit_hist, '(I0,",",F0.4,",",F0.4,",",F0.4,",",ES12.4)', advance='no') &
         istep, E_min, E_max, bin_width, n_per_volume
      do i = 1, diag%n_hist_bins
         write(unit_hist, '(",",I0)', advance='no') hist(i)
      end do
      write(unit_hist, *)

      ! ファイルクローズ
      close(unit_hist)

      ! 配列解放
      deallocate(hist)

   end subroutine output_energy_histogram

   !-----------------------------------------------------------------------------
   ! 最終出力
   !-----------------------------------------------------------------------------
   subroutine output_statistics(particles, n_particles, plasma, init_p, diag, time)
      type(particle_t), intent(in) :: particles(:)
      integer, intent(in) :: n_particles
      type(plasma_params), intent(in) :: plasma
      type(init_params), intent(in) :: init_p
      type(diag_params), intent(in) :: diag
      real(dp), intent(in) :: time

      real(dp) :: E_ev
      real(dp) :: E_mean, E_std, T_eff
      real(dp) :: E_sum, E2_sum, E_i
      integer  :: i, n_alive

      ! エネルギー統計量の計算
      E_mean = 0.0d0
      E_std = 0.0d0
      T_eff = 0.0d0
      n_alive = 0
      E_sum = 0.0d0
      E2_sum = 0.0d0

      do i = 1, n_particles
         if (particles(i)%alive) then
            n_alive = n_alive + 1
            ! 粒��のエネルギー[eV]
            E_i = 0.5d0 * M_D_kg * (particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2) * J_TO_EV
            E_sum = E_sum + E_i
            E2_sum = E2_sum + E_i**2
         end if
      end do

      if (n_alive > 0) then
         E_mean = E_sum / n_alive
         E_std = sqrt(max(E2_sum / n_alive - E_mean*E_mean, 0.0d0))
         ! ��効温度[eV]: 3/2kT_eff = E_mean
         T_eff = E_mean * 2.0d0 / 3.0d0
      end if

      write(*,'(A,F10.4,A)') '  Final E_mean = ', E_mean, ' eV'
      write(*,'(A,F10.4,A)') '  Final T_eff  = ', T_eff, ' eV'
      write(*,'(A,F10.4,A)') '  Target T_i   = ', plasma%T_i, ' eV'

   end subroutine output_statistics

end module io
