!===============================================================================
! Module: io
! 入出力: 入力ファイル読み込み、粒子初期化、出力
!===============================================================================
module io
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV
   use data_types, only: particle_t, sim_params, plasma_params, init_params, &
      diag_params, score_data
   use random_utils, only: sample_maxwell_velocity, set_beam_velocity, init_rng, random_double
   implicit none

   private
   public :: read_input_file, initialize_particles
   public :: output_energy_histogram, output_deltaE_histogram, output_statistics
   public :: output_ntscrg_header, output_ntscrg_step, output_ntscrg_final

contains

   !---------------------------------------------------------------------------
   ! 入力ファイル(input.nml)の読み込み
   !---------------------------------------------------------------------------
   subroutine read_input_file(sim, plasma, init_p, diag)
      type(sim_params), intent(out)    :: sim
      type(plasma_params), intent(out) :: plasma
      type(init_params), intent(out)   :: init_p
      type(diag_params), intent(out)   :: diag

      integer, parameter :: unit_input = 20
      logical :: file_exist
      integer :: ios

      !simulation namelist変数
      integer  :: n_particles, n_steps, seed, tl_el_inner_samples
      real(dp) :: dt, weight_min
      logical  :: enable_cx, enable_el, enable_ei, use_isotropic
      character(len=256) :: cdf_file, tl_el_table_file
      character(len=256) :: output_ntscrg, output_hist, output_deltaE_hist

      !plasma namelist変数
      real(dp) :: n_i, T_i, n_e, T_e, u_x, u_y, u_z

      !particle_init namelist変数
      real(dp) :: n_init, E_init, T_init
      integer  :: init_mode

      !diagnostics namelist変数
      integer  :: output_interval, n_hist_bins
      real(dp) :: E_hist_min, E_hist_max
      integer  :: hist_timing(6)
      integer  :: n_dE_bins, dE_collect_steps
      real(dp) :: dE_hist_min, dE_hist_max

      namelist /simulation/ n_particles, n_steps, dt, seed, tl_el_inner_samples, &
         enable_cx, enable_el, enable_ei, use_isotropic, &
         weight_min, cdf_file, tl_el_table_file, output_ntscrg, output_hist, &
         output_deltaE_hist
      namelist /plasma_nml/ n_i, T_i, n_e, T_e, u_x, u_y, u_z
      namelist /particle_init/ n_init, E_init, T_init, init_mode
      namelist /diagnostics/ output_interval, n_hist_bins, &
         E_hist_min, E_hist_max, hist_timing, &
         n_dE_bins, dE_hist_min, dE_hist_max, dE_collect_steps

      !デフォルト値
      n_particles = 10000
      n_steps     = 1000
      dt          = 1.0d-8
      seed        = 12345
      tl_el_inner_samples = 16
      enable_cx   = .true.
      enable_el   = .true.
      enable_ei   = .false.
      use_isotropic = .false.
      weight_min  = 1.0d-10
      cdf_file    = 'dd_00_elastic.cdf'
      tl_el_table_file = 'tl_el_table.dat'
      output_ntscrg = 'ntscrg.csv'
      output_hist = 'energy_hist.csv'
      output_deltaE_hist = 'deltaE_hist.csv'

      n_i = 1.0d21
      T_i = 2.0d0
      n_e = 1.0d21
      T_e = 2.0d0
      u_x = 0.0d0
      u_y = 0.0d0
      u_z = 0.0d0

      n_init    = 5.0d19
      E_init    = 3.0d0
      T_init    = 2.0d0
      init_mode = 0

      output_interval = 100
      n_hist_bins     = 400
      E_hist_min      = 0.0d0
      E_hist_max      = 150.0d0
      hist_timing     = (/10, 20, 50, 100, 500, 1000/)
      n_dE_bins       = 400
      dE_hist_min     = -100.0d0
      dE_hist_max     = 100.0d0
      dE_collect_steps = 0

      !ファイル読み込み
      inquire(file='input.nml', exist=file_exist)
      if (.not. file_exist) then
         write(*,*) 'Warning: input.nml not found, using defaults'
      else
         open(unit=unit_input, file='input.nml', status='old', iostat=ios)
         if (ios /= 0) then
            write(*,'(A,I0)') 'Error: Failed to open input.nml, iostat=', ios
            stop 1
         end if

         read(unit_input, nml=simulation, iostat=ios)
         if (ios /= 0) call abort_namelist_read('simulation', ios)
         rewind(unit_input)

         read(unit_input, nml=plasma_nml, iostat=ios)
         if (ios /= 0) call abort_namelist_read('plasma_nml', ios)
         rewind(unit_input)

         read(unit_input, nml=particle_init, iostat=ios)
         if (ios /= 0) call abort_namelist_read('particle_init', ios)
         rewind(unit_input)

         read(unit_input, nml=diagnostics, iostat=ios)
         if (ios /= 0) call abort_namelist_read('diagnostics', ios)
         close(unit_input)
      end if

      !構造体に格納
      sim%n_particles   = n_particles
      sim%n_steps       = n_steps
      sim%dt            = dt
      sim%seed          = seed
      sim%elastic_inner_samples = tl_el_inner_samples
      sim%enable_cx     = enable_cx
      sim%enable_el     = enable_el
      sim%enable_ionization = enable_ei
      sim%use_isotropic = use_isotropic
      sim%weight_min    = weight_min
      sim%cdf_file      = cdf_file
      sim%elastic_tl_table_file = tl_el_table_file
      sim%output_ntscrg = output_ntscrg
      sim%output_hist   = output_hist
      sim%output_delta_energy_hist = output_deltaE_hist

      plasma%ion_density = n_i
      plasma%ion_temperature_eV = T_i
      plasma%electron_density = n_e
      plasma%electron_temperature_eV = T_e
      plasma%ion_flow_vx = u_x
      plasma%ion_flow_vy = u_y
      plasma%ion_flow_vz = u_z

      init_p%initial_density = n_init
      init_p%initial_energy_eV = E_init
      init_p%initial_temperature_eV = T_init
      init_p%initial_velocity_mode = init_mode

      diag%output_interval = output_interval
      diag%energy_hist_bin_count = n_hist_bins
      diag%energy_hist_min_eV    = E_hist_min
      diag%energy_hist_max_eV    = E_hist_max
      diag%hist_timing     = hist_timing
      diag%delta_energy_hist_bin_count = n_dE_bins
      diag%delta_energy_hist_min_eV    = dE_hist_min
      diag%delta_energy_hist_max_eV    = dE_hist_max
      diag%delta_energy_collect_steps  = dE_collect_steps

      !パラメータ表示
      write(*,'(A)') '=========================================='
      write(*,'(A)') ' 3D3V Non-Analog Monte Carlo Simulation'
      write(*,'(A)') '=========================================='
      write(*,'(A,I10)')     '  n_particles  = ', n_particles
      write(*,'(A,I10)')     '  n_steps      = ', n_steps
      write(*,'(A,ES12.4)')  '  dt           = ', dt
      write(*,'(A,I10)')     '  tl_el_inner  = ', tl_el_inner_samples
      write(*,'(A,L1)')      '  enable_cx    = ', enable_cx
      write(*,'(A,L1)')      '  enable_el    = ', enable_el
      write(*,'(A,L1)')      '  enable_ei    = ', enable_ei
      write(*,'(A,ES12.4)')  '  weight_min   = ', weight_min
      write(*,'(A,A)')       '  tl_el_table  = ', trim(tl_el_table_file)
      write(*,'(A,ES12.4)')  '  n_i          = ', n_i
      write(*,'(A,F8.2,A)')  '  T_i          = ', T_i, ' eV'
      write(*,'(A,F8.2,A)')  '  T_e          = ', T_e, ' eV'
      if (init_mode == 0) then
         write(*,'(A,F8.2,A)')  '  E_init       = ', E_init, ' eV'
      else
         write(*,'(A,F8.2,A)')  '  T_init       = ', T_init, ' eV'
      end if
      write(*,'(A)')         '=========================================='

   end subroutine read_input_file

   subroutine abort_namelist_read(section_name, ios)
      character(len=*), intent(in) :: section_name
      integer, intent(in) :: ios

      write(*,'(A,A,A,I0)') 'Error: Failed to read namelist /', trim(section_name), &
         '/ from input.nml, iostat=', ios
      stop 1
   end subroutine abort_namelist_read

   !---------------------------------------------------------------------------
   ! 粒子の初期化（位置=原点、重み=1.0）
   !---------------------------------------------------------------------------
   subroutine initialize_particles(particles, n_particles, init_p, seed)
      type(particle_t), intent(out) :: particles(:)
      integer, intent(in)           :: n_particles
      type(init_params), intent(in) :: init_p
      integer, intent(in)           :: seed

      integer :: i
      integer(8) :: particle_seed

      do i = 1, n_particles
         !粒子固有のRNG初期化（seedとindexから一意のシードを生成）
         particle_seed = int(seed, 8) * 6364136223846793005_8 + int(i, 8)
         call init_rng(particles(i)%rng, particle_seed)

         particles(i)%x = 0.0d0
         particles(i)%y = 0.0d0
         particles(i)%z = 0.0d0
         particles(i)%weight = 1.0d0
         particles(i)%alive = .true.
         particles(i)%collision_clock_elapsed = 0.0d0
         particles(i)%pending_track_time = 0.0d0
         particles(i)%pending_effective_track_time = 0.0d0
         particles(i)%collision_clock_target = -log(max(random_double(particles(i)%rng), 1.0d-30))

         if (init_p%initial_velocity_mode == 0) then
            call set_beam_velocity(init_p%initial_energy_eV, &
               particles(i)%vx, particles(i)%vy, particles(i)%vz)
         else
            call sample_maxwell_velocity(particles(i)%rng, init_p%initial_temperature_eV, &
               particles(i)%vx, particles(i)%vy, particles(i)%vz)
         end if
      end do
   end subroutine initialize_particles

   !---------------------------------------------------------------------------
   ! エネルギーヒストグラム出力
   !---------------------------------------------------------------------------
   subroutine output_energy_histogram(filename, particles, n_particles, &
      n_bins, E_min, E_max, step, n_timing, &
      hist_timing, first_call)
      character(len=*), intent(in) :: filename
      type(particle_t), intent(in) :: particles(:)
      integer, intent(in)  :: n_particles, n_bins, step, n_timing
      real(dp), intent(in) :: E_min, E_max
      integer, intent(in)  :: hist_timing(:)
      logical, intent(in)  :: first_call

      real(dp) :: bin_width, E_particle, weight_sum
      real(dp), allocatable :: hist(:)
      integer :: i, ibin, unit_out
      logical :: should_output

      !出力タイミングチェック
      should_output = .false.
      do i = 1, n_timing
         if (step == hist_timing(i)) then
            should_output = .true.
            exit
         end if
      end do
      if (.not. should_output) return

      allocate(hist(n_bins))
      hist = 0.0d0
      bin_width = (E_max - E_min) / dble(n_bins)
      weight_sum = 0.0d0

      do i = 1, n_particles
         if (.not. particles(i)%alive) cycle
         E_particle = 0.5d0 * M_D_kg * &
            (particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2) &
            * J_TO_EV
         ibin = int((E_particle - E_min) / bin_width) + 1
         if (ibin >= 1 .and. ibin <= n_bins) then
            hist(ibin) = hist(ibin) + particles(i)%weight
            weight_sum = weight_sum + particles(i)%weight
         end if
      end do

      !CSV出力
      unit_out = 40
      if (first_call) then
         open(unit=unit_out, file=filename, status='replace')
         write(unit_out,'(A)', advance='no') 'step,E_min,E_max,bin_width,n_bins'
         do i = 1, n_bins
            write(unit_out,'(A,I0)', advance='no') ',bin_', i
         end do
         write(unit_out,*)
      else
         open(unit=unit_out, file=filename, status='old', position='append')
      end if

      write(unit_out,'(I0,A,ES14.6,A,ES14.6,A,ES14.6,A,I0)', advance='no') &
         step, ',', E_min, ',', E_max, ',', bin_width, ',', n_bins
      do i = 1, n_bins
         write(unit_out,'(A,ES14.6)', advance='no') ',', hist(i)
      end do
      write(unit_out,*)

      close(unit_out)
      deallocate(hist)

      write(*,'(A,I6,A,ES12.4)') '  Histogram output at step ', step, &
         ', weight_sum = ', weight_sum

   end subroutine output_energy_histogram

   !---------------------------------------------------------------------------
   ! エネルギー移行量（ΔE）ヒストグラム出力
   !---------------------------------------------------------------------------
   subroutine output_deltaE_histogram(filename, dE_hist_el, dE_hist_cx, &
      n_bins, dE_min, dE_max, &
      n_init, n_particles, collect_time)
      character(len=*), intent(in) :: filename
      real(dp), intent(in) :: dE_hist_el(:), dE_hist_cx(:)
      integer, intent(in)  :: n_bins
      real(dp), intent(in) :: dE_min, dE_max
      real(dp), intent(in) :: n_init
      integer, intent(in)  :: n_particles
      real(dp), intent(in) :: collect_time  !集計期間 [s]

      real(dp) :: bin_width, norm
      integer :: i, unit_out

      bin_width = (dE_max - dE_min) / dble(n_bins)

      ! 蓄積された重みを [events / (m^3 * s)] つまり反応率密度に変換する係数
      ! density_rate = sum_weight * n_init / (N_particles * collect_time)
      norm = n_init / (dble(n_particles) * collect_time)

      !CSV出力
      unit_out = 41
      open(unit=unit_out, file=filename, status='replace')
      write(unit_out,'(A)') 'step,time,dE_min,dE_max,bin_width,n_bins'
      write(unit_out,'(I0,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,I0)') &
         0, ',', collect_time, ',', dE_min, ',', dE_max, ',', bin_width, ',', n_bins
      write(unit_out,'(A)') 'bin_center_eV,rate_EL_m-3s-1,rate_CX_m-3s-1'

      do i = 1, n_bins
         write(unit_out,'(ES14.6,A,ES14.6,A,ES14.6)') &
            dE_min + (dble(i) - 0.5d0) * bin_width, ',', &
            dE_hist_el(i) * norm, ',', &
            dE_hist_cx(i) * norm
      end do

      close(unit_out)

      write(*,'(A,ES12.4,A)') '  dE Histogram output (collect_time=', collect_time, ' s)'

   end subroutine output_deltaE_histogram

   !---------------------------------------------------------------------------
   ! ntscrg.csv ヘッダー出力
   !---------------------------------------------------------------------------
   subroutine output_ntscrg_header(filename)
      character(len=*), intent(in) :: filename
      integer :: unit_out

      unit_out = 50
      open(unit=unit_out, file=filename, status='replace')
      write(unit_out,'(A)') 'time[s],n_alive,weight_sum,' // &
         'A_Q_el[W/m3],' // &
         'CL_Q_cx[W/m3],CL_Q_el[W/m3],CL_Q_ei[W/m3],CL_Q_total[W/m3],CLInner_Q_el[W/m3],' // &
         'TR_Q_cx[W/m3],TR_Q_el[W/m3],TR_Q_ei[W/m3],TR_Q_total[W/m3],TRInner_Q_el[W/m3],TRPretab_Q_el[W/m3]'
      close(unit_out)
   end subroutine output_ntscrg_header

   !---------------------------------------------------------------------------
   ! ntscrg.csv 1ステップ分のスコアリング出力 [W/m3]
   ! 変換: score[J] * n_init / (N_particles * dt) = [W/m3]
   !---------------------------------------------------------------------------
   subroutine output_ntscrg_step(filename, step_score, &
      n_particles, n_init, dt, &
      istep, n_alive, weight_sum)
      character(len=*), intent(in) :: filename
      type(score_data), intent(in) :: step_score
      integer, intent(in)  :: n_particles, istep, n_alive
      real(dp), intent(in) :: n_init, dt, weight_sum

      integer :: unit_out
      real(dp) :: norm, time
      real(dp) :: a_el
      real(dp) :: cl_cx, cl_el, cl_ei, cl_total
      real(dp) :: cl_inner_el
      real(dp) :: tr_cx, tr_el, tr_ei, tr_total
      real(dp) :: tr_inner_el, tr_pretab_el

      unit_out = 50

      ! 変換係数: [J] → [W/m3]
      ! score は1ステップの全粒子による合計[J]
      ! W/m3 = score[J] / dt[s] * n_init[m^-3] / N_particles
      norm = n_init / (dble(n_particles) * dt)
      time = dble(istep) * dt

      a_el = step_score%a_el * norm
      cl_cx = step_score%cl_cx * norm
      cl_el = step_score%cl_el * norm
      cl_ei = step_score%cl_ei * norm
      cl_total = cl_cx + cl_el + cl_ei
      cl_inner_el = step_score%cl_inner_el * norm

      tr_cx = step_score%tr_cx * norm
      tr_el = step_score%tr_el * norm
      tr_ei = step_score%tr_ei * norm
      tr_total = tr_cx + tr_el + tr_ei
      tr_inner_el = step_score%tr_inner_el * norm
      tr_pretab_el = step_score%tr_pretab_el * norm

      open(unit=unit_out, file=filename, status='old', position='append')
      write(unit_out,'(ES14.6,A,I8,A,ES14.6,12(A,ES14.6))') &
         time, ',', n_alive, ',', weight_sum, &
         ',', a_el, &
         ',', cl_cx, ',', cl_el, ',', cl_ei, ',', cl_total, ',', cl_inner_el, &
         ',', tr_cx, ',', tr_el, ',', tr_ei, ',', tr_total, ',', tr_inner_el, ',', tr_pretab_el
      close(unit_out)

   end subroutine output_ntscrg_step

   !---------------------------------------------------------------------------
   ! 最終スコアリングサマリー（コンソール出力）
   !---------------------------------------------------------------------------
   subroutine output_ntscrg_final(score, n_particles, n_init, n_steps, dt)
      type(score_data), intent(in) :: score
      integer, intent(in)  :: n_particles, n_steps
      real(dp), intent(in) :: n_init, dt

      real(dp) :: norm
      real(dp) :: cl_total, tr_total
      character(len=*), parameter :: score_sep = &
         '================================================================================'

      ! 全ステップ累積スコア → 平均 W/m3
      norm = n_init / (dble(n_particles) * dt * dble(n_steps))

      cl_total = (score%cl_ei + score%cl_cx + score%cl_el) * norm
      tr_total = (score%tr_ei + score%tr_cx + score%tr_el) * norm

      write(*,'(A)') ''
      write(*,'(A)') score_sep
      write(*,'(A)') ' Scoring Results (time-averaged) [W/m3]'
      write(*,'(A)') score_sep
      write(*,'(A16,4A16)') '', 'EI', 'CX', 'EL', 'Total'
      write(*,'(A16,2A16,ES16.6,A16)') ' A(EL):', '', '', score%a_el*norm, ''
      write(*,'(A16,4ES16.6)') ' CL:', &
         score%cl_ei*norm, score%cl_cx*norm, score%cl_el*norm, cl_total
      write(*,'(A16,2A16,ES16.6,A16)') ' CLInner(EL):', '', '', score%cl_inner_el*norm, ''
      write(*,'(A16,4ES16.6)') ' TR:', &
         score%tr_ei*norm, score%tr_cx*norm, score%tr_el*norm, tr_total
      write(*,'(A16,2A16,ES16.6,A16)') ' TRInner(EL):', '', '', score%tr_inner_el*norm, ''
      write(*,'(A16,2A16,ES16.6,A16)') ' TRPretab(EL):', '', '', score%tr_pretab_el*norm, ''
      write(*,'(A)') score_sep

   end subroutine output_ntscrg_final

   !---------------------------------------------------------------------------
   ! 最終統計情報の出力
   !---------------------------------------------------------------------------
   subroutine output_statistics(particles, n_particles)
      type(particle_t), intent(in) :: particles(:)
      integer, intent(in) :: n_particles

      integer :: i, n_alive
      real(dp) :: E_sum, E_sq, E_mean, E_var, v_sq
      real(dp) :: weight_sum

      n_alive = 0
      E_sum = 0.0d0
      E_sq  = 0.0d0
      weight_sum = 0.0d0

      do i = 1, n_particles
         if (.not. particles(i)%alive) cycle
         n_alive = n_alive + 1
         v_sq = particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2
         E_sum = E_sum + 0.5d0 * M_D_kg * v_sq * J_TO_EV * particles(i)%weight
         E_sq  = E_sq  + (0.5d0 * M_D_kg * v_sq * J_TO_EV)**2 * particles(i)%weight
         weight_sum = weight_sum + particles(i)%weight
      end do

      write(*,'(A)') ''
      write(*,'(A)') '=========================================='
      write(*,'(A)') ' Final Statistics'
      write(*,'(A)') '=========================================='
      write(*,'(A,I10)')    '  Alive particles  = ', n_alive
      write(*,'(A,F12.4)')  '  Total weight     = ', weight_sum

      if (weight_sum > 0.0d0) then
         E_mean = E_sum / weight_sum
         E_var  = E_sq / weight_sum - E_mean**2
         write(*,'(A,F12.4,A)')  '  Mean energy      = ', E_mean, ' eV'
         write(*,'(A,F12.4,A)')  '  Eff. temperature = ', E_mean * 2.0d0 / 3.0d0, ' eV'
      end if
      write(*,'(A)') '=========================================='

   end subroutine output_statistics

end module io
