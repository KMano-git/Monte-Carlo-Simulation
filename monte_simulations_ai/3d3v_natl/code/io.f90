!===============================================================================
! Module: io
! 入出力: 入力ファイル読み込み、粒子初期化、出力
!===============================================================================
module io
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV
   use data_types, only: particle_t, sim_params, plasma_params, init_params, &
      diag_params, score_data
   use random_utils, only: sample_maxwell_velocity, set_beam_velocity
   implicit none

   private
   public :: read_input_file, initialize_particles
   public :: output_energy_histogram, output_statistics
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
      integer  :: n_particles, n_steps, seed
      real(dp) :: dt, weight_min
      logical  :: enable_cx, enable_el, enable_ei, use_isotropic
      character(len=256) :: cdf_file, output_ntscrg, output_hist

      !plasma namelist変数
      real(dp) :: n_i, T_i, n_e, T_e, u_x, u_y, u_z

      !particle_init namelist変数
      real(dp) :: n_init, E_init, T_init
      integer  :: init_mode

      !diagnostics namelist変数
      integer  :: output_interval, n_hist_bins
      real(dp) :: E_hist_min, E_hist_max
      integer  :: hist_timing(6)

      namelist /simulation/ n_particles, n_steps, dt, seed, &
         enable_cx, enable_el, enable_ei, use_isotropic, &
         weight_min, cdf_file, output_ntscrg, output_hist
      namelist /plasma_nml/ n_i, T_i, n_e, T_e, u_x, u_y, u_z
      namelist /particle_init/ n_init, E_init, T_init, init_mode
      namelist /diagnostics/ output_interval, n_hist_bins, &
         E_hist_min, E_hist_max, hist_timing

      !デフォルト値
      n_particles = 10000
      n_steps     = 1000
      dt          = 1.0d-8
      seed        = 12345
      enable_cx   = .true.
      enable_el   = .true.
      enable_ei   = .false.
      use_isotropic = .false.
      weight_min  = 1.0d-10
      cdf_file    = 'dd_00_elastic.cdf'
      output_ntscrg = 'ntscrg.csv'
      output_hist = 'energy_hist.csv'

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

      !ファイル読み込み
      inquire(file='input.nml', exist=file_exist)
      if (.not. file_exist) then
         write(*,*) 'Warning: input.nml not found, using defaults'
      else
         open(unit=unit_input, file='input.nml', status='old', iostat=ios)
         if (ios == 0) then
            read(unit_input, nml=simulation, iostat=ios)
            rewind(unit_input)
            read(unit_input, nml=plasma_nml, iostat=ios)
            rewind(unit_input)
            read(unit_input, nml=particle_init, iostat=ios)
            rewind(unit_input)
            read(unit_input, nml=diagnostics, iostat=ios)
            close(unit_input)
         end if
      end if

      !構造体に格納
      sim%n_particles   = n_particles
      sim%n_steps       = n_steps
      sim%dt            = dt
      sim%seed          = seed
      sim%enable_cx     = enable_cx
      sim%enable_el     = enable_el
      sim%enable_ei     = enable_ei
      sim%use_isotropic = use_isotropic
      sim%weight_min    = weight_min
      sim%cdf_file      = cdf_file
      sim%output_ntscrg = output_ntscrg
      sim%output_hist   = output_hist

      plasma%n_i = n_i
      plasma%T_i = T_i
      plasma%n_e = n_e
      plasma%T_e = T_e
      plasma%u_x = u_x
      plasma%u_y = u_y
      plasma%u_z = u_z

      init_p%n_init    = n_init
      init_p%E_init    = E_init
      init_p%T_init    = T_init
      init_p%init_mode = init_mode

      diag%output_interval = output_interval
      diag%n_hist_bins     = n_hist_bins
      diag%E_hist_min      = E_hist_min
      diag%E_hist_max      = E_hist_max
      diag%hist_timing     = hist_timing

      !パラメータ表示
      write(*,'(A)') '=========================================='
      write(*,'(A)') ' 3D3V Non-Analog Monte Carlo Simulation'
      write(*,'(A)') '=========================================='
      write(*,'(A,I10)')     '  n_particles  = ', n_particles
      write(*,'(A,I10)')     '  n_steps      = ', n_steps
      write(*,'(A,ES12.4)')  '  dt           = ', dt
      write(*,'(A,L1)')      '  enable_cx    = ', enable_cx
      write(*,'(A,L1)')      '  enable_el    = ', enable_el
      write(*,'(A,L1)')      '  enable_ei    = ', enable_ei
      write(*,'(A,ES12.4)')  '  weight_min   = ', weight_min
      write(*,'(A,ES12.4)')  '  n_i          = ', n_i
      write(*,'(A,F8.2,A)')  '  T_i          = ', T_i, ' eV'
      write(*,'(A,F8.2,A)')  '  T_e          = ', T_e, ' eV'
      write(*,'(A,F8.2,A)')  '  E_init       = ', E_init, ' eV'
      write(*,'(A)')         '=========================================='

   end subroutine read_input_file

   !---------------------------------------------------------------------------
   ! 粒子の初期化（位置=原点、重み=1.0）
   !---------------------------------------------------------------------------
   subroutine initialize_particles(particles, n_particles, init_p)
      type(particle_t), intent(out) :: particles(:)
      integer, intent(in)           :: n_particles
      type(init_params), intent(in) :: init_p

      integer :: i

      do i = 1, n_particles
         particles(i)%x = 0.0d0
         particles(i)%y = 0.0d0
         particles(i)%z = 0.0d0
         particles(i)%weight = 1.0d0
         particles(i)%alive = .true.

         if (init_p%init_mode == 0) then
            call set_beam_velocity(init_p%E_init, &
               particles(i)%vx, particles(i)%vy, particles(i)%vz)
         else
            call sample_maxwell_velocity(init_p%T_init, &
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
   ! ntscrg.csv ヘッダー出力
   !---------------------------------------------------------------------------
   subroutine output_ntscrg_header(filename)
      character(len=*), intent(in) :: filename
      integer :: unit_out

      unit_out = 50
      open(unit=unit_out, file=filename, status='replace')
      write(unit_out,'(A)') 'time[s],n_alive,weight_sum,' // &
         'CL_Q_cx[W/m3],CL_Q_el[W/m3],CL_Q_ei[W/m3],CL_Q_total[W/m3],' // &
         'TL_Q_cx[W/m3],TL_Q_el[W/m3],TL_Q_ei[W/m3],TL_Q_total[W/m3]'
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
      real(dp) :: cl_cx, cl_el, cl_ei, cl_total
      real(dp) :: tl_cx, tl_el, tl_ei, tl_total

      unit_out = 50

      ! 変換係数: [J] → [W/m3]
      ! score は1ステップの全粒子による合計[J]
      ! W/m3 = score[J] / dt[s] * n_init[m^-3] / N_particles
      norm = n_init / (dble(n_particles) * dt)
      time = dble(istep) * dt

      cl_cx = step_score%cl_cx * norm
      cl_el = step_score%cl_el * norm
      cl_ei = step_score%cl_ei * norm
      cl_total = cl_cx + cl_el + cl_ei

      tl_cx = step_score%tl_cx * norm
      tl_el = step_score%tl_el * norm
      tl_ei = step_score%tl_ei * norm
      tl_total = tl_cx + tl_el + tl_ei

      open(unit=unit_out, file=filename, status='old', position='append')
      write(unit_out,'(ES14.6,A,I8,A,ES14.6,8(A,ES14.6))') &
         time, ',', n_alive, ',', weight_sum, &
         ',', cl_cx, ',', cl_el, ',', cl_ei, ',', cl_total, &
         ',', tl_cx, ',', tl_el, ',', tl_ei, ',', tl_total
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
      real(dp) :: cl_total, tl_total

      ! 全ステップ累積スコア → 平均 W/m3
      norm = n_init / (dble(n_particles) * dt * dble(n_steps))

      cl_total = (score%cl_ei + score%cl_cx + score%cl_el) * norm
      tl_total = (score%tl_ei + score%tl_cx + score%tl_el) * norm

      write(*,'(A)') ''
      write(*,'(A)') '=========================================='
      write(*,'(A)') ' Scoring Results (time-averaged) [W/m3]'
      write(*,'(A)') '=========================================='
      write(*,'(A)')         '          EI              CX              EL              Total'
      write(*,'(A,4ES16.6)') ' CL: ', &
         score%cl_ei*norm, score%cl_cx*norm, score%cl_el*norm, cl_total
      write(*,'(A,4ES16.6)') ' TL: ', &
         score%tl_ei*norm, score%tl_cx*norm, score%tl_el*norm, tl_total
      write(*,'(A)') '=========================================='

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
