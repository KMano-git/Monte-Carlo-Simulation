!===============================================================================
! Module: io
! Input, particle initialization, and diagnostics output.
!===============================================================================
module io
   use constants, only: dp, M_D_kg, J_TO_EV
   use data_types, only: particle_t, sim_params, plasma_params, init_params, &
      diag_params, score_data, score_stats_data, source_terms_t, &
      estimator_score_t, estimator_stats_t, source_stats_t, running_stats_t, &
      system_totals_t, REACT_EI, REACT_CX, REACT_EL, REACT_TOTAL, N_REACTIONS
   use random_utils, only: sample_maxwell_velocity, set_beam_velocity, init_rng, &
      random_double
   use statistics, only: running_stddev, running_stderr
   implicit none

   private
   public :: read_input_file, initialize_particles
   public :: output_energy_histogram, output_deltaE_histogram, output_statistics
   public :: output_sources_header, output_sources_step
   public :: output_ntscrg_header, output_ntscrg_step, output_ntscrg_final
   public :: output_source_stats, output_balance

contains

   subroutine read_input_file(sim, plasma, init_p, diag)
      type(sim_params), intent(out)    :: sim
      type(plasma_params), intent(out) :: plasma
      type(init_params), intent(out)   :: init_p
      type(diag_params), intent(out)   :: diag

      integer, parameter :: unit_input = 20
      logical :: file_exist
      integer :: ios

      integer  :: n_particles, n_steps, seed
      real(dp) :: dt, weight_min
      logical  :: enable_cx, enable_el, enable_ei, compare_estimator
      character(len=256) :: cdf_file, tl_el_table_file
      character(len=256) :: output_sources, output_stats, output_balance
      character(len=256) :: output_ntscrg, output_hist, output_deltaE_hist

      real(dp) :: n_i, T_i, n_e, T_e, u_x, u_y, u_z
      real(dp) :: n_init, E_init, T_init
      integer  :: init_mode

      integer  :: output_interval, n_hist_bins
      real(dp) :: E_hist_min, E_hist_max
      integer  :: hist_timing(6)
      integer  :: n_dE_bins, dE_collect_steps
      real(dp) :: dE_hist_min, dE_hist_max

      namelist /simulation/ n_particles, n_steps, dt, seed, &
         enable_cx, enable_el, enable_ei, compare_estimator, &
         weight_min, cdf_file, tl_el_table_file, &
         output_sources, output_stats, output_balance, output_ntscrg, &
         output_hist, output_deltaE_hist
      namelist /plasma_nml/ n_i, T_i, n_e, T_e, u_x, u_y, u_z
      namelist /particle_init/ n_init, E_init, T_init, init_mode
      namelist /diagnostics/ output_interval, n_hist_bins, &
         E_hist_min, E_hist_max, hist_timing, &
         n_dE_bins, dE_hist_min, dE_hist_max, dE_collect_steps

      n_particles = 10000
      n_steps = 100
      dt = 1.0d-9
      seed = 12345
      enable_cx = .true.
      enable_el = .true.
      enable_ei = .true.
      compare_estimator = .false.
      weight_min = 1.0d-10
      cdf_file = 'dd_00_elastic.cdf'
      tl_el_table_file = 'tl_el_table.dat'
      output_sources = 'sources.csv'
      output_stats = 'source_stats.csv'
      output_balance = 'balance.csv'
      output_ntscrg = 'ntscrg.csv'
      output_hist = 'energy_hist.csv'
      output_deltaE_hist = 'delta_E_hist.csv'

      n_i = 1.0d21
      T_i = 2.0d0
      n_e = 1.0d21
      T_e = 2.0d0
      u_x = 0.0d0
      u_y = 0.0d0
      u_z = 0.0d0

      n_init = 5.0d19
      E_init = 3.0d0
      T_init = 2.0d0
      init_mode = 1

      output_interval = 100
      n_hist_bins = 400
      E_hist_min = 0.0d0
      E_hist_max = 150.0d0
      hist_timing = (/1, 10, 50, 100, 500, 1000/)
      n_dE_bins = 400
      dE_hist_min = -20.0d0
      dE_hist_max = 20.0d0
      dE_collect_steps = 0

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

      sim%n_particles = n_particles
      sim%n_steps = n_steps
      sim%dt = dt
      sim%seed = seed
      sim%enable_cx = enable_cx
      sim%enable_el = enable_el
      sim%enable_ionization = enable_ei
      sim%compare_estimator = compare_estimator
      sim%weight_min = weight_min
      sim%cdf_file = cdf_file
      sim%elastic_tl_table_file = tl_el_table_file
      sim%output_sources = output_sources
      sim%output_stats = output_stats
      sim%output_balance = output_balance
      sim%output_ntscrg = output_ntscrg
      sim%output_hist = output_hist
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
      diag%energy_hist_min_eV = E_hist_min
      diag%energy_hist_max_eV = E_hist_max
      diag%hist_timing = hist_timing
      diag%delta_energy_hist_bin_count = n_dE_bins
      diag%delta_energy_hist_min_eV = dE_hist_min
      diag%delta_energy_hist_max_eV = dE_hist_max
      diag%delta_energy_collect_steps = dE_collect_steps

      write(*,'(A)') '=========================================='
      write(*,'(A)') ' 3D3V Monte Carlo Simulation ver 1.0.0'
      write(*,'(A)') '=========================================='
      write(*,'(A,I10)')     '  n_particles      = ', n_particles
      write(*,'(A,I10)')     '  n_steps          = ', n_steps
      write(*,'(A,ES12.4)')  '  dt               = ', dt
      write(*,'(A,L1)')      '  enable_cx        = ', enable_cx
      write(*,'(A,L1)')      '  enable_el        = ', enable_el
      write(*,'(A,L1)')      '  enable_ei        = ', enable_ei
      write(*,'(A,L1)')      '  compare_estimator= ', compare_estimator
      write(*,'(A,A)')       '  cdf_file         = ', trim(cdf_file)
      write(*,'(A,ES12.4)')  '  n_i              = ', n_i
      write(*,'(A,F8.2,A)')  '  T_i              = ', T_i, ' eV'
      write(*,'(A,F8.2,A)')  '  T_e              = ', T_e, ' eV'
      write(*,'(A)')         '=========================================='
   end subroutine read_input_file

   subroutine abort_namelist_read(section_name, ios)
      character(len=*), intent(in) :: section_name
      integer, intent(in) :: ios

      write(*,'(A,A,A,I0)') 'Error: Failed to read namelist /', trim(section_name), &
         '/ from input.nml, iostat=', ios
      stop 1
   end subroutine abort_namelist_read

   subroutine initialize_particles(particles, n_particles, init_p, seed)
      type(particle_t), intent(out) :: particles(:)
      integer, intent(in) :: n_particles
      type(init_params), intent(in) :: init_p
      integer, intent(in) :: seed

      integer :: i
      integer(8) :: particle_seed

      do i = 1, n_particles
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

   subroutine output_sources_header(filename)
      character(len=*), intent(in) :: filename
      integer :: unit_out

      unit_out = 50
      open(unit=unit_out, file=filename, status='replace')
      write(unit_out,'(A)') 'time[s],estimator,reaction,' // &
         'Sn_plus[m^-3 s^-1],Sn_minus[m^-3 s^-1],Sn_net[m^-3 s^-1],' // &
         'Spx[N m^-3],Spy[N m^-3],Spz[N m^-3],We[W m^-3],Wi[W m^-3]'
      close(unit_out)
   end subroutine output_sources_header

   subroutine output_sources_step(filename, step_score, n_particles, n_init, dt, &
      istep, include_pretab)
      character(len=*), intent(in) :: filename
      type(score_data), intent(in) :: step_score
      integer, intent(in) :: n_particles, istep
      real(dp), intent(in) :: n_init, dt
      logical, intent(in) :: include_pretab

      integer :: unit_out
      real(dp) :: norm, time

      norm = n_init / (real(n_particles, dp) * dt)
      time = real(istep, dp) * dt

      unit_out = 50
      open(unit=unit_out, file=filename, status='old', position='append')
      call write_estimator_sources(unit_out, time, 'CL', step_score%cl, norm)
      call write_estimator_sources(unit_out, time, 'CL(avg)', step_score%cl_avg, norm)
      call write_estimator_sources(unit_out, time, 'TR', step_score%tr, norm)
      if (include_pretab) then
         call write_estimator_sources(unit_out, time, 'TR(pretab)', step_score%tr_pretab, norm)
      end if
      close(unit_out)
   end subroutine output_sources_step

   subroutine write_estimator_sources(unit_out, time, estimator_name, est, norm)
      integer, intent(in) :: unit_out
      real(dp), intent(in) :: time, norm
      character(len=*), intent(in) :: estimator_name
      type(estimator_score_t), intent(in) :: est
      integer :: ir

      do ir = 1, N_REACTIONS
         call write_source_row(unit_out, time, estimator_name, reaction_name(ir), &
            est%reaction(ir), norm)
      end do
   end subroutine write_estimator_sources

   subroutine write_source_row(unit_out, time, estimator_name, reaction, src, norm)
      integer, intent(in) :: unit_out
      real(dp), intent(in) :: time, norm
      character(len=*), intent(in) :: estimator_name, reaction
      type(source_terms_t), intent(in) :: src

      write(unit_out,'(ES14.6,2(",",A),8(",",ES14.6))') &
         time, trim(estimator_name), trim(reaction), &
         src%sn_plus * norm, src%sn_minus * norm, src%sn_net * norm, &
         src%sp(1) * norm, src%sp(2) * norm, src%sp(3) * norm, &
         src%we * norm, src%wi * norm
   end subroutine write_source_row

   subroutine output_ntscrg_header(filename)
      character(len=*), intent(in) :: filename
      integer :: unit_out

      unit_out = 51
      open(unit=unit_out, file=filename, status='replace')
      write(unit_out,'(A)') 'time[s],n_alive,weight_sum,' // &
         'CL_Wi_ei[W/m3],CL_Wi_cx[W/m3],CL_Wi_el[W/m3],CL_Wi_total[W/m3],' // &
         'CLavg_Wi_ei[W/m3],CLavg_Wi_cx[W/m3],CLavg_Wi_el[W/m3],' // &
         'CLavg_Wi_total[W/m3],TR_Wi_ei[W/m3],TR_Wi_cx[W/m3],' // &
         'TR_Wi_el[W/m3],TR_Wi_total[W/m3]'
      close(unit_out)
   end subroutine output_ntscrg_header

   subroutine output_ntscrg_step(filename, step_score, n_particles, n_init, dt, &
      istep, n_alive, weight_sum)
      character(len=*), intent(in) :: filename
      type(score_data), intent(in) :: step_score
      integer, intent(in) :: n_particles, istep, n_alive
      real(dp), intent(in) :: n_init, dt, weight_sum

      integer :: unit_out
      real(dp) :: norm, time

      norm = n_init / (real(n_particles, dp) * dt)
      time = real(istep, dp) * dt

      unit_out = 51
      open(unit=unit_out, file=filename, status='old', position='append')
      write(unit_out,'(ES14.6,A,I8,A,ES14.6,12(A,ES14.6))') &
         time, ',', n_alive, ',', weight_sum, &
         ',', step_score%cl%reaction(REACT_EI)%wi * norm, &
         ',', step_score%cl%reaction(REACT_CX)%wi * norm, &
         ',', step_score%cl%reaction(REACT_EL)%wi * norm, &
         ',', step_score%cl%reaction(REACT_TOTAL)%wi * norm, &
         ',', step_score%cl_avg%reaction(REACT_EI)%wi * norm, &
         ',', step_score%cl_avg%reaction(REACT_CX)%wi * norm, &
         ',', step_score%cl_avg%reaction(REACT_EL)%wi * norm, &
         ',', step_score%cl_avg%reaction(REACT_TOTAL)%wi * norm, &
         ',', step_score%tr%reaction(REACT_EI)%wi * norm, &
         ',', step_score%tr%reaction(REACT_CX)%wi * norm, &
         ',', step_score%tr%reaction(REACT_EL)%wi * norm, &
         ',', step_score%tr%reaction(REACT_TOTAL)%wi * norm
      close(unit_out)
   end subroutine output_ntscrg_step

   subroutine output_ntscrg_final(score, n_particles, n_init, n_steps, dt)
      type(score_data), intent(in) :: score
      integer, intent(in) :: n_particles, n_steps
      real(dp), intent(in) :: n_init, dt
      real(dp) :: norm

      norm = n_init / (real(n_particles, dp) * dt * real(n_steps, dp))

      write(*,'(A)') ''
      write(*,'(A)') '================================================================================'
      write(*,'(A)') ' Ion Energy Source Summary (time-averaged) [W/m3]'
      write(*,'(A)') '================================================================================'
      write(*,'(A12,4A16)') '', 'EI', 'CX', 'EL', 'Total'
      call write_console_wi('CL', score%cl, norm)
      call write_console_wi('CL(avg)', score%cl_avg, norm)
      call write_console_wi('TR', score%tr, norm)
      write(*,'(A)') '================================================================================'
   end subroutine output_ntscrg_final

   subroutine write_console_wi(name, est, norm)
      character(len=*), intent(in) :: name
      type(estimator_score_t), intent(in) :: est
      real(dp), intent(in) :: norm

      write(*,'(A12,4ES16.6)') trim(name), &
         est%reaction(REACT_EI)%wi * norm, &
         est%reaction(REACT_CX)%wi * norm, &
         est%reaction(REACT_EL)%wi * norm, &
         est%reaction(REACT_TOTAL)%wi * norm
   end subroutine write_console_wi

   subroutine output_source_stats(filename, stats, include_pretab)
      character(len=*), intent(in) :: filename
      type(score_stats_data), intent(in) :: stats
      logical, intent(in) :: include_pretab
      integer :: unit_out

      unit_out = 52
      open(unit=unit_out, file=filename, status='replace')
      write(unit_out,'(A)') 'estimator,reaction,quantity,mean,stddev,stderr,' // &
         'ci95_low,ci95_high,n_blocks'
      call write_estimator_stats(unit_out, 'CL', stats%cl)
      call write_estimator_stats(unit_out, 'CL(avg)', stats%cl_avg)
      call write_estimator_stats(unit_out, 'TR', stats%tr)
      if (include_pretab) call write_estimator_stats(unit_out, 'TR(pretab)', stats%tr_pretab)
      close(unit_out)
   end subroutine output_source_stats

   subroutine write_estimator_stats(unit_out, estimator_name, stats)
      integer, intent(in) :: unit_out
      character(len=*), intent(in) :: estimator_name
      type(estimator_stats_t), intent(in) :: stats
      integer :: ir

      do ir = 1, N_REACTIONS
         call write_source_stats(unit_out, estimator_name, reaction_name(ir), &
            stats%reaction(ir))
      end do
   end subroutine write_estimator_stats

   subroutine write_source_stats(unit_out, estimator_name, reaction, stats)
      integer, intent(in) :: unit_out
      character(len=*), intent(in) :: estimator_name, reaction
      type(source_stats_t), intent(in) :: stats

      call write_stat_row(unit_out, estimator_name, reaction, 'Sn_plus', stats%sn_plus)
      call write_stat_row(unit_out, estimator_name, reaction, 'Sn_minus', stats%sn_minus)
      call write_stat_row(unit_out, estimator_name, reaction, 'Sn_net', stats%sn_net)
      call write_stat_row(unit_out, estimator_name, reaction, 'Spx', stats%sp(1))
      call write_stat_row(unit_out, estimator_name, reaction, 'Spy', stats%sp(2))
      call write_stat_row(unit_out, estimator_name, reaction, 'Spz', stats%sp(3))
      call write_stat_row(unit_out, estimator_name, reaction, 'We', stats%we)
      call write_stat_row(unit_out, estimator_name, reaction, 'Wi', stats%wi)
   end subroutine write_source_stats

   subroutine write_stat_row(unit_out, estimator_name, reaction, quantity, stats)
      integer, intent(in) :: unit_out
      character(len=*), intent(in) :: estimator_name, reaction, quantity
      type(running_stats_t), intent(in) :: stats

      real(dp) :: stddev, stderr, ci95

      stddev = running_stddev(stats)
      stderr = running_stderr(stats)
      ci95 = 1.96d0 * stderr

      write(unit_out,'(3(A,","),5(ES14.6,","),I0)') &
         trim(estimator_name), trim(reaction), trim(quantity), &
         stats%mean, stddev, stderr, stats%mean - ci95, stats%mean + ci95, stats%n
   end subroutine write_stat_row

   subroutine output_balance(filename, initial_totals, final_totals, score, &
      n_particles, n_init)
      character(len=*), intent(in) :: filename
      type(system_totals_t), intent(in) :: initial_totals, final_totals
      type(score_data), intent(in) :: score
      integer, intent(in) :: n_particles
      real(dp), intent(in) :: n_init

      integer :: unit_out
      real(dp) :: norm

      norm = n_init / real(n_particles, dp)
      unit_out = 53
      open(unit=unit_out, file=filename, status='replace')
      write(unit_out,'(A)') 'estimator,quantity,initial,final,source_integral,' // &
         'residual,relative_residual'
      call write_balance_estimator(unit_out, 'CL', initial_totals, final_totals, &
         score%cl, norm)
      call write_balance_estimator(unit_out, 'CL(avg)', initial_totals, final_totals, &
         score%cl_avg, norm)
      call write_balance_estimator(unit_out, 'TR', initial_totals, final_totals, &
         score%tr, norm)
      close(unit_out)
   end subroutine output_balance

   subroutine write_balance_estimator(unit_out, estimator_name, initial_totals, &
      final_totals, est, norm)
      integer, intent(in) :: unit_out
      character(len=*), intent(in) :: estimator_name
      type(system_totals_t), intent(in) :: initial_totals, final_totals
      type(estimator_score_t), intent(in) :: est
      real(dp), intent(in) :: norm

      call write_balance_row(unit_out, estimator_name, 'particle', &
         initial_totals%particles * norm, final_totals%particles * norm, &
         est%reaction(REACT_TOTAL)%sn_net * norm)
      call write_balance_row(unit_out, estimator_name, 'momentum_x', &
         initial_totals%momentum(1) * norm, final_totals%momentum(1) * norm, &
         est%reaction(REACT_TOTAL)%sp(1) * norm)
      call write_balance_row(unit_out, estimator_name, 'momentum_y', &
         initial_totals%momentum(2) * norm, final_totals%momentum(2) * norm, &
         est%reaction(REACT_TOTAL)%sp(2) * norm)
      call write_balance_row(unit_out, estimator_name, 'momentum_z', &
         initial_totals%momentum(3) * norm, final_totals%momentum(3) * norm, &
         est%reaction(REACT_TOTAL)%sp(3) * norm)
      call write_balance_row(unit_out, estimator_name, 'energy_neutral_plus_ion', &
         initial_totals%energy * norm, final_totals%energy * norm, &
         est%reaction(REACT_TOTAL)%wi * norm)
      call write_balance_row(unit_out, estimator_name, 'energy_including_electron', &
         initial_totals%energy * norm, final_totals%energy * norm, &
         (est%reaction(REACT_TOTAL)%wi + est%reaction(REACT_TOTAL)%we) * norm)
   end subroutine write_balance_estimator

   subroutine write_balance_row(unit_out, estimator_name, quantity, initial, final, source)
      integer, intent(in) :: unit_out
      character(len=*), intent(in) :: estimator_name, quantity
      real(dp), intent(in) :: initial, final, source

      real(dp) :: residual, relative

      residual = final + source - initial
      if (abs(initial) > 1.0d-300) then
         relative = residual / initial
      else
         relative = 0.0d0
      end if

      write(unit_out,'(2(A,","),5(ES14.6,:,","))') trim(estimator_name), &
         trim(quantity), initial, final, source, residual, relative
   end subroutine write_balance_row

   subroutine output_energy_histogram(filename, particles, n_particles, &
      n_bins, E_min, E_max, step, n_timing, hist_timing, first_call)
      character(len=*), intent(in) :: filename
      type(particle_t), intent(in) :: particles(:)
      integer, intent(in) :: n_particles, n_bins, step, n_timing
      real(dp), intent(in) :: E_min, E_max
      integer, intent(in) :: hist_timing(:)
      logical, intent(in) :: first_call

      integer :: i, ibin, unit_out
      logical :: should_output
      real(dp) :: bin_width, E_particle, weight_sum
      real(dp), allocatable :: hist(:)

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
      weight_sum = 0.0d0
      bin_width = (E_max - E_min) / real(n_bins, dp)

      do i = 1, n_particles
         if (.not. particles(i)%alive) cycle
         E_particle = 0.5d0 * M_D_kg * &
            (particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2) * J_TO_EV
         ibin = int((E_particle - E_min) / bin_width) + 1
         if (ibin >= 1 .and. ibin <= n_bins) then
            hist(ibin) = hist(ibin) + particles(i)%weight
            weight_sum = weight_sum + particles(i)%weight
         end if
      end do

      unit_out = 54
      if (first_call) then
         open(unit=unit_out, file=filename, status='replace')
         write(unit_out,'(A)', advance='no') 'step,energy_center_eV,total_weight'
         do i = 1, n_bins
            write(unit_out,'(A,I0)', advance='no') ',bin_', i
         end do
         write(unit_out,*)
      else
         open(unit=unit_out, file=filename, status='old', position='append')
      end if

      write(unit_out,'(I8,A,ES14.6,A,ES14.6)', advance='no') step, ',', &
         0.5d0 * (E_min + E_max), ',', weight_sum
      do i = 1, n_bins
         write(unit_out,'(A,ES14.6)', advance='no') ',', hist(i)
      end do
      write(unit_out,*)
      close(unit_out)
      deallocate(hist)
   end subroutine output_energy_histogram

   subroutine output_deltaE_histogram(filename, dE_hist_el, dE_hist_cx, &
      n_bins, E_min, E_max, n_init, n_particles, collect_time)
      character(len=*), intent(in) :: filename
      real(dp), intent(in) :: dE_hist_el(:), dE_hist_cx(:)
      integer, intent(in) :: n_bins, n_particles
      real(dp), intent(in) :: E_min, E_max, n_init, collect_time

      integer :: i, unit_out
      real(dp) :: bin_width, center, norm

      if (collect_time <= 0.0d0) return
      bin_width = (E_max - E_min) / real(n_bins, dp)
      norm = n_init / (real(n_particles, dp) * collect_time)

      unit_out = 55
      open(unit=unit_out, file=filename, status='replace')
      write(unit_out,'(A)') 'deltaE_center_eV,EL_rate[m^-3 s^-1],CX_rate[m^-3 s^-1]'
      do i = 1, n_bins
         center = E_min + (real(i, dp) - 0.5d0) * bin_width
         write(unit_out,'(ES14.6,A,ES14.6,A,ES14.6)') center, ',', &
            dE_hist_el(i) * norm, ',', dE_hist_cx(i) * norm
      end do
      close(unit_out)
   end subroutine output_deltaE_histogram

   subroutine output_statistics(particles, n_particles)
      type(particle_t), intent(in) :: particles(:)
      integer, intent(in) :: n_particles

      integer :: i, n_alive
      real(dp) :: weight_sum, E_sum, E_sq, E_mean, E_std, v_sq, E_particle
      real(dp) :: vx_sum, vy_sum, vz_sum, vx_mean, vy_mean, vz_mean
      real(dp) :: E_bulk, E_thermal, T_equiv

      n_alive = 0
      weight_sum = 0.0d0
      E_sum = 0.0d0
      E_sq = 0.0d0
      vx_sum = 0.0d0
      vy_sum = 0.0d0
      vz_sum = 0.0d0

      do i = 1, n_particles
         if (.not. particles(i)%alive) cycle
         n_alive = n_alive + 1
         v_sq = particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2
         E_particle = 0.5d0 * M_D_kg * v_sq * J_TO_EV
         weight_sum = weight_sum + particles(i)%weight
         E_sum = E_sum + E_particle * particles(i)%weight
         E_sq = E_sq + E_particle**2 * particles(i)%weight
         vx_sum = vx_sum + particles(i)%vx * particles(i)%weight
         vy_sum = vy_sum + particles(i)%vy * particles(i)%weight
         vz_sum = vz_sum + particles(i)%vz * particles(i)%weight
      end do

      if (weight_sum > 0.0d0) then
         E_mean = E_sum / weight_sum
         E_std = sqrt(max(0.0d0, E_sq / weight_sum - E_mean**2))
         vx_mean = vx_sum / weight_sum
         vy_mean = vy_sum / weight_sum
         vz_mean = vz_sum / weight_sum
         E_bulk = 0.5d0 * M_D_kg * &
            (vx_mean**2 + vy_mean**2 + vz_mean**2) * J_TO_EV
         E_thermal = max(0.0d0, E_mean - E_bulk)
         T_equiv = (2.0d0 / 3.0d0) * E_thermal
      else
         E_mean = 0.0d0
         E_std = 0.0d0
         E_bulk = 0.0d0
         E_thermal = 0.0d0
         T_equiv = 0.0d0
      end if

      write(*,'(A)') ''
      write(*,'(A)') 'Final particle statistics:'
      write(*,'(A,I10)') '  Alive particles  = ', n_alive
      write(*,'(A,ES12.4)') '  Weight sum       = ', weight_sum
      write(*,'(A,F10.4,A)') '  Mean kinetic E   = ', E_mean, ' eV'
      write(*,'(A,F10.4,A)') '  Bulk kinetic E   = ', E_bulk, ' eV'
      write(*,'(A,F10.4,A)') '  Thermal energy   = ', E_thermal, ' eV'
      write(*,'(A,F10.4,A)') '  Neutral T_equiv  = ', T_equiv, ' eV'
      write(*,'(A,F10.4,A)') '  Energy std       = ', E_std, ' eV'
   end subroutine output_statistics

   pure function reaction_name(ir) result(name)
      integer, intent(in) :: ir
      character(len=8) :: name

      select case (ir)
      case (REACT_EI)
         name = 'EI'
      case (REACT_CX)
         name = 'CX'
      case (REACT_EL)
         name = 'EL'
      case (REACT_TOTAL)
         name = 'TOTAL'
      case default
         name = 'UNKNOWN'
      end select
   end function reaction_name

end module io
