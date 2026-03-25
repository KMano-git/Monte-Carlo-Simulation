!===============================================================================
! Module: io
! 入出力、初期化、CSV 出力
!===============================================================================
module io
   use constants, only: dp, M_D_kg, J_TO_EV
   use data_types, only: particle_t, sim_params, plasma_params, init_params, &
      diag_params, score_data, profile_data
   use random_utils, only: sample_maxwell_velocity, set_beam_velocity, init_rng, random_double
   implicit none

   private
   public :: read_input_file, initialize_particles
   public :: output_energy_histogram, output_deltaE_histogram, output_statistics
   public :: output_ntscrg_header, output_ntscrg_step, output_ntscrg_final
   public :: output_profile_x

contains

   subroutine read_input_file(sim, plasma, init_p, diag)
      type(sim_params), intent(out)    :: sim
      type(plasma_params), intent(out) :: plasma
      type(init_params), intent(out)   :: init_p
      type(diag_params), intent(out)   :: diag

      integer, parameter :: unit_input = 20
      logical :: file_exist
      integer :: ios

      integer  :: n_particles, n_steps, n_x_bins, seed
      real(dp) :: dt, x_min, x_max, gamma_in, weight_min
      logical  :: enable_cx, enable_el, enable_ei, use_isotropic, enable_tl_lookup
      character(len=256) :: cdf_file, output_ntscrg, output_profile, output_hist, output_deltaE_hist

      real(dp) :: n_i, T_i, n_e, T_e, u_x, u_y, u_z

      real(dp) :: x_init, y_init, z_init, E_init, T_init
      integer  :: init_mode

      integer  :: output_interval, n_hist_bins
      real(dp) :: E_hist_min, E_hist_max
      integer  :: hist_timing(6)
      integer  :: n_dE_bins, dE_collect_steps
      real(dp) :: dE_hist_min, dE_hist_max

      namelist /simulation/ n_particles, n_steps, n_x_bins, dt, seed, &
         x_min, x_max, gamma_in, enable_cx, enable_el, enable_ei, use_isotropic, &
         enable_tl_lookup, weight_min, cdf_file, output_ntscrg, output_profile, &
         output_hist, output_deltaE_hist
      namelist /plasma_nml/ n_i, T_i, n_e, T_e, u_x, u_y, u_z
      namelist /particle_init/ x_init, y_init, z_init, E_init, T_init, init_mode
      namelist /diagnostics/ output_interval, n_hist_bins, E_hist_min, E_hist_max, &
         hist_timing, n_dE_bins, dE_hist_min, dE_hist_max, dE_collect_steps

      n_particles = 4000
      n_steps = 200
      n_x_bins = 100
      dt = 1.0d-8
      seed = 12345
      x_min = 0.0d0
      x_max = 0.1d0
      gamma_in = 1.0d23
      enable_cx = .true.
      enable_el = .true.
      enable_ei = .true.
      use_isotropic = .false.
      enable_tl_lookup = .false.
      weight_min = 1.0d-10
      cdf_file = 'dd_00_elastic.cdf'
      output_ntscrg = 'ntscrg.csv'
      output_profile = 'profile_x.csv'
      output_hist = 'energy_hist.csv'
      output_deltaE_hist = 'deltaE_hist.csv'

      n_i = 1.0d21
      T_i = 10.0d0
      n_e = 1.0d21
      T_e = 10.0d0
      u_x = 0.0d0
      u_y = 0.0d0
      u_z = 0.0d0

      x_init = 0.0d0
      y_init = 0.0d0
      z_init = 0.0d0
      E_init = 15.0d0
      T_init = 2.0d0
      init_mode = 0

      output_interval = 20
      n_hist_bins = 200
      E_hist_min = 0.0d0
      E_hist_max = 150.0d0
      hist_timing = (/0, 20, 50, 100, 150, 200/)
      n_dE_bins = 200
      dE_hist_min = -20.0d0
      dE_hist_max = 20.0d0
      dE_collect_steps = 0

      inquire(file='input.nml', exist=file_exist)
      if (file_exist) then
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
      else
         write(*,*) 'Warning: input.nml not found, using defaults'
      end if

      sim%n_particles = n_particles
      sim%n_steps = n_steps
      sim%n_x_bins = n_x_bins
      sim%seed = seed
      sim%dt = dt
      sim%x_min = x_min
      sim%x_max = x_max
      sim%dx = (x_max - x_min) / dble(max(1, n_x_bins))
      sim%gamma_in = gamma_in
      sim%enable_cx = enable_cx
      sim%enable_el = enable_el
      sim%enable_ei = enable_ei
      sim%use_isotropic = use_isotropic
      sim%enable_tl_lookup = enable_tl_lookup
      sim%weight_min = weight_min
      sim%cdf_file = cdf_file
      sim%output_ntscrg = output_ntscrg
      sim%output_profile = output_profile
      sim%output_hist = output_hist
      sim%output_deltaE_hist = output_deltaE_hist

      plasma%n_i = n_i
      plasma%T_i = T_i
      plasma%n_e = n_e
      plasma%T_e = T_e
      plasma%u_x = u_x
      plasma%u_y = u_y
      plasma%u_z = u_z

      init_p%x_init = x_init
      init_p%y_init = y_init
      init_p%z_init = z_init
      init_p%E_init = E_init
      init_p%T_init = T_init
      init_p%init_mode = init_mode

      diag%output_interval = output_interval
      diag%n_hist_bins = n_hist_bins
      diag%E_hist_min = E_hist_min
      diag%E_hist_max = E_hist_max
      diag%hist_timing = hist_timing
      diag%n_dE_bins = n_dE_bins
      diag%dE_hist_min = dE_hist_min
      diag%dE_hist_max = dE_hist_max
      diag%dE_collect_steps = dE_collect_steps

      write(*,'(A)') '=========================================='
      write(*,'(A)') ' 1D3V Monte Carlo Simulation'
      write(*,'(A)') '=========================================='
      write(*,'(A,I10)')     '  n_particles  = ', sim%n_particles
      write(*,'(A,I10)')     '  n_steps      = ', sim%n_steps
      write(*,'(A,I10)')     '  n_x_bins     = ', sim%n_x_bins
      write(*,'(A,ES12.4)')  '  dt           = ', sim%dt
      write(*,'(A,2ES12.4)') '  x-range      = ', sim%x_min, sim%x_max
      write(*,'(A,ES12.4)')  '  gamma_in     = ', sim%gamma_in
      write(*,'(A,L1)')      '  enable_cx    = ', sim%enable_cx
      write(*,'(A,L1)')      '  enable_el    = ', sim%enable_el
      write(*,'(A,L1)')      '  enable_ei    = ', sim%enable_ei
      write(*,'(A,L1)')      '  enable_tllkp = ', sim%enable_tl_lookup
      write(*,'(A,F8.2,A)')  '  T_i          = ', plasma%T_i, ' eV'
      write(*,'(A,F8.2,A)')  '  T_e          = ', plasma%T_e, ' eV'
      if (init_p%init_mode == 0) then
         write(*,'(A,F8.2,A)') '  E_init       = ', init_p%E_init, ' eV'
      else
         write(*,'(A,F8.2,A)') '  T_init       = ', init_p%T_init, ' eV'
      end if
      write(*,'(A)') '=========================================='
   end subroutine read_input_file

   subroutine initialize_particles(particles, sim, init_p)
      type(particle_t), intent(out) :: particles(:)
      type(sim_params), intent(in) :: sim
      type(init_params), intent(in) :: init_p

      integer :: i
      integer(8) :: particle_seed

      do i = 1, sim%n_particles
         particle_seed = int(sim%seed, 8) * 6364136223846793005_8 + int(i, 8)
         call init_rng(particles(i)%rng, particle_seed)

         particles(i)%x = init_p%x_init
         particles(i)%y = init_p%y_init
         particles(i)%z = init_p%z_init
         particles(i)%weight = 1.0d0
         particles(i)%alive = .true.
         particles(i)%zint1 = 0.0d0
         particles(i)%zincx = -log(max(random_double(particles(i)%rng), 1.0d-30))
         particles(i)%tl_pending_time = 0.0d0
         particles(i)%tl_pending_x_start = particles(i)%x
         particles(i)%tl_pending_weight_start = particles(i)%weight

         if (init_p%init_mode == 0) then
            call set_beam_velocity(init_p%E_init, particles(i)%vx, particles(i)%vy, particles(i)%vz)
         else
            call sample_maxwell_velocity(particles(i)%rng, init_p%T_init, &
               particles(i)%vx, particles(i)%vy, particles(i)%vz)
         end if
      end do
   end subroutine initialize_particles

   subroutine output_energy_histogram(filename, particles, n_particles, n_bins, E_min, E_max, &
      step, hist_timing, first_call)
      character(len=*), intent(in) :: filename
      type(particle_t), intent(in) :: particles(:)
      integer, intent(in)  :: n_particles, n_bins, step
      real(dp), intent(in) :: E_min, E_max
      integer, intent(in)  :: hist_timing(:)
      logical, intent(inout) :: first_call

      real(dp) :: bin_width, E_particle
      real(dp), allocatable :: hist(:)
      integer :: i, ibin, unit_out
      logical :: should_output

      should_output = (step == 0)
      if (.not. should_output) then
         do i = 1, size(hist_timing)
            if (step == hist_timing(i)) then
               should_output = .true.
               exit
            end if
         end do
      end if
      if (.not. should_output) return

      allocate(hist(n_bins))
      hist = 0.0d0
      bin_width = (E_max - E_min) / dble(n_bins)

      do i = 1, n_particles
         if (.not. particles(i)%alive) cycle
         E_particle = 0.5d0 * M_D_kg * &
            (particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2) * J_TO_EV
         ibin = int((E_particle - E_min) / bin_width) + 1
         if (ibin >= 1 .and. ibin <= n_bins) hist(ibin) = hist(ibin) + particles(i)%weight
      end do

      unit_out = 40
      if (first_call) then
         open(unit=unit_out, file=filename, status='replace')
         write(unit_out,'(A)', advance='no') 'step,E_min,E_max,bin_width,n_bins'
         do i = 1, n_bins
            write(unit_out,'(A,I0)', advance='no') ',bin_', i
         end do
         write(unit_out,*)
         first_call = .false.
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
   end subroutine output_energy_histogram

   subroutine output_deltaE_histogram(filename, dE_hist_el, dE_hist_cx, n_bins, dE_min, dE_max, &
      sim, collect_time)
      character(len=*), intent(in) :: filename
      real(dp), intent(in) :: dE_hist_el(:), dE_hist_cx(:)
      integer, intent(in)  :: n_bins
      real(dp), intent(in) :: dE_min, dE_max
      type(sim_params), intent(in) :: sim
      real(dp), intent(in) :: collect_time

      real(dp) :: bin_width, norm
      integer :: i, unit_out

      if (collect_time <= 0.0d0) return

      bin_width = (dE_max - dE_min) / dble(n_bins)
      norm = sim%gamma_in / (dble(sim%n_particles) * collect_time)

      unit_out = 41
      open(unit=unit_out, file=filename, status='replace')
      write(unit_out,'(A)') 'bin_center_eV,CX_rate[m-2s-1],EL_rate[m-2s-1]'
      do i = 1, n_bins
         write(unit_out,'(ES14.6,A,ES14.6,A,ES14.6)') &
            dE_min + (dble(i) - 0.5d0) * bin_width, ',', &
            dE_hist_cx(i) * norm, ',', &
            dE_hist_el(i) * norm
      end do
      close(unit_out)
   end subroutine output_deltaE_histogram

   subroutine output_ntscrg_header(filename, enable_tl_lookup)
      character(len=*), intent(in) :: filename
      logical, intent(in) :: enable_tl_lookup
      integer :: unit_out

      unit_out = 50
      open(unit=unit_out, file=filename, status='replace')
      if (enable_tl_lookup) then
         write(unit_out,'(A)') 'time[s],n_alive,weight_sum,' // &
            'CL_Q_cx[W/m3],CL_Q_el[W/m3],CL_Q_ei[W/m3],CL_Q_total[W/m3],' // &
            'TL_Q_cx[W/m3],TL_Q_el[W/m3],TL_Q_ei[W/m3],TL_Q_total[W/m3],' // &
            'TLLu_Q_cx[W/m3],TLLu_Q_el[W/m3],TLLu_Q_ei[W/m3],TLLu_Q_total[W/m3]'
      else
         write(unit_out,'(A)') 'time[s],n_alive,weight_sum,' // &
            'CL_Q_cx[W/m3],CL_Q_el[W/m3],CL_Q_ei[W/m3],CL_Q_total[W/m3],' // &
            'TL_Q_cx[W/m3],TL_Q_el[W/m3],TL_Q_ei[W/m3],TL_Q_total[W/m3]'
      end if
      close(unit_out)
   end subroutine output_ntscrg_header

   subroutine output_ntscrg_step(filename, step_score, sim, istep, n_alive, weight_sum)
      character(len=*), intent(in) :: filename
      type(score_data), intent(in) :: step_score
      type(sim_params), intent(in) :: sim
      integer, intent(in)  :: istep, n_alive
      real(dp), intent(in) :: weight_sum

      integer :: unit_out
      real(dp) :: norm, time, domain_L
      real(dp) :: cl_cx, cl_el, cl_ei, cl_total
      real(dp) :: tl_cx, tl_el, tl_ei, tl_total
      real(dp) :: tl_lookup_cx, tl_lookup_el, tl_lookup_ei, tl_lookup_total

      domain_L = max(sim%x_max - sim%x_min, sim%dx)
      norm = sim%gamma_in / (dble(sim%n_particles) * sim%dt * domain_L)
      time = dble(istep) * sim%dt

      cl_cx = step_score%cl_cx * norm
      cl_el = step_score%cl_el * norm
      cl_ei = step_score%cl_ei * norm
      cl_total = cl_cx + cl_el + cl_ei

      tl_cx = step_score%tl_cx * norm
      tl_el = step_score%tl_el * norm
      tl_ei = step_score%tl_ei * norm
      tl_total = tl_cx + tl_el + tl_ei

      tl_lookup_cx = step_score%tl_lookup_cx * norm
      tl_lookup_el = step_score%tl_lookup_el * norm
      tl_lookup_ei = step_score%tl_lookup_ei * norm
      tl_lookup_total = tl_lookup_cx + tl_lookup_el + tl_lookup_ei

      unit_out = 50
      open(unit=unit_out, file=filename, status='old', position='append')
      if (sim%enable_tl_lookup) then
         write(unit_out,'(ES14.6,A,I8,A,ES14.6,12(A,ES14.6))') &
            time, ',', n_alive, ',', weight_sum, &
            ',', cl_cx, ',', cl_el, ',', cl_ei, ',', cl_total, &
            ',', tl_cx, ',', tl_el, ',', tl_ei, ',', tl_total, &
            ',', tl_lookup_cx, ',', tl_lookup_el, ',', tl_lookup_ei, ',', tl_lookup_total
      else
         write(unit_out,'(ES14.6,A,I8,A,ES14.6,8(A,ES14.6))') &
            time, ',', n_alive, ',', weight_sum, &
            ',', cl_cx, ',', cl_el, ',', cl_ei, ',', cl_total, &
            ',', tl_cx, ',', tl_el, ',', tl_ei, ',', tl_total
      end if
      close(unit_out)
   end subroutine output_ntscrg_step

   subroutine output_ntscrg_final(score, sim, n_steps_done)
      type(score_data), intent(in) :: score
      type(sim_params), intent(in) :: sim
      integer, intent(in)  :: n_steps_done

      real(dp) :: norm, domain_L, total_time
      real(dp) :: cl_total, tl_total, tl_lookup_total

      total_time = max(1, n_steps_done) * sim%dt
      domain_L = max(sim%x_max - sim%x_min, sim%dx)
      norm = sim%gamma_in / (dble(sim%n_particles) * total_time * domain_L)

      cl_total = (score%cl_ei + score%cl_cx + score%cl_el) * norm
      tl_total = (score%tl_ei + score%tl_cx + score%tl_el) * norm
      tl_lookup_total = (score%tl_lookup_ei + score%tl_lookup_cx + score%tl_lookup_el) * norm

      write(*,'(A)') ''
      write(*,'(A)') '=========================================='
      write(*,'(A)') ' Time-Averaged Global Source [W/m3]'
      write(*,'(A)') '=========================================='
      write(*,'(A)')         '          EI              CX              EL              Total'
      write(*,'(A,4ES16.6)') ' CL: ', &
         score%cl_ei * norm, score%cl_cx * norm, score%cl_el * norm, cl_total
      write(*,'(A,4ES16.6)') ' TL: ', &
         score%tl_ei * norm, score%tl_cx * norm, score%tl_el * norm, tl_total
      if (sim%enable_tl_lookup) then
         write(*,'(A,4ES16.6)') ' TLLookup: ', &
            score%tl_lookup_ei * norm, score%tl_lookup_cx * norm, &
            score%tl_lookup_el * norm, tl_lookup_total
      end if
      write(*,'(A)') '=========================================='
   end subroutine output_ntscrg_final

   subroutine output_profile_x(filename, profile, sim, n_steps_done)
      character(len=*), intent(in) :: filename
      type(profile_data), intent(in) :: profile
      type(sim_params), intent(in) :: sim
      integer, intent(in) :: n_steps_done

      integer :: i, unit_out
      real(dp) :: total_time, x_center
      real(dp) :: norm_power, norm_density, norm_count
      real(dp) :: cl_total, tl_total, tllu_total

      total_time = max(1, n_steps_done) * sim%dt
      norm_power = sim%gamma_in / (dble(sim%n_particles) * total_time * sim%dx)
      norm_density = sim%gamma_in / (dble(sim%n_particles) * sim%dx)
      norm_count = 1.0d0 / dble(sim%n_particles)

      unit_out = 51
      open(unit=unit_out, file=filename, status='replace')
      if (sim%enable_tl_lookup) then
         write(unit_out,'(A)') 'x_center[m],residence_eff_per_history[s],neutral_density[m-3],' // &
            'CL_Q_cx[W/m3],CL_Q_el[W/m3],CL_Q_ei[W/m3],CL_Q_total[W/m3],' // &
            'TL_Q_cx[W/m3],TL_Q_el[W/m3],TL_Q_ei[W/m3],TL_Q_total[W/m3],' // &
            'TLLu_Q_cx[W/m3],TLLu_Q_el[W/m3],TLLu_Q_ei[W/m3],TLLu_Q_total[W/m3],' // &
            'collision_count_cx_per_history[-],collision_count_el_per_history[-],' // &
            'collision_count_total_per_history[-]'
      else
         write(unit_out,'(A)') 'x_center[m],residence_eff_per_history[s],neutral_density[m-3],' // &
            'CL_Q_cx[W/m3],CL_Q_el[W/m3],CL_Q_ei[W/m3],CL_Q_total[W/m3],' // &
            'TL_Q_cx[W/m3],TL_Q_el[W/m3],TL_Q_ei[W/m3],TL_Q_total[W/m3],' // &
            'collision_count_cx_per_history[-],collision_count_el_per_history[-],' // &
            'collision_count_total_per_history[-]'
      end if

      do i = 1, sim%n_x_bins
         x_center = sim%x_min + (dble(i) - 0.5d0) * sim%dx
         cl_total = (profile%cl_cx(i) + profile%cl_el(i) + profile%cl_ei(i)) * norm_power
         tl_total = (profile%tl_cx(i) + profile%tl_el(i) + profile%tl_ei(i)) * norm_power
         tllu_total = (profile%tl_lookup_cx(i) + profile%tl_lookup_el(i) + &
            profile%tl_lookup_ei(i)) * norm_power

         if (sim%enable_tl_lookup) then
            write(unit_out,'(ES14.6,17(A,ES14.6))') &
               x_center, ',', profile%residence_eff(i) / dble(sim%n_particles), ',', &
               profile%residence_eff(i) * norm_density, ',', &
               profile%cl_cx(i) * norm_power, ',', profile%cl_el(i) * norm_power, ',', &
               profile%cl_ei(i) * norm_power, ',', cl_total, ',', &
               profile%tl_cx(i) * norm_power, ',', profile%tl_el(i) * norm_power, ',', &
               profile%tl_ei(i) * norm_power, ',', tl_total, ',', &
               profile%tl_lookup_cx(i) * norm_power, ',', profile%tl_lookup_el(i) * norm_power, ',', &
               profile%tl_lookup_ei(i) * norm_power, ',', tllu_total, ',', &
               profile%collision_count_cx(i) * norm_count, ',', &
               profile%collision_count_el(i) * norm_count, ',', &
               profile%collision_count_total(i) * norm_count
         else
            write(unit_out,'(ES14.6,13(A,ES14.6))') &
               x_center, ',', profile%residence_eff(i) / dble(sim%n_particles), ',', &
               profile%residence_eff(i) * norm_density, ',', &
               profile%cl_cx(i) * norm_power, ',', profile%cl_el(i) * norm_power, ',', &
               profile%cl_ei(i) * norm_power, ',', cl_total, ',', &
               profile%tl_cx(i) * norm_power, ',', profile%tl_el(i) * norm_power, ',', &
               profile%tl_ei(i) * norm_power, ',', tl_total, ',', &
               profile%collision_count_cx(i) * norm_count, ',', &
               profile%collision_count_el(i) * norm_count, ',', &
               profile%collision_count_total(i) * norm_count
         end if
      end do

      close(unit_out)
   end subroutine output_profile_x

   subroutine output_statistics(particles, sim)
      type(particle_t), intent(in) :: particles(:)
      type(sim_params), intent(in) :: sim

      integer :: i, n_alive
      real(dp) :: E_sum, E_sq, E_mean, E_var, v_sq, weight_sum

      n_alive = 0
      E_sum = 0.0d0
      E_sq = 0.0d0
      weight_sum = 0.0d0

      do i = 1, sim%n_particles
         if (.not. particles(i)%alive) cycle
         n_alive = n_alive + 1
         v_sq = particles(i)%vx**2 + particles(i)%vy**2 + particles(i)%vz**2
         E_sum = E_sum + 0.5d0 * M_D_kg * v_sq * J_TO_EV * particles(i)%weight
         E_sq  = E_sq + (0.5d0 * M_D_kg * v_sq * J_TO_EV)**2 * particles(i)%weight
         weight_sum = weight_sum + particles(i)%weight
      end do

      write(*,'(A)') ''
      write(*,'(A)') '=========================================='
      write(*,'(A)') ' Final Statistics'
      write(*,'(A)') '=========================================='
      write(*,'(A,I10)')   '  Alive particles  = ', n_alive
      write(*,'(A,ES12.4)') '  Total weight     = ', weight_sum
      if (weight_sum > 0.0d0) then
         E_mean = E_sum / weight_sum
         E_var = max(E_sq / weight_sum - E_mean**2, 0.0d0)
         write(*,'(A,F12.4,A)') '  Mean energy      = ', E_mean, ' eV'
         write(*,'(A,F12.4,A)') '  Eff. temperature = ', E_mean * 2.0d0 / 3.0d0, ' eV'
         write(*,'(A,F12.4,A)') '  Energy std       = ', sqrt(E_var), ' eV'
      end if
      write(*,'(A)') '=========================================='
   end subroutine output_statistics

end module io
