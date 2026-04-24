!===============================================================================
! Program: monte_carlo_3d3v
! 3D3V non-analog Monte Carlo simulation with CL / CL(avg) / TR sources.
!===============================================================================
program monte_carlo_3d3v
   use constants
   use data_types
   use random_utils, only: init_random_seed, random_double
   use cdf_reader, only: load_elastic_cdf
   use tl_el_table_reader, only: load_tl_el_table
   use cross_sections, only: compute_nu_max, sigma_v_max, ionization_rate_coeff
   use dynamics, only: advance_particle, evaluate_collision_event, &
      collision_cx, collision_el, update_weight, COLL_NONE, COLL_CX, COLL_EL
   use scoring, only: score_collision_estimators, score_el_actual_collision, &
      score_track_length_estimator, score_pretabulated_track_length
   use source_terms, only: zero_score_data, add_score_data_active, &
      finalize_score_totals, compute_system_totals
   use statistics, only: reset_score_stats, update_score_stats
   use io
!$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
   implicit none

   type(sim_params) :: sim
   type(plasma_params) :: plasma
   type(init_params) :: init_p
   type(diag_params) :: diag
   type(score_data) :: total_score, step_score, particle_score, thread_score_local
   type(score_data), allocatable :: thread_scores(:)
   type(score_stats_data) :: score_stats
   type(system_totals_t) :: initial_totals, final_totals
   type(particle_t), allocatable :: particles(:)
   type(particle_t) :: p_before

   real(dp), allocatable :: dE_hist_el(:), dE_hist_cx(:)
   real(dp), allocatable :: thread_dE_hist_el(:,:), thread_dE_hist_cx(:,:)

   integer :: istep, ip, ierr, table_ierr
   integer :: coll_type, n_coll_cx, n_coll_el, n_coll_total
   integer :: n_alive, dE_start_step, dE_actual_steps
   integer :: count_start, count_end, count_rate
   integer :: n_threads_max, tid
   real(dp) :: vx_i, vy_i, vz_i, v_rel, E_rel, delta_E
   real(dp) :: time_remaining, step_dt, t_col, r, nu_max
   real(dp) :: ion_loss_rate, segment_eff_time, dE_bin_width
   real(dp) :: weight_sum, stats_norm, dE_collect_time, time_sec
   logical :: is_collision_event, collect_dE, first_hist

   call read_input_file(sim, plasma, init_p, diag)
   call init_random_seed(sim%seed)

   call load_elastic_cdf(sim%cdf_file, ierr)
   if (ierr /= 0) then
      write(*,*) 'Error: Failed to load CDF data'
      stop 1
   end if

   if (sim%compare_estimator) then
      call load_tl_el_table(sim%elastic_tl_table_file, table_ierr)
      if (table_ierr == 0) then
         write(*,'(A,A)') ' TL EL pretab table loaded: ', trim(sim%elastic_tl_table_file)
      else
         write(*,'(A,A)') ' TL EL pretab table not loaded: ', trim(sim%elastic_tl_table_file)
      end if
   end if

   call compute_nu_max(plasma%ion_density)
   nu_max = plasma%ion_density * sigma_v_max

   if (sim%enable_ionization) then
      ion_loss_rate = plasma%electron_density * &
         ionization_rate_coeff(plasma%electron_temperature_eV)
   else
      ion_loss_rate = 0.0d0
   end if

   allocate(particles(sim%n_particles))
   call initialize_particles(particles, sim%n_particles, init_p, sim%seed)
   call compute_system_totals(particles, initial_totals)

   call zero_score_data(total_score)
   call reset_score_stats(score_stats)
   call output_sources_header(sim%output_sources)
   call output_ntscrg_header(sim%output_ntscrg)

   allocate(dE_hist_el(diag%delta_energy_hist_bin_count), &
      dE_hist_cx(diag%delta_energy_hist_bin_count))
   dE_hist_el = 0.0d0
   dE_hist_cx = 0.0d0

   n_threads_max = 1
!$ n_threads_max = omp_get_max_threads()
   allocate(thread_scores(n_threads_max))
   allocate(thread_dE_hist_el(diag%delta_energy_hist_bin_count, n_threads_max), &
      thread_dE_hist_cx(diag%delta_energy_hist_bin_count, n_threads_max))

   dE_bin_width = (diag%delta_energy_hist_max_eV - diag%delta_energy_hist_min_eV) / &
      real(diag%delta_energy_hist_bin_count, dp)

   if (diag%delta_energy_collect_steps <= 0 .or. &
      diag%delta_energy_collect_steps >= sim%n_steps) then
      dE_start_step = 1
   else
      dE_start_step = sim%n_steps - diag%delta_energy_collect_steps + 1
   end if

   first_hist = .true.
   stats_norm = init_p%initial_density / (real(sim%n_particles, dp) * sim%dt)

   write(*,'(A)') ''
   write(*,'(A)') 'Starting simulation...'
!$ write(*,'(A,I0)') 'OpenMP max threads: ', n_threads_max
   call system_clock(count_start, count_rate)

   do istep = 1, sim%n_steps
      call zero_score_data(step_score)
      do tid = 1, n_threads_max
         call zero_score_data(thread_scores(tid))
      end do
      thread_dE_hist_el = 0.0d0
      thread_dE_hist_cx = 0.0d0

      collect_dE = (istep >= dE_start_step)
      n_coll_cx = 0
      n_coll_el = 0
      n_coll_total = 0
      n_alive = 0
      weight_sum = 0.0d0

!$omp parallel &
!$omp   private(ip, tid, particle_score, thread_score_local, p_before, coll_type, &
!$omp           vx_i, vy_i, vz_i, v_rel, E_rel, delta_E, time_remaining, &
!$omp           step_dt, t_col, r, segment_eff_time, is_collision_event)
      tid = 1
!$    tid = omp_get_thread_num() + 1
      call zero_score_data(thread_score_local)

!$omp do reduction(+:n_alive, weight_sum, n_coll_cx, n_coll_el, n_coll_total) &
!$omp    schedule(static)
      do ip = 1, sim%n_particles
         if (.not. particles(ip)%alive) cycle

         call zero_score_data(particle_score)
         time_remaining = sim%dt

         do while (time_remaining > 0.0d0 .and. particles(ip)%alive)
            if (nu_max > 1.0d-30) then
               t_col = max((particles(ip)%collision_clock_target - &
                  particles(ip)%collision_clock_elapsed) / nu_max, 0.0d0)
            else
               t_col = time_remaining * 2.0d0
            end if

            if (t_col <= time_remaining) then
               step_dt = t_col
               is_collision_event = .true.
            else
               step_dt = time_remaining
               is_collision_event = .false.
            end if

            particles(ip)%collision_clock_elapsed = particles(ip)%collision_clock_elapsed + &
               step_dt * nu_max

            if (step_dt > 0.0d0) then
               segment_eff_time = compute_effective_track_time( &
                  particles(ip)%weight, ion_loss_rate, step_dt)
               particles(ip)%pending_track_time = particles(ip)%pending_track_time + step_dt
               particles(ip)%pending_effective_track_time = &
                  particles(ip)%pending_effective_track_time + segment_eff_time
            end if

            if (sim%enable_ionization) then
               call update_weight(particles(ip), plasma, step_dt, sim%weight_min)
            end if

            call advance_particle(particles(ip), step_dt)
            time_remaining = max(0.0d0, time_remaining - step_dt)

            if (.not. particles(ip)%alive) then
               call flush_pending_track_scores(particles(ip), plasma, sim, particle_score)
               exit
            end if

            if (is_collision_event) then
               call evaluate_collision_event(particles(ip), plasma, sim, &
                  vx_i, vy_i, vz_i, v_rel, E_rel, coll_type)

               if (coll_type /= COLL_NONE) then
                  call flush_pending_track_scores(particles(ip), plasma, sim, particle_score)
                  n_coll_total = n_coll_total + 1

                  p_before = particles(ip)
                  call score_collision_estimators(p_before, plasma, vx_i, vy_i, vz_i, &
                     v_rel, E_rel, sim%enable_cx, sim%enable_el, &
                     sim%enable_ionization, particle_score)

                  if (coll_type == COLL_CX) then
                     call collision_cx(particles(ip), vx_i, vy_i, vz_i, delta_E)
                     n_coll_cx = n_coll_cx + 1
                     if (collect_dE) call add_deltaE_hist(thread_dE_hist_cx(:,tid), &
                        delta_E, particles(ip)%weight)
                  else if (coll_type == COLL_EL) then
                     call collision_el(particles(ip), vx_i, vy_i, vz_i, delta_E)
                     call score_el_actual_collision(p_before, particles(ip), particle_score)
                     n_coll_el = n_coll_el + 1
                     if (collect_dE) call add_deltaE_hist(thread_dE_hist_el(:,tid), &
                        delta_E, p_before%weight)
                  end if

                  r = random_double(particles(ip)%rng)
                  particles(ip)%collision_clock_target = -log(max(r, 1.0d-30))
                  particles(ip)%collision_clock_elapsed = 0.0d0
               else
                  r = random_double(particles(ip)%rng)
                  particles(ip)%collision_clock_target = -log(max(r, 1.0d-30))
                  particles(ip)%collision_clock_elapsed = 0.0d0
               end if
            end if
         end do

         if (particles(ip)%alive) then
            call flush_pending_track_scores(particles(ip), plasma, sim, particle_score)
            n_alive = n_alive + 1
            weight_sum = weight_sum + particles(ip)%weight
         end if

         call add_score_data_active(thread_score_local, particle_score, &
            sim%compare_estimator)
      end do
!$omp end do

      thread_scores(tid) = thread_score_local
!$omp end parallel

      do tid = 1, n_threads_max
         call add_score_data_active(step_score, thread_scores(tid), &
            sim%compare_estimator)
         dE_hist_el = dE_hist_el + thread_dE_hist_el(:,tid)
         dE_hist_cx = dE_hist_cx + thread_dE_hist_cx(:,tid)
      end do

      call finalize_score_totals(step_score, sim%compare_estimator)
      call add_score_data_active(total_score, step_score, sim%compare_estimator)
      call update_score_stats(score_stats, step_score, stats_norm, sim%compare_estimator)

      call output_sources_step(sim%output_sources, step_score, sim%n_particles, &
         init_p%initial_density, sim%dt, istep, sim%compare_estimator)
      call output_ntscrg_step(sim%output_ntscrg, step_score, sim%n_particles, &
         init_p%initial_density, sim%dt, istep, n_alive, weight_sum)

      if (mod(istep, diag%output_interval) == 0) then
         write(*,'(A,I6,A,I6,A,I8,A,I6,A,I6)') &
            '  Step ', istep, '/', sim%n_steps, &
            '  alive=', n_alive, '  CX=', n_coll_cx, '  EL=', n_coll_el
      end if

      call output_energy_histogram(sim%output_hist, particles, sim%n_particles, &
         diag%energy_hist_bin_count, diag%energy_hist_min_eV, diag%energy_hist_max_eV, &
         istep, 6, diag%hist_timing, first_hist)
      call update_first_hist_flag(first_hist, istep, diag%hist_timing)
   end do

   call system_clock(count_end)
   time_sec = real(count_end - count_start, dp) / real(count_rate, dp)
   write(*,'(A,F12.4,A)') 'Elapsed time: ', time_sec, ' seconds'
   write(*,'(A)') 'Simulation complete.'

   call output_ntscrg_final(total_score, sim%n_particles, init_p%initial_density, &
      sim%n_steps, sim%dt)
   call output_statistics(particles, sim%n_particles)
   call output_source_stats(sim%output_stats, score_stats, sim%compare_estimator)

   call compute_system_totals(particles, final_totals)
   call output_balance(sim%output_balance, initial_totals, final_totals, &
      total_score, sim%n_particles, init_p%initial_density)

   dE_actual_steps = sim%n_steps - dE_start_step + 1
   dE_collect_time = real(dE_actual_steps, dp) * sim%dt
   call output_deltaE_histogram(sim%output_delta_energy_hist, dE_hist_el, dE_hist_cx, &
      diag%delta_energy_hist_bin_count, diag%delta_energy_hist_min_eV, &
      diag%delta_energy_hist_max_eV, init_p%initial_density, sim%n_particles, &
      dE_collect_time)

   deallocate(particles)
   deallocate(dE_hist_el, dE_hist_cx)
   deallocate(thread_scores, thread_dE_hist_el, thread_dE_hist_cx)

contains

   pure function compute_effective_track_time(weight, loss_rate, dt) result(eff_time)
      real(dp), intent(in) :: weight, loss_rate, dt
      real(dp) :: eff_time

      if (loss_rate > 1.0d-30) then
         eff_time = (weight / loss_rate) * (1.0d0 - exp(-loss_rate * dt))
      else
         eff_time = weight * dt
      end if
   end function compute_effective_track_time

   subroutine flush_pending_track_scores(p, plasma, sim, local_score)
      type(particle_t), intent(inout) :: p
      type(plasma_params), intent(in) :: plasma
      type(sim_params), intent(in) :: sim
      type(score_data), intent(inout) :: local_score

      if (p%pending_effective_track_time <= 0.0d0) then
         p%pending_track_time = 0.0d0
         p%pending_effective_track_time = 0.0d0
         return
      end if

      call score_track_length_estimator(p, plasma, p%pending_effective_track_time, &
         sim%enable_cx, sim%enable_el, sim%enable_ionization, local_score)
      if (sim%compare_estimator) then
         call score_pretabulated_track_length(p, plasma, &
            p%pending_effective_track_time, sim%enable_el, local_score)
      end if

      p%pending_track_time = 0.0d0
      p%pending_effective_track_time = 0.0d0
   end subroutine flush_pending_track_scores

   subroutine add_deltaE_hist(hist, delta_E, weight)
      real(dp), intent(inout) :: hist(:)
      real(dp), intent(in) :: delta_E, weight
      integer :: ibin_local

      ibin_local = int((delta_E * J_TO_EV - diag%delta_energy_hist_min_eV) / &
         dE_bin_width) + 1
      if (ibin_local >= 1 .and. ibin_local <= size(hist)) then
         hist(ibin_local) = hist(ibin_local) + weight
      end if
   end subroutine add_deltaE_hist

   subroutine update_first_hist_flag(first_hist_flag, step, hist_timing)
      logical, intent(inout) :: first_hist_flag
      integer, intent(in) :: step
      integer, intent(in) :: hist_timing(:)
      integer :: i

      if (.not. first_hist_flag) return
      do i = 1, size(hist_timing)
         if (step == hist_timing(i)) then
            first_hist_flag = .false.
            return
         end if
      end do
   end subroutine update_first_hist_flag

end program monte_carlo_3d3v
