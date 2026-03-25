!===============================================================================
! Program: monte_carlo_1d3v
! 1D 空間・3D 速度/位置追跡のモンテカルロ輸送
!===============================================================================
program monte_carlo_1d3v
   use constants, only: dp, J_TO_EV
   use data_types
   use random_utils, only: init_random_seed, random_double
   use cdf_reader, only: load_elastic_cdf
   use cross_sections, only: compute_nu_max, sigma_v_max
   use dynamics, only: advance_particle, evaluate_collision_event, collision_cx, &
      collision_el, update_weight, COLL_NONE, COLL_CX, COLL_EL
   use scoring, only: zero_score, add_score, init_profile_data, free_profile_data, &
      score_collision_estimator, flush_pending_track_scores, accumulate_collision_profile
   use io
   implicit none

   type(sim_params)    :: sim
   type(plasma_params) :: plasma
   type(init_params)   :: init_p
   type(diag_params)   :: diag
   type(score_data)    :: score
   type(score_data)    :: step_score
   type(score_data)    :: event_score
   type(profile_data)  :: profile
   type(particle_t), allocatable :: particles(:)

   real(dp), allocatable :: dE_hist_el(:), dE_hist_cx(:)

   integer :: istep, ip, ierr
   integer :: n_steps_done
   integer :: n_alive
   integer :: n_coll_cx, n_coll_el
   real(dp) :: weight_sum
   logical :: first_hist, collect_dE, use_delta_hist
   integer :: dE_start_step, dE_actual_steps
   real(dp) :: dE_collect_time

   integer :: coll_type_l
   real(dp) :: vx_i_l, vy_i_l, vz_i_l
   real(dp) :: v_rel_l, E_rel_l, delta_E_l
   type(particle_t) :: p_copy
   real(dp) :: time_remaining, step_dt, t_col, t_boundary, nu_max
   integer :: event_kind
   real(dp) :: r, dE_bin_width
   integer :: ibin_l

   call read_input_file(sim, plasma, init_p, diag)
   call init_random_seed(sim%seed)

   call load_elastic_cdf(sim%cdf_file, ierr)
   if (ierr /= 0) then
      write(*,*) 'Error: Failed to load CDF data'
      stop 1
   end if

   call compute_nu_max(plasma%n_i)

   allocate(particles(sim%n_particles))
   call initialize_particles(particles, sim, init_p)
   call init_profile_data(profile, sim%n_x_bins)
   call zero_score(score)

   call output_ntscrg_header(sim%output_ntscrg, sim%enable_tl_lookup)

   first_hist = .true.
   call output_energy_histogram(sim%output_hist, particles, sim%n_particles, &
      diag%n_hist_bins, diag%E_hist_min, diag%E_hist_max, 0, diag%hist_timing, first_hist)

   use_delta_hist = sim%enable_tl_lookup .and. diag%n_dE_bins > 0
   if (use_delta_hist) then
      allocate(dE_hist_el(diag%n_dE_bins), dE_hist_cx(diag%n_dE_bins))
      dE_hist_el = 0.0d0
      dE_hist_cx = 0.0d0
      dE_bin_width = (diag%dE_hist_max - diag%dE_hist_min) / dble(diag%n_dE_bins)
      if (diag%dE_collect_steps <= 0 .or. diag%dE_collect_steps >= sim%n_steps) then
         dE_start_step = 1
      else
         dE_start_step = sim%n_steps - diag%dE_collect_steps + 1
      end if
   else
      dE_start_step = sim%n_steps + 1
      call delete_if_exists(sim%output_deltaE_hist)
   end if

   write(*,'(A)') ''
   write(*,'(A)') 'Starting simulation...'

   n_steps_done = 0
   do istep = 1, sim%n_steps
      call zero_score(step_score)
      n_coll_cx = 0
      n_coll_el = 0
      collect_dE = use_delta_hist .and. (istep >= dE_start_step)

      do ip = 1, sim%n_particles
         if (.not. particles(ip)%alive) cycle

         particles(ip)%tl_pending_time = 0.0d0
         particles(ip)%tl_pending_x_start = particles(ip)%x
         particles(ip)%tl_pending_weight_start = particles(ip)%weight

         time_remaining = sim%dt
         do while (time_remaining > 1.0d-14 .and. particles(ip)%alive)
            nu_max = plasma%n_i * sigma_v_max

            if (nu_max > 1.0d-30) then
               t_col = max((particles(ip)%zincx - particles(ip)%zint1) / nu_max, 0.0d0)
            else
               t_col = huge(1.0d0)
            end if

            t_boundary = time_to_x_boundary(particles(ip), sim)
            step_dt = time_remaining
            event_kind = 0

            if (t_boundary < step_dt - 1.0d-14) then
               step_dt = t_boundary
               event_kind = 2
            end if
            if (t_col < step_dt - 1.0d-14) then
               step_dt = t_col
               event_kind = 1
            end if

            if (step_dt < 0.0d0) step_dt = 0.0d0

            if (particles(ip)%tl_pending_time <= 0.0d0) then
               particles(ip)%tl_pending_x_start = particles(ip)%x
               particles(ip)%tl_pending_weight_start = particles(ip)%weight
            end if
            particles(ip)%tl_pending_time = particles(ip)%tl_pending_time + step_dt

            if (nu_max > 1.0d-30) particles(ip)%zint1 = particles(ip)%zint1 + step_dt * nu_max

            if (sim%enable_ei) then
               call update_weight(particles(ip), plasma, step_dt, sim%weight_min)
            end if

            call advance_particle(particles(ip), step_dt)
            time_remaining = max(0.0d0, time_remaining - step_dt)

            if (.not. particles(ip)%alive) then
               call flush_pending_track_scores(particles(ip), plasma, sim, profile, step_score)
               exit
            end if

            if (event_kind == 2) then
               call flush_pending_track_scores(particles(ip), plasma, sim, profile, step_score)
               particles(ip)%alive = .false.
               exit
            else if (event_kind == 1) then
               call evaluate_collision_event(particles(ip), plasma, sim, &
                  vx_i_l, vy_i_l, vz_i_l, v_rel_l, E_rel_l, coll_type_l)

               if (coll_type_l /= COLL_NONE) then
                  call flush_pending_track_scores(particles(ip), plasma, sim, profile, step_score)
                  call zero_score(event_score)

                  if (coll_type_l == COLL_CX) then
                     delta_E_l = 0.0d0
                     call score_collision_estimator(particles(ip), plasma, &
                        vx_i_l, vy_i_l, vz_i_l, v_rel_l, E_rel_l, coll_type_l, &
                        delta_E_l, sim%enable_ei, event_score)
                     call add_score(step_score, event_score)
                     call accumulate_collision_profile(profile, sim, particles(ip)%x, &
                        event_score, coll_type_l, particles(ip)%weight)
                     call collision_cx(particles(ip), vx_i_l, vy_i_l, vz_i_l, delta_E_l)
                     n_coll_cx = n_coll_cx + 1
                  else if (coll_type_l == COLL_EL) then
                     p_copy = particles(ip)
                     call collision_el(p_copy, vx_i_l, vy_i_l, vz_i_l, sim%use_isotropic, delta_E_l)
                     call score_collision_estimator(particles(ip), plasma, &
                        vx_i_l, vy_i_l, vz_i_l, v_rel_l, E_rel_l, coll_type_l, &
                        delta_E_l, sim%enable_ei, event_score)
                     call add_score(step_score, event_score)
                     call accumulate_collision_profile(profile, sim, particles(ip)%x, &
                        event_score, coll_type_l, particles(ip)%weight)
                     particles(ip) = p_copy
                     n_coll_el = n_coll_el + 1
                  end if

                  if (collect_dE) then
                     ibin_l = int((delta_E_l * J_TO_EV - diag%dE_hist_min) / dE_bin_width) + 1
                     if (ibin_l >= 1 .and. ibin_l <= diag%n_dE_bins) then
                        if (coll_type_l == COLL_CX) then
                           dE_hist_cx(ibin_l) = dE_hist_cx(ibin_l) + particles(ip)%weight
                        else if (coll_type_l == COLL_EL) then
                           dE_hist_el(ibin_l) = dE_hist_el(ibin_l) + particles(ip)%weight
                        end if
                     end if
                  end if
               end if

               r = random_double(particles(ip)%rng)
               particles(ip)%zincx = -log(max(r, 1.0d-30))
               particles(ip)%zint1 = 0.0d0
            end if
         end do

         if (particles(ip)%alive .and. particles(ip)%tl_pending_time > 0.0d0) then
            call flush_pending_track_scores(particles(ip), plasma, sim, profile, step_score)
         end if
      end do

      call add_score(score, step_score)
      call count_alive_particles(particles, sim%n_particles, n_alive, weight_sum)
      call output_ntscrg_step(sim%output_ntscrg, step_score, sim, istep, n_alive, weight_sum)
      call output_energy_histogram(sim%output_hist, particles, sim%n_particles, &
         diag%n_hist_bins, diag%E_hist_min, diag%E_hist_max, istep, diag%hist_timing, first_hist)

      if (mod(istep, max(1, diag%output_interval)) == 0 .or. istep == sim%n_steps) then
         write(*,'(A,I6,A,I6,A,I8,A,I6,A,I6)') &
            '  Step ', istep, '/', sim%n_steps, '  alive=', n_alive, &
            '  CX=', n_coll_cx, '  EL=', n_coll_el
      end if

      n_steps_done = istep
      if (n_alive == 0) exit
   end do

   write(*,'(A)') ''
   write(*,'(A)') 'Simulation complete.'

   call output_ntscrg_final(score, sim, n_steps_done)
   call output_profile_x(sim%output_profile, profile, sim, n_steps_done)
   call output_statistics(particles, sim)

   if (use_delta_hist) then
      dE_actual_steps = max(1, n_steps_done - dE_start_step + 1)
      dE_collect_time = dble(dE_actual_steps) * sim%dt
      call output_deltaE_histogram(sim%output_deltaE_hist, dE_hist_el, dE_hist_cx, &
         diag%n_dE_bins, diag%dE_hist_min, diag%dE_hist_max, sim, dE_collect_time)
      deallocate(dE_hist_el, dE_hist_cx)
   end if

   call free_profile_data(profile)
   deallocate(particles)

contains

   real(dp) function time_to_x_boundary(p, sim) result(t)
      type(particle_t), intent(in) :: p
      type(sim_params), intent(in) :: sim

      if (abs(p%vx) < 1.0d-20) then
         t = huge(1.0d0)
      else if (p%vx > 0.0d0) then
         if (p%x >= sim%x_max) then
            t = 0.0d0
         else
            t = max((sim%x_max - p%x) / p%vx, 0.0d0)
         end if
      else
         if (p%x <= sim%x_min) then
            t = 0.0d0
         else
            t = max((sim%x_min - p%x) / p%vx, 0.0d0)
         end if
      end if
   end function time_to_x_boundary

   subroutine count_alive_particles(particles, n_particles, n_alive, weight_sum)
      type(particle_t), intent(in) :: particles(:)
      integer, intent(in) :: n_particles
      integer, intent(out) :: n_alive
      real(dp), intent(out) :: weight_sum

      integer :: i

      n_alive = 0
      weight_sum = 0.0d0
      do i = 1, n_particles
         if (.not. particles(i)%alive) cycle
         n_alive = n_alive + 1
         weight_sum = weight_sum + particles(i)%weight
      end do
   end subroutine count_alive_particles

   subroutine delete_if_exists(filename)
      character(len=*), intent(in) :: filename

      logical :: exists
      integer :: unit, ios

      inquire(file=trim(filename), exist=exists)
      if (.not. exists) return

      open(newunit=unit, file=trim(filename), status='old', action='readwrite', iostat=ios)
      if (ios /= 0) return
      close(unit, status='delete')
   end subroutine delete_if_exists

end program monte_carlo_1d3v
