!===============================================================================
! Program: build_tl_el_table
! Offline Monte Carlo で EL pre-tabulated TL table を生成する
!===============================================================================
program build_tl_el_table
   use constants, only: dp, M_D_kg, EV_TO_J
   use data_types, only: sim_params, plasma_params, init_params, diag_params, &
      particle_t
   use io, only: read_input_file
   use cdf_reader, only: load_elastic_cdf
   use random_utils, only: init_rng
   use scoring, only: tl_el_sample_stats
   use tl_el_table_reader, only: build_default_tl_el_grids, write_tl_el_table
   implicit none

   type(sim_params) :: sim
   type(plasma_params) :: plasma
   type(init_params) :: init_p
   type(diag_params) :: diag
   type(particle_t) :: p
   real(dp), allocatable :: energy_grid(:), temp_grid(:)
   real(dp), allocatable :: rate_table(:,:)

   integer :: ierr, ios, i, j, k, unit_stats
   integer :: n_table_samples
   real(dp) :: ut, variance_rate, stddev_rate, stderr_rate, rel_stderr
   real(dp) :: max_rel_stderr, max_stderr
   real(dp) :: rel_stderr_thresholds(3)
   integer :: n_rel_gt_threshold(3), n_total_cells
   character(len=256) :: output_file, stats_output_file
   integer(8) :: seed_base
   logical :: file_exist

   namelist /tl_el_table_builder/ n_table_samples, output_file, stats_output_file, &
      rel_stderr_thresholds

   call read_input_file(sim, plasma, init_p, diag)

   n_table_samples = 2000
   output_file = trim(sim%elastic_tl_table_file)
   stats_output_file = 'tl_el_table_stats.csv'
   rel_stderr_thresholds = (/0.01d0, 0.05d0, 0.10d0/)

   inquire(file='input.nml', exist=file_exist)
   if (file_exist) then
      open(unit=71, file='input.nml', status='old', iostat=ios)
      if (ios == 0) then
         read(71, nml=tl_el_table_builder, iostat=ios)
         close(71)
      end if
   end if

   call load_elastic_cdf(sim%cdf_file, ierr)
   if (ierr /= 0) then
      write(*,*) 'Error: Failed to load CDF data'
      stop 1
   end if

      call build_default_tl_el_grids(energy_grid, temp_grid)
      allocate(rate_table(size(energy_grid), size(temp_grid)))

      p%x = 0.0d0
   p%y = 0.0d0
   p%z = 0.0d0
   p%weight = 1.0d0
   p%alive = .true.
   p%collision_clock_target = 0.0d0
      p%collision_clock_elapsed = 0.0d0
      p%pending_track_time = 0.0d0
      p%pending_effective_track_time = 0.0d0

      write(*,'(A,I8,A)') ' Building TL EL table with ', n_table_samples, ' samples/grid...'

      unit_stats = 72
      open(unit=unit_stats, file=stats_output_file, status='replace', action='write', iostat=ios)
      if (ios /= 0) then
         write(*,*) 'Error: Failed to open stats output file'
         stop 1
      end if
      write(unit_stats,'(A)') 'energy_eV_amu,temp_eV_amu,n_samples,mean_rate_Wm3,stddev_rate_Wm3,' // &
         'variance_rate_Wm3_2,stderr_rate_Wm3,rel_stderr_abs_mean'

   max_rel_stderr = 0.0d0
   max_stderr = 0.0d0
   n_rel_gt_threshold = 0
   n_total_cells = size(energy_grid) * size(temp_grid)

      do j = 1, size(temp_grid)
         plasma%ion_temperature_eV = 2.0d0 * temp_grid(j)
         do i = 1, size(energy_grid)
            ! energy_grid is specific_energy [eV/amu] from the CDF metadata.
            ut = sqrt(8.0d0 * energy_grid(i) * EV_TO_J / M_D_kg)
         p%vx = plasma%ion_flow_vx + ut
         p%vy = plasma%ion_flow_vy
         p%vz = plasma%ion_flow_vz

            seed_base = int(sim%seed, 8) + 104729_8 * int(j - 1, 8) + int(i, 8)
            call init_rng(p%rng, seed_base)
            call tl_el_sample_stats(p, plasma, sim%use_isotropic, n_table_samples, &
               rate_table(i, j), variance_rate, stddev_rate, stderr_rate)

            if (abs(rate_table(i, j)) > 1.0d-30) then
               rel_stderr = stderr_rate / abs(rate_table(i, j))
            else if (stderr_rate > 0.0d0) then
               rel_stderr = huge(1.0d0)
            else
               rel_stderr = 0.0d0
            end if

            do k = 1, size(rel_stderr_thresholds)
               if (rel_stderr > rel_stderr_thresholds(k)) then
                  n_rel_gt_threshold(k) = n_rel_gt_threshold(k) + 1
               end if
            end do
            max_rel_stderr = max(max_rel_stderr, rel_stderr)
            max_stderr = max(max_stderr, stderr_rate)

            write(unit_stats,'(ES14.6,A,ES14.6,A,I0,4(A,ES14.6),A,ES14.6)') &
               energy_grid(i), ',', temp_grid(j), ',', n_table_samples, &
               ',', rate_table(i, j), ',', stddev_rate, ',', variance_rate, &
               ',', stderr_rate, ',', rel_stderr
         end do

         if (mod(j, 10) == 0 .or. j == size(temp_grid)) then
            write(*,'(A,I4,A,I4)') '   finished temp row ', j, '/', size(temp_grid)
         end if
      end do

      close(unit_stats)

      call write_tl_el_table(output_file, energy_grid, temp_grid, rate_table, ierr)
      if (ierr /= 0) then
         write(*,*) 'Error: Failed to write TL EL table'
         stop 1
      end if

      write(*,'(A,A)') ' TL EL table written to: ', trim(output_file)
      write(*,'(A,A)') ' TL EL stats written to: ', trim(stats_output_file)
      do i = 1, size(rel_stderr_thresholds)
         write(*,'(A,F6.2,A,I6,A,I6)') ' rel stderr > ', 100.0d0 * rel_stderr_thresholds(i), &
            '%: ', n_rel_gt_threshold(i), ' / ', n_total_cells
      end do
      write(*,'(A,ES12.4)') ' max rel stderr: ', max_rel_stderr
      write(*,'(A,ES12.4)') ' max stderr: ', max_stderr
end program build_tl_el_table
