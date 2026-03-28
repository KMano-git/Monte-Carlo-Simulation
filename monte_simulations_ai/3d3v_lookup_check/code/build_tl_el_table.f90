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
   use scoring, only: tl_el_sample_avg
   use tl_el_table_reader, only: build_default_tl_el_grids, write_tl_el_table
   implicit none

   type(sim_params) :: sim
   type(plasma_params) :: plasma
   type(init_params) :: init_p
   type(diag_params) :: diag
   type(particle_t) :: p
   real(dp), allocatable :: energy_grid(:), temp_grid(:)
   real(dp), allocatable :: rate_table(:,:)

   integer :: ierr, ios, i, j
   integer :: n_table_samples
   real(dp) :: ut
   character(len=256) :: output_file
   integer(8) :: seed_base
   logical :: file_exist

   namelist /tl_el_table_builder/ n_table_samples, output_file

   call read_input_file(sim, plasma, init_p, diag)

   n_table_samples = 2000
   output_file = trim(sim%tl_el_table_file)

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
   p%zincx = 0.0d0
   p%zint1 = 0.0d0
   p%tl_pending_time = 0.0d0
   p%tl_pending_eff_time = 0.0d0

   write(*,'(A,I8,A)') ' Building TL EL table with ', n_table_samples, ' samples/grid...'

   do j = 1, size(temp_grid)
      plasma%T_i = 2.0d0 * temp_grid(j)
      do i = 1, size(energy_grid)
         ut = sqrt(4.0d0 * energy_grid(i) * EV_TO_J / M_D_kg)
         p%vx = plasma%u_x + ut
         p%vy = plasma%u_y
         p%vz = plasma%u_z

         seed_base = int(sim%seed, 8) + 104729_8 * int(j - 1, 8) + int(i, 8)
         call init_rng(p%rng, seed_base)
         rate_table(i, j) = tl_el_sample_avg(p, plasma, sim%use_isotropic, n_table_samples)
      end do

      if (mod(j, 10) == 0 .or. j == size(temp_grid)) then
         write(*,'(A,I4,A,I4)') '   finished temp row ', j, '/', size(temp_grid)
      end if
   end do

   call write_tl_el_table(output_file, energy_grid, temp_grid, rate_table, ierr)
   if (ierr /= 0) then
      write(*,*) 'Error: Failed to write TL EL table'
      stop 1
   end if

   write(*,'(A,A)') ' TL EL table written to: ', trim(output_file)
end program build_tl_el_table
