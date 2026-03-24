!===============================================================================
! Module: io
! input and output: input file reading, particle initialization, output
!===============================================================================
module io
   use constants, only: dp, M_D_kg, EV_TO_J, J_TO_EV
   use data_types, only: particle_t, sim_params, plasma_params, init_params, diag_params, score_data
   use random_utils, only: sample_maxwell_velocity, set_beam_velocity, init_rng, random_double
   implicit none

   private
   public :: read_input_file, initialize_particles
   public :: output_energy_histogram, output_deltaE_histogram, output_statistics
   public :: output_ntscrg_header, output_ntscrg_step, output_ntscrg_final

contains

   !---------------------------------------------------------------------------
   ! read input.nml
   !---------------------------------------------------------------------------
   subroutine read_input_file(sim, plasma, init_p, diag)
      type(sim_params), intent(out)    :: sim
      type(plasma_params), intent(out) :: plasma
      type(init_params), intent(out)   :: init_p
      type(diag_params), intent(out)   :: diag

      integer, parameter :: unit_input = 20
      logical :: file_exist
      integer :: ios

      ! simulation namelist variables
      integer  :: n_particles, n_steps, seed
      real(dp) :: dt, weight_min
      logical  :: enable_cx, enable_el, enable_ei, use_isotropic
      character(len=256) :: cdf_file, output_ntscrg, output_hist, output_deltaE_hist

      ! plasma namelist variables
      real(dp) :: n_i, T_i, n_e, T_e, u_x, u_y, u_z

      ! particle_init namelist variables
      real(dp) :: n_init, E_init, T_init
      integer  :: init_mode

      ! diagnostics namelist variables
      integer  :: output_interval, n_hist_bins
      real(dp) :: E_hist_min, E_hist_max
      integer  :: hist_timing(6)
      integer  :: n_dE_bins
      real(dp) :: dE_hist_min, dE_hist_max
      integer  :: dE_collect_steps

      ! namelist
      namelist /simulation/ n_particles, n_steps, dt, seed, &
         enable_cx, enable_el, enable_ei, use_isotropic, &
         weight_min, cdf_file, output_ntscrg, output_hist, &
         output_deltaE_hist
      namelist /plasma_nml/ n_i, T_i, n_e, T_e, u_x, u_y, u_z
      namelist /particle_init/ n_init, E_init, T_init, init_mode
      namelist /diagnostics/ output_interval, n_hist_bins, &
         E_hist_min, E_hist_max, hist_timing, &
         n_dE_bins, dE_hist_min, dE_hist_max, dE_collect_steps

      ! default values
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

      ! store values in struct
      sim%n_particles = n_particles
      sim%n_steps = n_steps
      sim%dt = dt
      sim%seed = seed
      sim%enable_cx = enable_cx
      sim%enable_el = enable_el
      sim%enable_ei = enable_ei
      sim%use_isotropic = use_isotropic
      sim%weight_min = weight_min
      sim%cdf_file = cdf_file
      sim%output_ntscrg = output_ntscrg
      sim%output_hist = output_hist
      sim%output_deltaE_hist = output_deltaE_hist

      plasma%n_i = n_i
      plasma%T_i = T_i
      plasma%n_e = n_e
      plasma%T_e = T_e
      plasma%u_x = u_x
      plasma%u_y = u_y
      plasma%u_z = u_z

      init_p%n_init = n_init
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

      ! print parameters
      write(*,'(A)') '=========================================='
      write(*,'(A)') '3D3V Non-Analog Monte Carlo Simulation'
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
      if (init_mode == 0) then
         write(*,'(A,F8.2,A)')  '  E_init       = ', E_init, ' eV'
      else
         write(*,'(A,F8.2,A)')  '  T_init       = ', T_init, ' eV'
      end if
      write(*,'(A)')         '=========================================='
   end subroutine read_input_file

   subroutine initialize_particles(particles, n_particles, init_p, seed)
      type(particle_t), intent(out) :: particles(:)
      integer, intent(in)           :: n_particles
      type(init_params), intent(in) :: init_p
      integer, intent(in)           :: seed

      integer(8) :: i
      integer(8) :: particle_seed

      do i = 1, n_particles
         ! initialize particle seed
         particle_seed = int(seed, 8) * 6364136223846793005_8 + int(i, 8)
         call init_rng(particles(i)%rng, particle_seed)

         ! initialize particle position
         particles(i)%x = 0.0d0
         particles(i)%y = 0.0d0
         particles(i)%z = 0.0d0

         ! initialize particle velocity
         if (init_p%init_mode == 0) then
            call set_beam_velocity(init_p%E_init, particles(i)%vx, particles(i)%vy, particles(i)%vz)
         else
            call sample_maxwell_velocity(particles(i)%rng, init_p%T_init, particles(i)%vx, particles(i)%vy, particles(i)%vz)
         end if

         ! initialize particle weight
         particles(i)%weight = 1.0d0
      end do
   end subroutine initialize_particles

end module io
