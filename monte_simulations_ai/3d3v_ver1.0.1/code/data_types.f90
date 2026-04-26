!===============================================================================
! Module: data_types
! Core data definitions for particles, inputs, source scores, and statistics.
!===============================================================================
module data_types
   use constants, only : dp
   implicit none

   integer, parameter :: N_REACTIONS = 4
   integer, parameter :: REACT_EI = 1
   integer, parameter :: REACT_CX = 2
   integer, parameter :: REACT_EL = 3
   integer, parameter :: REACT_TOTAL = 4

   !---------------------------------------------------------------------------
   ! RNG state for xoshiro256+
   !---------------------------------------------------------------------------
   type :: rng_state
      integer(8) :: s(4)
   end type rng_state

   !---------------------------------------------------------------------------
   ! Particle data
   !---------------------------------------------------------------------------
   type :: particle_t
      real(dp) :: x, y, z
      real(dp) :: vx, vy, vz
      real(dp) :: weight
      logical  :: alive
      real(dp) :: collision_clock_target
      real(dp) :: collision_clock_elapsed
      real(dp) :: pending_track_time
      real(dp) :: pending_effective_track_time
      type(rng_state) :: rng
   end type particle_t

   !---------------------------------------------------------------------------
   ! Simulation parameters
   !---------------------------------------------------------------------------
   type :: sim_params
      integer  :: n_particles
      integer  :: n_steps
      real(dp) :: dt
      integer  :: seed
      logical  :: enable_cx
      logical  :: enable_el
      logical  :: enable_ionization
      logical  :: compare_estimator
      real(dp) :: weight_min
      character(len=256) :: cdf_file
      character(len=256) :: elastic_tl_table_file
      character(len=256) :: output_sources
      character(len=256) :: output_stats
      character(len=256) :: output_balance
      character(len=256) :: output_ntscrg
      character(len=256) :: output_hist
      character(len=256) :: output_delta_energy_hist
   end type sim_params

   !---------------------------------------------------------------------------
   ! Plasma parameters
   !---------------------------------------------------------------------------
   type :: plasma_params
      real(dp) :: ion_density
      real(dp) :: ion_temperature_eV
      real(dp) :: electron_density
      real(dp) :: electron_temperature_eV
      real(dp) :: ion_flow_vx, ion_flow_vy, ion_flow_vz
   end type plasma_params

   !---------------------------------------------------------------------------
   ! Initial neutral distribution
   !---------------------------------------------------------------------------
   type :: init_params
      real(dp) :: initial_density
      real(dp) :: initial_energy_eV
      real(dp) :: initial_temperature_eV
      integer  :: initial_velocity_mode
   end type init_params

   !---------------------------------------------------------------------------
   ! Diagnostics
   !---------------------------------------------------------------------------
   type :: diag_params
      integer  :: output_interval
      integer  :: energy_hist_bin_count
      real(dp) :: energy_hist_min_eV
      real(dp) :: energy_hist_max_eV
      integer  :: hist_timing(6)
      integer  :: delta_energy_hist_bin_count
      real(dp) :: delta_energy_hist_min_eV
      real(dp) :: delta_energy_hist_max_eV
      integer  :: delta_energy_collect_steps
   end type diag_params

   !---------------------------------------------------------------------------
   ! Plasma-side source terms.
   ! sn_* are particle weights. sp is momentum [kg m/s]. we/wi are energy [J].
   !---------------------------------------------------------------------------
   type :: source_terms_t
      real(dp) :: sn_plus
      real(dp) :: sn_minus
      real(dp) :: sn_net
      real(dp) :: sp(3)
      real(dp) :: we
      real(dp) :: wi
   end type source_terms_t

   type :: estimator_score_t
      type(source_terms_t) :: reaction(N_REACTIONS)
   end type estimator_score_t

   type :: score_data
      type(estimator_score_t) :: cl
      type(estimator_score_t) :: cl_avg
      type(estimator_score_t) :: tr
      type(estimator_score_t) :: tr_pretab
   end type score_data

   !---------------------------------------------------------------------------
   ! Running statistics for source output.
   !---------------------------------------------------------------------------
   type :: running_stats_t
      integer :: n
      real(dp) :: mean
      real(dp) :: m2
   end type running_stats_t

   type :: source_stats_t
      type(running_stats_t) :: sn_plus
      type(running_stats_t) :: sn_minus
      type(running_stats_t) :: sn_net
      type(running_stats_t) :: sp(3)
      type(running_stats_t) :: we
      type(running_stats_t) :: wi
   end type source_stats_t

   type :: estimator_stats_t
      type(source_stats_t) :: reaction(N_REACTIONS)
   end type estimator_stats_t

   type :: score_stats_data
      type(estimator_stats_t) :: cl
      type(estimator_stats_t) :: cl_avg
      type(estimator_stats_t) :: tr
      type(estimator_stats_t) :: tr_pretab
   end type score_stats_data

   type :: system_totals_t
      real(dp) :: particles
      real(dp) :: momentum(3)
      real(dp) :: energy
   end type system_totals_t

end module data_types
