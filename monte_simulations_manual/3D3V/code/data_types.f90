!===============================================================================
! Module: data_types
! particle data, simulation parameters, and scoring data structures
!===============================================================================
module data_types
   use constants, only : dp
   implicit none

   !---------------------------------------------------------------------------
   ! random number generator state (xoshiro256+ for 64bit x 4)
   !---------------------------------------------------------------------------
   type :: rng_state
      integer(8) :: s(4)
   end type rng_state

   !---------------------------------------------------------------------------
   ! particle data structures (3D3V + weight + rng)
   !---------------------------------------------------------------------------
   type :: particle_t
      real(dp) :: x, y, z       ! position [m]
      real(dp) :: vx, vy, vz    ! velocity [m/s]
      real(dp) :: weight        ! particle weight (decreases due to ionization)
      logical  :: alive         ! alive flag
      real(dp) :: zincx         ! random distance to next collision
      real(dp) :: zint1         ! accumulated distance
      type(rng_state) :: rng    ! particle-specific random number state
   end type particle_t

   !---------------------------------------------------------------------------
   ! simulation parameters
   !---------------------------------------------------------------------------
   type :: sim_params
      integer  :: n_particles    ! number of particles
      integer  :: n_steps        ! number of steps
      real(dp) :: dt             ! time step [s]
      integer  :: seed           ! random seed
      logical  :: enable_cx      ! charge exchange On/Off
      logical  :: enable_el      ! elastic scattering On/Off
      logical  :: enable_ei      ! ionization On/Off
      logical  :: use_isotropic  ! isotropic scattering flag
      real(dp) :: weight_min     ! minimum weight threshold (Russian roulette)
      character(len=256) :: cdf_file       ! CDF file name
      character(len=256) :: output_ntscrg  ! scoring output file name
      character(len=256) :: output_hist    ! energy histogram output file name
      character(len=256) :: output_deltaE_hist ! energy transfer histogram output file name
   end type sim_params

   !---------------------------------------------------------------------------
   ! plasma parameters
   !---------------------------------------------------------------------------
   type :: plasma_params
      real(dp) :: n_i            ! background ion density [m^-3]
      real(dp) :: T_i            ! background ion temperature [eV]
      real(dp) :: n_e            ! background electron density [m^-3]
      real(dp) :: T_e            ! background electron temperature [eV]
      real(dp) :: u_x, u_y, u_z  ! background ion flow velocity [m/s]
   end type plasma_params

   !---------------------------------------------------------------------------
   ! particle initial condition parameters
   !---------------------------------------------------------------------------
   type :: init_params
      real(dp) :: n_init         ! initial particle density [m^-3]
      real(dp) :: E_init         ! initial particle energy [eV]
      real(dp) :: T_init         ! initial particle temperature [eV]
      integer  :: init_mode      ! initial particle distribution (0: mono, 1: maxwell)
   end type init_params

   !---------------------------------------------------------------------------
   ! diagnostic parameters
   !---------------------------------------------------------------------------
   type :: diag_params
      ! energy histogram
      integer  :: output_interval  ! output interval
      integer  :: n_hist_bins      ! number of histogram bins
      real(dp) :: E_hist_min       ! energy minimum [eV]
      real(dp) :: E_hist_max       ! energy maximum [eV]
      integer  :: hist_timing(6)   ! histogram output timing
      ! energy transfer histogram
      integer  :: n_dE_bins        ! number of energy transfer bins
      real(dp) :: dE_hist_min      ! energy transfer minimum [eV]
      real(dp) :: dE_hist_max      ! energy transfer maximum [eV]
      integer  :: dE_collect_steps ! energy transfer collection steps (0 = all steps)
   end type diag_params

   !---------------------------------------------------------------------------
   ! scoring data structures (CL/TL simultaneous storage)
   !---------------------------------------------------------------------------
   type :: score_data
      ! Collision Estimator [J]
      real(dp) :: cl_ei           ! ionization contribution
      real(dp) :: cl_cx           ! charge exchange contribution
      real(dp) :: cl_el           ! elastic scattering contribution
      ! Track-Length Estimator [J]
      real(dp) :: tl_ei           ! ionization contribution
      real(dp) :: tl_cx           ! charge exchange contribution
      real(dp) :: tl_el           ! elastic scattering contribution
      ! Table-Lookup Track-Length Estimator [J]
      real(dp) :: tl_lookup_ei    ! ionization contribution (table lookup)
      real(dp) :: tl_lookup_cx    ! charge exchange contribution (table lookup)
      real(dp) :: tl_lookup_el    ! elastic scattering contribution (table lookup)
   end type score_data

end module data_types
