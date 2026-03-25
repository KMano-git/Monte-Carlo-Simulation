!===============================================================================
! Module: data_types
! 1D3V 用データ型定義
!===============================================================================
module data_types
   use constants, only : dp
   implicit none

   type :: rng_state
      integer(8) :: s(4)
   end type rng_state

   type :: particle_t
      real(dp) :: x, y, z
      real(dp) :: vx, vy, vz
      real(dp) :: weight
      logical  :: alive
      real(dp) :: zincx
      real(dp) :: zint1
      real(dp) :: tl_pending_time
      real(dp) :: tl_pending_x_start
      real(dp) :: tl_pending_weight_start
      type(rng_state) :: rng
   end type particle_t

   type :: sim_params
      integer  :: n_particles
      integer  :: n_steps
      integer  :: n_x_bins
      integer  :: seed
      real(dp) :: dt
      real(dp) :: x_min
      real(dp) :: x_max
      real(dp) :: dx
      real(dp) :: gamma_in
      real(dp) :: weight_min
      logical  :: enable_cx
      logical  :: enable_el
      logical  :: enable_ei
      logical  :: use_isotropic
      logical  :: enable_tl_lookup
      character(len=256) :: cdf_file
      character(len=256) :: output_ntscrg
      character(len=256) :: output_profile
      character(len=256) :: output_hist
      character(len=256) :: output_deltaE_hist
   end type sim_params

   type :: plasma_params
      real(dp) :: n_i
      real(dp) :: T_i
      real(dp) :: n_e
      real(dp) :: T_e
      real(dp) :: u_x, u_y, u_z
   end type plasma_params

   type :: init_params
      real(dp) :: x_init
      real(dp) :: y_init
      real(dp) :: z_init
      real(dp) :: E_init
      real(dp) :: T_init
      integer  :: init_mode
   end type init_params

   type :: diag_params
      integer  :: output_interval
      integer  :: n_hist_bins
      real(dp) :: E_hist_min
      real(dp) :: E_hist_max
      integer  :: hist_timing(6)
      integer  :: n_dE_bins
      integer  :: dE_collect_steps
      real(dp) :: dE_hist_min
      real(dp) :: dE_hist_max
   end type diag_params

   type :: score_data
      real(dp) :: cl_ei
      real(dp) :: cl_cx
      real(dp) :: cl_el
      real(dp) :: tl_ei
      real(dp) :: tl_cx
      real(dp) :: tl_el
      real(dp) :: tl_lookup_ei
      real(dp) :: tl_lookup_cx
      real(dp) :: tl_lookup_el
   end type score_data

   type :: profile_data
      real(dp), allocatable :: residence_eff(:)
      real(dp), allocatable :: cl_ei(:)
      real(dp), allocatable :: cl_cx(:)
      real(dp), allocatable :: cl_el(:)
      real(dp), allocatable :: tl_ei(:)
      real(dp), allocatable :: tl_cx(:)
      real(dp), allocatable :: tl_el(:)
      real(dp), allocatable :: tl_lookup_ei(:)
      real(dp), allocatable :: tl_lookup_cx(:)
      real(dp), allocatable :: tl_lookup_el(:)
      real(dp), allocatable :: collision_count_cx(:)
      real(dp), allocatable :: collision_count_el(:)
      real(dp), allocatable :: collision_count_total(:)
   end type profile_data

end module data_types
