!===============================================================================
! Module: cdf_reader
! Reads dd_00_elastic.cdf and provides interpolation routines for elastic
! collision cross sections, 2D reaction-rate tables, and scattering-angle CDFs.
!===============================================================================
module cdf_reader
   use constants, only: dp, PI
   implicit none

   !---------------------------------------------------------------------------
   ! Grid size parameters
   !---------------------------------------------------------------------------
   integer, parameter :: N_ENERGY_SIGMA = 101  ! number of energy grid points for cross section
   integer, parameter :: N_ENERGY_ANGLE = 51   ! number of energy grid points for 2D tables / scattering angle
   integer, parameter :: N_TEMP_GRID    = 51   ! number of temperature grid points
   integer, parameter :: N_PROB         = 251  ! number of probability grid points for CDF

   !---------------------------------------------------------------------------
   ! 1D grids  (all in CDF-native units: eV for energy, eV/amu for temperature)
   !---------------------------------------------------------------------------
   real(dp) :: energy_grid_sigma(N_ENERGY_SIGMA)  ! energy grid for cross section [eV]
   real(dp) :: sigma_elastic(N_ENERGY_SIGMA)      ! elastic cross section [CGS, scaled by mult_sigma]
   real(dp) :: energy_grid_angle(N_ENERGY_ANGLE)  ! energy grid for 2D tables & scattering angle [eV]
   real(dp) :: temp_grid(N_TEMP_GRID)             ! temperature grid [eV/amu]
   real(dp) :: prob_grid(N_PROB)                  ! cumulative probability grid [dimensionless, 0-1]

   !---------------------------------------------------------------------------
   ! 2D table data  (N_ENERGY_ANGLE x N_TEMP_GRID, stored in CDF internal units)
   ! Unit conversion to SI is performed at query time via mult_* factors.
   !---------------------------------------------------------------------------
   real(dp) :: table_reaction_rate(N_ENERGY_ANGLE, N_TEMP_GRID)
   real(dp) :: table_I_1_0(N_ENERGY_ANGLE, N_TEMP_GRID)
   real(dp) :: table_I_1_1_up(N_ENERGY_ANGLE, N_TEMP_GRID)
   real(dp) :: table_I_1_2_up2(N_ENERGY_ANGLE, N_TEMP_GRID)

   !---------------------------------------------------------------------------
   ! Scattering-angle CDF: angle_cdf(i_prob, i_energy) [rad]
   !---------------------------------------------------------------------------
   real(dp) :: angle_cdf(N_PROB, N_ENERGY_ANGLE)

   !---------------------------------------------------------------------------
   ! Pre-computed log-scale grid parameters for fast index lookup
   !---------------------------------------------------------------------------
   real(dp) :: log_E_grid_sigma_min, log_E_grid_sigma_dlt
   real(dp) :: log_E_grid_2d_min,    log_E_grid_2d_dlt
   real(dp) :: log_E_grid_angle_min, log_E_grid_angle_dlt
   real(dp) :: log_T_grid_min, log_T_grid_dlt

   !---------------------------------------------------------------------------
   ! Unit conversion multipliers  (CDF internal units -> SI)
   !   mult_sigma          : cross section       1e-4   [cm^2] -> [m^2]
   !   mult_reaction_rate  : reaction rate        1e-6   [cm^3/s] -> [m^3/s]
   !   mult_I_1_0          : I_{1,0} integral     1e-6   [cm^3/s] -> [m^3/s]
   !   mult_I_1_1_up       : I_{1,1}*u'           1e-8   [cm^3*cm/s^2] -> [m^3*m/s^2]
   !   mult_I_1_2_up2      : I_{1,2}*u'^2         1e-10  [cm^3*cm^2/s^3] -> [m^3*m^2/s^3]
   !---------------------------------------------------------------------------
   real(dp) :: mult_sigma, mult_reaction_rate, mult_I_1_0, mult_I_1_1_up, mult_I_1_2_up2

   logical :: cdf_loaded = .false.

   public :: load_elastic_cdf, get_sigma_elastic, sample_scattering_angle
   public :: get_reaction_rate, get_I_1_0, get_I_1_1_up, get_I_1_2_up2

contains

   !============================================================================
   ! Load the CDF file and populate all grids and tables.
   !
   ! The CDF text file contains four named arrays:
   !   xs_min, xs_max  : grid boundary metadata
   !   xs_mult         : unit-conversion multipliers
   !   xs_data_tab     : flattened data (cross sections + 2D tables + angle CDF)
   !============================================================================
   subroutine load_elastic_cdf(filename, ierr)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ierr

      integer :: iunit, ios, i, j, k
      character(len=1000000) :: buffer

      ! Total expected data points in xs_data_tab:
      !   101 (sigma) + 4*51*51 (2D tables) + 251*51 (angle CDF) = 23,306
      real(dp) :: data_array(23320)
      real(dp) :: min_array(60), max_array(60), mult_array(80)
      integer :: data_idx, n_read_min, n_read_max, n_read_mult
      integer :: buf_len

      logical :: found_min, found_max, found_mult, found_data

      ierr = 0
      iunit = 21
      data_idx = 0
      found_min  = .false.
      found_max  = .false.
      found_mult = .false.
      found_data = .false.

      ! Open the CDF file
      open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         ierr = 1
         write(*,*) 'Error: Cannot open file ', trim(filename)
         return
      end if

      ! Scan file line by line, extracting each named array
      do while (.true.)
         read(iunit, '(A)', iostat=ios) buffer
         if (ios /= 0) exit

         buf_len = len_trim(buffer)
         if (buf_len == 0) cycle

         if (.not. found_min .and. &
            index(buffer, 'xs_min') > 0 .and. index(buffer, '=') > 0) then
            call extract_values(iunit, buffer, 60, min_array, n_read_min)
            found_min = .true.
            cycle
         end if

         if (.not. found_max .and. &
            index(buffer, 'xs_max') > 0 .and. index(buffer, '=') > 0) then
            call extract_values(iunit, buffer, 60, max_array, n_read_max)
            found_max = .true.
            cycle
         end if

         if (.not. found_mult .and. &
            index(buffer, 'xs_mult') > 0 .and. index(buffer, '=') > 0) then
            call extract_values(iunit, buffer, 80, mult_array, n_read_mult)
            found_mult = .true.
            cycle
         end if

         if (.not. found_data .and. &
            index(buffer, 'xs_data_tab') > 0 .and. index(buffer, '=') > 0) then
            call extract_values(iunit, buffer, 23320, data_array, data_idx)
            found_data = .true.
            cycle
         end if
      end do

      close(iunit)

      ! Validate that all required arrays were found
      if (.not. found_data .or. .not. found_min .or. &
         .not. found_max  .or. .not. found_mult) then
         ierr = 2
         write(*,*) 'Error: Required arrays (xs_min, xs_max, xs_mult, xs_data_tab) not fully found'
         return
      end if
      if (data_idx < 23306) then
         ierr = 3
         write(*,*) 'Error: Insufficient data points in xs_data_tab (expected >= 23306, got', data_idx, ')'
         return
      end if

      !------------------------------------------------------------------------
      ! Build logarithmic grids from metadata arrays
      !------------------------------------------------------------------------

      ! 1. Energy grid for elastic cross section  (Row 1 -> index 1)
      !    Grid range: [ min_array(1), max_array(1) ] in [eV]
      log_E_grid_sigma_min = log(min_array(1))
      log_E_grid_sigma_dlt = log(max_array(1) / min_array(1)) / (real(N_ENERGY_SIGMA, dp) - 1.0d0)
      do i = 1, N_ENERGY_SIGMA
         energy_grid_sigma(i) = exp(log_E_grid_sigma_min + (real(i, dp) - 1.0d0) * log_E_grid_sigma_dlt)
      end do

      ! 2. Energy grid for 2D tables  (Row 2 -> index 4)
      !    Grid range: [ min_array(4), max_array(4) ] in [eV]
      log_E_grid_2d_min = log(min_array(4))
      log_E_grid_2d_dlt = log(max_array(4) / min_array(4)) / (real(N_ENERGY_ANGLE, dp) - 1.0d0)

      ! 3. Temperature grid  (Row 3 -> index 5)
      !    Grid range: [ min_array(5), max_array(5) ] in [eV/amu]
      log_T_grid_min = log(min_array(5))
      log_T_grid_dlt = log(max_array(5) / min_array(5)) / (real(N_TEMP_GRID, dp) - 1.0d0)
      do i = 1, N_TEMP_GRID
         temp_grid(i) = exp(log_T_grid_min + (real(i, dp) - 1.0d0) * log_T_grid_dlt)
      end do

      ! 4. Energy grid for scattering-angle CDF  (Row 6 -> index 17)
      !    Grid range: [ min_array(17), max_array(17) ] in [eV]
      !    Note: this overwrites energy_grid_angle, which is shared with the 2D tables.
      log_E_grid_angle_min = log(min_array(17))
      log_E_grid_angle_dlt = log(max_array(17) / min_array(17)) / (real(N_ENERGY_ANGLE, dp) - 1.0d0)
      do i = 1, N_ENERGY_ANGLE
         energy_grid_angle(i) = exp(log_E_grid_angle_min + (real(i, dp) - 1.0d0) * log_E_grid_angle_dlt)
      end do

      ! 5. Probability grid for CDF (linear, Row 6 -> index 16)
      !    Grid range: [ min_array(16), max_array(16) ], dimensionless [0, 1]
      do i = 1, N_PROB
         prob_grid(i) = min_array(16) &
            + (real(i, dp) - 1.0d0) * (max_array(16) - min_array(16)) / (real(N_PROB, dp) - 1.0d0)
      end do

      !------------------------------------------------------------------------
      ! Store unit-conversion multipliers (CDF internal units -> SI)
      !------------------------------------------------------------------------
      mult_sigma         = mult_array(1)   ! 1e-4:   [cm^2] -> [m^2]
      mult_reaction_rate = mult_array(5)   ! 1e-6:   [cm^3/s] -> [m^3/s]
      mult_I_1_0         = mult_array(9)   ! 1e-6:   [cm^3/s] -> [m^3/s]
      mult_I_1_1_up      = mult_array(13)  ! 1e-8:   [cm^3*cm/s^2] -> [m^3*m/s^2]
      mult_I_1_2_up2     = mult_array(17)  ! 1e-10:  [cm^3*cm^2/s^3] -> [m^3*m^2/s^3]

      !------------------------------------------------------------------------
      ! Unpack flattened data array into structured tables
      !------------------------------------------------------------------------

      ! (a) Elastic cross section: indices 1..N_ENERGY_SIGMA
      !     Stored value is multiplied by mult_sigma to convert [cm^2] -> [m^2]
      do i = 1, N_ENERGY_SIGMA
         sigma_elastic(i) = data_array(i) * mult_sigma
      end do

      ! (b) 2D tables (each N_ENERGY_ANGLE * N_TEMP_GRID = 2601 values)
      !     Stored in CDF internal units; conversion applied at query time.
      k = N_ENERGY_SIGMA
      do j = 1, N_TEMP_GRID
         do i = 1, N_ENERGY_ANGLE
            k = k + 1
            table_reaction_rate(i, j) = data_array(k)
         end do
      end do

      do j = 1, N_TEMP_GRID
         do i = 1, N_ENERGY_ANGLE
            k = k + 1
            table_I_1_0(i, j) = data_array(k)
         end do
      end do

      do j = 1, N_TEMP_GRID
         do i = 1, N_ENERGY_ANGLE
            k = k + 1
            table_I_1_1_up(i, j) = data_array(k)
         end do
      end do

      do j = 1, N_TEMP_GRID
         do i = 1, N_ENERGY_ANGLE
            k = k + 1
            table_I_1_2_up2(i, j) = data_array(k)
         end do
      end do
      ! k should now be 10505 (= 101 + 4*51*51)

      ! (c) Scattering-angle CDF: starts at flat index 10506
      !     Values are scattering angles in [rad]
      do j = 1, N_ENERGY_ANGLE
         do i = 1, N_PROB
            k = k + 1
            angle_cdf(i, j) = data_array(k)
         end do
      end do

      cdf_loaded = .true.
      write(*,*) 'CDF file loaded successfully.'

   end subroutine load_elastic_cdf

   !============================================================================
   ! Extract numeric values from a CDF text line of the form "name = v1, v2, ..."
   ! Continues reading subsequent lines until a semicolon terminator is found.
   !============================================================================
   subroutine extract_values(iunit, buffer_in, max_vals, out_array, n_read)
      integer, intent(in)   :: iunit
      character(len=*), intent(in) :: buffer_in
      integer, intent(in)   :: max_vals
      real(dp), intent(out) :: out_array(max_vals)
      integer, intent(out)  :: n_read

      character(len=1000000) :: buffer
      character(len=64) :: token
      integer :: buf_pos, buf_len, start_pos, ios
      character(len=1) :: ch
      logical :: in_data

      n_read = 0

      ! Start parsing after the '=' sign
      buf_pos = index(buffer_in, '=')
      buffer = buffer_in(buf_pos + 1:)
      in_data = .true.

      do while (in_data)
         buf_len = len_trim(buffer)
         buf_pos = 1

         do while (buf_pos <= buf_len)
            ! Skip whitespace, commas, and tabs
            do while (buf_pos <= buf_len)
               ch = buffer(buf_pos:buf_pos)
               if (ch /= ' ' .and. ch /= ',' .and. ch /= char(9)) exit
               buf_pos = buf_pos + 1
            end do
            if (buf_pos > buf_len) exit

            ! Semicolon marks end of data
            if (buffer(buf_pos:buf_pos) == ';') then
               in_data = .false.
               buf_pos = buf_pos + 1
            end if

            ! Extract token (delimited by space, comma, tab, or semicolon)
            start_pos = buf_pos
            do while (buf_pos <= buf_len)
               ch = buffer(buf_pos:buf_pos)
               if (ch == ' ' .or. ch == ',' .or. ch == char(9) .or. ch == ';') exit
               buf_pos = buf_pos + 1
            end do

            if (buf_pos > start_pos) then
               token = buffer(start_pos:buf_pos-1)
               n_read = n_read + 1
               if (n_read <= max_vals) then
                  read(token, *, iostat=ios) out_array(n_read)
                  if (ios /= 0) n_read = n_read - 1  ! discard invalid token
               end if
            end if
         end do

         ! Read next line if data continues
         if (in_data) then
            read(iunit, '(A)', iostat=ios) buffer
            if (ios /= 0) exit
         end if
      end do
   end subroutine extract_values

   !============================================================================
   ! Interpolate the elastic cross section at a given energy.
   !
   ! Input:  energy [eV]
   ! Output: sigma  [m^2]  (already converted from [cm^2] during loading)
   !
   ! Uses log-log interpolation on the 1D energy grid.
   !============================================================================
   function get_sigma_elastic(energy) result(sigma)
      real(dp), intent(in) :: energy
      real(dp) :: sigma

      real(dp) :: E_clipped, log_E, t
      integer  :: idx

      if (.not. cdf_loaded) then
         sigma = 0.0d0
         return
      end if

      ! Clamp energy to grid bounds [eV]
      E_clipped = max(energy_grid_sigma(1), min(energy, energy_grid_sigma(N_ENERGY_SIGMA)))
      log_E = log(E_clipped)

      ! Find lower grid index
      idx = int((log_E - log_E_grid_sigma_min) / log_E_grid_sigma_dlt) + 1
      idx = max(1, min(idx, N_ENERGY_SIGMA - 1))

      ! Fractional position within the grid cell
      t = (log_E - (log_E_grid_sigma_min + (idx - 1) * log_E_grid_sigma_dlt)) / log_E_grid_sigma_dlt
      t = max(0.0d0, min(1.0d0, t))

      ! Log-log linear interpolation -> result in [m^2]
      sigma = exp((1.0d0 - t) * log(sigma_elastic(idx)) + t * log(sigma_elastic(idx + 1)))
   end function get_sigma_elastic

   !============================================================================
   ! Sample a scattering angle from the CDF via inverse-transform sampling.
   !
   ! Input:  energy [eV], rand_p [dimensionless, 0-1]
   ! Output: chi    [rad]
   !
   ! Performs bilinear interpolation in (energy, probability) space.
   !============================================================================
   function sample_scattering_angle(energy, rand_p) result(chi)
      real(dp), intent(in) :: energy
      real(dp), intent(in) :: rand_p
      real(dp) :: chi

      real(dp) :: E_clipped, log_E, t, s
      integer :: e_idx, p_idx
      real(dp) :: cdf_low(N_PROB), cdf_high(N_PROB), cdf_interp(N_PROB)

      if (.not. cdf_loaded) then
         chi = PI / 2.0d0  ! default to 90 degrees
         return
      end if

      ! Clamp energy to grid bounds [eV]
      E_clipped = max(energy_grid_angle(1), min(energy, energy_grid_angle(N_ENERGY_ANGLE)))
      log_E = log(E_clipped)

      ! Find lower energy index on the angle-CDF grid
      e_idx = int((log_E - log_E_grid_angle_min) / log_E_grid_angle_dlt) + 1
      e_idx = max(1, min(e_idx, N_ENERGY_ANGLE - 1))

      ! Fractional position within the energy grid cell
      t = (log_E - (log_E_grid_angle_min + (e_idx - 1) * log_E_grid_angle_dlt)) / log_E_grid_angle_dlt
      t = max(0.0d0, min(1.0d0, t))

      ! Interpolate angle CDF between adjacent energy grid points
      cdf_low  = angle_cdf(:, e_idx)
      cdf_high = angle_cdf(:, e_idx + 1)
      cdf_interp = (1.0d0 - t) * cdf_low + t * cdf_high

      ! Find lower probability index via linear search
      p_idx = 1
      do while (p_idx < N_PROB .and. prob_grid(p_idx + 1) < rand_p)
         p_idx = p_idx + 1
      end do
      p_idx = min(p_idx, N_PROB - 1)

      ! Linear interpolation within the probability cell -> result in [rad]
      s = (rand_p - prob_grid(p_idx)) / (prob_grid(p_idx + 1) - prob_grid(p_idx))
      s = max(0.0d0, min(1.0d0, s))

      chi = (1.0d0 - s) * cdf_interp(p_idx) + s * cdf_interp(p_idx + 1)
   end function sample_scattering_angle

   !============================================================================
   ! 2D bilinear interpolation on log-scale energy and temperature grids.
   !
   ! Input:  val_table(N_ENERGY_ANGLE, N_TEMP_GRID) - table in CDF internal units
   !         energy [eV], temp [eV/amu]
   ! Output: interpolated value in CDF internal units (caller applies mult_*)
   !============================================================================
   function interp_2d(val_table, energy, temp) result(val)
      real(dp), intent(in) :: val_table(N_ENERGY_ANGLE, N_TEMP_GRID)
      real(dp), intent(in) :: energy, temp
      real(dp) :: val

      real(dp) :: E_clipped, T_clipped, log_E, log_T
      real(dp) :: t_E, t_T, v1, v2, v3, v4
      integer :: e_idx, t_idx

      if (.not. cdf_loaded) then
         val = 0.0d0
         return
      end if

      ! Clamp energy [eV] and temperature [eV/amu] to their respective grid bounds
      E_clipped = max(energy_grid_angle(1), min(energy, energy_grid_angle(N_ENERGY_ANGLE)))
      T_clipped = max(temp_grid(1), min(temp, temp_grid(N_TEMP_GRID)))
      log_E = log(E_clipped)
      log_T = log(T_clipped)

      ! Locate lower energy index on the 2D energy grid
      e_idx = int((log_E - log_E_grid_2d_min) / log_E_grid_2d_dlt) + 1
      e_idx = max(1, min(e_idx, N_ENERGY_ANGLE - 1))
      t_E = (log_E - (log_E_grid_2d_min + (e_idx - 1) * log_E_grid_2d_dlt)) / log_E_grid_2d_dlt

      ! Locate lower temperature index
      t_idx = int((log_T - log_T_grid_min) / log_T_grid_dlt) + 1
      t_idx = max(1, min(t_idx, N_TEMP_GRID - 1))
      t_T = (log_T - (log_T_grid_min + (t_idx - 1) * log_T_grid_dlt)) / log_T_grid_dlt

      ! Four corner values for bilinear interpolation
      v1 = val_table(e_idx,     t_idx)
      v2 = val_table(e_idx + 1, t_idx)
      v3 = val_table(e_idx,     t_idx + 1)
      v4 = val_table(e_idx + 1, t_idx + 1)

      ! Bilinear interpolation
      val = (1.0d0 - t_E) * (1.0d0 - t_T) * v1 + &
         t_E  * (1.0d0 - t_T) * v2 + &
         (1.0d0 - t_E) *          t_T  * v3 + &
         t_E  *          t_T  * v4
   end function interp_2d

   !============================================================================
   ! Public query functions for 2D tables.
   ! Each interpolates in (energy [eV], temperature [eV/amu]) and applies the
   ! corresponding unit-conversion multiplier to return values in SI units.
   !============================================================================

   ! Reaction rate: interpolated value * mult -> [m^3/s]  (from CDF internal [cm^3/s])
   function get_reaction_rate(energy, temp) result(val)
      real(dp), intent(in) :: energy, temp
      real(dp) :: val
      val = interp_2d(table_reaction_rate, energy, temp) * mult_reaction_rate
   end function get_reaction_rate

   ! I_{1,0} integral: interpolated value * mult -> [m^3/s]  (from CDF internal [cm^3/s])
   function get_I_1_0(energy, temp) result(val)
      real(dp), intent(in) :: energy, temp
      real(dp) :: val
      val = interp_2d(table_I_1_0, energy, temp) * mult_I_1_0
   end function get_I_1_0

   ! I_{1,1}*u': interpolated value * mult -> [SI]  (from CDF internal [CGS])
   function get_I_1_1_up(energy, temp) result(val)
      real(dp), intent(in) :: energy, temp
      real(dp) :: val
      val = interp_2d(table_I_1_1_up, energy, temp) * mult_I_1_1_up
   end function get_I_1_1_up

   ! I_{1,2}*u'^2: interpolated value * mult -> [SI]  (from CDF internal [CGS])
   function get_I_1_2_up2(energy, temp) result(val)
      real(dp), intent(in) :: energy, temp
      real(dp) :: val
      val = interp_2d(table_I_1_2_up2, energy, temp) * mult_I_1_2_up2
   end function get_I_1_2_up2

end module cdf_reader
