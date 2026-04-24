!===============================================================================
! Module: tl_el_table_reader
! Pre-tabulated elastic TL table の入出力と補間
!===============================================================================
module tl_el_table_reader
   use constants, only: dp
   implicit none

   character(len=*), parameter :: TL_EL_TABLE_MAGIC = 'TL_EL_TABLE_V1'
   integer, parameter :: DEFAULT_N_ENERGY = 51
   integer, parameter :: DEFAULT_N_TEMP = 51
   real(dp), parameter :: DEFAULT_E_MIN = 0.001d0
   real(dp), parameter :: DEFAULT_E_MAX = 100.0d0
   real(dp), parameter :: DEFAULT_T_MIN = 0.00049665961148444d0
   real(dp), parameter :: DEFAULT_T_MAX = 49.665961148444d0

   real(dp), allocatable :: tl_el_energy_grid(:)
   real(dp), allocatable :: tl_el_temp_grid(:)
   real(dp), allocatable :: tl_el_rate_table(:,:)
   logical :: tl_el_loaded = .false.

   private
   public :: build_default_tl_el_grids
   public :: write_tl_el_table, load_tl_el_table
   public :: get_tl_el_rate_from_table, tl_el_table_is_loaded

contains

   subroutine build_default_tl_el_grids(energy_grid, temp_grid)
      real(dp), allocatable, intent(out) :: energy_grid(:), temp_grid(:)

      integer :: i
      real(dp) :: log_e_min, log_e_dlt
      real(dp) :: log_t_min, log_t_dlt

      allocate(energy_grid(DEFAULT_N_ENERGY))
      allocate(temp_grid(DEFAULT_N_TEMP))

      log_e_min = log(DEFAULT_E_MIN)
      log_e_dlt = log(DEFAULT_E_MAX / DEFAULT_E_MIN) / real(DEFAULT_N_ENERGY - 1, dp)
      do i = 1, DEFAULT_N_ENERGY
         energy_grid(i) = exp(log_e_min + real(i - 1, dp) * log_e_dlt)
      end do

      log_t_min = log(DEFAULT_T_MIN)
      log_t_dlt = log(DEFAULT_T_MAX / DEFAULT_T_MIN) / real(DEFAULT_N_TEMP - 1, dp)
      do i = 1, DEFAULT_N_TEMP
         temp_grid(i) = exp(log_t_min + real(i - 1, dp) * log_t_dlt)
      end do
   end subroutine build_default_tl_el_grids

   subroutine write_tl_el_table(filename, energy_grid, temp_grid, rate_table, ierr)
      character(len=*), intent(in) :: filename
      real(dp), intent(in) :: energy_grid(:), temp_grid(:)
      real(dp), intent(in) :: rate_table(size(energy_grid), size(temp_grid))
      integer, intent(out) :: ierr

      integer :: i, j, ios, unit_out

      ierr = 0
      unit_out = 61

      open(unit=unit_out, file=filename, status='replace', action='write', iostat=ios)
      if (ios /= 0) then
         ierr = 1
         return
      end if

      write(unit_out, '(A)') TL_EL_TABLE_MAGIC
      write(unit_out, *) size(energy_grid), size(temp_grid)

      do i = 1, size(energy_grid)
         write(unit_out, '(ES24.16)') energy_grid(i)
      end do

      do j = 1, size(temp_grid)
         write(unit_out, '(ES24.16)') temp_grid(j)
      end do

      do j = 1, size(temp_grid)
         do i = 1, size(energy_grid)
            write(unit_out, '(ES24.16)') rate_table(i, j)
         end do
      end do

      close(unit_out)
   end subroutine write_tl_el_table

   subroutine load_tl_el_table(filename, ierr)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ierr

      integer :: i, j, ios, unit_in
      integer :: n_energy, n_temp
      character(len=64) :: magic

      ierr = 0
      unit_in = 62
      tl_el_loaded = .false.

      open(unit=unit_in, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         ierr = 1
         return
      end if

      read(unit_in, '(A)', iostat=ios) magic
      if (ios /= 0 .or. trim(magic) /= TL_EL_TABLE_MAGIC) then
         ierr = 2
         close(unit_in)
         return
      end if

      read(unit_in, *, iostat=ios) n_energy, n_temp
      if (ios /= 0 .or. n_energy < 2 .or. n_temp < 2) then
         ierr = 3
         close(unit_in)
         return
      end if

      if (allocated(tl_el_energy_grid)) deallocate(tl_el_energy_grid)
      if (allocated(tl_el_temp_grid)) deallocate(tl_el_temp_grid)
      if (allocated(tl_el_rate_table)) deallocate(tl_el_rate_table)
      allocate(tl_el_energy_grid(n_energy))
      allocate(tl_el_temp_grid(n_temp))
      allocate(tl_el_rate_table(n_energy, n_temp))

      do i = 1, n_energy
         read(unit_in, *, iostat=ios) tl_el_energy_grid(i)
         if (ios /= 0) then
            ierr = 4
            close(unit_in)
            tl_el_loaded = .false.
            return
         end if
      end do

      do j = 1, n_temp
         read(unit_in, *, iostat=ios) tl_el_temp_grid(j)
         if (ios /= 0) then
            ierr = 5
            close(unit_in)
            tl_el_loaded = .false.
            return
         end if
      end do

      do j = 1, n_temp
         do i = 1, n_energy
            read(unit_in, *, iostat=ios) tl_el_rate_table(i, j)
            if (ios /= 0) then
               ierr = 6
               close(unit_in)
               tl_el_loaded = .false.
               return
            end if
         end do
      end do

      close(unit_in)
      tl_el_loaded = .true.
   end subroutine load_tl_el_table

   logical function tl_el_table_is_loaded()
      tl_el_table_is_loaded = tl_el_loaded
   end function tl_el_table_is_loaded

   function get_tl_el_rate_from_table(energy, temp) result(val)
      real(dp), intent(in) :: energy, temp
      real(dp) :: val

      integer :: i_energy, i_temp
      real(dp) :: e_clipped, t_clipped
      real(dp) :: s_energy, s_temp
      real(dp) :: v11, v21, v12, v22

      if (.not. tl_el_loaded) then
         val = 0.0d0
         return
      end if

      e_clipped = max(tl_el_energy_grid(1), min(energy, tl_el_energy_grid(size(tl_el_energy_grid))))
      t_clipped = max(tl_el_temp_grid(1), min(temp, tl_el_temp_grid(size(tl_el_temp_grid))))

      i_energy = find_lower_index(tl_el_energy_grid, e_clipped)
      i_temp = find_lower_index(tl_el_temp_grid, t_clipped)

      s_energy = safe_fraction(tl_el_energy_grid(i_energy), tl_el_energy_grid(i_energy + 1), e_clipped)
      s_temp = safe_fraction(tl_el_temp_grid(i_temp), tl_el_temp_grid(i_temp + 1), t_clipped)

      v11 = tl_el_rate_table(i_energy, i_temp)
      v21 = tl_el_rate_table(i_energy + 1, i_temp)
      v12 = tl_el_rate_table(i_energy, i_temp + 1)
      v22 = tl_el_rate_table(i_energy + 1, i_temp + 1)

      val = (1.0d0 - s_energy) * (1.0d0 - s_temp) * v11 + &
         s_energy * (1.0d0 - s_temp) * v21 + &
         (1.0d0 - s_energy) * s_temp * v12 + &
         s_energy * s_temp * v22
   end function get_tl_el_rate_from_table

   integer function find_lower_index(grid, x)
      real(dp), intent(in) :: grid(:)
      real(dp), intent(in) :: x

      integer :: i

      find_lower_index = 1
      do i = 1, size(grid) - 1
         if (x <= grid(i + 1)) then
            find_lower_index = i
            return
         end if
      end do
      find_lower_index = size(grid) - 1
   end function find_lower_index

   real(dp) function safe_fraction(x0, x1, x)
      real(dp), intent(in) :: x0, x1, x

      if (abs(x1 - x0) <= 1.0d-30) then
         safe_fraction = 0.0d0
      else
         safe_fraction = (x - x0) / (x1 - x0)
      end if
      safe_fraction = max(0.0d0, min(1.0d0, safe_fraction))
   end function safe_fraction

end module tl_el_table_reader
