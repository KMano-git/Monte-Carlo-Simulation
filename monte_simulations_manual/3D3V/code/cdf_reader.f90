!===============================================================================
! Module: cdf_reader
! CDFファイル読み込みと補間
!===============================================================================
module cdf_reader
   use constants, only: dp, PI, CM2_TO_M2
   implicit none

   ! 弾性散乱データ
   integer, parameter :: N_ENERGY_SIGMA = 101    ! 断面積エネルギー点数
   integer, parameter :: N_ENERGY_ANGLE = 51     ! 散乱角エネルギー点数
   integer, parameter :: N_TEMP_GRID    = 51     ! 温度点数
   integer, parameter :: N_PROB = 251            ! 確率グリッド点数

   real(dp) :: energy_grid_sigma(N_ENERGY_SIGMA)  ! 断面積用エネルギーグリッド [eV/amu]
   real(dp) :: sigma_elastic(N_ENERGY_SIGMA)      ! 弾性散乱断面積 [m²]
   real(dp) :: energy_grid_angle(N_ENERGY_ANGLE)  ! 散乱角用エネルギーグリッド [eV/amu]
   real(dp) :: temp_grid(N_TEMP_GRID)             ! 温度グリッド [eV/amu]
   real(dp) :: prob_grid(N_PROB)                  ! 確率グリッド

   ! 2D テーブルデータ (サイズ 51x51、データ元のCDF単位はCGS系)
   real(dp) :: table_reaction_rate(N_ENERGY_ANGLE, N_TEMP_GRID)
   real(dp) :: table_I_1_0(N_ENERGY_ANGLE, N_TEMP_GRID)
   real(dp) :: table_I_1_1_up(N_ENERGY_ANGLE, N_TEMP_GRID)
   real(dp) :: table_I_1_2_up2(N_ENERGY_ANGLE, N_TEMP_GRID)

   real(dp) :: angle_cdf(N_PROB, N_ENERGY_ANGLE)  ! 散乱角CDF [rad]

   ! 対数グリッドパラメータ
   real(dp) :: log_E_grid_sigma_min, log_E_grid_sigma_dlt
   real(dp) :: log_E_grid_angle_min, log_E_grid_angle_dlt
   real(dp) :: log_T_grid_min, log_T_grid_dlt

   real(dp) :: sigma_unit_mult = CM2_TO_M2
   real(dp) :: reaction_rate_mult = 1.0d-6
   real(dp) :: I_1_0_mult = 1.0d-6
   real(dp) :: I_1_1_mult = 1.0d-8
   real(dp) :: I_1_2_mult = 1.0d-10

   logical :: cdf_loaded = .false.

   public :: load_elastic_cdf, get_sigma_elastic, sample_scattering_angle
   public :: get_reaction_rate, get_I_1_0, get_I_1_1_up, get_I_1_2_up2

contains

   !============================================================================
   ! CDFファイル読み込み
   !============================================================================
   subroutine load_elastic_cdf(filename, ierr)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ierr

      integer :: iunit, ios, i, j, k
      real(dp) :: data_array(23320)
      real(dp) :: xs_mult_flat(80)
      integer :: data_idx, n_xs_mult
      integer :: file_size

      ierr  = 0
      iunit = 20

      ! ファイルサイズを確認
      open(unit=iunit, file=filename, status='old', action='read', iostat=ios, &
         access='stream', form='unformatted')
      if (ios /= 0) then
         ierr = 1
         write(*,*) 'Error: Cannot open file ', trim(filename)
         return
      end if

      inquire(unit=iunit, size=file_size)
      close(iunit)

      write(*,'(A,I10,A)') ' CDF file size: ', file_size, ' bytes'

      call read_named_real_array(filename, 'xs_mult', xs_mult_flat, n_xs_mult, ios)
      if (ios == 0 .and. n_xs_mult >= 20 * 4) then
         sigma_unit_mult = xs_mult_flat(1)
         reaction_rate_mult = xs_mult_flat(5)
         I_1_0_mult = xs_mult_flat(9)
         I_1_1_mult = xs_mult_flat(13)
         I_1_2_mult = xs_mult_flat(17)
      else
         write(*,*) 'Warning: xs_mult not found; using built-in unit multipliers'
         sigma_unit_mult = CM2_TO_M2
         reaction_rate_mult = 1.0d-6
         I_1_0_mult = 1.0d-6
         I_1_1_mult = 1.0d-8
         I_1_2_mult = 1.0d-10
      end if

      call read_named_real_array(filename, 'xs_data_tab', data_array, data_idx, ios)
      if (ios /= 0) then
         ierr = 2
         write(*,*) 'Error: xs_data_tab not found'
         return
      end if

      write(*,'(A,I8)') ' CDF data points read: ', data_idx

      if (data_idx < 23306) then
         ierr = 3
         write(*,*) 'Error: Insufficient data in CDF file'
         return
      end if

      !------------------------------------------------------------------------------
      ! 読み込んだ数値を配列に格納
      !------------------------------------------------------------------------------

      ! エネルギーグリッドの構築（対数等間隔: 0.001 - 100 eV）
      log_E_grid_sigma_min = log(0.001d0)
      log_E_grid_sigma_dlt = log(100.0d0/0.001d0) / (N_ENERGY_SIGMA - 1.0d0)
      do i = 1, N_ENERGY_SIGMA
         energy_grid_sigma(i) = exp(log_E_grid_sigma_min + (i-1)*log_E_grid_sigma_dlt)
      end do

      log_E_grid_angle_min = log(0.001d0)
      log_E_grid_angle_dlt = log(100.0d0/0.001d0) / (N_ENERGY_ANGLE - 1.0d0)
      do i = 1, N_ENERGY_ANGLE
         energy_grid_angle(i) = exp(log_E_grid_angle_min + (i-1)*log_E_grid_angle_dlt)
      end do

      ! 温度グリッドの構築（CDFから取得した範囲に基づく）
      log_T_grid_min = log(0.00049665961148444d0)
      log_T_grid_dlt = log(49.665961148444d0 / 0.00049665961148444d0) / (N_TEMP_GRID - 1.0d0)
      do i = 1, N_TEMP_GRID
         temp_grid(i) = exp(log_T_grid_min + (i-1)*log_T_grid_dlt)
      end do

      ! 確率グリッド（0 - 1）
      do i = 1, N_PROB
         prob_grid(i) = (i-1.0d0) / (N_PROB-1.0d0)
      end do

      ! 断面積データ（インデックス 1-101、単位: cm² → m²）
      do i = 1, N_ENERGY_SIGMA
         sigma_elastic(i) = data_array(i) * sigma_unit_mult
      end do

      ! 2Dテーブル群 (Col-major: Eが連続)
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

      ! 散乱角CDFデータ（インデックス 10506-23306、Column-major形式）
      k = 10505
      do j = 1, N_ENERGY_ANGLE
         do i = 1, N_PROB
            k = k + 1
            angle_cdf(i, j) = data_array(k)
         end do
      end do

      cdf_loaded = .true.
      write(*,*) 'CDF file loaded successfully.'
      write(*,'(A,5ES12.4)') '   xs_mult(y): ', sigma_unit_mult, reaction_rate_mult, &
         I_1_0_mult, I_1_1_mult, I_1_2_mult
      write(*,'(A,ES12.4,A,ES12.4)') '   Sigma range [m^2]: ', &
         minval(sigma_elastic), ' - ', maxval(sigma_elastic)

   end subroutine load_elastic_cdf

   !============================================================================
   ! named array の実数データを読み込む
   !============================================================================
   subroutine read_named_real_array(filename, var_name, values, n_read, ierr)
      character(len=*), intent(in) :: filename, var_name
      real(dp), intent(out) :: values(:)
      integer, intent(out) :: n_read
      integer, intent(out) :: ierr

      integer :: iunit, ios
      integer :: buf_len, buf_pos, start_pos
      character(len=1000000) :: buffer
      character(len=64) :: token
      character(len=1) :: ch
      logical :: in_data, found
      real(dp) :: parsed_value

      values = 0.0d0
      n_read = 0
      ierr = 0
      found = .false.
      in_data = .false.
      iunit = 21

      open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         ierr = 1
         return
      end if

      do
         read(iunit, '(A)', iostat=ios) buffer
         if (ios /= 0) exit

         if (index(buffer, trim(var_name)) <= 0) cycle

         buf_pos = index(buffer, '=')
         if (buf_pos <= 0) cycle

         buffer = buffer(buf_pos+1:)
         found = .true.
         in_data = .true.
         exit
      end do

      if (.not. found) then
         ierr = 2
         close(iunit)
         return
      end if

      do while (in_data)
         buf_len = len_trim(buffer)
         buf_pos = 1

         do while (buf_pos <= buf_len)
            do while (buf_pos <= buf_len)
               ch = buffer(buf_pos:buf_pos)
               if (ch /= ' ' .and. ch /= ',' .and. ch /= char(9)) exit
               buf_pos = buf_pos + 1
            end do
            if (buf_pos > buf_len) exit

            if (buffer(buf_pos:buf_pos) == ';') then
               in_data = .false.
               exit
            end if

            start_pos = buf_pos
            do while (buf_pos <= buf_len)
               ch = buffer(buf_pos:buf_pos)
               if (ch == ' ' .or. ch == ',' .or. ch == ';' .or. ch == char(9)) exit
               buf_pos = buf_pos + 1
            end do

            if (buf_pos > start_pos .and. n_read < size(values)) then
               token = buffer(start_pos:buf_pos-1)
               read(token, *, iostat=ios) parsed_value
               if (ios == 0) then
                  n_read = n_read + 1
                  values(n_read) = parsed_value
               end if
            end if
         end do

         if (in_data) then
            read(iunit, '(A)', iostat=ios) buffer
            if (ios /= 0) exit
         end if
      end do

      close(iunit)
   end subroutine read_named_real_array

   !============================================================================
   ! 散乱断面積の取得
   !============================================================================
   function get_sigma_elastic(energy) result(sigma)
      real(dp), intent(in) :: energy
      real(dp) :: sigma

      real(dp) :: E_clipped, log_E
      real(dp) :: t
      integer :: idx

      if (.not. cdf_loaded) then
         sigma = 0.0d0 ! 衝突を起こさない
         return
      end if

      ! エネルギー [eV/amu] をクリップ
      E_clipped = max(energy_grid_sigma(1), min(energy, energy_grid_sigma(N_ENERGY_SIGMA)))
      log_E = log(E_clipped)

      ! インデックス計算 (対数等間隔グリッド)
      idx = int((log_E - log_E_grid_sigma_min) / log_E_grid_sigma_dlt) + 1
      idx = max(1, min(idx, N_ENERGY_SIGMA - 1))

      ! 対数線形補間
      t = (log_E - (log_E_grid_sigma_min + (idx-1)*log_E_grid_sigma_dlt)) / log_E_grid_sigma_dlt
      t = max(0.0d0, min(1.0d0, t))

      sigma = exp((1.0d0 - t) * log(sigma_elastic(idx)) + t * log(sigma_elastic(idx+1)))

   end function get_sigma_elastic

   !-----------------------------------------------------------------------------
   ! 散乱角のサンプリング（CDFから逆変換）
   !-----------------------------------------------------------------------------
   function sample_scattering_angle(energy, rand_p) result(chi)
      real(dp), intent(in) :: energy   ! 衝突エネルギー [eV/amu]
      real(dp), intent(in) :: rand_p   ! 乱数 [0, 1]
      real(dp) :: chi                  ! 散乱角 [rad]

      real(dp) :: E_clipped, log_E
      real(dp) :: t, s
      integer :: e_idx, p_idx
      real(dp) :: cdf_low(N_PROB), cdf_high(N_PROB), cdf_interp(N_PROB)

      if (.not. cdf_loaded) then
         chi = PI / 2.0d0  ! デフォルト値
         return
      end if

      ! エネルギーをクリップ
      E_clipped = max(energy_grid_angle(1), min(energy, energy_grid_angle(N_ENERGY_ANGLE)))
      log_E = log(E_clipped)

      ! エネルギーインデックス計算
      e_idx = int((log_E - log_E_grid_angle_min) / log_E_grid_angle_dlt) + 1
      e_idx = max(1, min(e_idx, N_ENERGY_ANGLE - 1))

      ! エネルギー補間係数
      t = (log_E - (log_E_grid_angle_min + (e_idx-1)*log_E_grid_angle_dlt)) / log_E_grid_angle_dlt
      t = max(0.0d0, min(1.0d0, t))

      ! CDFを補間
      cdf_low = angle_cdf(:, e_idx)
      cdf_high = angle_cdf(:, e_idx+1)
      cdf_interp = (1.0d0 - t) * cdf_low + t * cdf_high

      ! 確率インデックス
      p_idx = 1
      do while (p_idx < N_PROB .and. prob_grid(p_idx+1) < rand_p)
         p_idx = p_idx + 1
      end do
      p_idx = min(p_idx, N_PROB - 1)

      ! 確率方向の線形補間
      s = (rand_p - prob_grid(p_idx)) / (prob_grid(p_idx+1) - prob_grid(p_idx))
      s = max(0.0d0, min(1.0d0, s))

      chi = (1.0d0 - s) * cdf_interp(p_idx) + s * cdf_interp(p_idx+1)

   end function sample_scattering_angle

   !============================================================================
   ! 2D 双一次補間用内部関数
   ! reaction_rate / I_1_x は CDF 側で y 値も log spacing なので、
   ! DEGAS の xdt_eval2 と同じく log(value) を補間してから exp に戻す。
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

      E_clipped = max(energy_grid_angle(1), min(energy, energy_grid_angle(N_ENERGY_ANGLE)))
      T_clipped = max(temp_grid(1), min(temp, temp_grid(N_TEMP_GRID)))
      log_E = log(E_clipped)
      log_T = log(T_clipped)

      e_idx = int((log_E - log_E_grid_angle_min) / log_E_grid_angle_dlt) + 1
      e_idx = max(1, min(e_idx, N_ENERGY_ANGLE - 1))
      t_E = (log_E - (log_E_grid_angle_min + (e_idx-1)*log_E_grid_angle_dlt)) / log_E_grid_angle_dlt

      t_idx = int((log_T - log_T_grid_min) / log_T_grid_dlt) + 1
      t_idx = max(1, min(t_idx, N_TEMP_GRID - 1))
      t_T = (log_T - (log_T_grid_min + (t_idx-1)*log_T_grid_dlt)) / log_T_grid_dlt

      v1 = val_table(e_idx, t_idx)
      v2 = val_table(e_idx+1, t_idx)
      v3 = val_table(e_idx, t_idx+1)
      v4 = val_table(e_idx+1, t_idx+1)

      if (min(v1, v2, v3, v4) > 0.0d0) then
         val = exp((1.0d0 - t_E) * (1.0d0 - t_T) * log(v1) + &
            t_E           * (1.0d0 - t_T) * log(v2) + &
            (1.0d0 - t_E) * t_T           * log(v3) + &
            t_E           * t_T           * log(v4))
      else
         val = (1.0d0 - t_E) * (1.0d0 - t_T) * v1 + &
            t_E           * (1.0d0 - t_T) * v2 + &
            (1.0d0 - t_E) * t_T           * v3 + &
            t_E           * t_T           * v4
      end if
   end function interp_2d

   !----------------------------------------------------------------------------
   ! MKS 単位でルックアップ変数を取得
   !----------------------------------------------------------------------------
   function get_reaction_rate(energy, temp) result(val)
      real(dp), intent(in) :: energy, temp
      real(dp) :: val
      val = interp_2d(table_reaction_rate, energy, temp) * reaction_rate_mult
   end function get_reaction_rate

   function get_I_1_0(energy, temp) result(val)
      real(dp), intent(in) :: energy, temp
      real(dp) :: val
      val = interp_2d(table_I_1_0, energy, temp) * I_1_0_mult
   end function get_I_1_0

   function get_I_1_1_up(energy, temp) result(val)
      real(dp), intent(in) :: energy, temp
      real(dp) :: val
      val = interp_2d(table_I_1_1_up, energy, temp) * I_1_1_mult
   end function get_I_1_1_up

   function get_I_1_2_up2(energy, temp) result(val)
      real(dp), intent(in) :: energy, temp
      real(dp) :: val
      val = interp_2d(table_I_1_2_up2, energy, temp) * I_1_2_mult
   end function get_I_1_2_up2

end module cdf_reader
