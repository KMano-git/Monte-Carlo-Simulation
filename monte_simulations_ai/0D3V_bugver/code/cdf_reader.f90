!===============================================================================
! Module: cdf_reader
! CDFファイル読み込みと補間（最適化版）
!===============================================================================
module cdf_reader
   use constants, only: dp, CM2_TO_M2, PI
   implicit none

   ! 弾性散乱データ
   integer, parameter :: N_ENERGY_SIGMA = 101    ! 断面積エネルギー点数
   integer, parameter :: N_ENERGY_ANGLE = 51     ! 散乱角エネルギー点数
   integer, parameter :: N_PROB = 251            ! 確率グリッド点数

   real(dp) :: energy_grid_sigma(N_ENERGY_SIGMA)   ! 断面積用エネルギーグリッド [eV]
   real(dp) :: sigma_elastic(N_ENERGY_SIGMA)       ! 弾性散乱断面積 [m²]
   real(dp) :: energy_grid_angle(N_ENERGY_ANGLE)   ! 散乱角用エネルギーグリッド [eV]
   real(dp) :: prob_grid(N_PROB)                   ! 確率グリッド
   real(dp) :: angle_cdf(N_PROB, N_ENERGY_ANGLE)   ! 散乱角CDF [rad]

   ! 対数グリッドパラメータ
   real(dp) :: log_E_grid_sigma_min, log_E_grid_sigma_dlt
   real(dp) :: log_E_grid_angle_min, log_E_grid_angle_dlt

   logical :: cdf_loaded = .false.

contains

   !-----------------------------------------------------------------------------
   ! CDFファイル読み込み（ストリーム読み込み版）
   !-----------------------------------------------------------------------------
   subroutine load_elastic_cdf(filename, ierr)
      character(len=*), intent(in)  :: filename
      integer, intent(out)          :: ierr

      integer :: iunit, ios, i, j, k
      character(len=1000000), allocatable :: buffer
      character(len=64) :: token
      real(dp), allocatable :: data_array(:)
      integer :: data_idx, buf_pos, buf_len
      integer :: start_pos
      logical :: in_data
      integer :: file_size
      character(len=1) :: ch

      ierr = 0
      iunit = 20
      data_idx = 0
      in_data = .false.

      allocate(buffer)
      allocate(data_array(23320))

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

      ! テキストモードで再度開く
      open(unit=iunit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         ierr = 1
         return
      end if

      ! xs_data_tab を探す
      do while (.true.)
         read(iunit, '(A)', iostat=ios) buffer
         if (ios /= 0) exit

         buf_len = len_trim(buffer)
         if (buf_len == 0) cycle

         ! xs_data_tab の開始を検出
         buf_pos = index(buffer, 'xs_data_tab')
         if (buf_pos > 0) then
            ! '=' の位置を探す
            buf_pos = index(buffer, '=')
            if (buf_pos > 0) then
               buffer = buffer(buf_pos+1:)
               in_data = .true.
               exit
            end if
         end if
      end do

      if (.not. in_data) then
         ierr = 2
         write(*,*) 'Error: xs_data_tab not found'
         close(iunit)
         return
      end if

      ! 数値データを読み込む
      do while (in_data)
         buf_len = len_trim(buffer)
         buf_pos = 1

         do while (buf_pos <= buf_len)
            ! 空白・カンマをスキップ
            do while (buf_pos <= buf_len)
               ch = buffer(buf_pos:buf_pos)
               if (ch /= ' ' .and. ch /= ',' .and. ch /= char(9)) exit
               buf_pos = buf_pos + 1
            end do
            if (buf_pos > buf_len) exit

            ! セミコロンで終了
            if (buffer(buf_pos:buf_pos) == ';') then
               in_data = .false.
               exit
            end if

            ! トークン終端を探す
            start_pos = buf_pos
            do while (buf_pos <= buf_len)
               ch = buffer(buf_pos:buf_pos)
               if (ch == ' ' .or. ch == ',' .or. ch == ';' .or. ch == char(9)) exit
               buf_pos = buf_pos + 1
            end do

            if (buf_pos > start_pos) then
               token = buffer(start_pos:buf_pos-1)
               data_idx = data_idx + 1
               if (data_idx <= 23320) then
                  read(token, *, iostat=ios) data_array(data_idx)
                  if (ios /= 0) data_idx = data_idx - 1
               end if
            end if
         end do

         ! 次の行を読む
         if (in_data) then
            read(iunit, '(A)', iostat=ios) buffer
            if (ios /= 0) exit
         end if
      end do

      close(iunit)

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

      ! 確率グリッド（0 - 1）
      do i = 1, N_PROB
         prob_grid(i) = (i-1.0d0) / (N_PROB-1.0d0)
      end do

      ! 断面積データ（インデックス 1-101、単位: cm² → m²）
      do i = 1, N_ENERGY_SIGMA
         sigma_elastic(i) = data_array(i) * CM2_TO_M2
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
      write(*,'(A,ES12.4,A,ES12.4)') '   Sigma range [m^2]: ', &
         minval(sigma_elastic), ' - ', maxval(sigma_elastic)

      deallocate(buffer)
      deallocate(data_array)

   end subroutine load_elastic_cdf

   !-----------------------------------------------------------------------------
   ! 弾性散乱断面積の補間（対数補間）
   !-----------------------------------------------------------------------------
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

      ! エネルギーをクリップ
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
      real(dp), intent(in) :: energy   ! 衝突エネルギー [eV]
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

end module cdf_reader
