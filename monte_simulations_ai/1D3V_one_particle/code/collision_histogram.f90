!===============================================================================
! Program: collision_histogram
! 衝突時のエネルギー交換量ヒストグラム生成
! 指定温度での CX/EL 反応によるエネルギー変化分布を計算
!===============================================================================
program collision_histogram
   use constants
   use types
   use cdf_reader
   use cross_sections
   implicit none

   ! パラメータ
   real(dp) :: T_neutral    ! 中性粒子温度 [eV] (E_neutral=0の時に使用)
   real(dp) :: E_neutral    ! 中性粒子固定エネルギー [eV] (>0で固定エネルギーモード)
   real(dp) :: T_ion        ! イオン温度 [eV]
   real(dp) :: n_bg         ! 背景密度 [m^-3]
   integer  :: n_samples    ! サンプル数
   integer  :: n_bins       ! ビン数

   ! 結果配列
   real(dp), allocatable :: delta_E_cx(:), delta_E_el(:)
   real(dp), allocatable :: weight_cx(:), weight_el(:)
   real(dp), allocatable :: hist_cx(:), hist_el(:)
   real(dp), allocatable :: bin_centers(:)

   ! 作業変数
   real(dp) :: E_min, E_max, bin_width
   real(dp) :: v_th_n, v_th_i
   real(dp) :: vn(3), vi(3), v_rel(3), v_rel_mag, v_cm(3)
   real(dp) :: E_n, E_i, E_rel
   real(dp) :: sig_cx, sig_el
   real(dp) :: chi, phi, rand_p, rand_phi
   real(dp) :: vn_new(3), E_n_new
   integer  :: i, bin_idx, ierr
   integer  :: seed

   ! 統計
   real(dp) :: total_rate_cx, total_rate_el
   real(dp) :: mean_dE_cx, mean_dE_el
   real(dp) :: simple_mean_cx, simple_mean_el

   character(len=256) :: cdf_file, output_file
   integer, parameter :: OUT_UNIT = 40

   namelist /histogram_params/ T_neutral, E_neutral, T_ion, n_bg, n_samples, n_bins, &
      seed, cdf_file, output_file

   !-----------------------------------------------------------------------------
   ! デフォルト値
   !-----------------------------------------------------------------------------
   T_neutral = 1.0d0
   E_neutral = 0.0d0   ! 0 = Maxwell分布モード, >0 = 固定エネルギーモード
   T_ion = 5.0d0
   n_bg = 1.0d19
   n_samples = 100000
   n_bins = 80
   seed = 54321
   cdf_file = 'dd_00_elastic.cdf'
   output_file = 'collision_histogram.dat'

   !-----------------------------------------------------------------------------
   ! Namelist読み込み
   !-----------------------------------------------------------------------------
   open(unit=10, file='histogram_input.nml', status='old', action='read', iostat=ierr)
   if (ierr == 0) then
      read(10, nml=histogram_params, iostat=ierr)
      close(10)
      write(*,*) 'Input file loaded: histogram_input.nml'
   else
      write(*,*) 'Using default parameters (histogram_input.nml not found)'
   end if

   !-----------------------------------------------------------------------------
   ! 初期化
   !-----------------------------------------------------------------------------
   call init_random_seed(seed)

   write(*,*) '=================================================='
   write(*,*) ' Collision Energy Transfer Histogram'
   write(*,*) '=================================================='
   if (E_neutral > 0.0d0) then
      write(*,'(A,F8.2,A)')   ' Neutral energy (FIXED): ', E_neutral, ' eV'
   else
      write(*,'(A,F8.2,A)')   ' Neutral temperature:    ', T_neutral, ' eV'
   end if
   write(*,'(A,F8.2,A)')   ' Ion temperature:     ', T_ion, ' eV'
   write(*,'(A,ES12.3,A)') ' Background density:  ', n_bg, ' m^-3'
   write(*,'(A,I10)')      ' Number of samples:   ', n_samples
   write(*,'(A,I6)')       ' Number of bins:      ', n_bins
   write(*,*) '=================================================='

   ! CDF読み込み
   call load_elastic_cdf(trim(cdf_file), ierr)
   if (ierr /= 0) then
      write(*,*) 'Error: Failed to load CDF file'
      stop 1
   end if

   ! 配列確保
   allocate(delta_E_cx(n_samples), delta_E_el(n_samples))
   allocate(weight_cx(n_samples), weight_el(n_samples))
   allocate(hist_cx(n_bins), hist_el(n_bins), bin_centers(n_bins))

   hist_cx = 0.0d0
   hist_el = 0.0d0

   ! 熱速度
   v_th_n = sqrt(T_neutral * EV_TO_J / M_D)
   v_th_i = sqrt(T_ion * EV_TO_J / M_D)

   !-----------------------------------------------------------------------------
   ! サンプリングループ
   !-----------------------------------------------------------------------------
   write(*,*) 'Computing energy transfers...'

   do i = 1, n_samples
      ! 中性粒子速度
      if (E_neutral > 0.0d0) then
         ! 固定エネルギーモード: +x方向に入射
         call set_fixed_energy_velocity(E_neutral, vn)
      else
         ! Maxwell分布モード
         call sample_maxwell(v_th_n, vn)
      end if

      ! イオン速度（Maxwell分布）
      call sample_maxwell(v_th_i, vi)

      ! 相対速度
      v_rel = vn - vi
      v_rel_mag = sqrt(v_rel(1)**2 + v_rel(2)**2 + v_rel(3)**2)
      if (v_rel_mag < 1.0d-30) v_rel_mag = 1.0d-30

      ! エネルギー計算
      E_n = 0.5d0 * M_D * (vn(1)**2 + vn(2)**2 + vn(3)**2) * J_TO_EV
      E_i = 0.5d0 * M_D * (vi(1)**2 + vi(2)**2 + vi(3)**2) * J_TO_EV
      E_rel = 0.5d0 * (M_D/2.0d0) * v_rel_mag**2 * J_TO_EV  ! 相対エネルギー（μ = M_D/2）

      ! 断面積
      sig_cx = sigma_cx(E_rel)
      sig_el = get_sigma_elastic(E_rel)

      ! 反応率重み → 衝突頻度 ν = n_bg × σ × v [s^-1]
      weight_cx(i) = n_bg * sig_cx * v_rel_mag
      weight_el(i) = n_bg * sig_el * v_rel_mag

      ! CX: 速度交換
      delta_E_cx(i) = E_i - E_n

      ! EL: 弾性散乱
      call random_number(rand_p)
      call random_number(rand_phi)
      chi = sample_scattering_angle(E_rel, rand_p)
      phi = 2.0d0 * PI * rand_phi

      ! 重心速度
      v_cm = 0.5d0 * (vn + vi)

      ! 散乱後の相対速度
      call rotate_velocity(v_rel, v_rel_mag, chi, phi, vn_new)

      ! 実験室系での中性粒子速度
      vn_new = v_cm + 0.5d0 * vn_new

      E_n_new = 0.5d0 * M_D * (vn_new(1)**2 + vn_new(2)**2 + vn_new(3)**2) * J_TO_EV
      delta_E_el(i) = E_n_new - E_n
   end do

   !-----------------------------------------------------------------------------
   ! 統計計算
   !-----------------------------------------------------------------------------
   total_rate_cx = sum(weight_cx)
   total_rate_el = sum(weight_el)

   mean_dE_cx = sum(weight_cx * delta_E_cx) / total_rate_cx
   mean_dE_el = sum(weight_el * delta_E_el) / total_rate_el

   simple_mean_cx = sum(delta_E_cx) / n_samples
   simple_mean_el = sum(delta_E_el) / n_samples

   write(*,*) ''
   write(*,*) '===== Results ====='
   write(*,*) 'Charge Exchange (CX):'
   write(*,'(A,F12.4,A)') '  Simple mean ΔE = ', simple_mean_cx, ' eV'
   write(*,'(A,F12.4,A)') '  Weighted mean ΔE = ', mean_dE_cx, ' eV'
   write(*,*) 'Elastic Scattering (EL):'
   write(*,'(A,F12.4,A)') '  Simple mean ΔE = ', simple_mean_el, ' eV'
   write(*,'(A,F12.4,A)') '  Weighted mean ΔE = ', mean_dE_el, ' eV'

   !-----------------------------------------------------------------------------
   ! ヒストグラム作成
   !-----------------------------------------------------------------------------
   ! 範囲設定（1-99パーセンタイル）
   call get_percentiles(delta_E_cx, delta_E_el, n_samples, E_min, E_max)
   bin_width = (E_max - E_min) / n_bins

   do i = 1, n_bins
      bin_centers(i) = E_min + (i - 0.5d0) * bin_width
   end do

   ! ヒストグラム積算
   do i = 1, n_samples
      bin_idx = int((delta_E_cx(i) - E_min) / bin_width) + 1
      if (bin_idx >= 1 .and. bin_idx <= n_bins) then
         hist_cx(bin_idx) = hist_cx(bin_idx) + weight_cx(i)
      end if

      bin_idx = int((delta_E_el(i) - E_min) / bin_width) + 1
      if (bin_idx >= 1 .and. bin_idx <= n_bins) then
         hist_el(bin_idx) = hist_el(bin_idx) + weight_el(i)
      end if
   end do

   ! 衍突頻度分布に正規化（n_samplesで割ってサンプル数非依存、単位: s⁻¹/eV）
   hist_cx = hist_cx / (dble(n_samples) * bin_width)
   hist_el = hist_el / (dble(n_samples) * bin_width)

   !-----------------------------------------------------------------------------
   ! 出力
   !-----------------------------------------------------------------------------
   open(unit=OUT_UNIT, file=trim(output_file), status='replace', action='write')

   write(OUT_UNIT, '(A)') '# Collision Energy Transfer Histogram'
   if (E_neutral > 0.0d0) then
      write(OUT_UNIT, '(A,F8.2,A,F8.2,A)') '# E_neutral = ', E_neutral, ' eV (FIXED), T_ion = ', T_ion, ' eV'
   else
      write(OUT_UNIT, '(A,F8.2,A,F8.2,A)') '# T_neutral = ', T_neutral, ' eV, T_ion = ', T_ion, ' eV'
   end if
   write(OUT_UNIT, '(A,ES10.3,A)') '# n_bg = ', n_bg, ' m^-3'
   write(OUT_UNIT, '(A,F12.4,A,F12.4,A)') '# CX weighted mean = ', mean_dE_cx, ' eV, EL weighted mean = ', mean_dE_el, ' eV'
   write(OUT_UNIT, '(A)') '# bin_center[eV]  CX_freq[s^-1/eV]  EL_freq[s^-1/eV]'

   do i = 1, n_bins
      write(OUT_UNIT, '(ES14.6,2(A,ES14.6))') bin_centers(i), '  ', hist_cx(i), '  ', hist_el(i)
   end do

   close(OUT_UNIT)

   write(*,*) ''
   write(*,'(A,A)') ' Histogram written to: ', trim(output_file)

   ! クリーンアップ
   deallocate(delta_E_cx, delta_E_el, weight_cx, weight_el)
   deallocate(hist_cx, hist_el, bin_centers)

contains

   !-----------------------------------------------------------------------------
   ! 乱数シード初期化
   ! seed_val < 0 の場合はシステム時刻からランダムなseedを生成
   !-----------------------------------------------------------------------------
   subroutine init_random_seed(seed_val)
      integer, intent(in) :: seed_val
      integer, allocatable :: seed(:)
      integer :: n, actual_seed
      integer :: time_values(8)

      call random_seed(size=n)
      allocate(seed(n))

      if (seed_val < 0) then
         ! システム時刻からseedを生成
         call date_and_time(values=time_values)
         actual_seed = time_values(5) * 3600000 + &  ! 時
            time_values(6) * 60000 + &    ! 分
            time_values(7) * 1000 + &     ! 秒
            time_values(8)                 ! ミリ秒
         write(*,'(A,I12)') ' Using random seed: ', actual_seed
      else
         actual_seed = seed_val
      end if

      seed = actual_seed
      call random_seed(put=seed)
   end subroutine init_random_seed

   !-----------------------------------------------------------------------------
   ! Maxwell分布から速度サンプリング
   !-----------------------------------------------------------------------------
   subroutine sample_maxwell(v_th, v)
      real(dp), intent(in) :: v_th
      real(dp), intent(out) :: v(3)
      real(dp) :: u1, u2, u3, u4

      call random_number(u1)
      call random_number(u2)
      call random_number(u3)
      call random_number(u4)

      v(1) = v_th * sqrt(-2.0d0 * log(max(u1, 1.0d-30))) * cos(2.0d0 * PI * u2)
      v(2) = v_th * sqrt(-2.0d0 * log(max(u1, 1.0d-30))) * sin(2.0d0 * PI * u2)
      v(3) = v_th * sqrt(-2.0d0 * log(max(u3, 1.0d-30))) * cos(2.0d0 * PI * u4)
   end subroutine sample_maxwell

   !-----------------------------------------------------------------------------
   ! 固定エネルギー速度設定（+x方向）
   !-----------------------------------------------------------------------------
   subroutine set_fixed_energy_velocity(E_eV, v)
      real(dp), intent(in) :: E_eV
      real(dp), intent(out) :: v(3)
      real(dp) :: v_mag

      v_mag = sqrt(2.0d0 * E_eV * EV_TO_J / M_D)
      v(1) = v_mag
      v(2) = 0.0d0
      v(3) = 0.0d0
   end subroutine set_fixed_energy_velocity

   !-----------------------------------------------------------------------------
   ! 速度ベクトルの回転（散乱）
   !-----------------------------------------------------------------------------
   subroutine rotate_velocity(v_rel, v_mag, chi, phi, v_new)
      real(dp), intent(in) :: v_rel(3), v_mag, chi, phi
      real(dp), intent(out) :: v_new(3)

      real(dp) :: ex(3), ey(3), ez(3)
      real(dp) :: sin_chi, cos_chi, sin_phi, cos_phi
      real(dp) :: e_perp

      sin_chi = sin(chi)
      cos_chi = cos(chi)
      sin_phi = sin(phi)
      cos_phi = cos(phi)

      ! 相対速度方向の単位ベクトル
      ez = v_rel / v_mag

      ! 垂直方向の基底ベクトル
      if (abs(ez(3)) > 0.9d0) then
         e_perp = sqrt(ez(1)**2 + ez(2)**2)
         if (e_perp > 1.0d-10) then
            ex(1) = -ez(2) / e_perp
            ex(2) = ez(1) / e_perp
            ex(3) = 0.0d0
         else
            ex = (/1.0d0, 0.0d0, 0.0d0/)
         end if
      else
         e_perp = sqrt(ez(1)**2 + ez(3)**2)
         if (e_perp > 1.0d-10) then
            ex(1) = -ez(3) / e_perp
            ex(2) = 0.0d0
            ex(3) = ez(1) / e_perp
         else
            ex = (/1.0d0, 0.0d0, 0.0d0/)
         end if
      end if

      ! 外積: ey = ez × ex
      ey(1) = ez(2) * ex(3) - ez(3) * ex(2)
      ey(2) = ez(3) * ex(1) - ez(1) * ex(3)
      ey(3) = ez(1) * ex(2) - ez(2) * ex(1)

      ! 回転後の速度
      v_new = v_mag * (sin_chi * cos_phi * ex + sin_chi * sin_phi * ey + cos_chi * ez)
   end subroutine rotate_velocity

   !-----------------------------------------------------------------------------
   ! パーセンタイル計算（簡易版）
   !-----------------------------------------------------------------------------
   subroutine get_percentiles(arr1, arr2, n, p_min, p_max)
      real(dp), intent(in) :: arr1(n), arr2(n)
      integer, intent(in) :: n
      real(dp), intent(out) :: p_min, p_max

      real(dp) :: all_min, all_max, margin

      all_min = min(minval(arr1), minval(arr2))
      all_max = max(maxval(arr1), maxval(arr2))
      margin = 0.05d0 * (all_max - all_min)

      p_min = all_min - margin
      p_max = all_max + margin
   end subroutine get_percentiles

end program collision_histogram
