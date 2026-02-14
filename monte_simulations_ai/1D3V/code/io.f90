!===============================================================================
! Module: io
! 入出力処理（結果出力、乱数初期化）
!===============================================================================
module io
   use constants, only: dp, EV_TO_J
   use types, only: score_grid, sim_params
   implicit none

contains

   !-----------------------------------------------------------------------------
   ! 乱数シードの初期化
   !-----------------------------------------------------------------------------
   subroutine init_random_seed(seed_val)
      integer, intent(in) :: seed_val
      integer, allocatable :: seed(:)
      integer :: n

      call random_seed(size=n)
      allocate(seed(n))
      seed = seed_val
      call random_seed(put=seed)

   end subroutine init_random_seed

   !-----------------------------------------------------------------------------
   ! 結果のCSV出力（物理単位 [W/m²]）
   ! 正規化: (結果/n_particles) * gamma_in = [W/m²]
   !-----------------------------------------------------------------------------
   subroutine write_results(filename, grid, params)
      character(len=*), intent(in) :: filename
      type(score_grid), intent(in) :: grid
      type(sim_params), intent(in) :: params

      integer :: i, iunit
      real(dp) :: x_center
      real(dp) :: total_time, ev_to_watt, norm_factor
      real(dp) :: power_cx, power_el, power_ei, power_total
      real(dp) :: power_cx_tl, power_el_tl, power_ei_tl, power_total_tl
      real(dp) :: sum_cx, sum_el, sum_ei, sum_total
      real(dp) :: sum_cx_tl, sum_el_tl, sum_ei_tl, sum_total_tl

      ! 変換係数
      total_time = params%n_steps * params%dt  ! 総シミュレーション時間 [s]
      ev_to_watt = EV_TO_J  ! 1 eV = 1.60218e-19 J
      ! 正規化: (1/n_particles) * gamma_in → 物理量 [/m²]
      norm_factor = params%gamma_in / dble(params%n_particles)

      iunit = 20
      open(unit=iunit, file=filename, status='replace', action='write')

      ! ヘッダー（物理単位 W/m²）
      write(iunit, '(A)') 'x_center[m],CL_Power_CX[W/m2],CL_Power_EL[W/m2],' // &
         'CL_Power_EI[W/m2],CL_Power_Total[W/m2],' // &
         'TL_Power_CX[W/m2],TL_Power_EL[W/m2],TL_Power_EI[W/m2],TL_Power_Total[W/m2]'

      ! 合計値の初期化
      sum_cx = 0.0d0
      sum_el = 0.0d0
      sum_ei = 0.0d0
      sum_total = 0.0d0
      sum_cx_tl = 0.0d0
      sum_el_tl = 0.0d0
      sum_ei_tl = 0.0d0
      sum_total_tl = 0.0d0

      ! データ
      do i = 1, params%n_grid
         x_center = params%x_min + (i - 0.5d0) * params%dx

         ! eV → W/m² (= eV/s * gamma_in/n_particles * eV→J)
         power_cx = grid%cl_energy_cx(i) / total_time * ev_to_watt * norm_factor
         power_el = grid%cl_energy_el(i) / total_time * ev_to_watt * norm_factor
         power_ei = grid%cl_energy_ei(i) / total_time * ev_to_watt * norm_factor
         power_total = grid%cl_energy(i) / total_time * ev_to_watt * norm_factor

         power_cx_tl = grid%tl_energy_cx(i) / total_time * ev_to_watt * norm_factor
         power_el_tl = grid%tl_energy_el(i) / total_time * ev_to_watt * norm_factor
         power_ei_tl = grid%tl_energy_ei(i) / total_time * ev_to_watt * norm_factor
         power_total_tl = grid%tl_energy(i) / total_time * ev_to_watt * norm_factor

         write(iunit, '(ES14.6,8(A,ES14.6))') &
            x_center, ',', power_cx, ',', power_el, ',', power_ei, ',', power_total, &
            ',', power_cx_tl, ',', power_el_tl, ',', power_ei_tl, ',', power_total_tl

         ! 合計
         sum_cx = sum_cx + power_cx
         sum_el = sum_el + power_el
         sum_ei = sum_ei + power_ei
         sum_total = sum_total + power_total
         sum_cx_tl = sum_cx_tl + power_cx_tl
         sum_el_tl = sum_el_tl + power_el_tl
         sum_ei_tl = sum_ei_tl + power_ei_tl
         sum_total_tl = sum_total_tl + power_total_tl
      end do

      close(iunit)
      write(*,'(A,A)') ' Results written to: ', filename

      ! 総計の表示
      write(*,*) ''
      write(*,*) ' --- Power Density (CL method) ---'
      write(*,'(A,ES12.4,A)') '   Charge Exchange: ', sum_cx, ' W/m²'
      write(*,'(A,ES12.4,A)') '   Elastic:         ', sum_el, ' W/m²'
      write(*,'(A,ES12.4,A)') '   Ionization:      ', sum_ei, ' W/m²'
      write(*,'(A,ES12.4,A)') '   Total:           ', sum_total, ' W/m²'
      write(*,*) ''
      write(*,*) ' --- Power Density (TL method) ---'
      write(*,'(A,ES12.4,A)') '   Charge Exchange: ', sum_cx_tl, ' W/m²'
      write(*,'(A,ES12.4,A)') '   Elastic:         ', sum_el_tl, ' W/m²'
      write(*,'(A,ES12.4,A)') '   Ionization:      ', sum_ei_tl, ' W/m²'
      write(*,'(A,ES12.4,A)') '   Total:           ', sum_total_tl, ' W/m²'

   end subroutine write_results

   !-----------------------------------------------------------------------------
   ! 密度のCSV出力（物理単位）
   ! 滞在時間と対応するフラックス密度を出力
   !-----------------------------------------------------------------------------
   subroutine write_density(filename, grid, params)
      character(len=*), intent(in) :: filename
      type(score_grid), intent(in) :: grid
      type(sim_params), intent(in) :: params

      integer :: i, iunit
      real(dp) :: x_center
      real(dp) :: total_time
      real(dp) :: residence_time, relative_density

      ! 正規化係数
      total_time = params%n_steps * params%dt  ! 総シミュレーション時間 [s]

      iunit = 21
      open(unit=iunit, file=filename, status='replace', action='write')

      ! ヘッダー
      write(iunit, '(A)') 'x_center[m],residence_time[s],flux_density[m-2s-1]'

      ! データ
      do i = 1, params%n_grid
         x_center = params%x_min + (i - 0.5d0) * params%dx

         ! 滞在時間: density(i) / n_particles は 1粒子あたりの滞在時間 [s]
         residence_time = grid%density(i) / dble(params%n_particles)

         ! フラックス密度: 滞在時間 / 総時間 * gamma_in = 各グリッドを通過するフラックス密度
         relative_density = (residence_time / total_time) * params%gamma_in

         write(iunit, '(ES14.6,A,ES14.6,A,ES14.6)') &
            x_center, ',', residence_time, ',', relative_density
      end do

      close(iunit)
      write(*,'(A,A)') ' Density written to: ', filename

   end subroutine write_density

end module io
