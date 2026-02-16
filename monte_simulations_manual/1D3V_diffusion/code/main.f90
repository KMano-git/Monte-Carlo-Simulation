!===============================================================================
!メインループ
!基本的にはサブルーチンの呼び出しのみで完結したい
!===============================================================================
program main
   use constants
   use particle_data
   use cdf_reader
   use random_utils
   use cross_sections
   use check_coll
   use collisions
   use io

   implicit none

   type(sim_params) :: sim
   type(plasma_params) :: plasma
   type(init_params) :: init_p
   type(diag_params) :: diag
   type(particle_t), allocatable :: particles(:)

   real(dp) :: vx_i, vy_i, vz_i
   integer :: coll_type

   ! 統計量
   real(dp) :: time
   real(dp) :: delta_E_total, delta_E_cx_total, delta_E_el_total
   real(dp) :: delta_E_avg, delta_E_cx_avg, delta_E_el_avg
   real(dp) :: delta_E_cx, delta_E_el

   integer :: i, istep, ihist
   integer :: ierr
   integer, parameter :: unit_ntscrg = 25
   integer, parameter :: unit_hist = 20

   !-----------------------------------------------------------------------------
   ! 初期化
   !-----------------------------------------------------------------------------
   write(*,'(A)') '============================================'
   write(*,'(A)') '  0D-3V Monte Carlo Simulation'
   write(*,'(A)') '============================================'
   write(*,*)

   ! パラメータ読み込み(input.nml)
   call read_input_file(sim, plasma, init_p, diag)

   ! CDFファイルの読み込み
   if (.not. sim%use_isotropic) then
      call load_elastic_cdf(trim(sim%cdf_file), ierr)
      if (ierr /= 0) then
         write(*,*) 'Error: Failed to load CDF file: ', trim(sim%cdf_file)
         stop
      end if
   end if

   ! 乱数シード初期化
   call init_random_seed(sim%seed)

   ! 粒子配列の確保
   allocate(particles(sim%n_particles))

   ! 粒子初期化
   call initialize_particles(particles, sim%n_particles, init_p)

   ! 出力ファイルオープン
   open(unit=unit_ntscrg, file=trim(sim%output_ntscrg), status='replace')
   write(unit_ntscrg, '(A)') 'time[s], delta_E_total[W/m^3], delta_E_cx[W/m^3], delta_E_el[W/m^3]'

   ! メインループ
   time  = 0.0d0
   ! 開始時点のエネルギー分布を出力、ファイルの初期化
   open(unit=unit_hist, file=trim(sim%output_hist), status='replace')
   close(unit_hist)
   call output_energy_histogram(particles, sim%n_particles, sim, diag, init_p, 0)

   ihist = 1
   do istep = 1, sim%n_steps
      time = time + sim%dt
      delta_E_total    = 0.0d0
      delta_E_cx_total = 0.0d0
      delta_E_el_total = 0.0d0

      ! 各粒子について衝突判定
      do i = 1, sim%n_particles
         if (.not. particles(i)%alive) cycle

         ! 衝突種判定
         call check_collision(particles(i), plasma, sim, sim%dt, vx_i, vy_i, vz_i, coll_type)

         ! 衝突処理
         if (coll_type == COLL_CX) then
            call collision_cx(particles(i), vx_i, vy_i, vz_i, delta_E_cx)
            delta_E_cx_total = delta_E_cx_total + delta_E_cx
            delta_E_total    = delta_E_total + delta_E_cx
         else if (coll_type == COLL_EL) then
            call collision_el(particles(i), vx_i, vy_i, vz_i, sim%use_isotropic, delta_E_el)
            delta_E_el_total = delta_E_el_total + delta_E_el
            delta_E_total    = delta_E_total + delta_E_el
         else if (coll_type == COLL_EI) then
            !call collision_ei(particles(i), sim%dt, vx_i, vy_i, vz_i, delta_E)
            !ここは未実装
         end if
      end do

      ! プラズマエネルギーへの寄与（1粒子あたり）
      delta_E_avg = delta_E_total / dble(sim%n_particles)
      delta_E_cx_avg = delta_E_cx_total / dble(sim%n_particles)
      delta_E_el_avg = delta_E_el_total / dble(sim%n_particles)

      ! 単位体積あたり、単位時間あたりのプラズマエネルギーへの寄与 [W/m^3]
      delta_E_avg = delta_E_avg * init_p%n_init / sim%dt
      delta_E_cx_avg = delta_E_cx_avg * init_p%n_init / sim%dt
      delta_E_el_avg = delta_E_el_avg * init_p%n_init / sim%dt

      ! ntscrg.csvへの書き込み
      write(unit_ntscrg,'(ES15.6E3,",",ES15.6E3,",",ES15.6E3,",",ES15.6E3)') &
         time, delta_E_avg, delta_E_cx_avg, delta_E_el_avg

      ! ヒストグラム出力、hist_timingを5個で固定しているので注意
      if (ihist <= 5) then
         if (istep == diag%hist_timing(ihist)) then
            call output_energy_histogram(particles, sim%n_particles, sim, diag, init_p, istep)
            ihist = ihist + 1
         end if
      end if

   end do

   ! 最終出力
   write(*,*)
   write(*,'(A)') '============================================'
   write(*,'(A)') '  Simulation completed!'
   write(*,'(A)') '============================================'
   write(*,*)
   call output_statistics(particles, sim%n_particles, plasma, init_p, diag, time)

   ! クリーンアップ
   close(unit_ntscrg)
   deallocate(particles)

end program main
