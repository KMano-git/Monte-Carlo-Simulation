!===============================================================================
! 1D3V Monte Carlo Simulation — Flux Injection Model
! 中性粒子が x=x_min から一定フラックスで入射し、
! 背景プラズマ（一様・定常）と衝突しながら輸送される
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
   real(dp) :: delta_E_cx, delta_E_el
   real(dp) :: weight  ! 1計算粒子あたりの重み

   integer :: i, istep, ihist, n_alive, n_injected
   integer :: ierr
   integer, parameter :: unit_ntscrg = 25

   !-----------------------------------------------------------------------------
   ! 初期化
   !-----------------------------------------------------------------------------
   write(*,'(A)') '============================================'
   write(*,'(A)') '  1D-3V Monte Carlo Simulation (Flux Model)'
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

   ! 粒子配列の確保（全て dead 状態で初期化）
   allocate(particles(sim%n_particles))
   call initialize_particles(particles, sim%n_particles)

   ! 重みの計算: w = Γ_in [m^-2 s^-1] * A [m^2] * dt [s] / n_inject_per_step
   ! A = 1.0 m^2（1D なので入射面の面積は 1）
   weight = init_p%gamma_in * 1.0d0 * sim%dt / dble(sim%n_inject_per_step)

   write(*,'(A,ES12.4)') '  Flux Gamma_in       = ', init_p%gamma_in
   write(*,'(A,I8)')     '  n_inject_per_step   = ', sim%n_inject_per_step
   write(*,'(A,ES12.4)') '  Weight per particle = ', weight
   write(*,'(A,I8)')     '  n_particles (max)   = ', sim%n_particles
   write(*,'(A,F8.4,A,F8.4)') '  Domain: [', sim%x_min, ', ', sim%x_max
   write(*,*)

   ! 出力ファイルオープン
   open(unit=unit_ntscrg, file=trim(sim%output_ntscrg), status='replace')
   write(unit_ntscrg, '(A)') 'time[s],n_alive,delta_E_total[J],delta_E_cx[J],delta_E_el[J],n_injected'

   ! 初期ヒストグラム出力（ファイル初期化）
   open(unit=20, file=trim(sim%output_hist), status='replace')
   close(20)
   open(unit=21, file=trim(sim%output_statx), status='replace')
   close(21)
   call output_energy_histogram(particles, sim, diag, 0)
   call output_position_histogram(particles, sim, 0)

   !-----------------------------------------------------------------------------
   ! メインループ
   !-----------------------------------------------------------------------------
   time  = 0.0d0
   ihist = 1

   do istep = 1, sim%n_steps
      time = time + sim%dt
      delta_E_total    = 0.0d0
      delta_E_cx_total = 0.0d0
      delta_E_el_total = 0.0d0

      !----------------------------------------------------------------------
      ! Step 1: 粒子注入（dead スロットに新しい入射粒子を配置）
      !----------------------------------------------------------------------
      call inject_particles(particles, sim%n_particles, sim, init_p, weight, n_injected)

      !----------------------------------------------------------------------
      ! Step 2: 位置更新 + 境界判定
      !----------------------------------------------------------------------
      do i = 1, sim%n_particles
         if (.not. particles(i)%alive) cycle

         ! 位置を更新
         particles(i)%x = particles(i)%x + particles(i)%vx * sim%dt

         ! 境界外に出たら消滅
         if (particles(i)%x < sim%x_min .or. particles(i)%x > sim%x_max) then
            particles(i)%alive = .false.
         end if
      end do

      !----------------------------------------------------------------------
      ! Step 3: 衝突判定 & 処理
      !----------------------------------------------------------------------
      do i = 1, sim%n_particles
         if (.not. particles(i)%alive) cycle

         ! 衝突種判定
         call check_collision(particles(i), plasma, sim, sim%dt, vx_i, vy_i, vz_i, coll_type)

         ! 衝突処理
         if (coll_type == COLL_CX) then
            call collision_cx(particles(i), vx_i, vy_i, vz_i, delta_E_cx)
            delta_E_cx_total = delta_E_cx_total + delta_E_cx * weight
            delta_E_total    = delta_E_total + delta_E_cx * weight
         else if (coll_type == COLL_EL) then
            call collision_el(particles(i), vx_i, vy_i, vz_i, sim%use_isotropic, delta_E_el)
            delta_E_el_total = delta_E_el_total + delta_E_el * weight
            delta_E_total    = delta_E_total + delta_E_el * weight
         else if (coll_type == COLL_EI) then
            ! 未実装
         end if
      end do

      !----------------------------------------------------------------------
      ! Step 4: 生存粒子数の計算
      !----------------------------------------------------------------------
      n_alive = 0
      do i = 1, sim%n_particles
         if (particles(i)%alive) n_alive = n_alive + 1
      end do

      !----------------------------------------------------------------------
      ! Step 5: ntscrg.csv 出力
      !----------------------------------------------------------------------
      write(unit_ntscrg,'(ES15.6E3,",",I0,",",ES15.6E3,",",ES15.6E3,",",ES15.6E3,",",I0)') &
         time, n_alive, delta_E_total, delta_E_cx_total, delta_E_el_total, n_injected

      !----------------------------------------------------------------------
      ! Step 6: ヒストグラム出力
      !----------------------------------------------------------------------
      if (ihist <= 5) then
         if (istep == diag%hist_timing(ihist)) then
            call output_energy_histogram(particles, sim, diag, istep)
            call output_position_histogram(particles, sim, istep)
            ihist = ihist + 1
         end if
      end if

      ! 進捗表示
      if (mod(istep, max(sim%n_steps/10, 1)) == 0) then
         write(*,'(A,I6,A,I6,A,I8)') '  Step ', istep, ' / ', sim%n_steps, &
            '  alive = ', n_alive
      end if

   end do

   !-----------------------------------------------------------------------------
   ! 最終出力
   !-----------------------------------------------------------------------------
   write(*,*)
   write(*,'(A)') '============================================'
   write(*,'(A)') '  Simulation completed!'
   write(*,'(A)') '============================================'
   write(*,*)
   call output_statistics(particles, sim%n_particles, plasma, time, weight)

   ! 最終ヒストグラム出力（hist_timingで既に出力済みでない場合のみ）
   if (diag%hist_timing(min(ihist-1, 5)) /= sim%n_steps) then
      call output_energy_histogram(particles, sim, diag, sim%n_steps)
      call output_position_histogram(particles, sim, sim%n_steps)
   end if

   ! クリーンアップ
   close(unit_ntscrg)
   deallocate(particles)

end program main
