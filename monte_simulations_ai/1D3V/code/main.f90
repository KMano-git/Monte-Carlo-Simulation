!===============================================================================
! Program: monte_carlo
! 1次元 モンテカルロ中性粒子シミュレーション
!===============================================================================
program monte_carlo
   use constants
   use types
   use cdf_reader
   use cross_sections
   use scoring
   use dynamics
   use io
   implicit none

   ! シミュレーションパラメータ
   type(sim_params)    :: sim
   type(plasma_params) :: plasma
   type(init_params)   :: init
   type(score_grid)    :: grid

   ! 粒子配列
   type(particle), allocatable :: particles(:)

   ! ループカウンタ・変数
   integer :: i, step, n_alive, total_collisions
   integer :: coll_type, n_cx, n_el, n_ei
   real(dp) :: delta_E
   real(dp) :: progress
   integer :: ierr

   ! Namelistパラメータ
   integer  :: n_particles
   integer  :: n_steps
   integer  :: n_grid
   real(dp) :: dt
   real(dp) :: x_min, x_max
   real(dp) :: n_bg, n_e, T_bg, T_e
   real(dp) :: x_init, E_init, T_init
   real(dp) :: gamma_in   ! 入射粒子フラックス [m^-2 s^-1]
   integer  :: init_mode  ! 0: 単一エネルギー, 1: Maxwell分布
   integer  :: seed
   character(len=256) :: cdf_file
   character(len=256) :: output_file
   character(len=256) :: density_file

   namelist /simulation/ n_particles, n_steps, n_grid, dt, x_min, x_max, &
      seed, cdf_file, output_file
   namelist /plasma_nml/ n_bg, n_e, T_bg, T_e
   namelist /particle_init/ x_init, E_init, T_init, init_mode, gamma_in

   !-----------------------------------------------------------------------------
   ! デフォルト値の設定
   !-----------------------------------------------------------------------------
   n_particles = 10000
   n_steps = 10000
   n_grid = 100
   dt = 1.0d-8
   x_min = 0.0d0
   x_max = 0.1d0
   seed = 12345
   cdf_file = 'dd_00_elastic.cdf'
   output_file = 'ntscrg.csv'
   density_file = 'ntden.csv'

   n_bg = 1.0d19
   n_e  = 1.0d19    ! 電子密度（デフォルトはイオン密度と同じ）
   T_bg = 2.0d0
   T_e  = 2.0d0     ! 電子温度

   x_init = 0.0d0
   E_init = 5.0d0
   T_init = 2.0d0
   init_mode = 0
   gamma_in = 1.0d23  ! デフォルト入射フラックス [m^-2 s^-1]

   !-----------------------------------------------------------------------------
   ! Namelistファイルの読み込み
   !-----------------------------------------------------------------------------
   open(unit=10, file='input.nml', status='old', action='read', iostat=ierr)
   if (ierr == 0) then
      read(10, nml=simulation, iostat=ierr)
      rewind(10)
      read(10, nml=plasma_nml, iostat=ierr)
      rewind(10)
      read(10, nml=particle_init, iostat=ierr)
      close(10)
      write(*,*) 'Input file loaded: input.nml'
   else
      write(*,*) 'Using default parameters (input.nml not found)'
   end if

   !-----------------------------------------------------------------------------
   ! パラメータの設定
   !-----------------------------------------------------------------------------
   sim%n_particles = n_particles
   sim%n_steps = n_steps
   sim%n_grid = n_grid
   sim%dt = dt
   sim%x_min = x_min
   sim%x_max = x_max
   sim%dx = (x_max - x_min) / n_grid
   sim%gamma_in = gamma_in

   plasma%n_bg = n_bg
   plasma%n_e  = n_e
   plasma%T_bg = T_bg
   plasma%T_e  = T_e

   init%x_init = x_init
   init%E_init = E_init
   init%T_init = T_init

   !-----------------------------------------------------------------------------
   ! 乱数の初期化
   !-----------------------------------------------------------------------------
   call init_random_seed(seed)

   !-----------------------------------------------------------------------------
   ! パラメータの表示
   !-----------------------------------------------------------------------------
   write(*,*) '=============================================='
   write(*,*) ' Monte Carlo Neutral Particle Simulation'
   write(*,*) '=============================================='
   write(*,'(A,I8)')      ' Particles:       ', sim%n_particles
   write(*,'(A,I8)')      ' Steps:           ', sim%n_steps
   write(*,'(A,I8)')      ' Grid points:     ', sim%n_grid
   write(*,'(A,ES12.3)')  ' Time step [s]:   ', sim%dt
   write(*,'(A,ES12.3)')  ' Domain [m]:      ', sim%x_max - sim%x_min
   write(*,'(A,ES12.3)')  ' Background n [m^-3]: ', plasma%n_bg
   write(*,'(A,F8.2)')    ' Background T [eV]:   ', plasma%T_bg
   write(*,'(A,F8.2)')    ' Initial E [eV]:      ', init%E_init
   write(*,*) '=============================================='

   !-----------------------------------------------------------------------------
   ! CDFファイルの読み込み
   !-----------------------------------------------------------------------------
   call load_elastic_cdf(trim(cdf_file), ierr)
   if (ierr /= 0) then
      write(*,*) 'Error: Failed to load CDF file'
      stop 1
   end if

   !-----------------------------------------------------------------------------
   ! スコアリンググリッドの初期化
   !-----------------------------------------------------------------------------
   call init_score_grid(grid, sim%n_grid)

   !-----------------------------------------------------------------------------
   ! 粒子の初期化
   !-----------------------------------------------------------------------------
   allocate(particles(sim%n_particles))

   do i = 1, sim%n_particles
      particles(i)%x = init%x_init
      particles(i)%weight = 1.0d0
      particles(i)%alive = .true.

      if (init_mode == 0) then
         ! 単一エネルギー（+x方向）
         call set_mono_velocity(init%E_init, 1, &
            particles(i)%vx, particles(i)%vy, particles(i)%vz)
      else
         ! Maxwell分布
         call sample_maxwell_velocity(init%T_init, &
            particles(i)%vx, particles(i)%vy, particles(i)%vz)
         ! +x方向にのみ進む（vx > 0にする）
         particles(i)%vx = abs(particles(i)%vx)
      end if
   end do

   !-----------------------------------------------------------------------------
   ! メインループ
   !-----------------------------------------------------------------------------
   total_collisions = 0
   n_cx = 0
   n_el = 0
   n_ei = 0

   write(*,*) 'Starting simulation...'

   do step = 1, sim%n_steps
      n_alive = 0

      do i = 1, sim%n_particles
         if (.not. particles(i)%alive) cycle

         ! Track-Length Estimator（毎ステップ）
         call score_tl(grid, particles(i), sim, plasma, sim%dt)

         ! 衝突判定
         call check_collision(particles(i), plasma, sim%dt, coll_type, delta_E)

         if (coll_type /= COLL_NONE) then
            ! Collision Estimator（衝突時のみ）- 衝突タイプも渡す
            call score_cl(grid, particles(i), delta_E, coll_type, sim)
            total_collisions = total_collisions + 1

            select case (coll_type)
             case (COLL_CX)
               n_cx = n_cx + 1
             case (COLL_EL)
               n_el = n_el + 1
             case (COLL_EI)
               n_ei = n_ei + 1
            end select
         end if

         ! 粒子の推進
         call advance_particle(particles(i), sim%dt, sim)

         if (particles(i)%alive) n_alive = n_alive + 1
      end do

      ! 進捗表示
      if (mod(step, sim%n_steps/10) == 0 .or. step == sim%n_steps) then
         progress = 100.0d0 * step / sim%n_steps
         write(*,'(A,F6.1,A,I8,A,I8)') ' Progress: ', progress, '%, Alive: ', n_alive, &
            ', Collisions: ', total_collisions
      end if

      ! 全粒子が消滅したら終了
      if (n_alive == 0) then
         write(*,*) 'All particles absorbed. Stopping simulation.'
         exit
      end if
   end do

   !-----------------------------------------------------------------------------
   ! 結果の出力
   !-----------------------------------------------------------------------------
   write(*,*) ''
   write(*,*) '=============================================='
   write(*,*) ' Simulation Results'
   write(*,*) '=============================================='
   write(*,'(A,I8)') ' Total collisions:  ', total_collisions
   write(*,'(A,I8)') '   Charge Exchange: ', n_cx
   write(*,'(A,I8)') '   Elastic:         ', n_el
   write(*,'(A,I8)') '   Ionization:      ', n_ei
   write(*,*) '=============================================='

   ! CSVファイルへの出力
   call write_results(trim(output_file), grid, sim)
   call write_density(trim(density_file), grid, sim)

   !-----------------------------------------------------------------------------
   ! クリーンアップ
   !-----------------------------------------------------------------------------
   call free_score_grid(grid)
   deallocate(particles)

   write(*,*) 'Simulation completed successfully!'

end program monte_carlo
