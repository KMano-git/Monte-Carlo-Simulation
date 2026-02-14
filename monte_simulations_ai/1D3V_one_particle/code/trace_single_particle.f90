!===============================================================================
! Program: trace_single_particle
! 1次元 モンテカルロ中性粒子シミュレーション - 単一粒子追跡版
! 1粒子の時間発展を追跡し、衝突イベントごとの物理状態を記録
!===============================================================================
program trace_single_particle
   use constants
   use types
   use cdf_reader
   use cross_sections
   use dynamics
   implicit none

   ! シミュレーションパラメータ
   type(sim_params)    :: sim
   type(plasma_params) :: plasma
   type(particle)      :: p

   ! ループカウンタ・変数
   integer :: step, max_collisions
   integer :: coll_type, n_total_coll, n_cx, n_el, n_ei
   real(dp) :: delta_E, time_elapsed
   real(dp) :: v, E_particle
   integer :: ierr

   ! Namelistパラメータ
   integer  :: n_steps
   real(dp) :: dt
   real(dp) :: x_min, x_max
   real(dp) :: n_bg, n_e, T_bg, T_e
   real(dp) :: x_init, E_init
   integer  :: seed
   character(len=256) :: cdf_file
   character(len=256) :: trace_file

   ! 出力ファイルユニット
   integer, parameter :: TRACE_UNIT = 30

   namelist /simulation/ n_steps, dt, x_min, x_max, seed, cdf_file
   namelist /plasma_nml/ n_bg, n_e, T_bg, T_e
   namelist /particle_init/ x_init, E_init
   namelist /trace_output/ trace_file, max_collisions

   !-----------------------------------------------------------------------------
   ! デフォルト値の設定
   !-----------------------------------------------------------------------------
   n_steps = 100000
   dt = 1.0d-8
   x_min = 0.0d0
   x_max = 0.1d0
   seed = 12345
   cdf_file = 'dd_00_elastic.cdf'

   n_bg = 1.0d19
   n_e  = 1.0d19
   T_bg = 2.0d0
   T_e  = 2.0d0

   x_init = 0.0d0
   E_init = 5.0d0

   trace_file = 'particle_trace.csv'
   max_collisions = 1000  ! 最大記録衝突数

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
      rewind(10)
      read(10, nml=trace_output, iostat=ierr)
      close(10)
      write(*,*) 'Input file loaded: input.nml'
   else
      write(*,*) 'Using default parameters (input.nml not found)'
   end if

   !-----------------------------------------------------------------------------
   ! パラメータの設定
   !-----------------------------------------------------------------------------
   sim%n_steps = n_steps
   sim%dt = dt
   sim%x_min = x_min
   sim%x_max = x_max
   sim%dx = (x_max - x_min) / 100  ! グリッド幅（参照用）

   plasma%n_bg = n_bg
   plasma%n_e  = n_e
   plasma%T_bg = T_bg
   plasma%T_e  = T_e

   !-----------------------------------------------------------------------------
   ! 乱数の初期化
   !-----------------------------------------------------------------------------
   call init_random_seed(seed)

   !-----------------------------------------------------------------------------
   ! パラメータの表示
   !-----------------------------------------------------------------------------
   write(*,*) '=================================================='
   write(*,*) ' Single Particle Trace Monte Carlo Simulation'
   write(*,*) '=================================================='
   write(*,'(A,I8)')      ' Max steps:        ', sim%n_steps
   write(*,'(A,ES12.3)')  ' Time step [s]:    ', sim%dt
   write(*,'(A,ES12.3)')  ' Domain [m]:       ', sim%x_max - sim%x_min
   write(*,'(A,ES12.3)')  ' Background n [m^-3]: ', plasma%n_bg
   write(*,'(A,F8.2)')    ' Background T [eV]:   ', plasma%T_bg
   write(*,'(A,F8.2)')    ' Initial E [eV]:      ', E_init
   write(*,*) '=================================================='

   !-----------------------------------------------------------------------------
   ! CDFファイルの読み込み
   !-----------------------------------------------------------------------------
   call load_elastic_cdf(trim(cdf_file), ierr)
   if (ierr /= 0) then
      write(*,*) 'Error: Failed to load CDF file'
      stop 1
   end if

   !-----------------------------------------------------------------------------
   ! 粒子の初期化
   !-----------------------------------------------------------------------------
   p%x = x_init
   p%weight = 1.0d0
   p%alive = .true.

   ! 単一エネルギー（+x方向）
   call set_mono_velocity(E_init, 1, p%vx, p%vy, p%vz)

   !-----------------------------------------------------------------------------
   ! 出力ファイルの準備
   !-----------------------------------------------------------------------------
   open(unit=TRACE_UNIT, file=trim(trace_file), status='replace', action='write')
   write(TRACE_UNIT, '(A)') 'step,time[s],x[m],vx[m/s],vy[m/s],vz[m/s],E[eV],coll_type,delta_E[eV]'

   ! 初期状態を記録
   v = sqrt(p%vx**2 + p%vy**2 + p%vz**2)
   E_particle = 0.5d0 * M_D * v**2 * J_TO_EV
   write(TRACE_UNIT, '(I8,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,I2,A,ES14.6)') &
      0, ',', 0.0d0, ',', p%x, ',', p%vx, ',', p%vy, ',', p%vz, ',', E_particle, ',', 0, ',', 0.0d0

   !-----------------------------------------------------------------------------
   ! メインループ - 単一粒子追跡
   !-----------------------------------------------------------------------------
   n_total_coll = 0
   n_cx = 0
   n_el = 0
   n_ei = 0
   time_elapsed = 0.0d0

   write(*,*) 'Starting particle trace...'

   do step = 1, sim%n_steps
      if (.not. p%alive) exit

      time_elapsed = step * sim%dt

      ! 衝突判定
      call check_collision(p, plasma, sim%dt, coll_type, delta_E)

      if (coll_type /= COLL_NONE) then
         n_total_coll = n_total_coll + 1

         select case (coll_type)
          case (COLL_CX)
            n_cx = n_cx + 1
          case (COLL_EL)
            n_el = n_el + 1
          case (COLL_EI)
            n_ei = n_ei + 1
         end select

         ! 衝突後の状態を記録
         v = sqrt(p%vx**2 + p%vy**2 + p%vz**2)
         E_particle = 0.5d0 * M_D * v**2 * J_TO_EV

         write(TRACE_UNIT, '(I8,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,I2,A,ES14.6)') &
            step, ',', time_elapsed, ',', p%x, ',', p%vx, ',', p%vy, ',', p%vz, ',', &
            E_particle, ',', coll_type, ',', delta_E
      end if

      ! 粒子の推進
      call advance_particle(p, sim%dt, sim)

      ! 進捗表示
      if (mod(step, n_steps/10) == 0) then
         write(*,'(A,I8,A,I8,A,I8)') ' Step: ', step, ', Collisions: ', n_total_coll, &
            ', Alive: ', merge(1, 0, p%alive)
      end if
   end do

   ! 境界超過で消滅した場合、最終状態を出力（衝突タイプ0）
   if (.not. p%alive .and. (p%x < sim%x_min .or. p%x > sim%x_max)) then
      v = sqrt(p%vx**2 + p%vy**2 + p%vz**2)
      E_particle = 0.5d0 * M_D * v**2 * J_TO_EV
      write(TRACE_UNIT, '(I8,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,I2,A,ES14.6)') &
         step, ',', time_elapsed, ',', p%x, ',', p%vx, ',', p%vy, ',', p%vz, ',', &
         E_particle, ',', 0, ',', 0.0d0
   end if

   close(TRACE_UNIT)

   !-----------------------------------------------------------------------------
   ! 結果の表示
   !-----------------------------------------------------------------------------
   write(*,*) ''
   write(*,*) '=================================================='
   write(*,*) ' Simulation Results'
   write(*,*) '=================================================='
   write(*,'(A,ES12.4,A)')  ' Simulation time:  ', time_elapsed, ' s'
   write(*,'(A,I8)')        ' Total collisions: ', n_total_coll
   write(*,'(A,I8)')        '   Charge Exchange:', n_cx
   write(*,'(A,I8)')        '   Elastic:        ', n_el
   write(*,'(A,I8)')        '   Ionization:     ', n_ei
   if (p%alive) then
      write(*,*) ' Final status: ALIVE'
   else
      write(*,*) ' Final status: TERMINATED (boundary or ionized)'
   end if
   write(*,*) '=================================================='
   write(*,'(A,A)') ' Trace data written to: ', trim(trace_file)

contains

   !-----------------------------------------------------------------------------
   ! 乱数シードの初期化
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

end program trace_single_particle
