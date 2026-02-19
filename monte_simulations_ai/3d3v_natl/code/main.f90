!===============================================================================
! Program: monte_carlo_3d3v_natl
! 3D3V Non-Analog Track-Length Monte Carlo Simulation
!
! フロー:
!   1. パラメータ読み込み, CDF読み込み
!   2. ν_max テーブル事前計算
!   3. 粒子初期化
!   4. ステップループ:
!      a. 重み更新: w = w * exp(-R_a * dt)
!      b. TLスコアリング（毎ステップ、仮想サンプリング付き）
!      c. Null-Collision判定 → Rejection判定
!      d. 実衝突時: CLスコアリング → CX/EL処理
!   5. 診断出力（ntscrg.csv: 各ステップ [W/m3]）
!
! OpenMP: -fopenmp付きでコンパイルすると粒子ループが並列化される
!         -fopenmpなしでも !$omp はコメントとして無視され逐次実行可能
!===============================================================================
program monte_carlo_3d3v_natl
   use constants
   use data_types
   use random_utils, only: init_random_seed
   use cdf_reader, only: load_elastic_cdf
   use cross_sections, only: compute_nu_max
   use dynamics, only: advance_particle, check_collision_nonanalog, &
      collision_cx, collision_el, update_weight, &
      COLL_NONE, COLL_CX, COLL_EL
   use scoring, only: score_collision_estimator, score_track_length_estimator
   use io
   implicit none

   !---------------------------------------------------------------------------
   ! 変数宣言
   !---------------------------------------------------------------------------
   type(sim_params)    :: sim
   type(plasma_params) :: plasma
   type(init_params)   :: init_p
   type(diag_params)   :: diag
   type(score_data)    :: score        !累積スコア
   type(score_data)    :: step_score   !1ステップ分のスコア
   type(particle_t), allocatable :: particles(:)

   integer :: istep, ip, ierr
   logical :: first_hist

   integer :: n_coll_cx, n_coll_el, n_coll_total
   integer :: n_alive
   real(dp) :: weight_sum

   !OpenMP reduction用のスカラー変数
   real(dp) :: s_cl_cx, s_cl_el, s_cl_ei
   real(dp) :: s_tl_cx, s_tl_el, s_tl_ei

   !スレッドプライベート変数（並列ループ内で使用）
   integer :: coll_type_l
   real(dp) :: vx_i_l, vy_i_l, vz_i_l
   real(dp) :: v_rel_l, E_rel_l
   real(dp) :: delta_E_l
   type(score_data) :: my_score

   !時間計測用
   integer :: count_start, count_end, count_rate
   real(dp) :: time_sec

   !---------------------------------------------------------------------------
   ! 初期化
   !---------------------------------------------------------------------------
   !入力ファイル読み込み
   call read_input_file(sim, plasma, init_p, diag)

   !乱数シード初期化
   call init_random_seed(sim%seed)

   !CDFデータ読み込み
   call load_elastic_cdf(sim%cdf_file, ierr)
   if (ierr /= 0) then
      write(*,*) 'Error: Failed to load CDF data'
      stop 1
   end if

   !ν_max 事前計算（σ_s * v_rel の最大値をスキャン）
   call compute_nu_max(plasma%n_i)

   !粒子配列の確保・初期化
   allocate(particles(sim%n_particles))
   call initialize_particles(particles, sim%n_particles, init_p, sim%seed)

   !スコアデータ初期化
   score%cl_ei = 0.0d0; score%cl_cx = 0.0d0; score%cl_el = 0.0d0
   score%tl_ei = 0.0d0; score%tl_cx = 0.0d0; score%tl_el = 0.0d0

   !ntscrg.csvヘッダー出力
   call output_ntscrg_header(sim%output_ntscrg)

   !初期ヒストグラム出力
   first_hist = .true.

   !---------------------------------------------------------------------------
   ! メインループ
   !---------------------------------------------------------------------------
   write(*,'(A)') ''
   write(*,'(A)') 'Starting simulation...'

   call system_clock(count_start, count_rate)

   do istep = 1, sim%n_steps

      n_coll_cx = 0
      n_coll_el = 0
      n_coll_total = 0
      n_alive = 0
      weight_sum = 0.0d0

      !1ステップ分のスコアをリセット（reduction用スカラー）
      s_cl_ei = 0.0d0; s_cl_cx = 0.0d0; s_cl_el = 0.0d0
      s_tl_ei = 0.0d0; s_tl_cx = 0.0d0; s_tl_el = 0.0d0

      !$omp parallel do &
      !$omp   private(ip, coll_type_l, vx_i_l, vy_i_l, vz_i_l, &
      !$omp           v_rel_l, E_rel_l, delta_E_l, my_score) &
      !$omp   reduction(+:n_alive, weight_sum, n_coll_cx, n_coll_el, &
      !$omp              n_coll_total, s_cl_cx, s_cl_el, s_cl_ei, &
      !$omp              s_tl_cx, s_tl_el, s_tl_ei) &
      !$omp   schedule(static)
      do ip = 1, sim%n_particles
         if (.not. particles(ip)%alive) cycle
         n_alive = n_alive + 1
         weight_sum = weight_sum + particles(ip)%weight

         !--- Step a: 重み更新 ---
         if (sim%enable_ei) then
            call update_weight(particles(ip), plasma, sim%dt, sim%weight_min)
            if (.not. particles(ip)%alive) cycle
         end if

         !--- Step b: Track-Length Estimator ---
         !スレッドローカルなスコア変数を使用
         my_score%cl_ei = 0.0d0; my_score%cl_cx = 0.0d0; my_score%cl_el = 0.0d0
         my_score%tl_ei = 0.0d0; my_score%tl_cx = 0.0d0; my_score%tl_el = 0.0d0

         call score_track_length_estimator(particles(ip), plasma, sim%dt, &
            sim%use_isotropic, sim%enable_ei, my_score)

         !--- Step c: Null-Collision判定 → Rejection判定 ---
         call check_collision_nonanalog(particles(ip), plasma, sim, sim%dt, &
            vx_i_l, vy_i_l, vz_i_l, &
            v_rel_l, E_rel_l, coll_type_l)

         !--- Step d: 実衝突時の処理 ---
         if (coll_type_l /= COLL_NONE) then
            n_coll_total = n_coll_total + 1

            if (coll_type_l == COLL_CX) then
               !荷電交換: CLスコアリング（衝突前の状態で）
               delta_E_l = 0.0d0
               call score_collision_estimator(particles(ip), plasma, &
                  v_rel_l, E_rel_l, coll_type_l, delta_E_l, &
                  sim%enable_ei, my_score)
               !衝突処理
               call collision_cx(particles(ip), vx_i_l, vy_i_l, vz_i_l, &
                  delta_E_l)
               n_coll_cx = n_coll_cx + 1

            else if (coll_type_l == COLL_EL) then
               !弾性散乱: 衝突処理（delta_Eを取得）
               call collision_el(particles(ip), vx_i_l, vy_i_l, vz_i_l, &
                  sim%use_isotropic, delta_E_l)
               !CLスコアリング（delta_Eを使用）
               call score_collision_estimator(particles(ip), plasma, &
                  v_rel_l, E_rel_l, coll_type_l, delta_E_l, &
                  sim%enable_ei, my_score)
               n_coll_el = n_coll_el + 1
            end if
         end if

         !--- スレッドローカルスコアをreduction変数に加算 ---
         s_cl_cx = s_cl_cx + my_score%cl_cx
         s_cl_el = s_cl_el + my_score%cl_el
         s_cl_ei = s_cl_ei + my_score%cl_ei
         s_tl_cx = s_tl_cx + my_score%tl_cx
         s_tl_el = s_tl_el + my_score%tl_el
         s_tl_ei = s_tl_ei + my_score%tl_ei

         !--- 位置更新 ---
         call advance_particle(particles(ip), sim%dt)

      end do
      !$omp end parallel do

      !--- reduction結果をstep_scoreに格納 ---
      step_score%cl_cx = s_cl_cx; step_score%cl_el = s_cl_el
      step_score%cl_ei = s_cl_ei
      step_score%tl_cx = s_tl_cx; step_score%tl_el = s_tl_el
      step_score%tl_ei = s_tl_ei

      !--- ステップスコアを累積スコアに加算 ---
      score%cl_ei = score%cl_ei + step_score%cl_ei
      score%cl_cx = score%cl_cx + step_score%cl_cx
      score%cl_el = score%cl_el + step_score%cl_el
      score%tl_ei = score%tl_ei + step_score%tl_ei
      score%tl_cx = score%tl_cx + step_score%tl_cx
      score%tl_el = score%tl_el + step_score%tl_el

      !--- ntscrg.csv出力（毎ステップ） ---
      call output_ntscrg_step(sim%output_ntscrg, step_score, &
         sim%n_particles, init_p%n_init, sim%dt, &
         istep, n_alive, weight_sum)

      !--- 進捗表示 ---
      if (mod(istep, diag%output_interval) == 0) then
         write(*,'(A,I6,A,I6,A,I8,A,I6,A,I6)') &
            '  Step ', istep, '/', sim%n_steps, &
            '  alive=', n_alive, &
            '  CX=', n_coll_cx, '  EL=', n_coll_el
      end if

      !--- エネルギーヒストグラム出力 ---
      call output_energy_histogram(sim%output_hist, particles, sim%n_particles, &
         diag%n_hist_bins, diag%E_hist_min, diag%E_hist_max, &
         istep, 6, diag%hist_timing, first_hist)
      if (first_hist) then
         do ip = 1, 6
            if (istep == diag%hist_timing(ip)) then
               first_hist = .false.
               exit
            end if
         end do
      end if

   end do

   call system_clock(count_end)
   time_sec = real(count_end - count_start, dp) / real(count_rate, dp)
   write(*,'(A,F12.4,A)') 'Elapsed time: ', time_sec, ' seconds'

   !---------------------------------------------------------------------------
   ! 最終出力
   !---------------------------------------------------------------------------
   write(*,'(A)') ''
   write(*,'(A)') 'Simulation complete.'

   !スコアリングサマリー（コンソール）
   call output_ntscrg_final(score, sim%n_particles, init_p%n_init, &
      sim%n_steps, sim%dt)

   !最終統計
   call output_statistics(particles, sim%n_particles)

   !後片付け
   deallocate(particles)

end program monte_carlo_3d3v_natl
