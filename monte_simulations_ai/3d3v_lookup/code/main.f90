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
   use random_utils, only: init_random_seed, random_double
   use cdf_reader, only: load_elastic_cdf
   use cross_sections, only: compute_nu_max, sigma_v_max, ionization_rate_coeff
   use dynamics, only: advance_particle, evaluate_collision_event, &
      collision_cx, collision_el, update_weight, &
      COLL_NONE, COLL_CX, COLL_EL
   use scoring, only: score_collision_estimator, score_track_length_estimator, score_table_lookup_track_length
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
   real(dp), allocatable :: dE_hist_el(:), dE_hist_cx(:)

   integer :: istep, ip, ierr
   logical :: first_hist, collect_dE
   integer :: dE_start_step, dE_actual_steps
   real(dp) :: dE_collect_time

   integer :: n_coll_cx, n_coll_el, n_coll_total
   integer :: n_alive
   real(dp) :: weight_sum

   !OpenMP reduction用のスカラー変数
   real(dp) :: s_cl_cx, s_cl_el, s_cl_ei
   real(dp) :: s_tl_cx, s_tl_el, s_tl_ei
   real(dp) :: s_tl_lookup_cx, s_tl_lookup_el, s_tl_lookup_ei

   !スレッドプライベート変数（並列ループ内で使用）
   integer :: coll_type_l
   real(dp) :: vx_i_l, vy_i_l, vz_i_l
   real(dp) :: v_rel_l, E_rel_l
   real(dp) :: delta_E_l
   type(particle_t) :: p_copy
   type(score_data) :: my_score
   real(dp), allocatable :: l_dE_hist_el(:), l_dE_hist_cx(:)
   integer :: ibin_l
   real(dp) :: dE_bin_width
   real(dp) :: time_remaining, step_dt, t_col, r, nu_max
   real(dp) :: ion_loss_rate, segment_eff_time
   logical  :: is_collision_event
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

   if (sim%enable_ei) then
      ion_loss_rate = plasma%n_e * ionization_rate_coeff(plasma%T_e)
   else
      ion_loss_rate = 0.0d0
   end if

   !粒子配列の確保・初期化
   allocate(particles(sim%n_particles))
   call initialize_particles(particles, sim%n_particles, init_p, sim%seed)

   !スコアデータ初期化
   score%cl_ei = 0.0d0; score%cl_cx = 0.0d0; score%cl_el = 0.0d0
   score%tl_ei = 0.0d0; score%tl_cx = 0.0d0; score%tl_el = 0.0d0
   score%tl_lookup_ei = 0.0d0; score%tl_lookup_cx = 0.0d0; score%tl_lookup_el = 0.0d0

   !ntscrg.csvヘッダー出力
   call output_ntscrg_header(sim%output_ntscrg)

   !ヒストグラム配列の宣言と初期化
   allocate(dE_hist_el(diag%n_dE_bins), dE_hist_cx(diag%n_dE_bins))
   dE_hist_el = 0.0d0
   dE_hist_cx = 0.0d0
   dE_bin_width = (diag%dE_hist_max - diag%dE_hist_min) / dble(diag%n_dE_bins)

   !初期ヒ���トグラム出力
   first_hist = .true.

   !deltaEヒストグラム集計開始ステップの決定
   if (diag%dE_collect_steps <= 0 .or. diag%dE_collect_steps >= sim%n_steps) then
      dE_start_step = 1  !全ステップで集計
   else
      dE_start_step = sim%n_steps - diag%dE_collect_steps + 1
   end if

   !---------------------------------------------------------------------------
   ! メインループ
   !---------------------------------------------------------------------------
   write(*,'(A)') ''
   write(*,'(A)') 'Starting simulation...'

   call system_clock(count_start, count_rate)

   do istep = 1, sim%n_steps
      collect_dE = (istep >= dE_start_step)

      n_coll_cx = 0
      n_coll_el = 0
      n_coll_total = 0
      n_alive = 0
      weight_sum = 0.0d0

      !1ステップ分のスコアをリセット（reduction用���������カラー）
      s_cl_ei = 0.0d0; s_cl_cx = 0.0d0; s_cl_el = 0.0d0
      s_tl_ei = 0.0d0; s_tl_cx = 0.0d0; s_tl_el = 0.0d0
      s_tl_lookup_ei = 0.0d0; s_tl_lookup_cx = 0.0d0; s_tl_lookup_el = 0.0d0

      !$omp parallel do &
      !$omp   private(ip, coll_type_l, vx_i_l, vy_i_l, vz_i_l, &
      !$omp           v_rel_l, E_rel_l, delta_E_l, my_score, p_copy, &
      !$omp           time_remaining, step_dt, t_col, r, nu_max, segment_eff_time, &
      !$omp           is_collision_event, &
      !$omp           l_dE_hist_el, l_dE_hist_cx, ibin_l) &
      !$omp   reduction(+:n_alive, weight_sum, n_coll_cx, n_coll_el, &
      !$omp              n_coll_total, s_cl_cx, s_cl_el, s_cl_ei, &
      !$omp              s_tl_cx, s_tl_el, s_tl_ei, &
      !$omp              s_tl_lookup_cx, s_tl_lookup_el, s_tl_lookup_ei) &
      !$omp   schedule(static)
      do ip = 1, sim%n_particles
         if (.not. particles(ip)%alive) cycle
         n_alive = n_alive + 1
         weight_sum = weight_sum + particles(ip)%weight

         !スレッドローカルなスコア変数を使用（1パーティクル、1マクロステップ分）
         my_score%cl_ei = 0.0d0; my_score%cl_cx = 0.0d0; my_score%cl_el = 0.0d0
         my_score%tl_ei = 0.0d0; my_score%tl_cx = 0.0d0; my_score%tl_el = 0.0d0
         my_score%tl_lookup_ei = 0.0d0; my_score%tl_lookup_cx = 0.0d0; my_score%tl_lookup_el = 0.0d0

         allocate(l_dE_hist_el(diag%n_dE_bins), l_dE_hist_cx(diag%n_dE_bins))
         l_dE_hist_el = 0.0d0
         l_dE_hist_cx = 0.0d0

         time_remaining = sim%dt
         nu_max = plasma%n_i * sigma_v_max

         ! 仮想衝突までの距離(zincx)と累積距離(zint1)は粒子構造体に保持されているためここでは初期化しない

         !--- イベント駆動の内側ループ ---
         do while (time_remaining > 0.0d0 .and. particles(ip)%alive)

            ! 次の仮想衝突までの時間
            if (nu_max > 1.0d-30) then
               t_col = max((particles(ip)%zincx - particles(ip)%zint1) / nu_max, 0.0d0)
            else
               t_col = time_remaining * 2.0d0 ! 衝突なし
            end if

            is_collision_event = .false.
            if (t_col <= time_remaining) then
               step_dt = t_col
               is_collision_event = .true.
            else
               step_dt = time_remaining
            end if

            ! 累積距離の更新
            particles(ip)%zint1 = particles(ip)%zint1 + step_dt * nu_max

            ! 実衝突までの飛跡長を蓄積し、受理衝突時にまとめて flush する
            if (step_dt > 0.0d0) then
               segment_eff_time = compute_effective_track_time( &
                  particles(ip)%weight, ion_loss_rate, step_dt)
               particles(ip)%tl_pending_time = particles(ip)%tl_pending_time + step_dt
               particles(ip)%tl_pending_eff_time = particles(ip)%tl_pending_eff_time + &
                  segment_eff_time
            end if

            !--- 重み更新 (step_dt に基づいて計算) ---
            if (sim%enable_ei) then
               call update_weight(particles(ip), plasma, step_dt, sim%weight_min)
            end if

            !--- 位置更新 ---
            call advance_particle(particles(ip), step_dt)

            ! 残り時間の更新
            time_remaining = max(0.0d0, time_remaining - step_dt)

            ! 粒子が死亡した場合はループ脱出
            if (.not. particles(ip)%alive) then
               call flush_pending_track_scores(particles(ip), plasma, sim, my_score)
               exit
            end if

            !--- 仮想衝突イベントの評価 ---
            if (is_collision_event) then
               call evaluate_collision_event(particles(ip), plasma, sim, &
                  vx_i_l, vy_i_l, vz_i_l, &
                  v_rel_l, E_rel_l, coll_type_l)

               if (coll_type_l /= COLL_NONE) then
                  call flush_pending_track_scores(particles(ip), plasma, sim, my_score)
                  n_coll_total = n_coll_total + 1

                  if (coll_type_l == COLL_CX) then
                     delta_E_l = 0.0d0
                     ! 速度・エネルギー更新前にスコアリング (衝突前の高いエネルギーを用いる)
                     call score_collision_estimator(particles(ip), plasma, &
                        vx_i_l, vy_i_l, vz_i_l, &
                        v_rel_l, E_rel_l, coll_type_l, delta_E_l, &
                        sim%enable_ei, my_score)
                     call collision_cx(particles(ip), vx_i_l, vy_i_l, vz_i_l, &
                        delta_E_l)
                     n_coll_cx = n_coll_cx + 1

                     ! エネルギー変化量[eV]をヒストグラムに記録
                     if (collect_dE) then
                        ibin_l = int((delta_E_l * J_TO_EV - diag%dE_hist_min) / dE_bin_width) + 1
                        if (ibin_l >= 1 .and. ibin_l <= diag%n_dE_bins) then
                           l_dE_hist_cx(ibin_l) = l_dE_hist_cx(ibin_l) + particles(ip)%weight
                        end if
                     end if

                  else if (coll_type_l == COLL_EL) then
                     ! 弾性散乱では delta_E_l を得る必要があるため、粒子のコピーを使用して処理する
                     p_copy = particles(ip)
                     call collision_el(p_copy, vx_i_l, vy_i_l, vz_i_l, &
                        sim%use_isotropic, delta_E_l)
                     ! 元の粒子(更新前)を使ってスコアを記録
                     call score_collision_estimator(particles(ip), plasma, &
                        vx_i_l, vy_i_l, vz_i_l, &
                        v_rel_l, E_rel_l, coll_type_l, delta_E_l, &
                        sim%enable_ei, my_score)
                     ! 更新後の状態(と進んだRNG状態)を書き戻し
                     particles(ip) = p_copy
                     n_coll_el = n_coll_el + 1

                     ! エネルギー変化量[eV]をヒストグラムに記録
                     if (collect_dE) then
                        ibin_l = int((delta_E_l * J_TO_EV - diag%dE_hist_min) / dE_bin_width) + 1
                        if (ibin_l >= 1 .and. ibin_l <= diag%n_dE_bins) then
                           l_dE_hist_el(ibin_l) = l_dE_hist_el(ibin_l) + particles(ip)%weight
                        end if
                     end if
                  end if

                  ! 実際の衝突が発生した場合、サンプリング距離を更新
                  r = random_double(particles(ip)%rng)
                  particles(ip)%zincx = -log(max(r, 1.0d-30))
                  particles(ip)%zint1 = 0.0d0
               else
                  ! 空衝突の場合: イベント距離を再サンプリングして継続
                  r = random_double(particles(ip)%rng)
                  particles(ip)%zincx = -log(max(r, 1.0d-30))
                  particles(ip)%zint1 = 0.0d0
               end if
            end if

         end do
         !--- イベント駆動の内側ループ終了 ---

         !--- スレッドローカルスコアをreduction変数に加算 ---
         s_cl_cx = s_cl_cx + my_score%cl_cx
         s_cl_el = s_cl_el + my_score%cl_el
         s_cl_ei = s_cl_ei + my_score%cl_ei
         s_tl_cx = s_tl_cx + my_score%tl_cx
         s_tl_el = s_tl_el + my_score%tl_el
         s_tl_ei = s_tl_ei + my_score%tl_ei
         s_tl_lookup_cx = s_tl_lookup_cx + my_score%tl_lookup_cx
         s_tl_lookup_el = s_tl_lookup_el + my_score%tl_lookup_el
         s_tl_lookup_ei = s_tl_lookup_ei + my_score%tl_lookup_ei

         !$omp critical
         dE_hist_el(:) = dE_hist_el(:) + l_dE_hist_el(:)
         dE_hist_cx(:) = dE_hist_cx(:) + l_dE_hist_cx(:)
         !$omp end critical

         deallocate(l_dE_hist_el, l_dE_hist_cx)

      end do
      !$omp end parallel do

      if (istep == sim%n_steps) then
         !$omp parallel do &
         !$omp   private(ip, my_score) &
         !$omp   reduction(+:s_tl_cx, s_tl_el, s_tl_ei, &
         !$omp              s_tl_lookup_cx, s_tl_lookup_el, s_tl_lookup_ei) &
         !$omp   schedule(static)
         do ip = 1, sim%n_particles
            my_score%cl_ei = 0.0d0; my_score%cl_cx = 0.0d0; my_score%cl_el = 0.0d0
            my_score%tl_ei = 0.0d0; my_score%tl_cx = 0.0d0; my_score%tl_el = 0.0d0
            my_score%tl_lookup_ei = 0.0d0; my_score%tl_lookup_cx = 0.0d0
            my_score%tl_lookup_el = 0.0d0

            call flush_pending_track_scores(particles(ip), plasma, sim, my_score)

            s_tl_cx = s_tl_cx + my_score%tl_cx
            s_tl_el = s_tl_el + my_score%tl_el
            s_tl_ei = s_tl_ei + my_score%tl_ei
            s_tl_lookup_cx = s_tl_lookup_cx + my_score%tl_lookup_cx
            s_tl_lookup_el = s_tl_lookup_el + my_score%tl_lookup_el
            s_tl_lookup_ei = s_tl_lookup_ei + my_score%tl_lookup_ei
         end do
         !$omp end parallel do
      end if

      !--- reduction結果をstep_scoreに格納 ---
      step_score%cl_cx = s_cl_cx; step_score%cl_el = s_cl_el
      step_score%cl_ei = s_cl_ei
      step_score%tl_cx = s_tl_cx; step_score%tl_el = s_tl_el
      step_score%tl_ei = s_tl_ei
      step_score%tl_lookup_cx = s_tl_lookup_cx; step_score%tl_lookup_el = s_tl_lookup_el
      step_score%tl_lookup_ei = s_tl_lookup_ei

      !--- ステップスコアを累積スコアに加算 ---
      score%cl_ei = score%cl_ei + step_score%cl_ei
      score%cl_cx = score%cl_cx + step_score%cl_cx
      score%cl_el = score%cl_el + step_score%cl_el
      score%tl_ei = score%tl_ei + step_score%tl_ei
      score%tl_cx = score%tl_cx + step_score%tl_cx
      score%tl_el = score%tl_el + step_score%tl_el
      score%tl_lookup_ei = score%tl_lookup_ei + step_score%tl_lookup_ei
      score%tl_lookup_cx = score%tl_lookup_cx + step_score%tl_lookup_cx
      score%tl_lookup_el = score%tl_lookup_el + step_score%tl_lookup_el

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

   !--- 移行エネルギーヒストグラム出力（シミュレーション終了後に1回） ---
   dE_actual_steps = sim%n_steps - dE_start_step + 1
   dE_collect_time = dble(dE_actual_steps) * sim%dt
   call output_deltaE_histogram(sim%output_deltaE_hist, dE_hist_el, dE_hist_cx, &
      diag%n_dE_bins, diag%dE_hist_min, diag%dE_hist_max, &
      init_p%n_init, sim%n_particles, dE_collect_time)
   !後片付け
   deallocate(particles)
   deallocate(dE_hist_el, dE_hist_cx)

contains

   pure function compute_effective_track_time(weight, loss_rate, dt) result(eff_time)
      real(dp), intent(in) :: weight, loss_rate, dt
      real(dp) :: eff_time

      if (loss_rate > 1.0d-30) then
         eff_time = (weight / loss_rate) * (1.0d0 - exp(-loss_rate * dt))
      else
         eff_time = weight * dt
      end if
   end function compute_effective_track_time

   subroutine flush_pending_track_scores(p, plasma, sim, local_score)
      type(particle_t), intent(inout) :: p
      type(plasma_params), intent(in) :: plasma
      type(sim_params), intent(in)    :: sim
      type(score_data), intent(inout) :: local_score

      if (p%tl_pending_eff_time <= 0.0d0) then
         p%tl_pending_time = 0.0d0
         p%tl_pending_eff_time = 0.0d0
         return
      end if

      call score_track_length_estimator(p, plasma, p%tl_pending_eff_time, &
         sim%use_isotropic, sim%enable_cx, sim%enable_el, sim%enable_ei, local_score)
      call score_table_lookup_track_length(p, plasma, p%tl_pending_eff_time, &
         sim%enable_cx, sim%enable_el, sim%enable_ei, local_score)

      p%tl_pending_time = 0.0d0
      p%tl_pending_eff_time = 0.0d0
   end subroutine flush_pending_track_scores

end program monte_carlo_3d3v_natl
