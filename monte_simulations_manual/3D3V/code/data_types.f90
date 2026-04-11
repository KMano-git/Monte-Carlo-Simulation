!===============================================================================
! Module: data_types
! 粒子データ・シミュレーションパラメータ・スコアデータ定義
!===============================================================================
module data_types
   use constants, only : dp
   implicit none

   !---------------------------------------------------------------------------
   ! 乱数生成器の状態（xoshiro256+用、64bit×4）
   !---------------------------------------------------------------------------
   type :: rng_state
      integer(8) :: s(4)
   end type rng_state

   !---------------------------------------------------------------------------
   ! 粒子データ構造体（3D3V + 重み + RNG）
   !---------------------------------------------------------------------------
   type :: particle_t
      real(dp) :: x, y, z       !位置 [m]
      real(dp) :: vx, vy, vz    !速度 [m/s]
      real(dp) :: weight         !粒子重み（Non-Analogで電離により減衰）
      logical  :: alive          !生存フラグ
      real(dp) :: zincx          !次の仮想衝突までのランダムな距離
      real(dp) :: zint1          !現在までの累積距離
      real(dp) :: tl_pending_time      !未flushの飛行時間 [s]
      real(dp) :: tl_pending_eff_time  !未flushの有効重み時間 [s]
      type(rng_state) :: rng    !粒子固有の乱数状態
   end type particle_t

   !---------------------------------------------------------------------------
   ! シミュレーションパラメータ型
   !---------------------------------------------------------------------------
   type :: sim_params
      integer  :: n_particles    !粒子数
      integer  :: n_steps        !ステップ数
      real(dp) :: dt             !時間刻み [s]
      integer  :: seed           !乱数シード
      integer  :: tl_el_inner_samples !EL inner multi-sampling のサンプル数
      logical  :: enable_cx      !荷電交換 On/Off
      logical  :: enable_el      !弾性散乱 On/Off
      logical  :: enable_ei      !電離 On/Off
      logical  :: use_isotropic  !等方散乱フラグ
      real(dp) :: weight_min     !最小重み閾値（ロシアンルーレット用）
      character(len=256) :: cdf_file       !CDFファイル名
      character(len=256) :: tl_el_table_file !EL pre-tabulated TL 用テーブル
      character(len=256) :: output_ntscrg  !スコアリング出力ファイル名
      character(len=256) :: output_hist    !エネルギーヒストグラム出力ファイル名
      character(len=256) :: output_deltaE_hist !エネルギー変化ヒストグラム出力ファイル名
   end type sim_params

   !---------------------------------------------------------------------------
   ! プラズマパラメータ型
   !---------------------------------------------------------------------------
   type :: plasma_params
      real(dp) :: n_i            !背景イオン密度 [m^-3]
      real(dp) :: T_i            !背景イオン温度 [eV]
      real(dp) :: n_e            !背景電子密度 [m^-3]
      real(dp) :: T_e            !背景電子温度 [eV]
      real(dp) :: u_x, u_y, u_z  !背景イオン流速 [m/s]
   end type plasma_params

   !---------------------------------------------------------------------------
   ! 粒子初期条件型
   !---------------------------------------------------------------------------
   type :: init_params
      real(dp) :: n_init         !粒子密度 [m^-3]
      real(dp) :: E_init         !初期エネルギー [eV]
      real(dp) :: T_init         !初期温度 [eV]
      integer  :: init_mode      !粒子初期分布(0: mono, 1: maxwell)
   end type init_params

   !---------------------------------------------------------------------------
   ! 診断パラメータ型
   !---------------------------------------------------------------------------
   type :: diag_params
      integer  :: output_interval  !出力間隔
      integer  :: n_hist_bins      !ヒストグラム分割数
      real(dp) :: E_hist_min       !エネルギー最小値 [eV]
      real(dp) :: E_hist_max       !エネルギー最大値 [eV]
      integer  :: hist_timing(6)   !ヒストグラム出力タイミング
      ! エネルギー移行量ヒストグラム用
      integer  :: n_dE_bins
      real(dp) :: dE_hist_min      ![eV]
      real(dp) :: dE_hist_max      ![eV]
      integer  :: dE_collect_steps !集計ステップ��(0=全ステップ)
   end type diag_params

   !---------------------------------------------------------------------------
   ! スコアリングデータ型（CL/TL同時保持）
   !---------------------------------------------------------------------------
   type :: score_data
      ! Analog Estimator [J]
      real(dp) :: a_el            !弾性散乱寄与
      ! Collision Estimator [J]
      real(dp) :: cl_ei           !電離寄与
      real(dp) :: cl_cx           !荷電交換寄与
      real(dp) :: cl_el           !弾性散乱寄与
      real(dp) :: cl_inner_el     !inner multi-sampling CL
      ! Track-Length Estimator [J]
      real(dp) :: tr_ei           !電離寄与
      real(dp) :: tr_cx           !荷電交換寄与
      real(dp) :: tr_el           !弾性散乱寄与 (lookup)
      ! Elastic-only comparison add-ons [J]
      real(dp) :: tr_inner_el     !inner multi-sampling TR
      real(dp) :: tr_pretab_el    !pre-tabulated TR
   end type score_data

end module data_types
