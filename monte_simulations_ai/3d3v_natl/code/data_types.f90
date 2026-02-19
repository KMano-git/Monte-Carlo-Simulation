!===============================================================================
! Module: data_types
! 粒子データ・シミュレーションパラメータ・スコアデータ定義
!===============================================================================
module data_types
   use constants, only : dp
   implicit none

   !---------------------------------------------------------------------------
   ! 粒子データ構造体（3D3V + 重み）
   !---------------------------------------------------------------------------
   type :: particle_t
      real(dp) :: x, y, z       !位置 [m]
      real(dp) :: vx, vy, vz    !速度 [m/s]
      real(dp) :: weight         !粒子重み（Non-Analogで電離により減衰）
      logical  :: alive          !生存フラグ
   end type particle_t

   !---------------------------------------------------------------------------
   ! シミュレーションパラメータ型
   !---------------------------------------------------------------------------
   type :: sim_params
      integer  :: n_particles    !粒子数
      integer  :: n_steps        !ステップ数
      real(dp) :: dt             !時間刻み [s]
      integer  :: seed           !乱数シード
      logical  :: enable_cx      !荷電交換 On/Off
      logical  :: enable_el      !弾性散乱 On/Off
      logical  :: enable_ei      !電離 On/Off
      logical  :: use_isotropic  !等方散乱フラグ
      real(dp) :: weight_min     !最小重み閾値（ロシアンルーレット用）
      character(len=256) :: cdf_file       !CDFファイル名
      character(len=256) :: output_ntscrg  !スコアリング出力ファイル名
      character(len=256) :: output_hist    !エネルギーヒストグラム出力ファイル名
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
   end type diag_params

   !---------------------------------------------------------------------------
   ! スコアリングデータ型（CL/TL同時保持）
   !---------------------------------------------------------------------------
   type :: score_data
      ! Collision Estimator [J]
      real(dp) :: cl_ei           !電離寄与
      real(dp) :: cl_cx           !荷電交換寄与
      real(dp) :: cl_el           !弾性散乱寄与
      ! Track-Length Estimator [J]
      real(dp) :: tl_ei           !電離寄与
      real(dp) :: tl_cx           !荷電交換寄与
      real(dp) :: tl_el           !弾性散乱寄与
   end type score_data

end module data_types
