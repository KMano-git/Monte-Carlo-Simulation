!===============================================================================
! Module: particle_data
! 粒子データ定義
!===============================================================================
module particle_data
   use constants, only : dp
   implicit none

   !粒子データ構造体
   type :: particle_t
      real(dp) :: vx, vy, vz !速度 [m/s]
      logical  :: alive      !生存フラグ（電離で消滅時にfalse）
   end type particle_t

   !シミュレーションパラメータ型
   type :: sim_params
      integer  :: n_particles !粒子数
      integer  :: n_steps     !ステップ数
      real(dp) :: dt          !時間刻み [s]
      integer  :: seed        !乱数シード
      logical  :: enable_cx   !荷電交換 On/Off
      logical  :: enable_el   !弾性散乱 On/Off
      logical  :: enable_ei   !電離 On/Off ! 当分はオフで固定
      logical  :: use_isotropic !等方散乱フラグ
      character(len=256) :: cdf_file !CDFファイル名
      character(len=256) :: output_ntscrg !エネルギー増減出力ファイル名
      character(len=256) :: output_hist !エネルギーヒストグラム出力ファイル名
      character(len=256) :: output_log !ログ出力ファイル名
   end type sim_params

   !プラズマパラメータ型
   type :: plasma_params
      real(dp) :: n_i !背景イオン密度 [m^-3]
      real(dp) :: T_i !背景イオン温度 [eV]
      real(dp) :: n_e !背景電子密度 [m^-3]
      real(dp) :: T_e !背景電子温度 [eV]
      real(dp) :: u_x, u_y, u_z !背景イオン流速
   end type plasma_params

   !粒子初期条件型
   type :: init_params
      real(dp) :: n_init !粒子密度 [m^-3]
      real(dp) :: E_init !初期エネルギー [eV]
      real(dp) :: T_init !初期温度 [eV]
      integer  :: init_mode !粒子初期分布(0: mono, 1: maxwell)
   end type init_params

   !診断パラメータ型
   type :: diag_params
      integer  :: output_interval !出力間隔
      integer  :: n_hist_bins !ヒストグラム分割数
      real(dp) :: E_hist_min !エネルギー最小値 [eV]
      real(dp) :: E_hist_max !エネルギー最大値 [eV]
      integer  :: hist_timing(5) !ヒストグラム出力タイミング
   end type diag_params

end module particle_data
