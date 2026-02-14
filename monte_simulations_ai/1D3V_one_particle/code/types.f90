!===============================================================================
! Module: types
! データ型定義（粒子、グリッド等）
!===============================================================================
module types
   use constants, only: dp
   implicit none

   !-----------------------------------------------------------------------------
   ! 粒子型
   !-----------------------------------------------------------------------------
   type :: particle
      real(dp) :: x           ! 位置 [m]
      real(dp) :: vx, vy, vz  ! 速度成分 [m/s]
      real(dp) :: weight      ! 統計的重み
      logical  :: alive       ! 生存フラグ
   end type particle

   !-----------------------------------------------------------------------------
   ! シミュレーションパラメータ型
   !-----------------------------------------------------------------------------
   type :: sim_params
      integer  :: n_particles     ! テスト粒子数
      integer  :: n_steps         ! ステップ数
      integer  :: n_grid          ! グリッド数
      real(dp) :: dt              ! 時間刻み [s]
      real(dp) :: x_min           ! 領域下限 [m]
      real(dp) :: x_max           ! 領域上限 [m]
      real(dp) :: dx              ! グリッド幅 [m]
      real(dp) :: gamma_in        ! 入射フラックス [m^-2 s^-1]
   end type sim_params

   !-----------------------------------------------------------------------------
   ! プラズマパラメータ型
   !-----------------------------------------------------------------------------
   type :: plasma_params
      real(dp) :: n_bg      ! 背景イオン密度 [m^-3]
      real(dp) :: n_e       ! 電子密度 [m^-3] (通常 n_bg と同じ)
      real(dp) :: T_bg      ! 背景イオン温度 [eV]
      real(dp) :: T_e       ! 電子温度 [eV]
   end type plasma_params

   !-----------------------------------------------------------------------------
   ! 粒子初期条件型
   !-----------------------------------------------------------------------------
   type :: init_params
      real(dp) :: x_init    ! 初期位置 [m]
      real(dp) :: E_init    ! 初期エネルギー [eV]
      real(dp) :: T_init    ! 初期温度 [eV] (Maxwell分布用)
   end type init_params

   !-----------------------------------------------------------------------------
   ! スコアリンググリッド型
   !-----------------------------------------------------------------------------
   type :: score_grid
      ! CL法（衝突発生時）
      real(dp), allocatable :: cl_energy(:)      ! 全エネルギー付与 [eV]
      real(dp), allocatable :: cl_energy_cx(:)   ! CXエネルギー付与 [eV]
      real(dp), allocatable :: cl_energy_el(:)   ! ELエネルギー付与 [eV]
      real(dp), allocatable :: cl_energy_ei(:)   ! EIエネルギー付与 [eV]
      real(dp), allocatable :: cl_count(:)       ! 衝突カウント
      real(dp), allocatable :: cl_count_cx(:)    ! CX衝突カウント
      real(dp), allocatable :: cl_count_el(:)    ! EL衝突カウント
      real(dp), allocatable :: cl_count_ei(:)    ! EI衝突カウント
      ! TL法（毎ステップ期待値）
      real(dp), allocatable :: tl_energy(:)      ! 全エネルギー付与 [eV]
      real(dp), allocatable :: tl_energy_cx(:)   ! CXエネルギー付与 [eV]
      real(dp), allocatable :: tl_energy_el(:)   ! ELエネルギー付与 [eV]
      real(dp), allocatable :: tl_energy_ei(:)   ! EIエネルギー付与 [eV]
      real(dp), allocatable :: tl_count(:)       ! カウント
      ! 密度追跡
      real(dp), allocatable :: density(:)        ! 中性粒子滞在時間（密度相当）
   end type score_grid

end module types
