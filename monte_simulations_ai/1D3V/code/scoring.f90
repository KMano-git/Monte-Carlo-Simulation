!===============================================================================
! Module: scoring
! スコアリング手法（CL法、TL法）- CX/EL/EI別エネルギー記録対応
! 電離は電子温度ベースのレート係数を使用
!===============================================================================
module scoring
   use constants, only: dp, EV_TO_J, M_D, J_TO_EV
   use types, only: particle, score_grid, sim_params, plasma_params
   use cross_sections, only: sigma_cx, sigma_el, ionization_rate_coeff
   implicit none

contains

   !-----------------------------------------------------------------------------
   ! スコアリンググリッドの初期化
   !-----------------------------------------------------------------------------
   subroutine init_score_grid(grid, n_grid)
      type(score_grid), intent(out) :: grid
      integer, intent(in) :: n_grid

      ! CL法
      allocate(grid%cl_energy(n_grid))
      allocate(grid%cl_energy_cx(n_grid))
      allocate(grid%cl_energy_el(n_grid))
      allocate(grid%cl_energy_ei(n_grid))
      allocate(grid%cl_count(n_grid))
      allocate(grid%cl_count_cx(n_grid))
      allocate(grid%cl_count_el(n_grid))
      allocate(grid%cl_count_ei(n_grid))
      ! TL法
      allocate(grid%tl_energy(n_grid))
      allocate(grid%tl_energy_cx(n_grid))
      allocate(grid%tl_energy_el(n_grid))
      allocate(grid%tl_energy_ei(n_grid))
      allocate(grid%tl_count(n_grid))
      ! 密度
      allocate(grid%density(n_grid))

      grid%cl_energy = 0.0d0
      grid%cl_energy_cx = 0.0d0
      grid%cl_energy_el = 0.0d0
      grid%cl_energy_ei = 0.0d0
      grid%cl_count = 0.0d0
      grid%cl_count_cx = 0.0d0
      grid%cl_count_el = 0.0d0
      grid%cl_count_ei = 0.0d0
      grid%tl_energy = 0.0d0
      grid%tl_energy_cx = 0.0d0
      grid%tl_energy_el = 0.0d0
      grid%tl_energy_ei = 0.0d0
      grid%tl_count = 0.0d0
      grid%density = 0.0d0

   end subroutine init_score_grid

   !-----------------------------------------------------------------------------
   ! スコアリンググリッドの解放
   !-----------------------------------------------------------------------------
   subroutine free_score_grid(grid)
      type(score_grid), intent(inout) :: grid

      if (allocated(grid%cl_energy)) deallocate(grid%cl_energy)
      if (allocated(grid%cl_energy_cx)) deallocate(grid%cl_energy_cx)
      if (allocated(grid%cl_energy_el)) deallocate(grid%cl_energy_el)
      if (allocated(grid%cl_energy_ei)) deallocate(grid%cl_energy_ei)
      if (allocated(grid%cl_count)) deallocate(grid%cl_count)
      if (allocated(grid%cl_count_cx)) deallocate(grid%cl_count_cx)
      if (allocated(grid%cl_count_el)) deallocate(grid%cl_count_el)
      if (allocated(grid%cl_count_ei)) deallocate(grid%cl_count_ei)
      if (allocated(grid%tl_energy)) deallocate(grid%tl_energy)
      if (allocated(grid%tl_energy_cx)) deallocate(grid%tl_energy_cx)
      if (allocated(grid%tl_energy_el)) deallocate(grid%tl_energy_el)
      if (allocated(grid%tl_energy_ei)) deallocate(grid%tl_energy_ei)
      if (allocated(grid%tl_count)) deallocate(grid%tl_count)
      if (allocated(grid%density)) deallocate(grid%density)

   end subroutine free_score_grid

   !-----------------------------------------------------------------------------
   ! グリッドインデックスの計算
   !-----------------------------------------------------------------------------
   function get_grid_index(x, params) result(ix)
      real(dp), intent(in) :: x
      type(sim_params), intent(in) :: params
      integer :: ix

      ix = int((x - params%x_min) / params%dx) + 1
      ix = max(1, min(ix, params%n_grid))

   end function get_grid_index

   !-----------------------------------------------------------------------------
   ! Collision Estimator (CL法)
   ! 衝突発生時に呼び出し、衝突タイプ別にエネルギー変化を記録
   !-----------------------------------------------------------------------------
   subroutine score_cl(grid, p, delta_E, coll_type, params)
      type(score_grid), intent(inout) :: grid
      type(particle), intent(in) :: p
      real(dp), intent(in) :: delta_E        ! エネルギー変化 [eV]
      integer, intent(in)  :: coll_type      ! 衝突タイプ (1:CX, 2:EL, 3:EI)
      type(sim_params), intent(in) :: params

      integer :: ix

      ix = get_grid_index(p%x, params)

      ! 全体
      grid%cl_energy(ix) = grid%cl_energy(ix) + p%weight * delta_E
      grid%cl_count(ix) = grid%cl_count(ix) + p%weight

      ! 衝突タイプ別
      select case (coll_type)
       case (1)  ! CX
         grid%cl_energy_cx(ix) = grid%cl_energy_cx(ix) + p%weight * delta_E
         grid%cl_count_cx(ix) = grid%cl_count_cx(ix) + p%weight
       case (2)  ! EL
         grid%cl_energy_el(ix) = grid%cl_energy_el(ix) + p%weight * delta_E
         grid%cl_count_el(ix) = grid%cl_count_el(ix) + p%weight
       case (3)  ! EI
         grid%cl_energy_ei(ix) = grid%cl_energy_ei(ix) + p%weight * delta_E
         grid%cl_count_ei(ix) = grid%cl_count_ei(ix) + p%weight
      end select

   end subroutine score_cl

   !-----------------------------------------------------------------------------
   ! Track-Length Estimator (TL法)
   ! 毎ステップ呼び出し、期待値を積分（CX/EL/EI別）
   ! 電離は電子温度ベースのレート係数を使用
   !-----------------------------------------------------------------------------
   subroutine score_tl(grid, p, params, plasma, dt)
      type(score_grid), intent(inout) :: grid
      type(particle), intent(in) :: p
      type(sim_params), intent(in) :: params
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in) :: dt

      integer :: ix
      real(dp) :: v, E_particle, E_bg
      real(dp) :: sig_cx, sig_el
      real(dp) :: rate_cx, rate_el, rate_ei
      real(dp) :: delta_E_cx, delta_E_el, delta_E_ei
      real(dp) :: contrib_cx, contrib_el, contrib_ei

      ! 粒子速度と運動エネルギー
      v = sqrt(p%vx**2 + p%vy**2 + p%vz**2)
      E_particle = 0.5d0 * M_D * v**2 * J_TO_EV  ! [eV]
      E_bg = 1.5d0 * plasma%T_bg  ! 背景イオンの平均エネルギー [eV]

      ! CX/EL断面積（イオン-中性粒子衝突）
      sig_cx = sigma_cx(E_particle)
      sig_el = sigma_el(E_particle)

      ! イオン衝突レート（CX, EL）
      rate_cx = plasma%n_bg * sig_cx * v
      rate_el = plasma%n_bg * sig_el * v

      ! 電離レート（電子衝突電離）
      ! <σv>_ei は電子温度で決まる。レート = n_e * <σv>_ei
      rate_ei = plasma%n_e * ionization_rate_coeff(plasma%T_e)

      ! 期待されるエネルギー変化
      delta_E_cx = E_bg - E_particle       ! CX: 背景とエネルギー交換
      delta_E_el = 0.5d0 * (E_bg - E_particle)  ! EL: 部分的エネルギー移行
      delta_E_ei = -13.6d0                 ! EI: 電離エネルギー損失

      ! 寄与の計算
      contrib_cx = rate_cx * delta_E_cx * dt
      contrib_el = rate_el * delta_E_el * dt
      contrib_ei = rate_ei * delta_E_ei * dt

      ix = get_grid_index(p%x, params)

      ! 全体
      grid%tl_energy(ix) = grid%tl_energy(ix) + p%weight * (contrib_cx + contrib_el + contrib_ei)
      grid%tl_count(ix) = grid%tl_count(ix) + p%weight * dt

      ! 反応別
      grid%tl_energy_cx(ix) = grid%tl_energy_cx(ix) + p%weight * contrib_cx
      grid%tl_energy_el(ix) = grid%tl_energy_el(ix) + p%weight * contrib_el
      grid%tl_energy_ei(ix) = grid%tl_energy_ei(ix) + p%weight * contrib_ei

      ! 密度（滞在時間積算）
      grid%density(ix) = grid%density(ix) + p%weight * dt

   end subroutine score_tl

end module scoring
