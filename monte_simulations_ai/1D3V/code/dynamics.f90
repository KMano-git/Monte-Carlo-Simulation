!===============================================================================
! Module: dynamics
! 粒子力学: null-collision 法による衝突評価と衝突処理
!===============================================================================
module dynamics
   use constants, only: dp, M_D_kg, J_TO_EV, PI
   use data_types, only: particle_t, plasma_params, sim_params
   use random_utils, only: sample_maxwell_velocity_ion, random_double
   use cross_sections, only: sigma_cx, sigma_el, sigma_v_max, ionization_rate_coeff
   use cdf_reader, only: sample_scattering_angle
   implicit none

   private
   public :: advance_particle, evaluate_collision_event
   public :: collision_cx, collision_el
   public :: update_weight
   public :: COLL_NONE, COLL_CX, COLL_EL

   integer, parameter :: COLL_NONE = 0
   integer, parameter :: COLL_CX   = 1
   integer, parameter :: COLL_EL   = 2

contains

   subroutine advance_particle(p, dt)
      type(particle_t), intent(inout) :: p
      real(dp), intent(in) :: dt

      p%x = p%x + p%vx * dt
      p%y = p%y + p%vy * dt
      p%z = p%z + p%vz * dt
   end subroutine advance_particle

   subroutine update_weight(p, plasma, dt, weight_min)
      type(particle_t), intent(inout) :: p
      type(plasma_params), intent(in) :: plasma
      real(dp), intent(in) :: dt
      real(dp), intent(in) :: weight_min

      real(dp) :: R_a, P_survive

      R_a = plasma%n_e * ionization_rate_coeff(plasma%T_e)
      if (R_a > 1.0d-30) p%weight = p%weight * exp(-R_a * dt)

      if (p%weight < weight_min) then
         P_survive = 0.1d0
         if (random_double(p%rng) <= P_survive) then
            p%weight = p%weight / P_survive
         else
            p%alive = .false.
         end if
      end if
   end subroutine update_weight

   subroutine evaluate_collision_event(p, plasma, sim, vx_i, vy_i, vz_i, v_rel, E_rel, coll_type)
      type(particle_t), intent(inout) :: p
      type(plasma_params), intent(in) :: plasma
      type(sim_params), intent(in) :: sim
      real(dp), intent(out) :: vx_i, vy_i, vz_i
      real(dp), intent(out) :: v_rel
      real(dp), intent(out) :: E_rel
      integer, intent(out) :: coll_type

      real(dp) :: sig_cx_val, sig_el_val
      real(dp) :: sig_cx_eff, sig_el_eff, sig_s_eff
      real(dp) :: P_accept
      real(dp) :: r_accept, r_type

      coll_type = COLL_NONE
      vx_i = 0.0d0
      vy_i = 0.0d0
      vz_i = 0.0d0
      v_rel = 0.0d0
      E_rel = 0.0d0

      if (sigma_v_max < 1.0d-30) return
      if (.not. sim%enable_cx .and. .not. sim%enable_el) return

      call sample_maxwell_velocity_ion(p%rng, plasma%T_i, plasma, vx_i, vy_i, vz_i)

      v_rel = sqrt((p%vx - vx_i)**2 + (p%vy - vy_i)**2 + (p%vz - vz_i)**2)
      E_rel = 0.25d0 * M_D_kg * v_rel * v_rel * J_TO_EV

      sig_cx_val = sigma_cx(E_rel)
      sig_el_val = sigma_el(E_rel)

      if (sim%enable_cx) then
         sig_cx_eff = sig_cx_val
      else
         sig_cx_eff = 0.0d0
      end if

      if (sim%enable_el) then
         sig_el_eff = sig_el_val
      else
         sig_el_eff = 0.0d0
      end if

      sig_s_eff = sig_cx_eff + sig_el_eff
      if (sig_s_eff <= 1.0d-30) return

      P_accept = sig_s_eff * v_rel / sigma_v_max
      if (P_accept > 1.0d0) P_accept = 1.0d0

      r_accept = random_double(p%rng)
      if (r_accept > P_accept) return

      r_type = random_double(p%rng)
      if (r_type < sig_cx_eff / sig_s_eff) then
         coll_type = COLL_CX
      else
         coll_type = COLL_EL
      end if
   end subroutine evaluate_collision_event

   subroutine collision_cx(p, vx_i, vy_i, vz_i, delta_E)
      type(particle_t), intent(inout) :: p
      real(dp), intent(in)  :: vx_i, vy_i, vz_i
      real(dp), intent(out) :: delta_E

      real(dp) :: E_old, E_new

      E_old = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)

      p%vx = vx_i
      p%vy = vy_i
      p%vz = vz_i

      E_new = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)
      delta_E = E_new - E_old
   end subroutine collision_cx

   subroutine collision_el(p, vx_i, vy_i, vz_i, use_isotropic, delta_E)
      type(particle_t), intent(inout) :: p
      real(dp), intent(in)  :: vx_i, vy_i, vz_i
      logical, intent(in)   :: use_isotropic
      real(dp), intent(out) :: delta_E

      real(dp) :: E_old, E_new
      real(dp) :: vx_g, vy_g, vz_g
      real(dp) :: ux, uy, uz, u_mag
      real(dp) :: ux_n, uy_n, uz_n
      real(dp) :: chi, phi
      real(dp) :: r_chi, E_rel

      E_old = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)

      vx_g = 0.5d0 * (p%vx + vx_i)
      vy_g = 0.5d0 * (p%vy + vy_i)
      vz_g = 0.5d0 * (p%vz + vz_i)

      ux = p%vx - vx_i
      uy = p%vy - vy_i
      uz = p%vz - vz_i
      u_mag = sqrt(ux * ux + uy * uy + uz * uz)

      if (u_mag < 1.0d-30) then
         delta_E = 0.0d0
         return
      end if

      if (use_isotropic) then
         r_chi = random_double(p%rng)
         chi = acos(1.0d0 - 2.0d0 * r_chi)
      else
         E_rel = 0.25d0 * M_D_kg * u_mag * u_mag * J_TO_EV
         r_chi = random_double(p%rng)
         chi = sample_scattering_angle(0.5d0 * E_rel, r_chi)
      end if

      phi = 2.0d0 * PI * random_double(p%rng)
      call rotate_vector(ux, uy, uz, chi, phi, ux_n, uy_n, uz_n)

      p%vx = vx_g + 0.5d0 * ux_n
      p%vy = vy_g + 0.5d0 * uy_n
      p%vz = vz_g + 0.5d0 * uz_n

      E_new = 0.5d0 * M_D_kg * (p%vx**2 + p%vy**2 + p%vz**2)
      delta_E = E_new - E_old
   end subroutine collision_el

   subroutine rotate_vector(ux, uy, uz, chi, phi, ux_n, uy_n, uz_n)
      real(dp), intent(in)  :: ux, uy, uz
      real(dp), intent(in)  :: chi, phi
      real(dp), intent(out) :: ux_n, uy_n, uz_n

      real(dp) :: u_mag, u_perp
      real(dp) :: cos_chi, sin_chi, cos_phi, sin_phi
      real(dp) :: nx, ny, nz

      u_mag = sqrt(ux * ux + uy * uy + uz * uz)
      if (u_mag < 1.0d-30) then
         ux_n = ux
         uy_n = uy
         uz_n = uz
         return
      end if

      cos_chi = cos(chi)
      sin_chi = sin(chi)
      cos_phi = cos(phi)
      sin_phi = sin(phi)

      nx = ux / u_mag
      ny = uy / u_mag
      nz = uz / u_mag

      u_perp = sqrt(nx * nx + ny * ny)
      if (u_perp > 1.0d-10) then
         ux_n = u_mag * (nx * cos_chi + &
            (nx * nz * cos_phi - ny * sin_phi) * sin_chi / u_perp)
         uy_n = u_mag * (ny * cos_chi + &
            (ny * nz * cos_phi + nx * sin_phi) * sin_chi / u_perp)
         uz_n = u_mag * (nz * cos_chi - u_perp * cos_phi * sin_chi)
      else
         ux_n = u_mag * sin_chi * cos_phi
         uy_n = u_mag * sin_chi * sin_phi
         uz_n = u_mag * cos_chi * sign(1.0d0, nz)
      end if
   end subroutine rotate_vector

end module dynamics
