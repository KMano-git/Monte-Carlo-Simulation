!***********************************************************************
      subroutine thermal_force_para_by_q_para(coulomb_logarithm,
     >     mi, mz, ei, ez, ti, vz_para, vz_perp, ui, q_para, fth_para)
!  We use the assumption [vz// ~ vz] to evalute Fth_//.
! If vz_perp >> vz_//, the coef_block_two could change...
! (we suppose that this change may not be so big...)
!***********************************************************************
      use cphcns, only : c_eps_zero, c_ten_pi_sqrt_pi, cev
      implicit none
!::In/Out
      real*8, intent(in):: coulomb_logarithm
      real*8, intent(in):: mi, mz     !  [kg]
      real*8, intent(in):: ei, ez     !  [C]
      real*8, intent(in):: ti         !  [eV]
      real*8, intent(in):: vz_para, vz_perp, ui     ! [m/s]
      real*8, intent(in):: q_para     ! [joule/(m2.s)]
      real*8, intent(out):: fth_para  ! [N]

!::local variables
      integer file_n6
      real*8 sub_block_one, sub_block_two
      real*8 sub_block_three, sub_block_four
      real*8 coef_block_one, coef_block_two
      real*8  ui_para
      real*8 vz_tilde_para, vz_tilde_squared
      real*8 sqrt_mi_on_ti, inv_ti_joule
!::clear
!
!::start
!     zero-check
      if(ti.le.0.d0) then
         fth_para = 0.d0

         file_n6 = 6
!
         write(file_n6,*)
         write(file_n6,*) "Subroutine thermal_force_para_by_q_para,"
         write(file_n6,*) "  warning: Ti = 0[eV]"
         write(file_n6,*) "If this line appears many times, check."
         return
      endif

      sqrt_mi_on_ti = sqrt(mi/(2.d0*cev*ti) )
      inv_ti_joule = 1.d0 / (cev*ti)

      ui_para = ui    ! Assumption or Approximation
      vz_tilde_para = sqrt_mi_on_ti * (vz_para - ui_para)

      vz_tilde_squared = (vz_perp**2.d0 + (vz_para-ui_para)**2.d0)
     >     *(sqrt_mi_on_ti **2.d0)

      sub_block_one = coulomb_logarithm / c_ten_pi_sqrt_pi
      sub_block_two = ( ez * ei / c_eps_zero  )**2.d0
      sub_block_three = 1.d0 + ( mi/mz  )
      sub_block_four = sqrt_mi_on_ti * inv_ti_joule * inv_ti_joule

      coef_block_one = sub_block_one * sub_block_two
     >     * sub_block_three * sub_block_four

      coef_block_two = exp(-vz_tilde_squared)
     >     * (1.d0 - 2.d0 * vz_tilde_para * vz_tilde_para )

      fth_para = -coef_block_one * coef_block_two * q_para
!
!     cf. note page imp8  Fth_Reiser = Fth_IMPMC_Shimizu_ver
!                                    = Fth_JCP / 1.248
!     fth_para_reiser = fth_para / 1.248d0
!
! test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$$$      write(6,*)
!$$$      write(6,*) "test; subroutine thermal_force_para..."
!$$$      write(6,*) "coulomb_logarithm, mi, mz, ei, ez"
!$$$      write(6,'(10(E16.4e2))') coulomb_logarithm, mi, mz, ei, ez
!$$$      write(6,*) "ti, vz_para, vz_perp, ui, q_para"
!$$$      write(6,'(10(E16.4e2))') ti, vz_para, vz_perp, ui, q_para
!$$$      write(6,*) "fth_para"
!$$$      write(6,'(10(E16.4e2))') fth_para
!$$$      write(6,*)
!$$$      write(6,*) "vth_i, vz_tilde_para,  vz_tilde_squared"
!$$$      write(6,'(10(E16.4e2))') 1.d0/sqrt_mi_on_ti, vz_tilde_para,
!$$$     >     vz_tilde_squared
!$$$      write(6,*) "coef_block_one, two"
!$$$      write(6,'(10(E16.4e2))') coef_block_one,
!$$$     >     coef_block_two
!$$$
!$$$      write(6,*) "sub_block_one, two, three, four"
!$$$      write(6,'(10(E16.4e2))') sub_block_one,
!$$$     >  sub_block_two, sub_block_three, sub_block_four
!$$$
!$$$      write(6,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
      end subroutine
!!!
!!!
!!!
!***********************************************************************
      subroutine friction_force_theoretical_value(coulomb_logarithm,
     >     mi, mz, ei, ez, ti, ni, vz_para, vz_perp, ui, friction_para)
!***********************************************************************
      use cphcns, only : c_eps_zero, cev, cpi
      implicit none
!::In/Out
      real*8, intent(in):: coulomb_logarithm
      real*8, intent(in):: mi, mz     !  [kg]
      real*8, intent(in):: ei, ez     !  [C]
      real*8, intent(in):: ti         !  [eV]
      real*8, intent(in):: ni         !  [m^(-3)]
      real*8, intent(in):: vz_para, vz_perp, ui     ! [m/s]
      real*8, intent(out):: friction_para  ! [N]

!::local variables
      integer file_n6
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 sub_block_one, sub_block_two
!ik   real*8 sub_block_three, sub_block_four
      real*8 sub_block_two
      real*8 sub_block_three
      real*8 friction_coef_one, friction_coef_two
      real*8  ui_para
      real*8 vz_tilde_para, vz_tilde_squared
      real*8 vz_tilde_abs
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 sqrt_mi_on_ti, inv_ti_joule
      real*8 sqrt_mi_on_ti
      real*8 phi, phi_prime
!::clear
!
!::start
!     zero-check
      if(ti.le.0.d0) then
         friction_para = 0.d0

         file_n6 = 6
!
         write(file_n6,*)
         write(file_n6,*) "Subroutine thermal_force_para_by_q_para,"
         write(file_n6,*) "  warning: Ti = 0[eV]"
         write(file_n6,*) "If this line appears many times, check."
         return
      endif

      sqrt_mi_on_ti = sqrt(mi/(2.d0*cev*ti) )
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   inv_ti_joule = 1.d0 / (cev*ti)

      ui_para = ui    ! Assumption or Approximation
      vz_tilde_para = sqrt_mi_on_ti * (vz_para - ui_para)

      vz_tilde_squared = (vz_perp**2.d0 + (vz_para-ui_para)**2.d0)
     >     *(sqrt_mi_on_ti **2.d0)

! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   sub_block_one = coulomb_logarithm / c_ten_pi_sqrt_pi
      sub_block_two = ( ez * ei / c_eps_zero  )**2.d0
      sub_block_three = 1.d0 + ( mi/mz  )
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   sub_block_four = sqrt_mi_on_ti * inv_ti_joule * inv_ti_joule

      friction_coef_one = (coulomb_logarithm / (4.d0*cpi))
     >     * sub_block_two * sub_block_three
     >     * ni / (2.d0 * cev * ti)

      if (vz_tilde_squared.eq.0.d0) then
         friction_coef_two = 0.d0
      else
         vz_tilde_abs = sqrt(vz_tilde_squared)
         phi = derf(vz_tilde_abs)
         phi_prime = 2.d0*dexp(-vz_tilde_squared)/sqrt(cpi)
         friction_coef_two = (phi - vz_tilde_abs*phi_prime)/
     >        (vz_tilde_abs**3.d0)
      endif

      friction_para = - friction_coef_one * friction_coef_two
     >     * vz_tilde_para      ! FRICTION FORCE F0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
      end subroutine
