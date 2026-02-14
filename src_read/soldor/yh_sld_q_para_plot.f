!***********************************************************************
      subroutine yh_sld_q_para_plot
!***********************************************************************
      use cgdcom, only : grdx, grdy
      use cphcns, only : cev, cmp, cpi
      use cplcom, only : aion, aza, heat_flux_para_by_kappa_para
     >    , flmxi, vna, vne, vte, vti, vzf
      use cplmet, only : hdsp, hdxm, hdxp, hpit, hvsb, icel, itmax
     >    , itmpe, itmps, jcel, jtmax, jtmin, kce, kcw, kgdx, kgdy
      use csize,  only : ndx
      implicit none
!
!::local variables
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jt, jts, jte, j, jp, i, irg, ia, iy
!ik   integer  jchsl, icosl, icisl, i6
!ik   real*8   fcenh, fcslw
      integer  it, jt, jts, jte, j, jp, i, ia
!
!  connecting length
      real*8 l_connect(ndx)

      real*8 permittivity_eps_zero
      real*8 m_bulk
      real*8 q_para_helander, q_para_helander_limited
      real*8 q_para_helander_limited_sld_apprx
      real*8 tau_ii_soldor_file, tau_ii_helander
      real*8 tau_ii_soldor_file_num, tau_ii_soldor_file_denom
      real*8 tau_ii_helander_num, tau_ii_helander_denom
      real*8 kappa_ii_para_helander
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   real*8 kappa
!      real*8 kappa_ii_helander, kappa_ii_soldor
      real*8 temp_zf, temp_ne, temp_te, temp_ti, temp_ni
      real*8 temp_mass_number, temp_abs_charge_number
      real*8 coulog, coulomb_logarithm
      real*8 clambda_wesson
      real*8 alpha_helander, flux_limiter_helander
      real*8 q_para_upper_limit
      real*8 dl_connect, grad_ti, grad_ti_sld_apprx
      real*8 q_para_helander_w_gt_sld_apprx
      real*8 flux_limiter_helander_w_gt_sld_apprx


!::clear
!
!     Constants
      permittivity_eps_zero = 8.854187817d-12
!       Bulk plamsa ion species is D+
!       cmp = 1.67252e-27[kg] @ phcnst.f
      ia = 1
      m_bulk = aion(ia) * cmp

      alpha_helander = 0.5d0

!:: pm3d plot q_para
      !file open
      open(3020,file='sld_heat_flux_para_pm3d.dat')

      do it = 1, itmax
         jts = jtmin(it)
         jte = jtmax(it)
         if( it.ge.itmps .and. it.le.itmpe ) then
            jts = jts + 1
            jte = jte - 1
         endif
!
         do jt = jts, jte
            j  = jcel(jt,it)
            i  = icel(jt,it)
            jp  = jcel(jt+1,it)
!            irg = kreg(j,i)
!
!       q// helander evaluate (at cell center) for check
!       ***Attention: just rough check,
!           kappa --> cell center j, grad Ti --> cell boundary j+1/2
!           SONIC poloidal direction (outer to inner div) is defined to be POSITIVE.
!                                             2018/05/22 Y.Homma
            temp_zf = vzf(j,i) ; temp_ne = vne(j,i)
            temp_ni = vna(j,i,ia)
            temp_te = vte(j,i) ; temp_ti = vti(j,i)
!                       temp_te [eV], temp_ti[eV]
            temp_mass_number = aion(ia)
            temp_abs_charge_number = abs( aza(ia) )
            coulomb_logarithm = coulog( temp_zf, temp_ne,
     >           temp_te, temp_ti,
     >           temp_mass_number, temp_abs_charge_number )

!     Wesson Sec. 14.5 clambda for ion-ion collisions T<10keV
            clambda_wesson = 17.3d0 -0.5d0*log(temp_ne/1.0d20)
     >           + 1.5d0*log(temp_ti/1000.d0)

!     Tau_ii
            tau_ii_soldor_file_num = 2.0853d13 * (temp_ti**1.5d0)
     >           * sqrt(m_bulk / cmp)
            tau_ii_soldor_file_denom = coulomb_logarithm * temp_ni
     >           * (aza(ia)**4.d0)
            tau_ii_soldor_file
     >           = tau_ii_soldor_file_num / tau_ii_soldor_file_denom

!            tau_ii_helander_num = (12.d0/sqrt(2.d0)) * (cpi**1.5d0)
            tau_ii_helander_num = (12.d0) * (cpi**1.5d0)
     >           * (permittivity_eps_zero**2.d0)
     >           * sqrt(m_bulk) * ( (cev*temp_ti)**1.5d0 )
            tau_ii_helander_denom = coulomb_logarithm * temp_ni
     >           * ( (aza(ia)*cev )**4.d0)

            tau_ii_helander
     >           = tau_ii_helander_num / tau_ii_helander_denom

!     Kappa // classical definition Helander @ cell center (j,i)
!                           SOLDOR file definition is identical
            kappa_ii_para_helander = 3.9d0 * temp_ni * (cev*temp_ti)
     >           * tau_ii_helander / m_bulk

!     Grad Ti evaluation @ cell boundary (j+1/2,i)
!     SONIC poloidal direction (outer to inner div) is defined to be POSITIVE.
            dl_connect = ( hdxm(j,i) + hdxp(j,i) ) / hpit(j,i)
            grad_ti = ( vti(jp,i) - vti(j,i) )*cev / dl_connect
!                                               [joule/m]
            q_para_helander = - kappa_ii_para_helander * grad_ti
!                                               [joule/(m2.s)]

            q_para_upper_limit = alpha_helander*temp_ni*
     >           (temp_ti*cev)*sqrt(temp_ti*cev/m_bulk)

            flux_limiter_helander = 1.0d0
     >           /(1.0d0+dabs(q_para_helander/q_para_upper_limit))

            q_para_helander_limited =
     >           flux_limiter_helander*q_para_helander

!
            write(3020, '(14(E16.4e2))')
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)),
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           heat_flux_para_by_kappa_para(j,i),
     >           flmxi(j,i), coulomb_logarithm, clambda_wesson,
     >           tau_ii_soldor_file, tau_ii_helander,
     >           kappa_ii_para_helander, q_para_helander,
     >           flux_limiter_helander, q_para_helander_limited
            write(3020, '(14(E16.4e2))')
     >           grdx(kgdx(j,i,2), kgdy(j,i,2)),
     >           grdy(kgdx(j,i,2), kgdy(j,i,2)),
     >           heat_flux_para_by_kappa_para(j,i),
     >           flmxi(j,i), coulomb_logarithm, clambda_wesson,
     >           tau_ii_soldor_file, tau_ii_helander,
     >           kappa_ii_para_helander, q_para_helander,
     >           flux_limiter_helander, q_para_helander_limited
            write(3020, '(14(E16.4e2))')
     >           grdx(kgdx(j,i,3), kgdy(j,i,3)),
     >           grdy(kgdx(j,i,3), kgdy(j,i,3)),
     >           heat_flux_para_by_kappa_para(j,i),
     >           flmxi(j,i), coulomb_logarithm, clambda_wesson,
     >           tau_ii_soldor_file, tau_ii_helander,
     >           kappa_ii_para_helander, q_para_helander,
     >           flux_limiter_helander, q_para_helander_limited
            write(3020, '(14(E16.4e2))')
     >           grdx(kgdx(j,i,4), kgdy(j,i,4)),
     >           grdy(kgdx(j,i,4), kgdy(j,i,4)),
     >           heat_flux_para_by_kappa_para(j,i),
     >           flmxi(j,i), coulomb_logarithm, clambda_wesson,
     >           tau_ii_soldor_file, tau_ii_helander,
     >           kappa_ii_para_helander, q_para_helander,
     >           flux_limiter_helander, q_para_helander_limited
            write(3020, '(14(E16.4e2))')
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)),
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           heat_flux_para_by_kappa_para(j,i),
     >           flmxi(j,i), coulomb_logarithm, clambda_wesson,
     >           tau_ii_soldor_file, tau_ii_helander,
     >           kappa_ii_para_helander, q_para_helander,
     >           flux_limiter_helander, q_para_helander_limited

            write(3020, '(14(E16.4e2))')
         enddo

         write(3020, '(14(E16.4e2))')
      enddo

      !file close
      close(3020)
!
!     1D tube plot q_para
      l_connect = 0.d0

      open(3110,file='sld_heat_flux_para_pol_edge_42.dat')
      open(3120,file='sld_heat_flux_para_pol_sepx_24.dat')
      open(3130,file='sld_heat_flux_para_pol_sol_15.dat')
      open(3140,file='sld_heat_flux_para_pol_sol_05.dat')
      open(3150,file='sld_heat_flux_para_pol_outermost_02.dat')

      do it = 1, itmax
         jts = jtmin(it)
         jte = jtmax(it)
         if( it.ge.itmps .and. it.le.itmpe ) then
            jts = jts + 1
            jte = jte - 1
         endif

         call slenb(it,l_connect)
!
         do jt = jts, jte
            j  = jcel(jt,it)
            i  = icel(jt,it)
            jp  = jcel(jt+1,it)
!            irg = kreg(j,i)
!
!       q// helander evaluate (at cell center) for check
!       ***Attention: just rough check,
!           kappa --> cell center j, grad Ti --> cell boundary j+1/2
!           SONIC poloidal direction (outer to inner div) is defined to be POSITIVE.
!                                             2018/05/22 Y.Homma
            temp_zf = vzf(j,i) ; temp_ne = vne(j,i)
            temp_ni = vna(j,i,ia)
            temp_te = vte(j,i) ; temp_ti = vti(j,i)
!                       temp_te [eV], temp_ti[eV]
            temp_mass_number = aion(ia)
            temp_abs_charge_number = abs( aza(ia) )
            coulomb_logarithm = coulog( temp_zf, temp_ne,
     >           temp_te, temp_ti,
     >           temp_mass_number, temp_abs_charge_number )

!     Wesson Sec. 14.5 clambda for ion-ion collisions T<10keV
            clambda_wesson = 17.3d0 -0.5d0*log(temp_ne/1.0d20)
     >           + 1.5d0*log(temp_ti/1000.d0)

!     Tau_ii
            tau_ii_soldor_file_num = 2.0853d13 * (temp_ti**1.5d0)
     >           * sqrt(m_bulk / cmp)
            tau_ii_soldor_file_denom = coulomb_logarithm * temp_ni
     >           * (aza(ia)**4.d0)
            tau_ii_soldor_file
     >           = tau_ii_soldor_file_num / tau_ii_soldor_file_denom

!            tau_ii_helander_num = (12.d0/sqrt(2.d0)) * (cpi**1.5d0)
            tau_ii_helander_num = (12.d0) * (cpi**1.5d0)
     >           * (permittivity_eps_zero**2.d0)
     >           * sqrt(m_bulk) * ( (cev*temp_ti)**1.5d0 )
            tau_ii_helander_denom = coulomb_logarithm * temp_ni
     >           * ( (aza(ia)*cev )**4.d0)

            tau_ii_helander
     >           = tau_ii_helander_num / tau_ii_helander_denom

!     Kappa // classical definition Helander @ cell center (j,i)
!                           SOLDOR file definition is identical
            kappa_ii_para_helander = 3.9d0 * temp_ni * (cev*temp_ti)
     >           * tau_ii_helander / m_bulk

!     Grad Ti evaluation @ cell boundary (j+1/2,i)
!     SONIC poloidal direction (outer to inner div) is defined to be POSITIVE.
            dl_connect = ( hdxm(j,i) + hdxp(j,i) ) / hpit(j,i)
            grad_ti = ( vti(jp,i) - vti(j,i) )*cev / dl_connect
!                                               [joule/m]
            q_para_helander = - kappa_ii_para_helander * grad_ti
!                                               [joule/(m2.s)]
            q_para_upper_limit = alpha_helander*temp_ni*
     >           (temp_ti*cev)*sqrt(temp_ti*cev/m_bulk)

            flux_limiter_helander = 1.0d0
     >           /(1.0d0+dabs(q_para_helander/q_para_upper_limit))

            q_para_helander_limited =
     >           flux_limiter_helander*q_para_helander

!     q_para helander with grad Ti sld approximation
            grad_ti_sld_apprx =
     >         ( hdsp(j,i,kce)*vti(jp,i) - hdsp(j,i,kcw)*vti(j,i) )*cev
     >           / hvsb(j,i)
!                                               [joule/m]
            q_para_helander_w_gt_sld_apprx =
     >           -kappa_ii_para_helander * grad_ti_sld_apprx
            flux_limiter_helander_w_gt_sld_apprx = 1.d0
     >         /(1.0d0+
     >         dabs(q_para_helander_w_gt_sld_apprx/q_para_upper_limit))

            q_para_helander_limited_sld_apprx   =
     >           flux_limiter_helander_w_gt_sld_apprx *
     >           q_para_helander_w_gt_sld_apprx


            if(it.eq.42) then
               write(3110, '(2i6,15(E16.4e2))')
     >              j, i, l_connect(jt),
     >              heat_flux_para_by_kappa_para(j,i),
     >              flmxi(j,i), coulomb_logarithm, clambda_wesson,
     >              tau_ii_soldor_file, tau_ii_helander,
     >              kappa_ii_para_helander, q_para_helander,
     >              flux_limiter_helander, q_para_helander_limited,
     >              flux_limiter_helander_w_gt_sld_apprx,
     >              q_para_helander_limited_sld_apprx,
     >           heat_flux_para_by_kappa_para(j,i) / q_para_helander
            endif

            if(it.eq.24) then
               write(3120, '(2i6,15(E16.4e2))')
     >              j, i, l_connect(jt),
     >              heat_flux_para_by_kappa_para(j,i),
     >              flmxi(j,i), coulomb_logarithm, clambda_wesson,
     >              tau_ii_soldor_file, tau_ii_helander,
     >              kappa_ii_para_helander, q_para_helander,
     >              flux_limiter_helander, q_para_helander_limited,
     >              flux_limiter_helander_w_gt_sld_apprx,
     >              q_para_helander_limited_sld_apprx,
     >           heat_flux_para_by_kappa_para(j,i) / q_para_helander
            endif

            if(it.eq.15) then
               write(3130, '(2i6,15(E16.4e2))')
     >              j, i, l_connect(jt),
     >              heat_flux_para_by_kappa_para(j,i),
     >              flmxi(j,i), coulomb_logarithm, clambda_wesson,
     >              tau_ii_soldor_file, tau_ii_helander,
     >              kappa_ii_para_helander, q_para_helander,
     >              flux_limiter_helander, q_para_helander_limited,
     >              flux_limiter_helander_w_gt_sld_apprx,
     >              q_para_helander_limited_sld_apprx,
     >           heat_flux_para_by_kappa_para(j,i) / q_para_helander
            endif

            if(it.eq.5) then
               write(3140, '(2i6,15(E16.4e2))')
     >              j, i, l_connect(jt),
     >              heat_flux_para_by_kappa_para(j,i),
     >              flmxi(j,i), coulomb_logarithm, clambda_wesson,
     >              tau_ii_soldor_file, tau_ii_helander,
     >              kappa_ii_para_helander, q_para_helander,
     >              flux_limiter_helander, q_para_helander_limited,
     >              flux_limiter_helander_w_gt_sld_apprx,
     >              q_para_helander_limited_sld_apprx,
     >           heat_flux_para_by_kappa_para(j,i) / q_para_helander
            endif

            if(it.eq.2) then
               write(3150, '(2i6,15(E16.4e2))')
     >              j, i, l_connect(jt),
     >              heat_flux_para_by_kappa_para(j,i),
     >              flmxi(j,i), coulomb_logarithm, clambda_wesson,
     >              tau_ii_soldor_file, tau_ii_helander,
     >              kappa_ii_para_helander, q_para_helander,
     >              flux_limiter_helander, q_para_helander_limited,
     >              flux_limiter_helander_w_gt_sld_apprx,
     >              q_para_helander_limited_sld_apprx,
     >           heat_flux_para_by_kappa_para(j,i) / q_para_helander
            endif

         enddo
      enddo

      close(3110)
      close(3120)
      close(3130)
      close(3140)
      close(3150)

!
!====================

!
!::debug write

      write(3000,*) "it, jtmin(it), jtmax(it)"
      do it = 1, itmax
         write(3000,*) it, jtmin(it), jtmax(it)
      enddo

      return
      end subroutine
!
!============================
!============================
!============================
!***********************************************************************
      subroutine yh_sld_plasma_parameters_1_plot
!***********************************************************************
      use cgdcom, only : grdx, grdy, hbr, hbt, hbz
      use cimcom, only : ami, amz
      use cmeffz, only : vdnz
      use cntcom, only : mcel
      use cphcns, only : cev, cme, cmp, cpi
      use cplcom, only : aion, aza, flimi, flmxe, flmxi
     >    , heat_flux_para_by_kappa_para, heat_flux_para_elec, vna, vne
     >    , vnezef, vni, vte, vti, vva, vzf
      use cplimp, only : ismaxl
      use cplmet, only : hdxm, hdxp, hpit, icel, itmax, jcel, jtmax
     >    , jtmin, kgdx, kgdy
      use cplwrd, only : wime
      use csize,  only : ndx, ndy
      implicit none

!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jt, jts, jte, j, i, irg, ia, iy
      integer  it, jt, jts, jte, j, i, ia
      integer  jp
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  i_imp_species, iz, ic
!ik   integer  jchsl, icosl, icisl, i6
!ik   real*8   fcenh, fcslw
      integer  i_imp_species, ic
      integer mx, my
      real*8 br, bz, bt
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 b_pol(ndx,ndy), b_abs(ndx,ndy)
!ik   real*8 coulog, clambda, eps_zero
      real*8 b_abs(ndx,ndy)
      real*8 clambda, eps_zero
      real*8 mi, ni, ni_p, ti_ev, ti_p_ev, vth_i
      real*8 te_ev, te_p_ev, charge_imp
      real*8 ui, ui_p, dl_connect, sound_speed_ion
      real*8 field_elec_para(ndx,ndy)
      real*8 tau_ii_helander_num, tau_ii_helander_denom
      real*8 tau_ii_helander(ndx,ndy)
      real*8 mfp_ii(ndx,ndy) , mach_number(ndx,ndy)
      real*8 grad_ti_ev(ndx,ndy), grad_te_ev(ndx,ndy)
      real*8 grad_pe(ndx,ndy)
      real*8 grad_ni(ndx,ndy)
      real*8 grad_ui(ndx,ndy), grad_b_abs(ndx,ndy)
      real*8 ch_length_ti_ev(ndx,ndy), ch_length_ni(ndx,ndy)
      real*8 ch_length_ui(ndx,ndy), ch_length_b_abs(ndx,ndy)
      real*8 knudsen_ti_ev(ndx,ndy), knudsen_ni(ndx,ndy)
      real*8 knudsen_ti_ev_effective(ndx,ndy)
      real*8 knudsen_ui(ndx,ndy), knudsen_b_abs(ndx,ndy)
      real*8 knudsen_connect(ndx,ndy), abs_knudsen_max(ndx,ndy)
      real*8 pressure_e(ndx,ndy), pressure_i(ndx,ndy)
      real*8 nz_tot_after_wfac(ndx,ndy), nz_on_ne(ndx,ndy)
      real*8 tot_rad_rate(ndx,ndy)
      real*8 f_thi(ndx,ndy), f_the(ndx,ndy),  f_0(ndx,ndy)
      real*8 f_elec(ndx,ndy), f_tot_simple(ndx,ndy), f_tot_full(ndx,ndy)
      real*8 energy_dens_total_para(ndx,ndy)
      real*8 energy_dens_convective_para(ndx,ndy)

      real*8 mfp_ion_effective(ndx,ndy)
      real*8 omega_gm(ndx,ndy), flmxi_qgm(ndx,ndy)
      real*8 heat_flux_para_gm(ndx,ndy)

!
!  connecting length
      real*8 l_connect(ndx)
! added 2 lines organize local variables and include files by kamata 2021/05/31
! function
      real(8)    coulog

!::clear
      ch_length_ti_ev=0.d0; ch_length_ni=0.d0; ch_length_ui=0.d0
      ch_length_b_abs =0.d0
      knudsen_ti_ev=0.d0; knudsen_ni=0.d0; knudsen_ui=0.d0
      knudsen_b_abs=0.d0; knudsen_connect=0.d0
      knudsen_ti_ev_effective=0.d0

      nz_tot_after_wfac=0.d0; nz_on_ne=0.d0
      tot_rad_rate=0.d0
      field_elec_para=0.d0

      f_tot_simple=0.d0; f_tot_full=0.d0
      f_0=0.d0; f_thi=0.d0; f_the=0.d0; f_elec=0.d0

      mfp_ion_effective=0.d0; omega_gm=0.d0; flmxi_qgm=0.d0
      heat_flux_para_gm=0.d0
!
!:: pm3d
      !species specify
      ia = 1
      eps_zero = 8.854187817d-12
      mi = aion(ia) * cmp

      i_imp_species = 1
!     index to specify impurity species(NOT Zimp!)
      charge_imp = 4.d0 * cev
!     Impurity charge q_z[C]  !!! Ar 4+ supposed. !!!)
      write(2000,*) "test: ami/cmp, amz/cmp"
      write(2000,*)  ami/cmp, amz/cmp


      !Physics quantities evaluation
      do it = 1, itmax
         jts = jtmin(it)
         jte = jtmax(it)
!         if( it.ge.itmps .and. it.le.itmpe ) then
!          jts = jts + 1
!          jte = jte - 1
!         endif
!
         do jt = jts, jte
            j  = jcel(jt,it)
            i  = icel(jt,it)
            ic = mcel(j,i)  ! MC mesh number
            if(jt .eq. jte) then
               jp = j
            else
               jp  = jcel(jt+1,it)
            endif
!     Coulomb logarithm value
            clambda = coulog(
     >           vzf(j,i), vne(j,i), vte(j,i), vti(j,i),
     >           aion(ia), abs( aza(ia) )  )
!     Magnetic field  (*** at the corner 2, boundary j+1/2 ***)
            mx = kgdx(j,i,2); my = kgdy(j,i,2)
            br = hbr(mx,my); bz = hbz(mx,my)
            bt = hbt(mx,my)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik         b_pol(j,i) = sqrt(br*br + bz*bz)
            b_abs(j,i) = sqrt(br*br + bz*bz + bt*bt)
!     tau_ii
            ti_ev = vti(j,i);  ni = vna(j,i,ia)
!            tau_ii_helander_num = (12.d0/sqrt(2.d0)) * (cpi**1.5d0)
            tau_ii_helander_num = (12.d0) * (cpi**1.5d0)
     >           * (eps_zero**2.d0)
     >           * sqrt(mi) * ( (cev*ti_ev)**1.5d0 )
            tau_ii_helander_denom = clambda * ni
     >           * ( (aza(ia)*cev )**4.d0)
            tau_ii_helander(j,i)
     >           = tau_ii_helander_num / tau_ii_helander_denom
!     MFP i-i collision
            vth_i = sqrt(cev*ti_ev/mi)
            mfp_ii(j,i) = vth_i * tau_ii_helander(j,i)
!     pressures
            pressure_e(j,i) = vnezef(j,i) * vte(j,i) * cev
            pressure_i(j,i) = vna(j,i,ia) * vti(j,i) *cev
!     Impurity density and radiation
            nz_tot_after_wfac(j,i)
     >       = sum( vdnz(j,i, 0:ismaxL(i_imp_species), i_imp_species) )
            nz_on_ne(j,i) = nz_tot_after_wfac(j,i) / vnezef(j,i)
!                           ne includes elec. from ionized impurities.
            tot_rad_rate(j,i) = -wime(ic)
!     Mach number, same sign as ui//=vva is retained.
            sound_speed_ion = sqrt(cev*( vte(j,i) + vti(j,i)  ) / mi)
            mach_number(j,i) = vva(j,i,ia) / sound_speed_ion
!     Convective & Total energy flux density
            energy_dens_convective_para(j,i) =
     >           2.5d0*vna(j,i,ia)*(cev*vti(j,i))*vva(j,i,ia) +
     >           0.5d0*mi*vna(j,i,ia)*(vva(j,i,ia)**3.d0)
            energy_dens_total_para(j,i) =
     >           energy_dens_convective_para(j,i) +
     >           heat_flux_para_by_kappa_para(j,i)
         enddo                  !            do jt = jts, jte
      enddo     !      do it = 1, itmax


!     EVALUATE Grad Ti-ni-ui-B
!              Collisionalities
!              SONIC's static-elec. field
!              Forces
      do it = 1, itmax
         l_connect = 0.d0
         call slenb(it,l_connect)

         jts = jtmin(it); jte = jtmax(it)
!         if( it.ge.itmps .and. it.le.itmpe ) then
!            jts = jts + 1
!            jte = jte - 1
!         endif
!
         do jt = jts, jte
            j  = jcel(jt,it)
            i  = icel(jt,it)
            if(jt .eq. jte) then
               jp = j
            else
               jp  = jcel(jt+1,it)
            endif
!     Coulomb logarithm value
            clambda = coulog(
     >           vzf(j,i), vne(j,i), vte(j,i), vti(j,i),
     >           aion(ia), abs( aza(ia) )  )

!           Grad Ti evaluation @ cell boundary (j+1/2,i), same as q// helander.
!           SONIC poloidal direction (outer to inner div) is defined to be POSITIVE.
!            grad_ti = ( vti(jp,i) - vti(j,i) )*cev / dl_connect
!                                               [joule/m]
            ti_ev = vti(j,i);  ti_p_ev = vti(jp,i)
            te_ev = vte(j,i);  te_p_ev = vte(jp,i)
            ni = vna(j,i,ia);  ni_p = vna(jp,i,ia)
            ui = vva(j,i,ia);  ui_p = vva(jp,i,ia)
            dl_connect = ( hdxm(j,i) + hdxp(j,i) ) / hpit(j,i)
            if(dl_connect.le.0.d0) cycle

            grad_ti_ev(j,i) =  ( ti_p_ev - ti_ev ) / dl_connect
            grad_te_ev(j,i) =  ( te_p_ev - te_ev ) / dl_connect
!                                         [eV/m]
            grad_ni(j,i) =  ( ni_p - ni ) / dl_connect
            grad_ui(j,i) =  ( ui_p - ui ) / dl_connect
            grad_pe(j,i) =  ( pressure_e(jp,i) - pressure_e(j,i) ) /
     >           dl_connect
            grad_b_abs(j,i) = ( b_abs(jp,i) - b_abs(j,i) ) / dl_connect

!     Electric field (static) (cf. Note p.imp6 OR soniv doc. Eq.(2.2-11))
            field_elec_para(j,i) = -grad_pe(j,i)/(vnezef(j,i)*cev) -
     >           0.71d0*( grad_te_ev(j,i)*cev )/cev

!     Characteristic length evaluation
            if(abs(grad_ti_ev(j,i) ).gt.0.d0)  ch_length_ti_ev(j,i)
     >           = ti_ev / grad_ti_ev(j,i)
            if(abs(grad_ni(j,i) ).gt.0.d0)  ch_length_ni(j,i)
     >           = ni / grad_ni(j,i)
            if(abs(grad_ui(j,i) ).gt.0.d0)  ch_length_ui(j,i)
     >           = ui / grad_ui(j,i)
            if(abs(grad_b_abs(j,i) ).gt.0.d0)  ch_length_b_abs(j,i)
     >           = b_abs(j,i) / grad_b_abs(j,i)

!     Knudsen numbers
            if(abs(ch_length_ti_ev(j,i)).gt.0.d0) knudsen_ti_ev(j,i)
     >           = mfp_ii(j,i) / ch_length_ti_ev(j,i)
            if(abs(ch_length_ni(j,i)).gt.0.d0) knudsen_ni(j,i)
     >           = mfp_ii(j,i) / ch_length_ni(j,i)
            if(abs(ch_length_ui(j,i)).gt.0.d0) knudsen_ui(j,i)
     >           = mfp_ii(j,i) / ch_length_ui(j,i)
            if(abs(ch_length_b_abs(j,i)).gt.0.d0) knudsen_b_abs(j,i)
     >           = mfp_ii(j,i) / ch_length_b_abs(j,i)
            knudsen_connect(j,i) = mfp_ii(j,i) / l_connect(jte-1)
            abs_knudsen_max(j,i) = max( abs(knudsen_ti_ev(j,i) ),
     >           abs(knudsen_ni(j,i) ), abs(knudsen_ui(j,i) ),
     >        abs(knudsen_b_abs(j,i) ), abs(knudsen_connect(j,i) ) )

!           Effective Kn_Ti
            if(abs(flmxi(j,i) ).gt.0.d0) knudsen_ti_ev_effective(j,i)
     >           = ( (1.d0/flmxi(j,i)) -1.d0) * flimi / 3.9d0

!     qGM rough estimation
!         ion_MFP_effective(i-i & i-Z colls included, evaluated from SOLDOR info.)
            mfp_ion_effective(j,i) =
     >       knudsen_ti_ev_effective(j,i) * abs(ch_length_ti_ev(j,i) )
!         Omega_GM
            if(abs(ch_length_ui(j,i)).gt.0.d0) then
               if(abs(ch_length_b_abs(j,i)).gt.0.d0) then
                  omega_gm(j,i) =
     >              5.88d0 * mach_number(j,i) * mfp_ion_effective(j,i) *
     >     (1.d0 / ch_length_ui(j,i)  - 0.25d0 / ch_length_b_abs(j,i) )
               endif
            endif
!         flmxi_qGM
            if(abs(1.d0+omega_gm(j,i) ).gt.0.d0) then
               flmxi_qgm(j,i) =
     >           (1.d0 - 1.2d0*mach_number(j,i)*mach_number(j,i) ) /
     >           (1.d0 + omega_gm(j,i) )
            endif
!         qGM (heat_flux_para_GM)
            if(abs(flmxi(j,i)).gt.0.d0) then
               heat_flux_para_gm(j,i) =
     >           (heat_flux_para_by_kappa_para(j,i)/flmxi(j,i) ) *
     >           flmxi_qgm(j,i)
            endif

!     Force balance
            call thermal_force_para_by_q_para(clambda,
     >           ami, amz, 1.d0*cev, charge_imp, ti_ev,
     >           0.d0, sqrt(2.d0 * cev * ti_ev /amz), ui,
     >           heat_flux_para_by_kappa_para(j,i), f_thi(j,i) )
!!!
!!!     JUST TEST 2018/11/29  analysis for Fthi W/O collisionality effect
!            if(abs(flmxi(j,i) ).gt.0.d0) then
!               f_thi(j,i) = f_thi(j,i) / flmxi(j,i)
!            else
!               f_thi(j,i) = 0.d0
!            endif
!!!
            call thermal_force_para_by_q_para(clambda,
     >           cme, amz, (-cev), charge_imp, te_ev,
     >           0.d0, sqrt(2.d0 * cev * ti_ev /amz), ui,
     >           heat_flux_para_elec(j,i), f_the(j,i) )

!      subroutine thermal_force_para_by_q_para(coulomb_logarithm,
!     >     mi, mz, ei[C], ez[C], ti[eV], vz_para, vz_perp, ui, q_para, fth_para)
            call friction_force_theoretical_value(clambda,
     >           ami, amz, 1.d0*cev, charge_imp, ti_ev,
     >           ni,
     >           0.d0, sqrt(2.d0 * cev * ti_ev /amz), ui, f_0(j,i) )
!      subroutine friction_force_theoretical_value(coulomb_logarithm,
!     >     mi, mz, ei[C], ez[C], ti[eV], ni[m^(-3)], vz_para, vz_perp, ui, friction_para)

!     Electric field (static) force
            f_elec(j,i) = charge_imp  * field_elec_para(j,i)

!     Resultant force
            f_tot_simple(j,i) = f_0(j,i) + f_thi(j,i)
            f_tot_full(j,i) = f_0(j,i) + f_thi(j,i) + f_the(j,i) +
     >           f_elec(j,i)

!HERE ONCE REVIEW
         enddo                  !            do jt = jts, jte
      enddo                     !      do it = 1, itmax


!     FILE WRITING
      !file open 1
      open(3000,file='sld_plasma_parameters_1_pm3d.dat')
      open(3100,file='sld_plasma_parameters_2_pm3d.dat')
      open(3200,file='sld_plasma_parameters_3_pm3d.dat')
      open(3300,file='sld_plasma_parameters_4_pm3d.dat')
      do it = 1, itmax
         jts = jtmin(it)
         jte = jtmax(it)
!         if( it.ge.itmps .and. it.le.itmpe ) then
!            jts = jts + 1
!            jte = jte - 1
!         endif
!
!         do jt = jts, jte
         do jt = jts+1, jte-1    ! to remove extreme values at the ends
            j  = jcel(jt,it)
            i  = icel(jt,it)
!            irg = kreg(j,i)
!
!     PARAMETERS 1
            write(3000, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)),
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           vte(j,i), vti(j,i), vva(j,i,ia),
     >           b_abs(j,i), vne(j,i), vni(j,i), vna(j,i,ia),
     >           pressure_e(j,i), pressure_i(j,i),
     >           heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >           grad_ti_ev(j,i)
            write(3000, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,2), kgdy(j,i,2)),
     >           grdy(kgdx(j,i,2), kgdy(j,i,2)),
     >           vte(j,i), vti(j,i), vva(j,i,ia),
     >           b_abs(j,i), vne(j,i), vni(j,i), vna(j,i,ia),
     >           pressure_e(j,i), pressure_i(j,i),
     >           heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >           grad_ti_ev(j,i)
            write(3000, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,3), kgdy(j,i,3)),
     >           grdy(kgdx(j,i,3), kgdy(j,i,3)),
     >           vte(j,i), vti(j,i), vva(j,i,ia),
     >           b_abs(j,i), vne(j,i), vni(j,i), vna(j,i,ia),
     >           pressure_e(j,i), pressure_i(j,i),
     >           heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >           grad_ti_ev(j,i)
            write(3000, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,4), kgdy(j,i,4)),
     >           grdy(kgdx(j,i,4), kgdy(j,i,4)),
     >           vte(j,i), vti(j,i), vva(j,i,ia),
     >           b_abs(j,i), vne(j,i), vni(j,i), vna(j,i,ia),
     >           pressure_e(j,i), pressure_i(j,i),
     >           heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >           grad_ti_ev(j,i)
            write(3000, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)),
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           vte(j,i), vti(j,i), vva(j,i,ia),
     >           b_abs(j,i), vne(j,i), vni(j,i), vna(j,i,ia),
     >           pressure_e(j,i), pressure_i(j,i),
     >           heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >           grad_ti_ev(j,i)
            write(3000, '(15(E16.4e2))')
!
!     PARAMETERS 2
!
            write(3100, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)),
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           tau_ii_helander(j,i), mfp_ii(j,i),
     >           grad_ti_ev(j,i), grad_ni(j,i), grad_ui(j,i),
     >           grad_b_abs(j,i), ch_length_ti_ev(j,i),
     >           ch_length_ni(j,i), ch_length_ui(j,i),
     >           ch_length_b_abs(j,i)
            write(3100, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,2), kgdy(j,i,2)),
     >           grdy(kgdx(j,i,2), kgdy(j,i,2)),
     >           tau_ii_helander(j,i), mfp_ii(j,i),
     >           grad_ti_ev(j,i), grad_ni(j,i), grad_ui(j,i),
     >           grad_b_abs(j,i), ch_length_ti_ev(j,i),
     >           ch_length_ni(j,i), ch_length_ui(j,i),
     >           ch_length_b_abs(j,i)
            write(3100, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,3), kgdy(j,i,3)),
     >           grdy(kgdx(j,i,3), kgdy(j,i,3)),
     >           tau_ii_helander(j,i), mfp_ii(j,i),
     >           grad_ti_ev(j,i), grad_ni(j,i), grad_ui(j,i),
     >           grad_b_abs(j,i), ch_length_ti_ev(j,i),
     >           ch_length_ni(j,i), ch_length_ui(j,i),
     >           ch_length_b_abs(j,i)
            write(3100, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,4), kgdy(j,i,4)),
     >           grdy(kgdx(j,i,4), kgdy(j,i,4)),
     >           tau_ii_helander(j,i), mfp_ii(j,i),
     >           grad_ti_ev(j,i), grad_ni(j,i), grad_ui(j,i),
     >           grad_b_abs(j,i), ch_length_ti_ev(j,i),
     >           ch_length_ni(j,i), ch_length_ui(j,i),
     >           ch_length_b_abs(j,i)
            write(3100, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)),
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           tau_ii_helander(j,i), mfp_ii(j,i),
     >           grad_ti_ev(j,i), grad_ni(j,i), grad_ui(j,i),
     >           grad_b_abs(j,i), ch_length_ti_ev(j,i),
     >           ch_length_ni(j,i), ch_length_ui(j,i),
     >           ch_length_b_abs(j,i)
            write(3100, '(15(E16.4e2))')
!
!     PARAMETERS 3
!
            write(3200, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)),
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           knudsen_ti_ev(j,i), knudsen_ni(j,i), knudsen_ui(j,i),
     >           knudsen_b_abs(j,i), knudsen_connect(j,i),
     >           abs_knudsen_max(j,i)
            write(3200, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,2), kgdy(j,i,2)),
     >           grdy(kgdx(j,i,2), kgdy(j,i,2)),
     >           knudsen_ti_ev(j,i), knudsen_ni(j,i), knudsen_ui(j,i),
     >           knudsen_b_abs(j,i), knudsen_connect(j,i),
     >           abs_knudsen_max(j,i)
            write(3200, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,3), kgdy(j,i,3)),
     >           grdy(kgdx(j,i,3), kgdy(j,i,3)),
     >           knudsen_ti_ev(j,i), knudsen_ni(j,i), knudsen_ui(j,i),
     >           knudsen_b_abs(j,i), knudsen_connect(j,i),
     >           abs_knudsen_max(j,i)
            write(3200, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,4), kgdy(j,i,4)),
     >           grdy(kgdx(j,i,4), kgdy(j,i,4)),
     >           knudsen_ti_ev(j,i), knudsen_ni(j,i), knudsen_ui(j,i),
     >           knudsen_b_abs(j,i), knudsen_connect(j,i),
     >           abs_knudsen_max(j,i)
            write(3200, '(15(E16.4e2))')
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)),
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           knudsen_ti_ev(j,i), knudsen_ni(j,i), knudsen_ui(j,i),
     >           knudsen_b_abs(j,i), knudsen_connect(j,i),
     >           abs_knudsen_max(j,i)
            write(3200, '(15(E16.4e2))')
!
!     PARAMETERS 4
            write(3300, '(25(E16.4e2))')
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)),
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >           tot_rad_rate(j,i),
!     >           f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >           f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >           vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >           vva(j,i,ia),
     >           heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >           mach_number(j,i),
     >           energy_dens_convective_para(j,i),
     >           energy_dens_total_para(j,i),
     >           knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >           heat_flux_para_gm(j,i)
            write(3300, '(25(E16.4e2))')
     >           grdx(kgdx(j,i,2), kgdy(j,i,2)),
     >           grdy(kgdx(j,i,2), kgdy(j,i,2)),
     >           nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >           tot_rad_rate(j,i),
!     >           f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >           f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >           vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >           vva(j,i,ia),
     >           heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >           mach_number(j,i),
     >           energy_dens_convective_para(j,i),
     >           energy_dens_total_para(j,i),
     >           knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >           heat_flux_para_gm(j,i)
            write(3300, '(25(E16.4e2))')
     >           grdx(kgdx(j,i,3), kgdy(j,i,3)),
     >           grdy(kgdx(j,i,3), kgdy(j,i,3)),
     >           nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >           tot_rad_rate(j,i),
!     >           f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >           f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >           vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >           vva(j,i,ia),
     >           heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >           mach_number(j,i),
     >           energy_dens_convective_para(j,i),
     >           energy_dens_total_para(j,i),
     >           knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >           heat_flux_para_gm(j,i)
            write(3300, '(25(E16.4e2))')
     >           grdx(kgdx(j,i,4), kgdy(j,i,4)),
     >           grdy(kgdx(j,i,4), kgdy(j,i,4)),
     >           nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >           tot_rad_rate(j,i),
!     >           f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >           f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >           vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >           vva(j,i,ia),
     >           heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >           mach_number(j,i),
     >           energy_dens_convective_para(j,i),
     >           energy_dens_total_para(j,i),
     >           knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >           heat_flux_para_gm(j,i)
            write(3300, '(25(E16.4e2))')
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)),
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >           tot_rad_rate(j,i),
!     >           f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >           f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >           vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >           vva(j,i,ia),
     >           heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >           mach_number(j,i),
     >           energy_dens_convective_para(j,i),
     >           energy_dens_total_para(j,i),
     >           knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >           heat_flux_para_gm(j,i)
            write(3300, '(25(E16.4e2))')

         enddo
         write(3000, '(15(E16.4e2))')
         write(3100, '(15(E16.4e2))')
         write(3200, '(15(E16.4e2))')
         write(3300, '(25(E16.4e2))')
      enddo

      !file close
      close(3000);  close(3100); close(3200); close(3300)

!
!     1D tube plot q_para
      open(3010,file='sld_plasma_parameters_1_pol_edge_42.dat')
      open(3020,file='sld_plasma_parameters_1_pol_sepx_24.dat')
      open(3030,file='sld_plasma_parameters_1_pol_sol_15.dat')
      open(3040,file='sld_plasma_parameters_1_pol_sol_05.dat')
      open(3050,file='sld_plasma_parameters_1_pol_outermost_02.dat')

      open(3110,file='sld_plasma_parameters_2_pol_edge_42.dat')
      open(3120,file='sld_plasma_parameters_2_pol_sepx_24.dat')
      open(3130,file='sld_plasma_parameters_2_pol_sol_15.dat')
      open(3140,file='sld_plasma_parameters_2_pol_sol_05.dat')
      open(3150,file='sld_plasma_parameters_2_pol_outermost_02.dat')

      open(3210,file='sld_plasma_parameters_3_pol_edge_42.dat')
      open(3220,file='sld_plasma_parameters_3_pol_sepx_24.dat')
      open(3230,file='sld_plasma_parameters_3_pol_sol_15.dat')
      open(3240,file='sld_plasma_parameters_3_pol_sol_05.dat')
      open(3250,file='sld_plasma_parameters_3_pol_outermost_02.dat')

      open(3310,file='sld_plasma_parameters_4_pol_edge_42.dat')
      open(3320,file='sld_plasma_parameters_4_pol_sepx_24.dat')
      open(3330,file='sld_plasma_parameters_4_pol_sol_15.dat')
      open(3340,file='sld_plasma_parameters_4_pol_sol_05.dat')
      open(3350,file='sld_plasma_parameters_4_pol_outermost_02.dat')

      open(3321,file='sld_plasma_parameters_4_pol_sepx_23.dat')
      open(3322,file='sld_plasma_parameters_4_pol_sepx_22.dat')
      open(3323,file='sld_plasma_parameters_4_pol_sepx_21.dat')
      open(3324,file='sld_plasma_parameters_4_pol_sepx_20.dat')

      open(4000,file='sld_5_forces_pol_sepx_24.dat',status='replace')

      do it = 1, itmax
         l_connect = 0.d0
         jts = jtmin(it)
         jte = jtmax(it)
!         if( it.ge.itmps .and. it.le.itmpe ) then
!            jts = jts + 1
!            jte = jte - 1
!         endif

         call slenb(it,l_connect)
!
!         do jt = jts, jte
         do jt = jts+1, jte-1    ! to remove extreme values at the ends
            j  = jcel(jt,it)
            i  = icel(jt,it)
!            irg = kreg(j,i)
!
!           PARAMETERS 1
            if(it.eq.42) then
               write(3010, '(2i6,15(E16.4e2))')
     >              j, i, l_connect(jt),
     >              vte(j,i), vti(j,i), vva(j,i,ia),
     >              b_abs(j,i), vne(j,i), vni(j,i), vna(j,i,ia),
     >              pressure_e(j,i), pressure_i(j,i),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              grad_ti_ev(j,i)
            elseif(it.eq.24) then
               write(3020, '(2i6,15(E16.4e2))')
     >              j, i, l_connect(jt),
     >              vte(j,i), vti(j,i), vva(j,i,ia),
     >              b_abs(j,i), vne(j,i), vni(j,i), vna(j,i,ia),
     >              pressure_e(j,i), pressure_i(j,i),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              grad_ti_ev(j,i)
            elseif(it.eq.15) then
               write(3030, '(2i6,15(E16.4e2))')
     >              j, i, l_connect(jt),
     >              vte(j,i), vti(j,i), vva(j,i,ia),
     >              b_abs(j,i), vne(j,i), vni(j,i), vna(j,i,ia),
     >              pressure_e(j,i), pressure_i(j,i),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              grad_ti_ev(j,i)
            elseif(it.eq.05) then
               write(3040, '(2i6,15(E16.4e2))')
     >              j, i, l_connect(jt),
     >              vte(j,i), vti(j,i), vva(j,i,ia),
     >              b_abs(j,i), vne(j,i), vni(j,i), vna(j,i,ia),
     >              pressure_e(j,i), pressure_i(j,i),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              grad_ti_ev(j,i)
            elseif(it.eq.02) then
               write(3050, '(2i6,15(E16.4e2))')
     >              j, i, l_connect(jt),
     >              vte(j,i), vti(j,i), vva(j,i,ia),
     >              b_abs(j,i), vne(j,i), vni(j,i), vna(j,i,ia),
     >              pressure_e(j,i), pressure_i(j,i),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              grad_ti_ev(j,i)
            endif
!           PARAMETERS 2
            if(it.eq.42) then
               write(3110, '(2i6,12(E16.4e2))')
     >              j, i, l_connect(jt),
     >              tau_ii_helander(j,i), mfp_ii(j,i), grad_ti_ev(j,i),
     >              grad_ni(j,i), grad_ui(j,i), grad_b_abs(j,i),
     >              ch_length_ti_ev(j,i), ch_length_ni(j,i),
     >              ch_length_ui(j,i), ch_length_b_abs(j,i)
            elseif(it.eq.24) then
               write(3120, '(2i6,12(E16.4e2))')
     >              j, i, l_connect(jt),
     >              tau_ii_helander(j,i), mfp_ii(j,i), grad_ti_ev(j,i),
     >              grad_ni(j,i), grad_ui(j,i), grad_b_abs(j,i),
     >              ch_length_ti_ev(j,i), ch_length_ni(j,i),
     >              ch_length_ui(j,i), ch_length_b_abs(j,i)
            elseif(it.eq.15) then
               write(3130, '(2i6,12(E16.4e2))')
     >              j, i, l_connect(jt),
     >              tau_ii_helander(j,i), mfp_ii(j,i), grad_ti_ev(j,i),
     >              grad_ni(j,i), grad_ui(j,i), grad_b_abs(j,i),
     >              ch_length_ti_ev(j,i), ch_length_ni(j,i),
     >              ch_length_ui(j,i), ch_length_b_abs(j,i)
            elseif(it.eq.05) then
               write(3140, '(2i6,12(E16.4e2))')
     >              j, i, l_connect(jt),
     >              tau_ii_helander(j,i), mfp_ii(j,i), grad_ti_ev(j,i),
     >              grad_ni(j,i), grad_ui(j,i), grad_b_abs(j,i),
     >              ch_length_ti_ev(j,i), ch_length_ni(j,i),
     >              ch_length_ui(j,i), ch_length_b_abs(j,i)
            elseif(it.eq.02) then
               write(3150, '(2i6,12(E16.4e2))')
     >              j, i, l_connect(jt),
     >              tau_ii_helander(j,i), mfp_ii(j,i), grad_ti_ev(j,i),
     >              grad_ni(j,i), grad_ui(j,i), grad_b_abs(j,i),
     >              ch_length_ti_ev(j,i), ch_length_ni(j,i),
     >              ch_length_ui(j,i), ch_length_b_abs(j,i)
            endif
!           PARAMETERS 3
            if(it.eq.42) then
               write(3210, '(2i6,12(E16.4e2))')
     >              j, i, l_connect(jt),
     >           knudsen_ti_ev(j,i), knudsen_ni(j,i), knudsen_ui(j,i),
     >              knudsen_b_abs(j,i), knudsen_connect(j,i),
     >              abs_knudsen_max(j,i)
            elseif(it.eq.24) then
               write(3220, '(2i6,12(E16.4e2))')
     >              j, i, l_connect(jt),
     >           knudsen_ti_ev(j,i), knudsen_ni(j,i), knudsen_ui(j,i),
     >              knudsen_b_abs(j,i), knudsen_connect(j,i),
     >              abs_knudsen_max(j,i)
            elseif(it.eq.15) then
               write(3230, '(2i6,12(E16.4e2))')
     >              j, i, l_connect(jt),
     >           knudsen_ti_ev(j,i), knudsen_ni(j,i), knudsen_ui(j,i),
     >              knudsen_b_abs(j,i), knudsen_connect(j,i),
     >              abs_knudsen_max(j,i)
            elseif(it.eq.05) then
               write(3240, '(2i6,12(E16.4e2))')
     >              j, i, l_connect(jt),
     >           knudsen_ti_ev(j,i), knudsen_ni(j,i), knudsen_ui(j,i),
     >              knudsen_b_abs(j,i), knudsen_connect(j,i),
     >              abs_knudsen_max(j,i)
            elseif(it.eq.02) then
               write(3250, '(2i6,12(E16.4e2))')
     >              j, i, l_connect(jt),
     >           knudsen_ti_ev(j,i), knudsen_ni(j,i), knudsen_ui(j,i),
     >              knudsen_b_abs(j,i), knudsen_connect(j,i),
     >              abs_knudsen_max(j,i)
            endif

!           PARAMETERS 4
            if(it.eq.42) then
               write(3310, '(2i6,23(E16.4e2))')
     >              j, i, l_connect(jt),
     >              nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >              tot_rad_rate(j,i),
!     >              f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >              f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >              vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >              vva(j,i,ia),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              mach_number(j,i),
     >              energy_dens_convective_para(j,i),
     >              energy_dens_total_para(j,i),
     >              knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >              heat_flux_para_gm(j,i)
            elseif(it.eq.24) then
               write(3320, '(2i6,23(E16.4e2))')
     >              j, i, l_connect(jt),
     >              nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >              tot_rad_rate(j,i),
!     >              f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >              f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >              vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >              vva(j,i,ia),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              mach_number(j,i),
     >              energy_dens_convective_para(j,i),
     >              energy_dens_total_para(j,i),
     >              knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >              heat_flux_para_gm(j,i)
            elseif(it.eq.23) then
               write(3321, '(2i6,23(E16.4e2))')
     >              j, i, l_connect(jt),
     >              nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >              tot_rad_rate(j,i),
!     >              f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >              f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >              vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >              vva(j,i,ia),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              mach_number(j,i),
     >              energy_dens_convective_para(j,i),
     >              energy_dens_total_para(j,i),
     >              knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >              heat_flux_para_gm(j,i)
            elseif(it.eq.22) then
               write(3322, '(2i6,23(E16.4e2))')
     >              j, i, l_connect(jt),
     >              nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >              tot_rad_rate(j,i),
!     >              f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >              f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >              vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >              vva(j,i,ia),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              mach_number(j,i),
     >              energy_dens_convective_para(j,i),
     >              energy_dens_total_para(j,i),
     >              knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >              heat_flux_para_gm(j,i)
            elseif(it.eq.21) then
               write(3323, '(2i6,23(E16.4e2))')
     >              j, i, l_connect(jt),
     >              nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >              tot_rad_rate(j,i),
!     >              f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >              f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >              vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >              vva(j,i,ia),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              mach_number(j,i),
     >              energy_dens_convective_para(j,i),
     >              energy_dens_total_para(j,i),
     >              knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >              heat_flux_para_gm(j,i)
            elseif(it.eq.20) then
               write(3324, '(2i6,23(E16.4e2))')
     >              j, i, l_connect(jt),
     >              nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >              tot_rad_rate(j,i),
!     >              f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >              f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >              vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >              vva(j,i,ia),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              mach_number(j,i),
     >              energy_dens_convective_para(j,i),
     >              energy_dens_total_para(j,i),
     >              knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >              heat_flux_para_gm(j,i)

            elseif(it.eq.15) then
               write(3330, '(2i6,23(E16.4e2))')
     >              j, i, l_connect(jt),
     >              nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >              tot_rad_rate(j,i),
!     >              f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >              f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >              vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >              vva(j,i,ia),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              mach_number(j,i),
     >              energy_dens_convective_para(j,i),
     >              energy_dens_total_para(j,i),
     >              knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >              heat_flux_para_gm(j,i)
            elseif(it.eq.05) then
               write(3340, '(2i6,23(E16.4e2))')
     >              j, i, l_connect(jt),
     >              nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >              tot_rad_rate(j,i),
!     >              f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >              f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >              vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >              vva(j,i,ia),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              mach_number(j,i),
     >              energy_dens_convective_para(j,i),
     >              energy_dens_total_para(j,i),
     >              knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >              heat_flux_para_gm(j,i)
            elseif(it.eq.02) then
               write(3350, '(2i6,23(E16.4e2))')
     >              j, i, l_connect(jt),
     >              nz_tot_after_wfac(j,i), nz_on_ne(j,i),
     >              tot_rad_rate(j,i),
!     >              f_0(j,i), f_thi(j,i), f_tot_full(j,i),
     >              f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
     >              vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
     >              vva(j,i,ia),
     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
     >              mach_number(j,i),
     >              energy_dens_convective_para(j,i),
     >              energy_dens_total_para(j,i),
     >              knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
     >              heat_flux_para_gm(j,i)
            endif


! qGM rough check under qFL-solution
!    just test and comment out
            if(it.eq.24) then
!$$$               write(8888,'(2i,23(E16.4e2))')
!$$$     >              j, i, l_connect(jt),
!$$$     >              mfp_ion_effective(j,i), ch_length_ui(j,i),
!$$$     >              mfp_ion_effective(j,i)/ ch_length_ui(j,i),
!$$$     >              omega_gm(j,i),
!$$$     >              1.d0-1.2d0*(mach_number(j,i)**2.d0)

               write(4000,'(2i6,23(E16.4e2))')
     >              j, i, l_connect(jt),
     >              f_0(j,i), f_thi(j,i), f_the(j,i), f_elec(j,i),
     >              f_tot_simple(j,i), f_tot_full(j,i),
     >              flmxi(j,i), heat_flux_para_by_kappa_para(j,i),
     >              flmxe(j,i), heat_flux_para_elec(j,i)

            endif


         enddo
      enddo

      close(3010); close(3020); close(3030); close(3040); close(3050)
      close(3110); close(3120); close(3130); close(3140); close(3150)
      close(3210); close(3220); close(3230); close(3240); close(3250)
      close(3310); close(3320); close(3330); close(3340); close(3350)

      close(3321); close(3322); close(3323); close(3324)

      close(4000)
!

!just test
!      j=65; i=24
!      write(888,'(2i, 22(E16.4e2))')
!     >   j, i, vdnz(j,i,0,1), vdnz(j,i,1,1),vdnz(j,i,2,1),vdnz(j,i,3,1),
!     >     vdnz(j,i,4,1), vdnz(j,i,5,1),vdnz(j,i,6,1),vdnz(j,i,7,1),
!     >     vdnz(j,i,8,1), vdnz(j,i,9,1),vdnz(j,i,10,1),vdnz(j,i,11,1),
!     >   vdnz(j,i,12,1), vdnz(j,i,13,1),vdnz(j,i,14,1),vdnz(j,i,15,1),
!     >   vdnz(j,i,16,1), vdnz(j,i,17,1),vdnz(j,i,18,1),
!     > sum( vdnz(j,i, 0:ismaxL(i_imp_species), i_imp_species) ),
!     >     nz_tot_after_wfac(j,i)

      return
      end subroutine
!============================
!============================

!***********************************************************************
      subroutine yh_sld_kappa_interpolation_test
!***********************************************************************
      use cphcns, only : cev, cmp, cpi
      use cplcom, only : aion, aza, vna, vne, vte, vti, vzf
      use cplmet, only : hwtm, hwtp, icel, itmax, itmpe, itmps, jcel
     >    , jtmax, jtmin
      use csize,  only : ndx, ndy
      implicit none
!
!::local variables
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jt, jts, jte, j, jp, i, irg, ia, iy
!ik   integer  jchsl, icosl, icisl, i6
!ik   real*8   fcenh, fcslw
      integer  it, jt, jts, jte, j, jp, i, ia
!
!  connecting length
      real*8 l_connect(ndx)

      real*8 permittivity_eps_zero
      real*8 m_bulk
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 q_para_helander, q_para_helander_limited
!ik   real*8 tau_ii_soldor_file, tau_ii_helander
!ik   real*8 tau_ii_soldor_file_num, tau_ii_soldor_file_denom
      real*8 tau_ii_helander
      real*8 tau_ii_helander_num, tau_ii_helander_denom
      real*8 kappa_ii_para_helander
      real*8 kappa_helander_center(ndx,ndy)
      real*8 kappa_helander_jplushalf(ndx,ndy)

!      real*8 kappa_ii_helander, kappa_ii_soldor
      real*8 temp_zf, temp_ne, temp_te, temp_ti, temp_ni
      real*8 temp_mass_number, temp_abs_charge_number
      real*8 coulog, coulomb_logarithm
! deleted 4 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 clambda_wesson
!ik   real*8 alpha_helander, flux_limiter_helander
!ik   real*8 q_para_upper_limit
!ik   Real*8 dl_connect, grad_ti

!::interpolation method for diffusion at j+1/2
      real*8 zxicl
      real*8  fundav, fdfc, fdwm, fdwp, fddm, fddp
      fundav(fdfc,fdwm,fdwp,fddm,fddp) = fdfc/(fdwm/fddm+fdwp/fddp)

!::clear
!
!     Constants
      permittivity_eps_zero = 8.854187817d-12
!       Bulk plamsa ion species is D+
!       cmp = 1.67252e-27[kg] @ phcnst.f
      ia = 1
      m_bulk = aion(ia) * cmp

! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   alpha_helander = 0.5d0

!:: q para helander at cell center
      do it = 1, itmax
         jts = jtmin(it)
         jte = jtmax(it)
         if( it.ge.itmps .and. it.le.itmpe ) then
            jts = jts + 1
            jte = jte - 1
         endif
!
         do jt = jts, jte
            j  = jcel(jt,it)
            i  = icel(jt,it)
            jp  = jcel(jt+1,it)
!            irg = kreg(j,i)
!
!       q// helander evaluate (at cell center) for check
!       ***Attention: just rough check,
!           kappa --> cell center j, grad Ti --> cell boundary j+1/2
!           SONIC poloidal direction (outer to inner div) is defined to be POSITIVE.
!                                             2018/05/22 Y.Homma
            temp_zf = vzf(j,i) ; temp_ne = vne(j,i)
            temp_ni = vna(j,i,ia)
            temp_te = vte(j,i) ; temp_ti = vti(j,i)
!                       temp_te [eV], temp_ti[eV]
            temp_mass_number = aion(ia)
            temp_abs_charge_number = abs( aza(ia) )
            coulomb_logarithm = coulog( temp_zf, temp_ne,
     >           temp_te, temp_ti,
     >           temp_mass_number, temp_abs_charge_number )

!     Wesson Sec. 14.5 clambda for ion-ion collisions T<10keV
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik         clambda_wesson = 17.3d0 -0.5d0*log(temp_ne/1.0d20)
!ik  >           + 1.5d0*log(temp_ti/1000.d0)

!     Tau_ii
! deleted 6 lines organize local variables and include files by kamata 2021/05/31
!ik         tau_ii_soldor_file_num = 2.0853d13 * (temp_ti**1.5d0)
!ik  >           * sqrt(m_bulk / cmp)
!ik         tau_ii_soldor_file_denom = coulomb_logarithm * temp_ni
!ik  >           * (aza(ia)**4.d0)
!ik         tau_ii_soldor_file
!ik  >           = tau_ii_soldor_file_num / tau_ii_soldor_file_denom

!            tau_ii_helander_num = (12.d0/sqrt(2.d0)) * (cpi**1.5d0)
            tau_ii_helander_num = (12.d0) * (cpi**1.5d0)
     >           * (permittivity_eps_zero**2.d0)
     >           * sqrt(m_bulk) * ( (cev*temp_ti)**1.5d0 )
            tau_ii_helander_denom = coulomb_logarithm * temp_ni
     >           * ( (aza(ia)*cev )**4.d0)

            tau_ii_helander
     >           = tau_ii_helander_num / tau_ii_helander_denom

!     Kappa // classical definition Helander @ cell center (j,i)
!                           SOLDOR file definition is identical
            kappa_ii_para_helander = 3.9d0 * temp_ni * (cev*temp_ti)
     >           * tau_ii_helander / m_bulk

            kappa_helander_center(j,i) = kappa_ii_para_helander
         enddo
      enddo


!   Inter polation to cell boundary j+1/2
!        & data file writing
!     1D tube plot q_para


      open(3010,file='sld_kappa_interpolation_test_edge_42.dat')
      open(3020,file='sld_kappa_interpolation_test_sepx_24.dat')
      open(3030,file='sld_kappa_interpolation_test_sol_15.dat')
      open(3040,file='sld_kappa_interpolation_test_sol_05.dat')
      open(3050,file='sld_kappa_interpolation_test_outermost_02.dat')
!
      do it = 1, itmax
         jts = jtmin(it)
         jte = jtmax(it)
         if( it.ge.itmps .and. it.le.itmpe ) then
            jts = jts + 1
            jte = jte - 1
         endif
!
         l_connect = 0.d0
         call slenb(it,l_connect)
!
         do jt = jts, jte
            j  = jcel(jt,it)
            i  = icel(jt,it)
            jp  = jcel(jt+1,it)

!::d-coef Xi  at cell boundary (j+1/2,i)
!x    zxicl = 1.0d0/(hwtm(j,i)/xicl(j)+hwtp(j,i)/xicl(jp))
            zxicl =fundav(1.d0,hwtm(j,i),hwtp(j,i),
     >           kappa_helander_center(j,i),
     >           kappa_helander_center(jp,i) )

            kappa_helander_jplushalf(j,i) = zxicl

!     write
            if(it.eq.42) then
               write(3010, '(2i6,9(E16.4e2))')
     >              j, i, l_connect(jt),
     >              kappa_helander_center(j,i),
     >              kappa_helander_jplushalf(j,i)
            endif

            if(it.eq.24) then
               write(3020, '(2i6,9(E16.4e2))')
     >              j, i, l_connect(jt),
     >              kappa_helander_center(j,i),
     >              kappa_helander_jplushalf(j,i)
            endif

            if(it.eq.15) then
               write(3030, '(2i6,9(E16.4e2))')
     >              j, i, l_connect(jt),
     >              kappa_helander_center(j,i),
     >              kappa_helander_jplushalf(j,i)
            endif

            if(it.eq.05) then
               write(3040, '(2i6,9(E16.4e2))')
     >              j, i, l_connect(jt),
     >              kappa_helander_center(j,i),
     >              kappa_helander_jplushalf(j,i)
            endif

            if(it.eq.02) then
               write(3050, '(2i6,9(E16.4e2))')
     >              j, i, l_connect(jt),
     >              kappa_helander_center(j,i),
     >              kappa_helander_jplushalf(j,i)
            endif

         enddo
      enddo


      close(3010)
      close(3020)
      close(3030)
      close(3040)
      close(3050)

      return
      end subroutine
!============================
!============================

!***********************************************************************
      subroutine kappa_para_helander(temp_ni, temp_ti,
     >     coulomb_logarithm, kappa_ii_para_helander )
!                     temp_ni [m^{-3}], temp_ti[eV]
!***********************************************************************
      use cphcns, only : cev, cmp, cpi
      use cplcom, only : aion, aza
      implicit none
!
!:: in/out variables
      real*8, intent(in):: temp_ni, temp_ti
!                     temp_ni [m^{-3}], temp_ti[eV]
      real*8, intent(in)::  coulomb_logarithm
      real*8, intent(out):: kappa_ii_para_helander

!::local variables
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jt, jts, jte, j, jp, i, irg, ia, iy
!ik   integer  jchsl, icosl, icisl, i6
!ik   real*8   fcenh, fcslw
      integer  ia
!
      real*8 permittivity_eps_zero
      real*8 m_bulk
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 tau_ii_soldor_file, tau_ii_helander
!ik   real*8 tau_ii_soldor_file_num, tau_ii_soldor_file_denom
      real*8 tau_ii_helander
      real*8 tau_ii_helander_num, tau_ii_helander_denom

!::clear
!
!     Constants
      permittivity_eps_zero = 8.854187817d-12
!       Bulk plamsa ion species is D+
!       cmp = 1.67252e-27[kg] @ phcnst.f
      ia = 1
      m_bulk = aion(ia) * cmp

!     Tau_ii
!$$$            tau_ii_soldor_file_num = 2.0853d13 * (temp_ti**1.5d0)
!$$$     >           * sqrt(m_bulk / cmp)
!$$$            tau_ii_soldor_file_denom = coulomb_logarithm * temp_ni
!$$$     >           * (aza(ia)**4.d0)
!$$$            tau_ii_soldor_file
!$$$     >           = tau_ii_soldor_file_num / tau_ii_soldor_file_denom

!            tau_ii_helander_num = (12.d0/sqrt(2.d0)) * (cpi**1.5d0)
      tau_ii_helander_num = (12.d0) * (cpi**1.5d0)
     >     * (permittivity_eps_zero**2.d0)
     >     * sqrt(m_bulk) * ( (cev*temp_ti)**1.5d0 )
      tau_ii_helander_denom = coulomb_logarithm * temp_ni
     >     * ( (aza(ia)*cev )**4.d0)

      tau_ii_helander
     >     = tau_ii_helander_num / tau_ii_helander_denom

!     Kappa // classical definition Helander @ cell center (j,i)
!                           SOLDOR file definition is identical
      kappa_ii_para_helander = 3.9d0 * temp_ni * (cev*temp_ti)
     >     * tau_ii_helander / m_bulk


      return
      end subroutine
!============================
!============================


!***********************************************************************
      subroutine kappa_para_write
!***********************************************************************
      use cplcom, only : kappa_helander_test, xicl_test
      use cplmet, only : icel, itmax, itmpe, itmps, jcel,jtmax, jtmin
      use csize,  only : ndx
      implicit none

!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jt, jts, jte, j, jp, i, irg, ia, iy
      integer  it, jt, jts, jte, j, i
      real*8 l_connect(ndx)

      open(3010,file='sld_kappa_comparison_test_edge_42.dat')
      open(3020,file='sld_kappa_comparison_test_sepx_24.dat')
      open(3030,file='sld_kappa_comparison_test_sol_15.dat')
      open(3040,file='sld_kappa_comparison_test_sol_05.dat')
      open(3050,file='sld_kappa_comparison_test_outermost_02.dat')
!
      do it = 1, itmax
         jts = jtmin(it)
         jte = jtmax(it)
         if( it.ge.itmps .and. it.le.itmpe ) then
            jts = jts + 1
            jte = jte - 1
         endif
!
         l_connect = 0.d0
         call slenb(it,l_connect)
!
         do jt = jts, jte
            j  = jcel(jt,it)
            i  = icel(jt,it)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik         jp  = jcel(jt+1,it)

            if(it.eq.42) then
               write(3010, '(2i6,9(E16.4e2))')
     >              j, i, l_connect(jt),
     >              xicl_test(j,i), kappa_helander_test(j,i),
     >              xicl_test(j,i) / kappa_helander_test(j,i)
            endif

            if(it.eq.24) then
               write(3020, '(2i6,9(E16.4e2))')
     >              j, i, l_connect(jt), 
     >              xicl_test(j,i), kappa_helander_test(j,i),
     >              xicl_test(j,i) / kappa_helander_test(j,i)
            endif

            if(it.eq.15) then
               write(3030, '(2i6,9(E16.4e2))')
     >              j, i, l_connect(jt), 
     >              xicl_test(j,i), kappa_helander_test(j,i) ,
     >              xicl_test(j,i) / kappa_helander_test(j,i)
            endif

            if(it.eq.05) then
               write(3040, '(2i6,9(E16.4e2))')
     >              j, i, l_connect(jt), 
     >              xicl_test(j,i), kappa_helander_test(j,i) ,
     >              xicl_test(j,i) / kappa_helander_test(j,i)
            endif

            if(it.eq.02) then
               write(3050, '(2i6,9(E16.4e2))')
     >              j, i, l_connect(jt), 
     >              xicl_test(j,i), kappa_helander_test(j,i) ,
     >              xicl_test(j,i) / kappa_helander_test(j,i)
            endif

         enddo
      enddo


      close(3010)
      close(3020)
      close(3030)
      close(3040)
      close(3050)


      return
      end subroutine
!============================
!============================
