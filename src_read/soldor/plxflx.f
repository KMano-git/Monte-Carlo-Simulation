!***********************************************************************
      subroutine plxflx(crg,ix,kdr)
!***********************************************************************
!
!       creg  "odv", "idv", "sol", "man"
!         ix    poloidal mesh
!         kdr    +1 or -1
!
!           *--------*  j+1/2   kdr = +1
!           |        |
!           |        |   (j)
!           |        |
!           *--------*  j-1/2   kdr = -1
!
!
!-----------------------------------------------------------------------
      use cgdcom, only : grdx, grdy
      use cphcns, only : cpi
      use cplcom, only : ama, nion, vna, vne, vte, vti, vva
      use cplmet, only : hpit, hvsb, icmpe, icmps, icspx, icwl1, icwl2
      use cplpst, only : rhf_idp, rhf_imd, rhf_odp, rhf_omd
      use cplqcn, only : qfx_cd, qfx_cv, qfx_vh
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   character  crg*(*)
!ik   integer    ix, kdr
      character, intent(in) :: crg*(*)
      integer, intent(in) :: ix, kdr
!
!::local variables
      integer  jc, ics, ice, jsf, i6
      integer  ic, j, i, mx1, mx2, mx3, mx4, my1, my2, my3, my4
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  ia, m1a, m2a, m3, m4
!ik   real*8   smf1, smf3, smf4, x0, y0, x1, y1, x2, y2
      integer  ia, m1a, m3, m4
      real*8   smf1, smf3, smf4, x0, x1, y1, x2, y2
      real*8   dl, zsv, zsp, hcsw, hh, fcv3, fcv4, fcd3, fcd4
      real*8   f1, f3, f4, f34
      character  csgn*1, cdsn*20
!
!::mesh number (jc,ic)   KCSFUJI
      jc = ix
      if( crg.eq."idv" .or. crg.eq."odv" ) then
        ics = icwl1
        ice = icwl2
      elseif( crg.eq."sol" ) then
        ics = icwl1
        ice = icspx
      elseif( crg.eq."man" ) then
        ics = icmps
        ice = icmpe
      else
        call wexit("plxflx","invarid crg "//crg)
      endif
!
!::file name
      csgn = "p"
      if( kdr.lt.0 ) csgn = "m"
      write(cdsn,'(a,i3.3,a,a)') "Xflx_",ix,csgn,".txt"
      i6 = 21
      open(unit=i6,file=cdsn)
!
      write(i6,'(/2x,"*** plxflx ***  ",a,2x,"ix =",i3,"  kdr =",
     >  i2,2x,a)') crg, ix, kdr, trim(cdsn)
!
!::ion species
       ia  = 1
       m1a = 2*ia - 1
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik    m2a = 2*ia
       m3  = 2*nion + 1
       m4  = 2*nion + 2
!
!::jsf
      jsf = jc
      if( kdr.lt.0 ) jsf = jc-1
!
      write(i6,'(2x,"ix =",i5,"  jc =",i4,"  jsf =",i4,"  ics,ice =",
     >   2i4)') ix, jc, jsf, ics, ice
!
      smf1 = 0.0d0
      smf3 = 0.0d0
      smf4 = 0.0d0
!
      write(i6,'(3x,"ix",2x,"jsf",4x,"i",2x,"x1",7x,"y1",7x,
     >  "hcsw",5x,"hvsb",5x, "pitch",4x,
     >  "r_odp",4x,"r_idp",4x,"r_omd",4x,"r_imd",4x,
     >  "Ne",9x,"Ni",9x,"Te",9x,"Ti",9x,"Va",9x,
     >  "f1p",8x,"f3p",8x,"f4p",8x,"ftp",8x,
     >  "fcd3",7x,"fcv3",7x,"fcd4",7x,"fcv4")')
!
      do ic = ics, ice
      j = jc
      i = ic
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)  ! KSFUJI
      if( kdr.gt.0 ) then
      x1 = grdx(mx2,my2)
      y1 = grdy(mx2,my2)
      x2 = grdx(mx3,my3)
      y2 = grdy(mx3,my3)
      else
      x1 = grdx(mx1,my1)
      y1 = grdy(mx1,my1)
      x2 = grdx(mx4,my4)
      y2 = grdy(mx4,my4)
      endif
      x0 = 0.5d0*(x1+x2)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   y0 = 0.5d0*(y1+y2)
      dl = sqrt((x2-x1)**2+(y2-y1)**2)
!
      zsv = 2.0d0*cpi*x0*dl  ! S^theta
      zsp = hvsb(jsf,i)      ! S^theta*h*cosw
      if( zsv.eq.0.0d0 ) then
      hcsw = 1.0d0
      zsp  = 1.0d0
      else
      hcsw = zsp/zsv
      endif
      hh   = hpit(jsf,i)
!
      fcv3 = qfx_cv(jsf,i,m3)
      fcv4 = qfx_cv(jsf,i,m4)
      fcd3 = qfx_cd(jsf,i,m3) + qfx_vh(jsf,i)
      fcd4 = qfx_cd(jsf,i,m4)
      f1   = qfx_cv(jsf,i,m1a)/ama(ia)
      f3   = fcv3 + fcd3
      f4   = fcv4 + fcd4
      f34  = f3 + f4
      smf1 = smf1 + f1
      smf3 = smf3 + f3
      smf4 = smf4 + f4
!
      write(i6,'(2x,i3,2i5,0p9f9.5,1p15e11.3)')
     >  ix, jsf, i, x1, y1, hcsw, zsp, hh,
     >  rhf_odp(i), rhf_idp(i), rhf_omd(i), rhf_imd(i),
     >  vne(j,i), vna(j,i,ia), vte(j,i), vti(j,i), vva(j,i,ia),
     >  f1/zsp, f3/zsp, f4/zsp, f34/zsp,
     >  fcd3/zsp, fcv3/zsp, fcd4/zsp, fcv4/zsp
      enddo
!
      write(i6,'(2x)')
      write(i6,'(2x,"sum_f1p =",1pe11.3,"  sum_f3p =",1pe11.3,
     >  "  sum_f4p =",1pe11.3)') smf1, smf3, smf4
      close(i6)
!
!::summary
      write(n6,'(2x,a,2x,i3,2x,i2,2x,1p4e12.3)')
     >  "Xflx_"//crg, ix, kdr, smf1, smf3+smf4, smf3, smf4
!
      return
      end


!***********************************************************************
!     Heat flux & physics parameters along a magnetic flux tube  "it_flux_1"
      subroutine heat_flux_along_flux_tube(it_flux_1)
!***********************************************************************
      use cmeffz, only : vdnz
      use cntcom, only : mcel
      use cntmnt, only : sn0
      use cphcns, only : cev, cme, cmp
      use cplcom, only : aion, heat_flux_para_by_kappa_para
     >    , heat_flux_para_elec, flmxe, flmxi, vna, vnezef, vte, vti
     >    , vva
      use cplimp, only : ismaxl
      use cplmet, only : hpit, icel, itmax, jcel, jtmax, jtmin
      use cplwrd, only : wime
      use csize,  only :  ndx, ndy
      implicit none

!::Main variables
      integer, intent(in):: it_flux_1

!::local variables
      character(40) file_name_1
! modified 17/4 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jt, jts, jte, j, i, irg, ia, iy
!ik   integer  jp
!ik   integer  i_imp_species, iz, ic
!ik   integer  jchsl, icosl, icisl, i6
!ik   real*8   fcenh, fcslw
!ik   real*8 mx, my
!ik   real*8 br, bz, bt
!ik   real*8 b_pol(ndx,ndy), b_abs(ndx,ndy)
!ik   real*8 coulog, clambda, eps_zero
!ik   real*8 mi, ni, ni_p, ti_ev, ti_p_ev, vth_i
!ik   real*8 te_ev, te_p_ev, charge_imp
!ik   real*8 ui, ui_p, dl_connect, sound_speed_ion
!ik   real*8 field_elec_para(ndx,ndy)
!ik   real*8 mi, ni, ti_ev
!ik   real*8 tau_ii_helander_num, tau_ii_helander_denom
!ik   real*8 tau_ii_helander(ndx,ndy)
!ik   real*8 mfp_ii(ndx,ndy)
      integer  it, jt, jts, jte, j, i, ia
      integer  i_imp_species, ic
      real*8 mi
      real*8 sound_speed_ion
      real*8 mach_number(ndx,ndy), cs_ion(ndx,ndy)
! deleted 10 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 grad_ti_ev(ndx,ndy), grad_te_ev(ndx,ndy)
!ik   real*8 grad_pe(ndx,ndy)
!ik   real*8 grad_ni(ndx,ndy)
!ik   real*8 grad_ui(ndx,ndy), grad_b_abs(ndx,ndy)
!ik   real*8 ch_length_ti_ev(ndx,ndy), ch_length_ni(ndx,ndy)
!ik   real*8 ch_length_ui(ndx,ndy), ch_length_b_abs(ndx,ndy)
!ik   real*8 knudsen_ti_ev(ndx,ndy), knudsen_ni(ndx,ndy)
!ik   real*8 knudsen_ti_ev_effective(ndx,ndy)
!ik   real*8 knudsen_ui(ndx,ndy), knudsen_b_abs(ndx,ndy)
!ik   real*8 knudsen_connect(ndx,ndy), abs_knudsen_max(ndx,ndy)
      real*8 pressure_e(ndx,ndy), pressure_i(ndx,ndy)
      real*8 pressure_dyn_e(ndx,ndy), pressure_dyn_i(ndx,ndy)
      real*8 nz_tot_after_wfac(ndx,ndy), nz_on_ne(ndx,ndy)
      real*8 nz_ion_after_wfac(ndx,ndy)
      real*8 tot_rad_rate(ndx,ndy)
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 f_thi(ndx,ndy), f_the(ndx,ndy),  f_0(ndx,ndy)
!ik   real*8 f_elec(ndx,ndy), f_tot_simple(ndx,ndy), f_tot_full(ndx,ndy)
      real*8 energy_dens_total_para(ndx,ndy)
      real*8 energy_dens_convective_para(ndx,ndy)
      real*8 energy_dens_elec_conv_para(ndx,ndy)
      real*8 energy_dens_elec_para_tot(ndx,ndy)

      real*8 particle_flux_dens_para_i(ndx,ndy)

! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 mfp_ion_effective(ndx,ndy)
!ik   real*8 omega_gm(ndx,ndy), flmxi_qgm(ndx,ndy)
!ik   real*8 heat_flux_para_gm(ndx,ndy)

!
!  connecting length / poloidal lentgth
      real*8 l_connect(ndx), l_poloidal(ndx)

!Initail inputs/specifications
!:: pm3d
      !background bulk species specify
      ia = 1
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   eps_zero = 8.854187817d-12
      mi = aion(ia) * cmp

      i_imp_species = 1 !  Index to specify impurity species(NOT Zimp!)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   charge_imp = 4.d0 * cev !  Impurity charge supposed q_z[C] )
!      write(2000,*) "test: ami/cmp, amz/cmp"
!      write(2000,*)  ami/cmp, amz/cmp


!::zero clear
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik   ch_length_ti_ev=0.d0; ch_length_ni=0.d0; ch_length_ui=0.d0
!ik   ch_length_b_abs =0.d0
!ik   knudsen_ti_ev=0.d0; knudsen_ni=0.d0; knudsen_ui=0.d0
!ik   knudsen_b_abs=0.d0; knudsen_connect=0.d0
!ik   knudsen_ti_ev_effective=0.d0

      nz_tot_after_wfac=0.d0; nz_on_ne=0.d0
      tot_rad_rate=0.d0
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   field_elec_para=0.d0

! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   f_tot_simple=0.d0; f_tot_full=0.d0
!ik   f_0=0.d0; f_thi=0.d0; f_the=0.d0; f_elec=0.d0

! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   mfp_ion_effective=0.d0; omega_gm=0.d0; flmxi_qgm=0.d0
!ik   heat_flux_para_gm=0.d0


      !Physics quantities evaluation
!      do it = 1, itmax
      do it = 2, itmax-1
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
! added 1 line replace all include files with module files by kamata 2021/08/18
! delete 1 line by yamamoto 2023/02/24 
!            if( ic <= 0 ) cycle
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik         if(jt .eq. jte) then
!ik            jp = j
!ik         else
!ik            jp  = jcel(jt+1,it)
!ik         endif
!     Coulomb logarithm value
! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik         clambda = coulog(
!ik  >           vzf(j,i), vne(j,i), vte(j,i), vti(j,i),
!ik  >           aion(ia), abs( aza(ia) )  )
!     Magnetic field  (*** at the corner 2, boundary j+1/2 ***)
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik         mx = kgdx(j,i,2); my = kgdy(j,i,2)
!ik         br = hbr(mx,my); bz = hbz(mx,my)
!ik         bt = hbt(mx,my)
!ik         b_pol(j,i) = sqrt(br*br + bz*bz)
!ik         b_abs(j,i) = sqrt(br*br + bz*bz + bt*bt)
!     tau_ii
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik         ti_ev = vti(j,i);  ni = vna(j,i,ia)
!            tau_ii_helander_num = (12.d0/sqrt(2.d0)) * (cpi**1.5d0)
! deleted 7 lines organize local variables and include files by kamata 2021/05/31
!ik         tau_ii_helander_num = (12.d0) * (cpi**1.5d0)
!ik  >           * (eps_zero**2.d0)
!ik  >           * sqrt(mi) * ( (cev*ti_ev)**1.5d0 )
!ik         tau_ii_helander_denom = clambda * ni
!ik  >           * ( (aza(ia)*cev )**4.d0)
!ik         tau_ii_helander(j,i)
!ik  >           = tau_ii_helander_num / tau_ii_helander_denom
!     MFP i-i collision
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik         vth_i = sqrt(cev*ti_ev/mi)
!ik         mfp_ii(j,i) = vth_i * tau_ii_helander(j,i)
!     pressures
            pressure_e(j,i) = vnezef(j,i) * vte(j,i) * cev
            pressure_i(j,i) = vna(j,i,ia) * vti(j,i) *cev
            pressure_dyn_e(j,i) = 0.5d0 * cme * vnezef(j,i)
     >           * (vva(j,i,ia)**2.d0) ! dynamic elec. pressure
            pressure_dyn_i(j,i) = 0.5d0 * mi * vna(j,i,ia)
     >           * (vva(j,i,ia)**2.d0) ! dynamic ion pressure
!     Impurity density and radiation
            nz_tot_after_wfac(j,i)
     >       = sum( vdnz(j,i, 0:ismaxL(i_imp_species), i_imp_species) )
            nz_ion_after_wfac(j,i)
     >       = sum( vdnz(j,i, 1:ismaxL(i_imp_species), i_imp_species) )
            nz_on_ne(j,i) = nz_tot_after_wfac(j,i) / vnezef(j,i)
!               vnezef :=  ne includes elec. from ionized impurities.
            if(ic>0)then
              tot_rad_rate(j,i) = -wime(ic)
            endif
!     Mach number, same sign as ui//=vva is retained.
            sound_speed_ion = sqrt(cev*( vte(j,i) + vti(j,i)  ) / mi)
            mach_number(j,i) = vva(j,i,ia) / sound_speed_ion
            cs_ion(j,i) = sound_speed_ion
!     ION: Particle para flux density (physics definition along B,
!          the mesh shapes and metrics are not taken into account.)
            particle_flux_dens_para_i(j,i) = vna(j,i,ia)*vva(j,i,ia)
!     ION: Convective & Total energy flux density
            energy_dens_convective_para(j,i) =
     >           2.5d0*vna(j,i,ia)*(cev*vti(j,i))*vva(j,i,ia) +
     >           0.5d0*mi*vna(j,i,ia)*(vva(j,i,ia)**3.d0)
            energy_dens_total_para(j,i) =
     >           energy_dens_convective_para(j,i) +
     >           heat_flux_para_by_kappa_para(j,i)
!     ELECTRON: Convective & Total energy flux density
!               ui ~ ue supposed in SOLDOR as zero current.
            energy_dens_elec_conv_para(j,i) =
     >           2.5d0*pressure_e(j,i)*vva(j,i,ia) +
     >           0.5d0*cme*vnezef(j,i)*(vva(j,i,ia)**3.d0)
            energy_dens_elec_para_tot(j,i) =
     >           energy_dens_elec_conv_para(j,i) +
     >           heat_flux_para_elec(j,i)
         enddo                  !            do jt = jts, jte
      enddo     !      do it = 1, itmax


!     EVALUATE Grad Ti-ni-ui-B
!              Collisionalities
!              SONIC's static-elec. field
!              Forces
!      do it = 1, itmax
! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik   do it = 2, itmax-1
!ik           l_connect = 0.d0
!ik      call slenb(it,l_connect)

! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik      jts = jtmin(it); jte = jtmax(it)
!         if( it.ge.itmps .and. it.le.itmpe ) then
!            jts = jts + 1
!            jte = jte - 1
!         endif
!
! deleted 8 lines organize local variables and include files by kamata 2021/05/31
!ik      do jt = jts, jte
!ik         j  = jcel(jt,it)
!ik         i  = icel(jt,it)
!ik         if(jt .eq. jte) then
!ik            jp = j
!ik         else
!ik            jp  = jcel(jt+1,it)
!ik         endif
!     Coulomb logarithm value
! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik         clambda = coulog(
!ik  >           vzf(j,i), vne(j,i), vte(j,i), vti(j,i),
!ik  >           aion(ia), abs( aza(ia) )  )

!           Grad Ti evaluation @ cell boundary (j+1/2,i), same as q// helander.
!           SONIC poloidal direction (outer to inner div) is defined to be POSITIVE.
!            grad_ti = ( vti(jp,i) - vti(j,i) )*cev / dl_connect
!                                               [joule/m]
! deleted 7 lines organize local variables and include files by kamata 2021/05/31
!ik         ti_ev = vti(j,i);  ti_p_ev = vti(jp,i)
!ik         te_ev = vte(j,i);  te_p_ev = vte(jp,i)
!ik         ni = vna(j,i,ia);  ni_p = vna(jp,i,ia)
!ik         ui = vva(j,i,ia);  ui_p = vva(jp,i,ia)
!ik         if(hpit(j,i).le.0.d0) cycle
!ik         dl_connect = ( hdxm(j,i) + hdxp(j,i) ) / hpit(j,i)
!ik         if(dl_connect.le.0.d0) cycle

! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik         grad_ti_ev(j,i) =  ( ti_p_ev - ti_ev ) / dl_connect
!ik         grad_te_ev(j,i) =  ( te_p_ev - te_ev ) / dl_connect
!                                         [eV/m]
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik         grad_ni(j,i) =  ( ni_p - ni ) / dl_connect
!ik         grad_ui(j,i) =  ( ui_p - ui ) / dl_connect
!ik         grad_pe(j,i) =  ( pressure_e(jp,i) - pressure_e(j,i) ) /
!ik  >           dl_connect
!ik         grad_b_abs(j,i) = ( b_abs(jp,i) - b_abs(j,i) ) / dl_connect

!     Electric field (static) (cf. Note p.imp6 OR soniv doc. Eq.(2.2-11))
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik         field_elec_para(j,i) = -grad_pe(j,i)/(vnezef(j,i)*cev) -
!ik  >           0.71d0*( grad_te_ev(j,i)*cev )/cev

!     Characteristic length evaluation
! deleted 8 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(grad_ti_ev(j,i) ).gt.0.d0)  ch_length_ti_ev(j,i)
!ik  >           = ti_ev / grad_ti_ev(j,i)
!ik         if(abs(grad_ni(j,i) ).gt.0.d0)  ch_length_ni(j,i)
!ik  >           = ni / grad_ni(j,i)
!ik         if(abs(grad_ui(j,i) ).gt.0.d0)  ch_length_ui(j,i)
!ik  >           = ui / grad_ui(j,i)
!ik         if(abs(grad_b_abs(j,i) ).gt.0.d0)  ch_length_b_abs(j,i)
!ik  >           = b_abs(j,i) / grad_b_abs(j,i)

!     Knudsen numbers
! deleted 12 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(ch_length_ti_ev(j,i)).gt.0.d0) knudsen_ti_ev(j,i)
!ik  >           = mfp_ii(j,i) / ch_length_ti_ev(j,i)
!ik         if(abs(ch_length_ni(j,i)).gt.0.d0) knudsen_ni(j,i)
!ik  >           = mfp_ii(j,i) / ch_length_ni(j,i)
!ik         if(abs(ch_length_ui(j,i)).gt.0.d0) knudsen_ui(j,i)
!ik  >           = mfp_ii(j,i) / ch_length_ui(j,i)
!ik         if(abs(ch_length_b_abs(j,i)).gt.0.d0) knudsen_b_abs(j,i)
!ik  >           = mfp_ii(j,i) / ch_length_b_abs(j,i)
!ik         knudsen_connect(j,i) = mfp_ii(j,i) / l_connect(jte-1)
!ik         abs_knudsen_max(j,i) = max( abs(knudsen_ti_ev(j,i) ),
!ik  >           abs(knudsen_ni(j,i) ), abs(knudsen_ui(j,i) ),
!ik  >        abs(knudsen_b_abs(j,i) ), abs(knudsen_connect(j,i) ) )

!           Effective Kn_Ti
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(flmxi(j,i) ).gt.0.d0) knudsen_ti_ev_effective(j,i)
!ik  >           = ( (1.d0/flmxi(j,i)) -1.d0) * flimi / 3.9d0

!     qGM rough estimation
!         ion_MFP_effective(i-i & i-Z colls included, evaluated from SOLDOR info.)
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik         mfp_ion_effective(j,i) =
!ik  >       knudsen_ti_ev_effective(j,i) * abs(ch_length_ti_ev(j,i) )
!         Omega_GM
! deleted 7 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(ch_length_ui(j,i)).gt.0.d0) then
!ik            if(abs(ch_length_b_abs(j,i)).gt.0.d0) then
!ik               omega_gm(j,i) =
!ik  >              5.88d0 * mach_number(j,i) * mfp_ion_effective(j,i) *
!ik  >     (1.d0 / ch_length_ui(j,i)  - 0.25d0 / ch_length_b_abs(j,i) )
!ik            endif
!ik         endif
!         flmxi_qGM
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(1.d0+omega_gm(j,i) ).gt.0.d0) then
!ik            flmxi_qgm(j,i) =
!ik  >           (1.d0 - 1.2d0*mach_number(j,i)*mach_number(j,i) ) /
!ik  >           (1.d0 + omega_gm(j,i) )
!ik         endif
!         qGM (heat_flux_para_GM)
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(flmxi(j,i)).gt.0.d0) then
!ik            heat_flux_para_gm(j,i) =
!ik  >           (heat_flux_para_by_kappa_para(j,i)/flmxi(j,i) ) *
!ik  >           flmxi_qgm(j,i)
!ik         endif

! deleted 1 line organize local variables and include files by kamata 2021/11/18
!ik         if(limp==3)then
!     Force balance
! deleted 4 lines organize local variables and include files by kamata 2021/05/31
!ik         call thermal_force_para_by_q_para(clambda,
!ik  >           ami, amz, 1.d0*cev, charge_imp, ti_ev,
!ik  >           0.d0, sqrt(2.d0 * cev * ti_ev /amz), ui,
!ik  >           heat_flux_para_by_kappa_para(j,i), f_thi(j,i) )
!!!
!!!     JUST TEST 2018/11/29  analysis for Fthi W/O collisionality effect
!            if(abs(flmxi(j,i) ).gt.0.d0) then
!               f_thi(j,i) = f_thi(j,i) / flmxi(j,i)
!            else
!               f_thi(j,i) = 0.d0
!            endif
!!!
! deleted 4 lines organize local variables and include files by kamata 2021/05/31
!ik         call thermal_force_para_by_q_para(clambda,
!ik  >           cme, amz, (-cev), charge_imp, te_ev,
!ik  >           0.d0, sqrt(2.d0 * cev * ti_ev /amz), ui,
!ik  >           heat_flux_para_elec(j,i), f_the(j,i) )

!      subroutine thermal_force_para_by_q_para(coulomb_logarithm,
!     >     mi, mz, ei[C], ez[C], ti[eV], vz_para, vz_perp, ui, q_para, fth_para)
! deleted 4 lines organize local variables and include files by kamata 2021/05/31
!ik         call friction_force_theoretical_value(clambda,
!ik  >           ami, amz, 1.d0*cev, charge_imp, ti_ev,
!ik  >           ni,
!ik  >           0.d0, sqrt(2.d0 * cev * ti_ev /amz), ui, f_0(j,i) )
!      subroutine friction_force_theoretical_value(coulomb_logarithm,
!     >     mi, mz, ei[C], ez[C], ti[eV], ni[m^(-3)], vz_para, vz_perp, ui, friction_para)
! deleted 1 line organize local variables and include files by kamata 2021/11/18
!ik         endif !limp

!     Electric field (static) force
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik         f_elec(j,i) = charge_imp  * field_elec_para(j,i)

!     Resultant force
! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik         f_tot_simple(j,i) = f_0(j,i) + f_thi(j,i)
!ik         f_tot_full(j,i) = f_0(j,i) + f_thi(j,i) + f_the(j,i) +
!ik  >           f_elec(j,i)

! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik      enddo                  !            do jt = jts, jte
!ik   enddo                     !      do it = 1, itmax
!
!
!===========================================
!     FILE WRITING
!===========================================
!     1D tube plot q_para

      write(file_name_1,'(a,i3.3,a)')
     >     'heat_flux_para_it_', it_flux_1, '.dat'

      open(3000,file=trim(file_name_1), status='replace')

!      it_count: do it = 1, itmax
      it_count: do it = 2, itmax-1
         l_connect = 0.d0
         l_poloidal = 0.d0
         jts = jtmin(it)
         jte = jtmax(it)
!         if( it.ge.itmps .and. it.le.itmpe ) then
!            jts = jts + 1
!            jte = jte - 1
!         endif

         if(it.eq.it_flux_1) then
            call slenb(it,l_connect)
            call plenc(it,l_poloidal)
!
            jt_count: do jt = jts, jte
!            jt_count: do jt = jts+1, jte-1 ! to remove extreme values at the ends
              j  = jcel(jt,it)
              i  = icel(jt,it)
!             irg = kreg(j,i)
!
              write(3000, '(4i6,45(E16.4e2))')
     >             jt, it, j, i,
     >             l_connect(jt), l_connect(jte-1) - l_connect(jt),
     >             l_poloidal(jt), l_poloidal(jte-1) - l_poloidal(jt),
     >             hpit(j,i), vnezef(j,i), vna(j,i,ia),
     >             vte(j,i), vti(j,i), vva(j,i,ia),
     >             sn0(j,i,ia),  ! neutral density [m-3]
     >             nz_ion_after_wfac(j,i), nz_tot_after_wfac(j,i),
     >          nz_ion_after_wfac(j,i)/vnezef(j,i), nz_on_ne(j,i),
     >             energy_dens_elec_conv_para(j,i),  ! q//_e_convective[W/m2]
     >             heat_flux_para_elec(j,i), ! q//_e_conductive[W/m2]
     >             energy_dens_convective_para(j,i), ! q//_i_convective[W/m2]
     >             heat_flux_para_by_kappa_para(j,i),  ! q//_i_conduction [W/m2]
     >             energy_dens_elec_para_tot(j,i), ! q//_e_cond+conv [W/m2]
     >             energy_dens_total_para(j,i), ! q//_i_cond+conv [W/m2]
     >   energy_dens_total_para(j,i) + energy_dens_elec_para_tot(j,i),
     >             particle_flux_dens_para_i(j,i),  ! Gamma_para_ion [s-1m-2]
     >             pressure_e(j,i), pressure_i(j,i),  ! Static pressure e & i [Pa]
     >             pressure_e(j,i)+pressure_i(j,i),
     >             cs_ion(j,i), mach_number(j,i),
     >             pressure_dyn_e(j,i)+pressure_dyn_i(j,i), ! dynamic pressure i+e [Pa]
     >             pressure_dyn_e(j,i)+pressure_dyn_i(j,i)
     >             +pressure_e(j,i)+pressure_i(j,i), ! Total plasma pressure [Pa]
     >             flmxe(j,i), flmxi(j,i),  ! heat flux limiter elec. & ion
     >             tot_rad_rate(j,i)        ! impurity radiation density[W/m3]

!$$$     >              nz_tot_after_wfac(j,i), nz_on_ne(j,i),
!$$$     >              tot_rad_rate(j,i),
!$$$!     >              f_0(j,i), f_thi(j,i), f_tot_full(j,i),
!$$$     >              f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
!$$$     >              vnezef(j,i), vna(j,i,ia),
!$$$     >              vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
!$$$     >              vva(j,i,ia),
!$$$     >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
!$$$     >              mach_number(j,i),
!$$$     >              energy_dens_convective_para(j,i),
!$$$     >              energy_dens_total_para(j,i),
!$$$     >              knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
!$$$     >              heat_flux_para_gm(j,i)

           enddo jt_count
        endif
      enddo it_count

      close(3000)

!
      return
      end subroutine
!***********************************************************************
!***********************************************************************
!***********************************************************************
!     Heat flux & physics parameters radial plot data output at "jcel(jt)"
      subroutine heat_flux_radial_output(flg_radial, flg_write_all,
     >     jcel_input)
!***********************************************************************
!     If flg_radial>0, radial profile data file on jcel(jt_input) is outputted.
!     If flg_write_all>0, physics parameters over all (j,i) points are outputted.
!     (flg_raidal, flg_write_all) must be eigher (1,0) or (0,1). Otherwise, error.
!
!     As flag (& flag checking), we use explicitely flg_radial and flg_write_all,
!     and no default function is supposed.
!     We can realize same function without flg_write_all, if a default fuction is supposed.
!     But it seems complicated, we did not choose the latter idea.
!
      use cgdcom, only : grdx, grdy
      use cmeffz, only : vdnz
      use cntcom, only : mcel
      use cntmnt, only : sn0
      use cphcns, only : cev, cme, cmp
      use cplcom, only : aion, flmxe, flmxi
     >    , heat_flux_para_by_kappa_para, heat_flux_para_elec, vna
     >    , vnezef, vte, vti, vva
      use cplimp, only : ismaxl
      use cplmet, only : hpit, icel, itmax, jcel, jtmax, jtmin, kgdx
     >    , kgdy
      use cplpst, only : rhf_idp, rhf_imd, rhf_odp, rhf_omd ! radial projected position
      use cplwrd, only : wime
      use csize,  only : ndx, ndy
      use cunit,  only : n6
      implicit none

!::Main variables
      integer, intent(in):: flg_radial, flg_write_all, jcel_input

!::local variables
      character(40) file_name_1
! modified 16/4 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, jt, jts, jte, j, i, irg, ia, iy
!ik   integer  jp
!ik   integer  i_imp_species, iz, ic
!ik   integer  jchsl, icosl, icisl, i6
!ik   real*8   fcenh, fcslw
!ik   real*8 mx, my
!ik   real*8 br, bz, bt
!ik   real*8 b_pol(ndx,ndy), b_abs(ndx,ndy)
!ik   real*8 coulog, clambda, eps_zero
!ik   real*8 mi, ni, ni_p, ti_ev, ti_p_ev, vth_i
!ik   real*8 te_ev, te_p_ev, charge_imp
!ik   real*8 ui, ui_p, dl_connect, sound_speed_ion
!ik   real*8 field_elec_para(ndx,ndy)
!ik   real*8 tau_ii_helander_num, tau_ii_helander_denom
!ik   real*8 tau_ii_helander(ndx,ndy)
!ik   real*8 mfp_ii(ndx,ndy)
      integer  it, jt, jts, jte, j, i, ia
      integer  i_imp_species, ic
      real*8 mi
      real*8 sound_speed_ion
      real*8 mach_number(ndx,ndy), cs_ion(ndx,ndy)
! deleted 10 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 grad_ti_ev(ndx,ndy), grad_te_ev(ndx,ndy)
!ik   real*8 grad_pe(ndx,ndy)
!ik   real*8 grad_ni(ndx,ndy)
!ik   real*8 grad_ui(ndx,ndy), grad_b_abs(ndx,ndy)
!ik   real*8 ch_length_ti_ev(ndx,ndy), ch_length_ni(ndx,ndy)
!ik   real*8 ch_length_ui(ndx,ndy), ch_length_b_abs(ndx,ndy)
!ik   real*8 knudsen_ti_ev(ndx,ndy), knudsen_ni(ndx,ndy)
!ik   real*8 knudsen_ti_ev_effective(ndx,ndy)
!ik   real*8 knudsen_ui(ndx,ndy), knudsen_b_abs(ndx,ndy)
!ik   real*8 knudsen_connect(ndx,ndy), abs_knudsen_max(ndx,ndy)
      real*8 pressure_e(ndx,ndy), pressure_i(ndx,ndy)
      real*8 pressure_dyn_e(ndx,ndy), pressure_dyn_i(ndx,ndy)
      real*8 nz_tot_after_wfac(ndx,ndy), nz_on_ne(ndx,ndy)
      real*8 nz_ion_after_wfac(ndx,ndy)
      real*8 tot_rad_rate(ndx,ndy)
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 f_thi(ndx,ndy), f_the(ndx,ndy),  f_0(ndx,ndy)
!ik   real*8 f_elec(ndx,ndy), f_tot_simple(ndx,ndy), f_tot_full(ndx,ndy)
      real*8 energy_dens_total_para(ndx,ndy)
      real*8 energy_dens_convective_para(ndx,ndy)
      real*8 energy_dens_elec_conv_para(ndx,ndy)
      real*8 energy_dens_elec_para_tot(ndx,ndy)

      real*8 particle_flux_dens_para_i(ndx,ndy)

! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8 mfp_ion_effective(ndx,ndy)
!ik   real*8 omega_gm(ndx,ndy), flmxi_qgm(ndx,ndy)
!ik   real*8 heat_flux_para_gm(ndx,ndy)

!
!  connecting length / poloidal lentgth
      real*8 l_connect(ndx), l_poloidal(ndx)

!Initail inputs/specifications
!:: pm3d
      !background bulk species specify
      ia = 1
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   eps_zero = 8.854187817d-12
      mi = aion(ia) * cmp

      i_imp_species = 1 !  Index to specify impurity species(NOT Zimp!)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   charge_imp = 4.d0 * cev !  Impurity charge supposed q_z[C] )
!      write(2000,*) "test: ami/cmp, amz/cmp"
!      write(2000,*)  ami/cmp, amz/cmp


!::zero clear
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik   ch_length_ti_ev=0.d0; ch_length_ni=0.d0; ch_length_ui=0.d0
!ik   ch_length_b_abs =0.d0
!ik   knudsen_ti_ev=0.d0; knudsen_ni=0.d0; knudsen_ui=0.d0
!ik   knudsen_b_abs=0.d0; knudsen_connect=0.d0
!ik   knudsen_ti_ev_effective=0.d0

      nz_tot_after_wfac=0.d0; nz_on_ne=0.d0
      tot_rad_rate=0.d0
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   field_elec_para=0.d0

! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   f_tot_simple=0.d0; f_tot_full=0.d0
!ik   f_0=0.d0; f_thi=0.d0; f_the=0.d0; f_elec=0.d0

! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   mfp_ion_effective=0.d0; omega_gm=0.d0; flmxi_qgm=0.d0
!ik   heat_flux_para_gm=0.d0


      !Physics quantities evaluation
!      do it = 1, itmax
      do it = 2, itmax-1
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
! added 1 line replace all include files with module files by kamata 2021/08/18
! delete 1 line by yamamoto 2023/02/24
!            if( ic <= 0 ) cycle
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik         if(jt .eq. jte) then
!ik            jp = j
!ik         else
!ik            jp  = jcel(jt+1,it)
!ik         endif
!     Coulomb logarithm value
! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik         clambda = coulog(
!ik  >           vzf(j,i), vne(j,i), vte(j,i), vti(j,i),
!ik  >           aion(ia), abs( aza(ia) )  )
!     Magnetic field  (*** at the corner 2, boundary j+1/2 ***)
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik         mx = kgdx(j,i,2); my = kgdy(j,i,2)
!ik         br = hbr(mx,my); bz = hbz(mx,my)
!ik         bt = hbt(mx,my)
!ik         b_pol(j,i) = sqrt(br*br + bz*bz)
!ik         b_abs(j,i) = sqrt(br*br + bz*bz + bt*bt)
!     tau_ii
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik         ti_ev = vti(j,i);  ni = vna(j,i,ia)
!            tau_ii_helander_num = (12.d0/sqrt(2.d0)) * (cpi**1.5d0)
! deleted 7 lines organize local variables and include files by kamata 2021/05/31
!ik         tau_ii_helander_num = (12.d0) * (cpi**1.5d0)
!ik  >           * (eps_zero**2.d0)
!ik  >           * sqrt(mi) * ( (cev*ti_ev)**1.5d0 )
!ik         tau_ii_helander_denom = clambda * ni
!ik  >           * ( (aza(ia)*cev )**4.d0)
!ik         tau_ii_helander(j,i)
!ik  >           = tau_ii_helander_num / tau_ii_helander_denom
!     MFP i-i collision
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik         vth_i = sqrt(cev*ti_ev/mi)
!ik         mfp_ii(j,i) = vth_i * tau_ii_helander(j,i)
!     pressures
            pressure_e(j,i) = vnezef(j,i) * vte(j,i) * cev
            pressure_i(j,i) = vna(j,i,ia) * vti(j,i) *cev
            pressure_dyn_e(j,i) = 0.5d0 * cme * vnezef(j,i)
     >           * (vva(j,i,ia)**2.d0) ! dynamic elec. pressure
            pressure_dyn_i(j,i) = 0.5d0 * mi * vna(j,i,ia)
     >           * (vva(j,i,ia)**2.d0) ! dynamic ion pressure
!     Impurity density and radiation
            nz_tot_after_wfac(j,i)
     >       = sum( vdnz(j,i, 0:ismaxL(i_imp_species), i_imp_species) )
            nz_ion_after_wfac(j,i)
     >       = sum( vdnz(j,i, 1:ismaxL(i_imp_species), i_imp_species) )
            nz_on_ne(j,i) = nz_tot_after_wfac(j,i) / vnezef(j,i)
!               vnezef :=  ne includes elec. from ionized impurities.
            if(ic>0) then
              tot_rad_rate(j,i) = -wime(ic)
            endif
!     Mach number, same sign as ui//=vva is retained.
            sound_speed_ion = sqrt(cev*( vte(j,i) + vti(j,i)  ) / mi)
            mach_number(j,i) = vva(j,i,ia) / sound_speed_ion
            cs_ion(j,i) = sound_speed_ion
!     ION: Particle para flux density (physics definition along B,
!          the mesh shapes and metrics are not taken into account.)
            particle_flux_dens_para_i(j,i) = vna(j,i,ia)*vva(j,i,ia)
!     ION: Convective & Total energy flux density
            energy_dens_convective_para(j,i) =
     >           2.5d0*vna(j,i,ia)*(cev*vti(j,i))*vva(j,i,ia) +
     >           0.5d0*mi*vna(j,i,ia)*(vva(j,i,ia)**3.d0)
            energy_dens_total_para(j,i) =
     >           energy_dens_convective_para(j,i) +
     >           heat_flux_para_by_kappa_para(j,i)
!     ELECTRON: Convective & Total energy flux density
!               ui ~ ue supposed in SOLDOR as zero current.
            energy_dens_elec_conv_para(j,i) =
     >           2.5d0*pressure_e(j,i)*vva(j,i,ia) +
     >           0.5d0*cme*vnezef(j,i)*(vva(j,i,ia)**3.d0)
            energy_dens_elec_para_tot(j,i) =
     >           energy_dens_elec_conv_para(j,i) +
     >           heat_flux_para_elec(j,i)
         enddo                  !            do jt = jts, jte
      enddo     !      do it = 1, itmax


!     EVALUATE Grad Ti-ni-ui-B
!              Collisionalities
!              SONIC's static-elec. field
!              Forces
!      do it = 1, itmax
! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik   do it = 2, itmax-1
!ik           l_connect = 0.d0
!ik      call slenb(it,l_connect) !FIXME:: 0div reported SYamoto

! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik      jts = jtmin(it); jte = jtmax(it)
!         if( it.ge.itmps .and. it.le.itmpe ) then
!            jts = jts + 1
!            jte = jte - 1
!         endif
!
! deleted 8 lines organize local variables and include files by kamata 2021/05/31
!ik      do jt = jts, jte
!ik         j  = jcel(jt,it)
!ik         i  = icel(jt,it)
!ik         if(jt .eq. jte) then
!ik            jp = j
!ik         else
!ik            jp  = jcel(jt+1,it)
!ik         endif
!     Coulomb logarithm value
! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik         clambda = coulog(
!ik  >           vzf(j,i), vne(j,i), vte(j,i), vti(j,i),
!ik  >           aion(ia), abs( aza(ia) )  )

!           Grad Ti evaluation @ cell boundary (j+1/2,i), same as q// helander.
!           SONIC poloidal direction (outer to inner div) is defined to be POSITIVE.
!            grad_ti = ( vti(jp,i) - vti(j,i) )*cev / dl_connect
!                                               [joule/m]
! deleted 7 lines organize local variables and include files by kamata 2021/05/31
!ik         ti_ev = vti(j,i);  ti_p_ev = vti(jp,i)
!ik         te_ev = vte(j,i);  te_p_ev = vte(jp,i)
!ik         ni = vna(j,i,ia);  ni_p = vna(jp,i,ia)
!ik         ui = vva(j,i,ia);  ui_p = vva(jp,i,ia)
!ik         if(hpit(j,i).le.0.d0) cycle
!ik         dl_connect = ( hdxm(j,i) + hdxp(j,i) ) / hpit(j,i)
!ik         if(dl_connect.le.0.d0) cycle

! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik         grad_ti_ev(j,i) =  ( ti_p_ev - ti_ev ) / dl_connect
!ik         grad_te_ev(j,i) =  ( te_p_ev - te_ev ) / dl_connect
!                                         [eV/m]
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik         grad_ni(j,i) =  ( ni_p - ni ) / dl_connect
!ik         grad_ui(j,i) =  ( ui_p - ui ) / dl_connect
!ik         grad_pe(j,i) =  ( pressure_e(jp,i) - pressure_e(j,i) ) /
!ik  >           dl_connect
!ik         grad_b_abs(j,i) = ( b_abs(jp,i) - b_abs(j,i) ) / dl_connect

!     Electric field (static) (cf. Note p.imp6 OR soniv doc. Eq.(2.2-11))
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik         field_elec_para(j,i) = -grad_pe(j,i)/(vnezef(j,i)*cev) -
!ik  >           0.71d0*( grad_te_ev(j,i)*cev )/cev

!     Characteristic length evaluation
! deleted 8 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(grad_ti_ev(j,i) ).gt.0.d0)  ch_length_ti_ev(j,i)
!ik  >           = ti_ev / grad_ti_ev(j,i)
!ik         if(abs(grad_ni(j,i) ).gt.0.d0)  ch_length_ni(j,i)
!ik  >           = ni / grad_ni(j,i)
!ik         if(abs(grad_ui(j,i) ).gt.0.d0)  ch_length_ui(j,i)
!ik  >           = ui / grad_ui(j,i)
!ik         if(abs(grad_b_abs(j,i) ).gt.0.d0)  ch_length_b_abs(j,i)
!ik  >           = b_abs(j,i) / grad_b_abs(j,i)

!     Knudsen numbers
! deleted 12 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(ch_length_ti_ev(j,i)).gt.0.d0) knudsen_ti_ev(j,i)
!ik  >           = mfp_ii(j,i) / ch_length_ti_ev(j,i)
!ik         if(abs(ch_length_ni(j,i)).gt.0.d0) knudsen_ni(j,i)
!ik  >           = mfp_ii(j,i) / ch_length_ni(j,i)
!ik         if(abs(ch_length_ui(j,i)).gt.0.d0) knudsen_ui(j,i)
!ik  >           = mfp_ii(j,i) / ch_length_ui(j,i)
!ik         if(abs(ch_length_b_abs(j,i)).gt.0.d0) knudsen_b_abs(j,i)
!ik  >           = mfp_ii(j,i) / ch_length_b_abs(j,i)
!ik         knudsen_connect(j,i) = mfp_ii(j,i) / l_connect(jte-1)
!ik         abs_knudsen_max(j,i) = max( abs(knudsen_ti_ev(j,i) ),
!ik  >           abs(knudsen_ni(j,i) ), abs(knudsen_ui(j,i) ),
!ik  >        abs(knudsen_b_abs(j,i) ), abs(knudsen_connect(j,i) ) )

!           Effective Kn_Ti
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(flmxi(j,i) ).gt.0.d0) knudsen_ti_ev_effective(j,i)
!ik  >           = ( (1.d0/flmxi(j,i)) -1.d0) * flimi / 3.9d0

!     qGM rough estimation
!         ion_MFP_effective(i-i & i-Z colls included, evaluated from SOLDOR info.)
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik         mfp_ion_effective(j,i) =
!ik  >       knudsen_ti_ev_effective(j,i) * abs(ch_length_ti_ev(j,i) )
!         Omega_GM
! deleted 7 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(ch_length_ui(j,i)).gt.0.d0) then
!ik            if(abs(ch_length_b_abs(j,i)).gt.0.d0) then
!ik               omega_gm(j,i) =
!ik  >              5.88d0 * mach_number(j,i) * mfp_ion_effective(j,i) *
!ik  >     (1.d0 / ch_length_ui(j,i)  - 0.25d0 / ch_length_b_abs(j,i) )
!ik            endif
!ik         endif
!         flmxi_qGM
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(1.d0+omega_gm(j,i) ).gt.0.d0) then
!ik            flmxi_qgm(j,i) =
!ik  >           (1.d0 - 1.2d0*mach_number(j,i)*mach_number(j,i) ) /
!ik  >           (1.d0 + omega_gm(j,i) )
!ik         endif
!         qGM (heat_flux_para_GM)
! deleted 5 lines organize local variables and include files by kamata 2021/05/31
!ik         if(abs(flmxi(j,i)).gt.0.d0) then
!ik            heat_flux_para_gm(j,i) =
!ik  >           (heat_flux_para_by_kappa_para(j,i)/flmxi(j,i) ) *
!ik  >           flmxi_qgm(j,i)
!ik         endif

! deleted 1 line organize local variables and include files by kamata 2021/11/18
!ik         if(limp==3)then
!     Force balance
! deleted 4 lines organize local variables and include files by kamata 2021/05/31
!ik         call thermal_force_para_by_q_para(clambda,
!ik  >           ami, amz, 1.d0*cev, charge_imp, ti_ev,
!ik  >           0.d0, sqrt(2.d0 * cev * ti_ev /amz), ui,
!ik  >           heat_flux_para_by_kappa_para(j,i), f_thi(j,i) )
!!!
!!!     JUST TEST 2018/11/29  analysis for Fthi W/O collisionality effect
!            if(abs(flmxi(j,i) ).gt.0.d0) then
!               f_thi(j,i) = f_thi(j,i) / flmxi(j,i)
!            else
!               f_thi(j,i) = 0.d0
!            endif
!!!
! deleted 4 lines organize local variables and include files by kamata 2021/05/31
!ik         call thermal_force_para_by_q_para(clambda,
!ik  >           cme, amz, (-cev), charge_imp, te_ev,
!ik  >           0.d0, sqrt(2.d0 * cev * ti_ev /amz), ui,
!ik  >           heat_flux_para_elec(j,i), f_the(j,i) )

!      subroutine thermal_force_para_by_q_para(coulomb_logarithm,
!     >     mi, mz, ei[C], ez[C], ti[eV], vz_para, vz_perp, ui, q_para, fth_para)
! deleted 4 lines organize local variables and include files by kamata 2021/05/31
!ik         call friction_force_theoretical_value(clambda,
!ik  >           ami, amz, 1.d0*cev, charge_imp, ti_ev,
!ik  >           ni,
!ik  >           0.d0, sqrt(2.d0 * cev * ti_ev /amz), ui, f_0(j,i) )
!      subroutine friction_force_theoretical_value(coulomb_logarithm,
!     >     mi, mz, ei[C], ez[C], ti[eV], ni[m^(-3)], vz_para, vz_perp, ui, friction_para)
! deleted 1 line organize local variables and include files by kamata 2021/11/18
!ik         endif !limp

!     Electric field (static) force
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik         f_elec(j,i) = charge_imp  * field_elec_para(j,i)

!     Resultant force
! deleted 3 lines organize local variables and include files by kamata 2021/05/31
!ik         f_tot_simple(j,i) = f_0(j,i) + f_thi(j,i)
!ik         f_tot_full(j,i) = f_0(j,i) + f_thi(j,i) + f_the(j,i) +
!ik  >           f_elec(j,i)

! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik      enddo                  !            do jt = jts, jte
!ik   enddo                     !      do it = 1, itmax
!
!
!===========================================
!     FILE WRITING
!===========================================
!     Flag check (refer to comments on top)
!   and file name definition
      if( (flg_radial>0).and.(flg_write_all==0) ) then
         write(file_name_1,'(a,i3.3,a)')
     >        'heat_flux_radial_jt_', jcel_input, '.dat'

      else if( (flg_radial==0).and.(flg_write_all>0) ) then
         write(file_name_1,'(a,i3.3,a)') 'heat_flux_all_cells.dat'

      else
       write(n6, *) "Subroutine 'heat_flux_radial_output' did not work."
       write(n6, *) "Flags must be eigher "
       write(n6, *) "(flg_radial, flg_write_all) = (1,0) or (0,1)."
         return
      endif


!     Radial or All cell profile of q_para and main physics parameters.
      open(3000,file=trim(file_name_1), status='replace')

!     write the common parameters necessary for reading at post-treatment.
      write(3000,*) "# common parameters for reading at post treatment:"
      write(3000,*) "# ndx = ", ndx
      write(3000,*) "# ndy = ", ndy

!      it_count: do it = 1, itmax
      it_count: do it = 2, itmax-1
         l_connect = 0.d0
         l_poloidal = 0.d0
         jts = jtmin(it)
         jte = jtmax(it)
!         ic = mcel(j,i)
!         if( it.ge.itmps .and. it.le.itmpe ) then
!            jts = jts + 1
!            jte = jte - 1
!         endif

         call slenb(it,l_connect)
         call plenc(it,l_poloidal)
!
         jt_count: do jt = jts, jte
!     jt_count: do jt = jts+1, jte-1 ! to remove extreme values at the ends
            j  = jcel(jt,it)
            i  = icel(jt,it)
!     irg = kreg(j,i)
!
!     Write decision, according to the flags:
            if( (flg_radial>0).and.(flg_write_all==0).and.
     >           (j.ne.jcel_input) ) then
               cycle

            else if( (flg_radial>0).and.(flg_write_all==0).and.
     >              (j==jcel_input) )  then
               continue

            else if( (flg_radial==0).and.(flg_write_all>0) ) then
               continue
            else
               return
            endif

            write(3000, '(4i6,50(E16.4e2))')
     >           jt, it, j, i,
     >           l_connect(jt), l_connect(jte-1) - l_connect(jt),
     >           l_poloidal(jt), l_poloidal(jte-1) - l_poloidal(jt),
     >           rhf_odp(i), rhf_idp(i), rhf_omd(i), rhf_imd(i),
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)),
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           hpit(j,i), vnezef(j,i), vna(j,i,ia),
     >           vte(j,i), vti(j,i), vva(j,i,ia),
     >           sn0(j,i,ia),   ! neutral density [m-3]
     >           nz_ion_after_wfac(j,i), nz_tot_after_wfac(j,i),
     >           nz_ion_after_wfac(j,i)/vnezef(j,i), nz_on_ne(j,i),
     >           energy_dens_elec_conv_para(j,i), ! q//_e_convective[W/m2]
     >           heat_flux_para_elec(j,i), ! q//_e_conductive[W/m2]
     >           energy_dens_convective_para(j,i), ! q//_i_convective[W/m2]
     >           heat_flux_para_by_kappa_para(j,i), ! q//_i_conduction [W/m2]
     >           energy_dens_elec_para_tot(j,i), ! q//_e_cond+conv [W/m2]
     >           energy_dens_total_para(j,i), ! q//_i_cond+conv [W/m2]
     >           energy_dens_total_para(j,i) +
     >           energy_dens_elec_para_tot(j,i), ! q//i tot + q//e tot [W/m2]
     >           particle_flux_dens_para_i(j,i), ! Gamma_para_ion [s-1m-2]
     >           pressure_e(j,i), pressure_i(j,i), ! Static pressure e & i [Pa]
     >           pressure_e(j,i)+pressure_i(j,i),
     >           cs_ion(j,i), mach_number(j,i),
     >           pressure_dyn_e(j,i)+pressure_dyn_i(j,i), ! dynamic pressure i+e [Pa]
     >           pressure_dyn_e(j,i)+pressure_dyn_i(j,i)
     >           +pressure_e(j,i)+pressure_i(j,i), ! Total plasma pressure [Pa]
     >           flmxe(j,i), flmxi(j,i), ! heat flux limiter elec. & ion
     >           tot_rad_rate(j,i) ! impurity radiation density[W/m3]

!$$$  >              nz_tot_after_wfac(j,i), nz_on_ne(j,i),
!$$$  >              tot_rad_rate(j,i),
!$$$  !     >              f_0(j,i), f_thi(j,i), f_tot_full(j,i),
!$$$  >              f_0(j,i), f_thi(j,i), f_tot_simple(j,i),
!$$$  >              vnezef(j,i), vna(j,i,ia),
!$$$  >              vte(j,i), vti(j,i), pressure_e(j,i), pressure_i(j,i),
!$$$     >              vva(j,i,ia),
!$$$  >              heat_flux_para_by_kappa_para(j,i), flmxi(j,i),
!$$$  >              mach_number(j,i),
!$$$  >              energy_dens_convective_para(j,i),
!$$$  >              energy_dens_total_para(j,i),
!$$$  >              knudsen_ti_ev_effective(j,i), flmxi_qgm(j,i),
!$$$  >              heat_flux_para_gm(j,i)

         enddo jt_count
      enddo it_count

      close(3000)

!
      return
      end subroutine
