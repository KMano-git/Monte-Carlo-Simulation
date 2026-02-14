!***********************************************************************
      subroutine imorbit(ip,dtunt)
!***********************************************************************
!
!   kinetic thermal force
!          D. Reiser at al., Nucl. Fusion 38 (1998) 165
!
!   comt   ip  ic   ix iy is   tt          Ez        Ti
!   gdti      glti      grdti     gdte      glte      grdte     ftot
!   ftot1     x         bet0      bet       vz0      dvz
!
!-----------------------------------------------------------------------
      use cimcom,    only : aimas, ami, amz, azmas, cfez, cfgte, cfgti
     >    , flmxi_at_impmc, gdez, gdte, gdti, glte, glti
     >    , heat_flux_para_at_ic_center, heat_flux_para_e_ic_center, ir
     >    , is, lfgti, lfgtzf, mdl_felec, mdl_fth_q_base
     >    , mdl_fthe_qe_base, mdl_fthi_wo_limiter, rr, vr, v, vv, vvr
     >    , vvz, vz, zz
      use cimp_loc,  only : dlbp
      use cntpls,    only : dene, teme, temi, vflw, zefm
      use com_eqdat, only : dr, dz, nr, rg, ubx, uby, zg
      use cphcns,    only : cev, cme
      use cunit,     only : lmspe, lmype
      implicit none
!
!::argument
      integer, intent(in) :: ip
      real(8), intent(in) :: dtunt
!
!::local variables
      integer  i, j
      real*8   wdtv
      integer  ic, iz
      real*8   r0, z0, vz0, r1, z1, vz1
      real*8   tx, ty, tc0, tc1, tc2, txy, tbx, tby
      integer  itx, ity, itx1, ity1
      real*8   dsp, zro
      real*8   grdti, grdte
!
!::local variables (ftot)
      real(8) :: ftot
      real(8) :: za, zfca, zcti, zcte
      real(8) :: zti, zvf, zvi2, xx2, xxp2, zbet, ckp0
      real(8) :: zcti1
      integer   is_deuterium
      real*8   zlne, ei, ez
      real*8   vz_para, vz_perp, q_para
      real*8   hzeff, hAi, hZi, zne, zte
      integer ia
      real*8   fth_i_para, fth_i_para_divided_cev
      real*8   f_elec_divided_cev, fth_e_divided_cev, fth_e_para
!
! function
      real(8)   coulog, funroh
      integer   imox, imoy
      integer iino
!
      iino(i,j) = nr*(j-1) + i
!
      ic  = ir(ip)
      iz  = is(ip)
      if( iz.le.0 ) return
!
      wdtv = dtunt*cev/amz
!
      r0  = rr(ip)
      z0  = zz(ip)
      vz0 = vz(ip)
!
!---------------------------------------------------------------------
!::Fz = electric field + e-thermal + i-thermal force
      ftot = 0.0d0
      fth_i_para_divided_cev = 0.d0
      fth_e_divided_cev = 0.d0
      f_elec_divided_cev = 0.d0

      if( temi(ic).gt.0.0d0 ) then
!
         if( mdl_fth_q_base.eq.1 ) then  ! Y.Homma 2018/06/01
!::Thermal force evaluated based on bulk-ion heat flux (Homma JCP13)
!  If mdl_fthe_qe_base = 1, kinetic thermal force contribution from electron is ON.
!  If mdl_felec = 1,  static-electric field force is ON.
            ia  = 1  ! Single bulk ion D+ (ia = 1)supposed
            is_deuterium = 1
!            hAi = aion(ia)  ! *** to improve?
!             hZi = aza(ia)  ! *** to improve?
            hAi = 2.d0 ! Single bulk ion D+ (ia = 1)supposed
            hZi = dfloat(is_deuterium)

            zvf = vflw(ic,ia)
            za  = iz
            zne = dene(ic)
            zte = teme(ic)
            zti = temi(ic)
            hzeff = zefm(ic)
!     hzeff = wzf(j) = vzf(j,i) = zefm(ic)@monte/mcplas.f
            zlne = coulog( hzeff, zne, zte, zti, hAi, hZi )

            ei = dfloat(is_deuterium)*cev ; ez = za*cev
            vz_para = vz(ip) ; vz_perp = vr(ip)
            q_para = heat_flux_para_at_ic_center(ic)
!                                                [joule/(m2.s)]
            call thermal_force_para_by_q_para(zlne,
     >           ami, amz, ei, ez, zti, vz_para, vz_perp, zvf, q_para,
     >           fth_i_para)
            fth_i_para_divided_cev = fth_i_para/cev
!
            if(mdl_fthi_wo_limiter.eq.1) then
               if( abs(flmxi_at_impmc(ic) ).gt.0.d0) then
                  fth_i_para_divided_cev =
     >                 (fth_i_para/cev) / flmxi_at_impmc(ic)
!               write(555,*) "mdl_fthi_wo_limiter ON: arrived@imorbit.f"
               endif
            endif

            if(mdl_fthe_qe_base.eq.1) then
!              Assumption: zlne:= Coulomb lambda is supposed
!                   to be common with Fth ion.
               call thermal_force_para_by_q_para(zlne,
     >              cme, amz, (-cev), ez, zte, vz_para, vz_perp, zvf,
     >              heat_flux_para_e_ic_center(ic),
     >              fth_e_para)

               fth_e_divided_cev = fth_e_para / cev
            endif

            if(mdl_felec.eq.1) then
               f_elec_divided_cev = cfez(iz)*gdez(ic)
            endif
            ftot  = f_elec_divided_cev + fth_e_divided_cev +
     >           fth_i_para_divided_cev
         else
            ckp0 = 1.0d0/0.5657d0
            ia  = 1   ! Single bulk ion (ia = 1)supposed
            zvf = vflw(ic,ia)
            za  = iz
            zti = temi(ic)
            zvi2 = 2.0d0*cev*zti/ami
            xx2  = (vvr(ip)+(vz(ip)-zvf)**2)/zvi2
            xxp2 = (vz(ip)-zvf)**2/zvi2
            zbet = 1.5d0*ckp0*(1.0d0-2.0d0*xxp2)*exp(-xx2)
            zcti1 = (1.0d0+aimas/azmas)*za**2*zbet
            zfca  = 1.0d0
            if( lfgtzf.eq.1 ) zfca  = 1.0d0/za**2*(za/zefm(ic)-1.0d0)*za
            zcti  = cfgti(iz)*zfca
            zcte  = cfgte(iz)*zfca
            if( lfgti.eq.1 ) zcti  = zcti1*zfca
!
!::   new limit of temperature gradient (up to 2018/05/31)
            grdti = dmin1( glti(ic), dabs(gdti(ic)) )
            grdti = dsign( grdti, gdti(ic) )
            grdte = dmin1( glte(ic), dabs(gdte(ic)) )
            grdte = dsign( grdte, gdte(ic) )

            fth_i_para_divided_cev = zcti*grdti
            fth_e_divided_cev = zcte*grdte
            f_elec_divided_cev = cfez(iz)*gdez(ic)

!     ftot  = cfez(iz)*gdez(ic)+zcte*grdte+zcti*grdti
!     ftot  = cfez(iz)*gdez(ic)+zcte*grdte+ fth_i_para_divided_cev

            ftot  = f_elec_divided_cev + fth_e_divided_cev +
     >           fth_i_para_divided_cev
         endif                  ! if( mdl_fth_q_base.eq.1 ) then

      endif

!
!::predict
      tx   = (r0-rg(1))/dr + 1.0d0
      ty   = (z0-zg(1))/dz + 1.0d0
      itx  = int(tx)
      ity  = int(ty)
      itx1 = itx + 1
      ity1 = ity + 1
      tx   = tx - dfloat(itx)
      ty   = ty - dfloat(ity)
      txy  = tx * ty
      tc0  = 1.0d0 - tx - ty + txy
      tc1  = tx - txy
      tc2  = ty - txy
      tbx  = tc0*ubx(iino(itx,ity))  + tc1*ubx(iino(itx1,ity))
     >     + tc2*ubx(iino(itx,ity1)) + txy*ubx(iino(itx1,ity1))
      tby  = tc0*uby(iino(itx,ity))  + tc1*uby(iino(itx1,ity))
     >     + tc2*uby(iino(itx,ity1)) + txy*uby(iino(itx1,ity1))
      r1  = r0  + 0.5d0*dtunt*vz0*tbx
      z1  = z0  + 0.5d0*dtunt*vz0*tby
      vz1 = vz0 + 0.5d0*wdtv*ftot
!
!::correct
      tx   = (r1-rg(1))/dr + 1.0d0
      ty   = (z1-zg(1))/dz + 1.0d0
      itx  = int(tx)
      ity  = int(ty)
      itx1 = itx + 1
      ity1 = ity + 1
      tx   = tx - dfloat(itx)
      ty   = ty - dfloat(ity)
      txy  = tx * ty
      tc0  = 1.0d0 - tx - ty + txy
      tc1  = tx - txy
      tc2  = ty - txy
      tbx  = tc0*ubx(iino(itx,ity))  + tc1*ubx(iino(itx1,ity))
     >     + tc2*ubx(iino(itx,ity1)) + txy*ubx(iino(itx1,ity1))
      tby  = tc0*uby(iino(itx,ity))  + tc1*uby(iino(itx1,ity))
     >     + tc2*uby(iino(itx,ity1)) + txy*uby(iino(itx1,ity1))
      rr(ip)  = r0  + dtunt*vz1*tbx
      zz(ip)  = z0  + dtunt*vz1*tby
      vz(ip)  = vz0 + wdtv*ftot
!-----
      vvz(ip) = vz(ip)*vz(ip)
      vv(ip)  = vvr(ip) + vvz(ip)
      vv(ip)  = dmax1( vv(ip), vvz(ip) )
      vvr(ip) = vv(ip) - vvz(ip)
      vr(ip)  = sqrt(vvr(ip))
      v(ip)   = sqrt(vv(ip))
!-----
!
!::limit  ! <== 05.06 11.18  08/8.17
      if( lmype.eq.lmspe ) then
      zro = funroh(ic,r0,z0)
      if( zro.gt.0.90d0 ) then
      dsp = sqrt((rr(ip)-r0)**2+(zz(ip)-z0)**2)
      if( dlbp(ic).gt.2.0d-2 .and. dsp.gt.1.5d0*dlbp(ic) ) then
      write(975,'(2x,"Warning imorbit",2x,i6,i8,2i5,1pe12.3,
     >  0p2f9.5, 1pe12.3, 0p2f9.5, 1p7e12.3)')
     >  ip, ic, imox(ic), imoy(ic), dtunt, rr(ip), zz(ip), vz(ip),
     >  r0, z0, vz0, zro, wdtv*ftot, tbx, tby, dsp, dlbp(ic)
      endif
      endif
      endif
!
      return
      end
