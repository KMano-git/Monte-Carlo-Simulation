!***********************************************************************
      subroutine impldif(ip,dtunt,taus)
!***********************************************************************
!
!   Langevan equation
!
!     dv/dt = -(v-Vf)/taus + Fm + 1/dlt*{<dv^2>*dlt}^1/2*rG
!     ds/dt = v
!
!     v(n+1)-v(n) = -(v(n)-Vf)*dlt/taus + Fm*dlt + {<dv^2>*dlt}^1/2*rG
!     s(n+1)-s(n) = v(n)*dlt
!
!     Vst = Vf + Fm*taus
!
!     M[v(n)] = v0*exp(-t/taus) + Vst*(1-exp(-t/taus))
!     D[v(n)] = 1/2*{1-exp(-2*t/taus)}*taus*<dv^2>
!
!     M[s(n)] = Vst*t + (v0-Vst)*taus*{1-exp(-t/taus)}
!     D[s(n)] = 1/2*{ 2*t/taus - 3 + 4*exp(-t/taus) - exp(-2*t/taus) }
!                     * taus^3*<dv^2>
!
!       Note.  sgm = <dv^2>   m^2/s^2/s
!
!                           2004/9/21    K. Shimizu
!-----------------------------------------------------------------------
      use cimcns,    only : fgau, gaus0, gaus1 
      use cimcom,    only : aimas, ami, amz, azmas, cfez, cfgte, cfgti
     >    , dfm, flmxi_at_impmc, fxm2, gdez, gdte, gdti, glte, glti
     >    , heat_flux_para_at_ic_center, heat_flux_para_e_ic_center, ir
     >    , is, lfgti, lfgtzf, mdl_felec, mdl_fth_q_base
     >    , mdl_fthe_qe_base, mdl_fthi_wo_limiter, rr, slnv, slw1, v, vr
     >    , vv, vvr, vvz, vz, zz
      use cntpls,    only : dene, teme, temi, vflw, zefm
      use com_eqdat, only : dr, dz, nr, rg, ubx, uby, zg 
      use cphcns,    only : cev, cme
      implicit none
!
!::argument
      integer, intent(in) :: ip
      real(8), intent(in) :: dtunt, taus
!
!::local variables
      integer  i, j
      integer  ic, iz
      real*8   r0, z0, ftot, r1, z1
      real*8   tx, ty, tc0, tc1, tc2, txy, tbx, tby
      integer  itx, ity, itx1, ity1
!
!::local variables for kinetic thermal force
      real*8  zti, zvf, zvi2, xx2, xxp2, zbet, ckp0
      integer ia
      real(8) :: zcti1, grdti, grdte
      real(8) :: za, zfca, zcti, zcte

      !Y.Homma
      integer  is_deuterium
      real*8   zlne, ei, ez
      real*8   vz_para, vz_perp, q_para
      real*8   hzeff, hAi, hZi, zne, zte
      real*8   fth_i_para, fth_i_para_divided_cev
      real*8   f_elec_divided_cev, fth_e_divided_cev, fth_e_para
!
!::local varibles for new model
      integer irn
      real*8  t, ex1, ex2, spmn, spdv, vpmn, vpdv, rng, sp, vp
      real*8  zvx, zvy, zvz
      real*8  Vst, v0, sgm, vf, Fm
!
      integer ms, iif, ig0, ig1
      real*8  dt, dts, slx1, uz, u, fx2, sl0, zvt, zvr, dlt
! function
      real(8)    coulog, random
      integer iino
!
      iino(i,j) = nr*(j-1) + i
!
      ic  = ir(ip)
      iz  = is(ip)
      za  = iz
      if( iz.le.0 ) return
!
      r0  = rr(ip)
      z0  = zz(ip)
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

!     ftot  = cfez(iz)*gdez(ic)+zcte*grdte+   zcti*grdti
!     ftot  = cfez(iz)*gdez(ic)+zcte*grdte+ fth_i_para_divided_cev
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
!
            zfca  = 1.0d0
            if( lfgtzf.eq.1 ) zfca  = 1.0d0/za**2*(za/zefm(ic)-1.0d0)*za
            zcti  = cfgti(iz)*zfca
            zcte  = cfgte(iz)*zfca
            if( lfgti.eq.1 ) zcti  = zcti1*zfca
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
!---------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!::sgm = {<dv^2>}   see imscat.f
!-----------------------------------------------------------------------
      ms  = 2
      dlt = taus/20.0d0
      dt  = dlt*is(ip)**2
      dts = sqrt(dt)
      slx1 = dts*slw1(ic,ms)
!-----
      zvt  = sqrt(temi(ic)*cev/amz)
      ig0  = int(2000.0*random(0)+1.0)
      zvx  = zvt*gaus0(ig0)
      ig0  = int(2000.0*random(0)+1.0)
      zvy  = zvt*gaus0(ig0)
      ig0  = int(2000.0*random(0)+1.0)
      zvz  = zvt*gaus0(ig0)
!
      v0  = vz(ip)
      vf  = vflw(ic,1)
      Fm  = ftot*cev/amz
      Vst = vf + Fm*taus
      zvz = zvz + Vst
!
      uz  = zvz - vf
      u   = sqrt(zvx**2+zvy**2+uz**2)
      u   = sqrt(3.0d0*temi(ic)*cev/amz)
!-----
      iif = int(dfm*slnv(ic,ms)*u+1.0)
      if( iif.ge.1000 ) iif = 1000
      fx2 = fxm2(iif,ms)
      sl0 = slx1*fx2        ! {<dv^2>*dlt}^1/2
      sgm = sl0**2          ! {<dv^2>*dlt}*(uz/u)**2
      sgm = sgm/dlt         ! {<dv^2>}
      sgm = sgm*2.0d0       !  factor
!
!-----------------------------------------------------------------------
!::(sp,vp) based on parallel diffusion model
!-----------------------------------------------------------------------
      t   = dtunt/taus
      ex1 = exp(-t)
      ex2 = exp(-2.0d0*t)
!
      spmn = Vst*taus*t+(v0-Vst)*taus*(1.0d0-ex1)
      spdv = 0.5d0*(2.0d0*t-3.0d0+4.0d0*ex1-ex2)*taus**3*sgm
!
      vpmn = v0*ex1+Vst*(1.0d0-ex1)
      vpdv = 0.5d0*(1.0d0-ex2)*taus*sgm
!
      irn = int(fgau*random(0) + 1.0d0)
      rng = gaus0(irn)
      sp  = spmn + sqrt(spdv)*rng
      irn = int(fgau*random(0) + 1.0d0)
      rng = gaus0(irn)
      vp  = vpmn + sqrt(vpdv)*rng
!-----------------------------------------------------------------------
!
!::predict
      tx   = (r0-rg(1))/dr + 1.0d0
      ty   = (z0-zg(1))/dz + 1.0d0
      itx  = int(tx)
      ity  = int(ty)
      itx1 = itx + 1
      ity1 = ity + 1
      tx   = tx - itx
      ty   = ty - ity
      txy  = tx * ty
      tc0  = 1.0d0 - tx - ty + txy
      tc1  = tx - txy
      tc2  = ty - txy
      tbx  = tc0*ubx(iino(itx,ity))  + tc1*ubx(iino(itx1,ity))
     >     + tc2*ubx(iino(itx,ity1)) + txy*ubx(iino(itx1,ity1))
      tby  = tc0*uby(iino(itx,ity))  + tc1*uby(iino(itx1,ity))
     >     + tc2*uby(iino(itx,ity1)) + txy*uby(iino(itx1,ity1))
!xx   tbz  = tc0*ubz(iino(itx,ity))  + tc1*ubz(iino(itx1,ity))
!xx  >     + tc2*ubz(iino(itx,ity1)) + txy*ubz(iino(itx1,ity1))
!
      r1  = r0  + 0.5d0*sp*tbx
      z1  = z0  + 0.5d0*sp*tby
!
!::correct
      tx   = (r1-rg(1))/dr + 1.0d0
      ty   = (z1-zg(1))/dz + 1.0d0
      itx  = int(tx)
      ity  = int(ty)
      itx1 = itx + 1
      ity1 = ity + 1
      tx   = tx - itx
      ty   = ty - ity
      txy  = tx * ty
      tc0  = 1.0d0 - tx - ty + txy
      tc1  = tx - txy
      tc2  = ty - txy
      tbx  = tc0*ubx(iino(itx,ity))  + tc1*ubx(iino(itx1,ity))
     >     + tc2*ubx(iino(itx,ity1)) + txy*ubx(iino(itx1,ity1))
      tby  = tc0*uby(iino(itx,ity))  + tc1*uby(iino(itx1,ity))
     >     + tc2*uby(iino(itx,ity1)) + txy*uby(iino(itx1,ity1))
!xx   tbz  = tc0*ubz(iino(itx,ity))  + tc1*ubz(iino(itx1,ity))
!xx  >     + tc2*ubz(iino(itx,ity1)) + txy*ubz(iino(itx1,ity1))
!
      rr(ip) = r0 + sp*tbx
      zz(ip) = z0 + sp*tby
      vz(ip) = vp
!
      zvt = sqrt(2.0d0*temi(ic)*cev/amz)
      ig1 = int(1000.0*random(0)+1.0)
      zvr = zvt*gaus1(ig1)
      vr(ip) = zvr
!-----
      vvz(ip) = vz(ip)*vz(ip)
      vvr(ip) = vr(ip)*vr(ip)
      vv(ip)  = vvr(ip) + vvz(ip)
      v(ip)   = sqrt(vv(ip))
!-----
!
      return
      end
