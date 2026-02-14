!***********************************************************************
      subroutine pxvisc(it)
!***********************************************************************
!
!  -------------------------
!   flux = -Diff*X*grad(U)
!  -------------------------
!
!                   Diff * X            *  grad(U)
!     -----------------------------------------------------------
!      Admt, dacf   Da                  *  grad(roa)   qfx_df(j,i,m1a)
!                   Da * va             *  grad(roa)   qfx_df(j,i,m2a)
!                   Da * (ea+pa)/roa    *  grad(roa)   qfx_df(j,i,m3)
!                   Da * 5/2*Te*za/ma   *  grad(roa)   qfx_df(j,i,m4)
!      Atmt, etcf   eta                 *  grad(va)    qfx_cd(j,i,m2a)
!      Aimt, xicf    xi                 *  grad(Ti)    qfx_cd(j,i,m3)
!      Aemt, xecf    xe                 *  grad(Te)    qfx_cd(j,i,m4)
!
!                   eta             *  1/2*grad(va^2)  qfx_vh(j,i)
!
!
!    Ad_(E,W,N,S) = hdsp*D^para + hdsv*D^ano
!
!    dacfp(j,i) = AdE + 1/2*(AdN(+)-AdS(-))
!    dacfm(j,i) = AdW + 1/2*(AdS(+)-AdN(-))
!
!         Dano : harmonic mean
!         X    : linear interpolation
!
!     plasma parameter at grid points                       00/10/10
!     bounday condition at divertor plates   D^ano = 0.0    00/12/12
!     new X interpolation     zwtm*q(j)+zwtp*q(j+1)         00/12/14
!     delta-X  (va)                                         00/12/19
!     add D^para in prv-region                              01/09/02
!     qcon                                                  02/10/25
!     Dif|j+1/2   interpolation method                      03/07/26
!
!-----------------------------------------------------------------------
!
!     G1a = - Da*grad(roa)              [zf11]  qfx_df(m1a)
!     G2a = - va*Da*grad(roa)               [zf21]  qfx_df(m2a)
!           - eta*grad(va)                  [zf22]  qfx_cd(m2a)
!     G3  = - sum[(ea+pa)/roa*Da*grad(roa)]    [zf31]  qfx_df(m3)
!           - sum[eta*grad(1/2*va**2)]         [zf32]  qfx_vh
!           - Ki*grad(Ti)                      [zf33]  qfx_cd(m3)
!     G4  = - sum[(ee+pe)/ne*za/ma*Da*grad(roa)]  [zf41]  qfx_df(m4)
!           - Ke*grad(Te)                         [zf43]  qfx_cd(m4)
!
!              note.  (ee+pe)/ne = 5/2*Te   roa = Q1a
!
!     flux limiter for ion parallel heat diffusivity   2010/11/30
!
!-----------------------------------------------------------------------
      use cphcns, only : cev, cme, cmp
      use cplcom, only : aion, ama, aza, c13, c23, c52, cc, daan, dacl
     >    , dacfe, dacfw, dd, ee, ehcfe, ehcfw, etan, etcfe, etcfw, etcl
     >    , ff, flime, flimi, flimv, flmet, flmxe, flmxi, hea
     >    , heat_flux_para_by_kappa_para, heat_flux_para_elec, hia
     >    , nion, q1a, q2a, r1ae, r1ai, r1av, r1aw, r2ai, r2av, r2aw
     >    , r3i,  r4e, vea, vna, vnag, vne, vni, vte, vteg, vti, vtig
     >    , vva, vvag, xean, xecfe, xecfw, xecl, xian, xicfe, xicfw
     >    , xicl, xste, xsti, xsva
      use cplmet, only : hdsp, hdsv, hvsb, hwtm, hwtp, icel, jcel, jtmax
     >    , kce, kcn, kcs, kcw, set_hdsv
      use cplqcn, only : qfx_cd, qfx_df, qfx_vh
      use csize,  only : ndeq, ndx
      use cunit,  only : n6
      implicit none
!
!::argument
      integer, intent(in) :: it
!
!::local variables
      integer  k1, k2, nequ
      integer  i, ia, jmax, jmax1, jw, j, jp, jm, jw1
      integer  m1a, m2a, m3, m4
      real*8   zwtm, zwtp, fcan, fccl
      real*8   zran, zras, zvan, zvas, ztin, ztis, zten, ztes
      real*8   zdacl, zdaan, zetcl, zetan, zxicl, zxian, zxecl, zxean
      real*8   zqcl, znao, ztio, zqfl, zfet
      real*8   zhai, zhae, zhaip, zhaep
      real*8   zneo, zteo, zfxe
      real*8   zf11, zf21
      real*8   zf22, zf31
      real*8   zf32, zf33
      real*8   zf41, zf43
      real*8   Ade, Adw, Adn, Ads, Adnp, Adnm, Adsp, Adsm
      real*8   Ate, Atw, Atn, Ats, Atnp, Atnm, Atsp, Atsm
      real*8   Ahe, Ahw, Ahn, Ahs
      real*8   Aie, Aiw, Ain, Ais, Ainp, Ainm, Aisp, Aism
      real*8   Aee, Aew, Aen, Aes, Aenp, Aenm, Aesp, Aesm
      real*8   fvsh
      integer m, n
!::to send heat flux to IMPMC
      real*8 Aie_limited, Aiw_limited
      real*8 heat_flux_density_para_temp_1,heat_flux_density_para_temp_2
      real*8 Aee_limited, Aew_limited
      real*8 heat_flux_para_elec_temp_1,heat_flux_para_elec_temp_2
!
!::separation cc,dd,ee  !  %%% 2002/11/8
      integer jwb
      real*8 zdp(ndeq,ndeq),zdm(ndeq,ndeq)
      real*8 zcc(ndeq,ndeq),zee(ndeq,ndeq)
      real*8 wcc(ndeq,ndeq),wdp(ndeq,ndeq),wdm(ndeq,ndeq),wee(ndeq,ndeq)
      real*8 fx_df(ndx,ndeq),fx_cd(ndx,ndeq),fx_vh(ndx)
!
!::option of flux limiter
      real*8   znio, zmsi, zfxi
!
!::interpolation method for diffusion at j+1/2
      real*8  fundav, fdfc, fdwm, fdwp, fddm, fddp
! function
      fundav(fdfc,fdwm,fdwp,fddm,fddp) = fdfc/(fdwm/fddm+fdwp/fddp)
!
      i  = icel(1,it)
      jmax  = jtmax(it)
      jmax1 = jmax - 1
      m3 = 2*nion + 1
      m4 = 2*nion + 2
!
!-----------------------------------------------------------------------
!::flux clear
!-----------------------------------------------------------------------
      nequ = 2*nion + 2
      do k2=1,nequ
      do k1=1,ndx
      fx_df(k1,k2) = 0.0d0
      fx_cd(k1,k2) = 0.0d0
      enddo
      enddo
      do k1=1,ndx
      fx_vh(k1) = 0.0d0
      enddo
!
!-----------------------------------------------------------------------
!::plasma parameter & trans-function at cell center (j,i)
!-----------------------------------------------------------------------
      do 110 jw = 1, jmax
      j = jcel(jw,it)
      do 115 ia = 1, nion
      hia(j,ia)  = (vea(j,i,ia)+vna(j,i,ia)*vti(j,i)*cev)/q1a(j,i,ia)
      hea(j,ia)  =  c52*aza(ia)/ama(ia)*vte(j,i)*cev
!
      r1av(j,ia) = -q2a(j,i,ia)/q1a(j,i,ia)**2
      r2av(j,ia) = 1.0d0/q1a(j,i,ia)
      r1aw(j,ia) = vva(j,i,ia)*r1av(j,ia)
      r2aw(j,ia) = vva(j,i,ia)*r2av(j,ia)
      r1ai(j,ia) = (c13*vva(j,i,ia)**2-vti(j,i)*cev/ama(ia))/vni(j,i)
      r2ai(j,ia) = -c23*vva(j,i,ia)/vni(j,i)
      r1ae(j,ia) = -aza(ia)/ama(ia)*vte(j,i)*cev/vne(j,i)
 115  continue
      r3i(j) = c23/vni(j,i)
      r4e(j) = c23/vne(j,i)
 110  continue
!
!-----------------------------------------------------------------------
!::diff. coef. & flux at cell boundary (j+1/2,i)
!-----------------------------------------------------------------------
      do 210 jw = 1, jmax1
      j   = jcel(jw,it)
      jp  = jcel(jw+1,it)
!
!::D^ano = 0.0 at d-plates
      fcan = 1.0d0              ! diff.ne.0
      fccl = 1.0d0
!
      zf31 = 0.0d0
      zf32 = 0.0d0
      zf41 = 0.0d0
!
      do 220 ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::parameter (roa,va) at N(j+1/2,i+1/2), S(j+1/2,i-1/2)
!x    zran= c14*(q1a(j,i,ia)+q1a(j,i+1,ia)+q1a(jp,i,ia)+q1a(jp,i+1,ia))
!x    zras= c14*(q1a(j,i,ia)+q1a(j,i-1,ia)+q1a(jp,i,ia)+q1a(jp,i-1,ia))
!x    zvan= c14*(vva(j,i,ia)+vva(j,i+1,ia)+vva(jp,i,ia)+vva(jp,i+1,ia))
!x    zvas= c14*(vva(j,i,ia)+vva(j,i-1,ia)+vva(jp,i,ia)+vva(jp,i-1,ia))
      zran = ama(ia)*vnag(j,i,ia)
      zras = ama(ia)*vnag(j,i-1,ia)
      zvan = vvag(j,i,ia)
      zvas = vvag(j,i-1,ia)
!
!::d-coef (Da) at cell boundary (j+1/2,i)
!x    zdaan =fcan/(hwtm(j,i)/daan(j,ia)+hwtp(j,i)/daan(jp,ia))
      zdaan =fundav(fcan,hwtm(j,i),hwtp(j,i),daan(j,ia),daan(jp,ia))
      zdacl = dacl(j,ia)   ! add in prv-region
!-----
      if(set_hdsv .ne. 0) then
          Ade = hdsp(j,i,kce)*zdacl + hdsv(j,i,kce)*zdaan
          Adw = hdsp(j,i,kcw)*zdacl + hdsv(j,i,kcw)*zdaan
          Adn = hdsp(j,i,kcn)*zdacl + hdsv(j,i,kcn)*zdaan
          Ads = hdsp(j,i,kcs)*zdacl + hdsv(j,i,kcs)*zdaan
      else ! as default
          Ade = hdsp(j,i,kce)*zdacl
          Adw = hdsp(j,i,kcw)*zdacl
          Adn = hdsp(j,i,kcn)*zdacl
          Ads = hdsp(j,i,kcs)*zdacl
      endif
!
      Adnp = dmax1( Adn, 0.0d0 )
      Adnm = Adn - Adnp
      Adsp = dmax1( Ads, 0.0d0 )
      Adsm = Ads - Adsp
      dacfE(jw,ia) = Ade + 0.5d0*(Adnp-Adsm)
      dacfW(jw,ia) = Adw + 0.5d0*(Adsp-Adnm)
!
!---------------------------------------------------------------
!::up-wind
!---------------------------------------------------------------
!
!::velocity across cell boundary
      zf11 = - Ade*q1a(jp,i,ia) + Adw*q1a(j,i,ia) - Adn*zran + Ads*zras
!
!::energy
      zhai  = (vea(j,i,ia)+vna(j,i,ia)*vti(j,i)*cev)/q1a(j,i,ia)
      zhae  =  c52*aza(ia)/ama(ia)*vte(j,i)*cev
      zhaip = (vea(jp,i,ia)+vna(jp,i,ia)*vti(jp,i)*cev)/q1a(jp,i,ia)
      zhaep =  c52*aza(ia)/ama(ia)*vte(jp,i)*cev
!
!::weit
      zwtm = hwtm(j,i)
      zwtp = hwtp(j,i)
!
!::up-wind
      zwtm = 1.0d0
      zwtp = 0.0d0
      if( zf11.lt.0.0d0 ) then
      zwtm = 0.0d0
      zwtp = 1.0d0
      endif
!
!::X-interpolation (X) at O(j+1/2,i)
      xsva(jw,ia) = zwtm*vva(j,i,ia)+zwtp*vva(jp,i,ia)
      xsti(jw,ia) = zwtm*zhai       +zwtp*zhaip
      xste(jw,ia) = zwtm*zhae       +zwtp*zhaep
!
!
!::d-coef Eta  at cell boundary (j-1/2,i)
!x    zetcl =1.0d0/(hwtm(j,i)/etcl(j,ia)+hwtp(j,i)/etcl(jp,ia))
!x    zetan = fcan/(hwtm(j,i)/etan(j,ia)+hwtp(j,i)/etan(jp,ia))
      zetcl =fundav(fccl,hwtm(j,i),hwtp(j,i),etcl(j,ia),etcl(jp,ia))
      zetan =fundav(fcan,hwtm(j,i),hwtp(j,i),etan(j,ia),etan(jp,ia))
!
!---------------------------------------------------------------
!::flux limit of viscosity
!---------------------------------------------------------------
      Ate = hdsp(j,i,kce)*zetcl
      Atw = hdsp(j,i,kcw)*zetcl
      zqcl = - Ate*vva(jp,i,ia) + Atw*vva(j,i,ia)
!
!--new intepolation of flux at cell boundary
      znao = hwtm(j,i)*vna(j,i,ia)+hwtp(j,i)*vna(jp,i,ia)
      ztio = hwtm(j,i)*vti(j,i)+hwtp(j,i)*vti(jp,i)
      zqfl = flimv*znao*ztio*cev*hvsb(j,i)
!
      zfet = 1.0d0
      if( flimv.gt.0.0d0 ) zfet = 1.0d0/(1.0d0+dabs(zqcl/zqfl))
      zetcl = zetcl*zfet
      flmet(j,i) = zfet
!
      if(set_hdsv .ne. 0) then
          Ate = hdsp(j,i,kce)*zetcl+hdsv(j,i,kce)*zetan
          Atw = hdsp(j,i,kcw)*zetcl+hdsv(j,i,kcw)*zetan
          Atn = hdsp(j,i,kcn)*zetcl+hdsv(j,i,kcn)*zetan
          Ats = hdsp(j,i,kcs)*zetcl+hdsv(j,i,kcs)*zetan
      else ! as default
          Ate = hdsp(j,i,kce)*zetcl
          Atw = hdsp(j,i,kcw)*zetcl
          Atn = hdsp(j,i,kcn)*zetcl
          Ats = hdsp(j,i,kcs)*zetcl
      endif
      Atnp = dmax1( Atn, 0.0d0 )
      Atnm = Atn - Atnp
      Atsp = dmax1( Ats, 0.0d0 )
      Atsm = Ats - Atsp
      etcfE(jw,ia) = Ate + 0.5d0*(Atnp-Atsm)
      etcfW(jw,ia) = Atw + 0.5d0*(Atsp-Atnm)
!
!::viscous heating G32
      fvsh = 1.0d0
      Ahe = Ate*fvsh
      Ahw = Atw*fvsh
      Ahn = Atn*fvsh
      Ahs = Ats*fvsh
      ehcfE(jw,ia) = etcfE(jw,ia)*fvsh
      ehcfW(jw,ia) = etcfW(jw,ia)*fvsh
!
!::flux
      zf11 = - Ade*q1a(jp,i,ia) + Adw*q1a(j,i,ia) - Adn*zran + Ads*zras
      zf21 = zf11*xsva(jw,ia)
      zf31 = zf11*xsti(jw,ia) + zf31
      zf41 = zf11*xste(jw,ia) + zf41
      zf22 = - Ate*vva(jp,i,ia) + Atw*vva(j,i,ia) - Atn*zvan + Ats*zvas
      zf32 =(- Ahe*vva(jp,i,ia)**2 + Ahw*vva(j,i,ia)**2
     >       - Ahn*zvan**2 + Ahs*zvas**2)*0.5d0 + zf32
      fx_df(jw,m1a) = zf11
      fx_df(jw,m2a) = zf21
      fx_cd(jw,m2a) = zf22
 220  continue  !  loop (ia)
!
      fx_df(jw,m3) = zf31
      fx_df(jw,m4) = zf41
      fx_vh(jw)    = zf32
!
!::parameter (Ti,Te) at N(j+1/2,i+1/2), S(j+1/2,i-1/2)
!x      ztin= c14*(vti(j,i)+vti(j,i+1)+vti(jp,i)+vti(jp,i+1))
!x      ztis= c14*(vti(j,i)+vti(j,i-1)+vti(jp,i)+vti(jp,i-1))
!x      zten= c14*(vte(j,i)+vte(j,i+1)+vte(jp,i)+vte(jp,i+1))
!x      ztes= c14*(vte(j,i)+vte(j,i-1)+vte(jp,i)+vte(jp,i-1))
      ztin= vtig(j,i)
      ztis= vtig(j,i-1)
      zten= vteg(j,i)
      ztes= vteg(j,i-1)
!
!::d-coef Xi  at cell boundary (j+1/2,i)
!x    zxicl = 1.0d0/(hwtm(j,i)/xicl(j)+hwtp(j,i)/xicl(jp))
!x    zxian =  fcan/(hwtm(j,i)/xian(j)+hwtp(j,i)/xian(jp))
      zxicl =fundav(fccl,hwtm(j,i),hwtp(j,i),xicl(j),xicl(jp))
      zxian =fundav(fcan,hwtm(j,i),hwtp(j,i),xian(j),xian(jp))
!
!================================================================
!::flux limit of ion heat diffusivity
      Aie = hdsp(j,i,kce)*zxicl
      Aiw = hdsp(j,i,kcw)*zxicl
      zqcl = - Aie*vti(jp,i) + Aiw*vti(j,i)
      zqcl = zqcl*cev
!
!--new interpolation of flux at cell boundary
      znio = hwtm(j,i)*vni(j,i)+hwtp(j,i)*vni(jp,i)
      ztio = hwtm(j,i)*vti(j,i)+hwtp(j,i)*vti(jp,i)
      zmsi = aion(1)*cmp
      zqfl = flimi*znio*ztio*cev*dsqrt(ztio*cev/zmsi)*hvsb(j,i)
!
      zfxi = 1.0d0
      if( flimi.gt.0.0d0 ) zfxi = 1.0d0/(1.0d0+dabs(zqcl/zqfl))
      zxicl = zxicl*zfxi
      flmxi(j,i) = zfxi
!
!::ion heat flux density qi//
!   containing only //diffusive heat transfer
!   at cell boundary  j-jp, i.e. j+1/2. Unit [joule/(m2*s)]
!   Cf. Note page "sld mtc 21 & 24"             2018/05/17 Y.Homma
      Aie_limited = hdsp(j,i,kce)*zxicl
      Aiw_limited = hdsp(j,i,kcw)*zxicl
      heat_flux_density_para_temp_1
     &     = - Aie_limited*vti(jp,i) + Aiw_limited*vti(j,i)

!  Non-zero check: hvsb
      if(hvsb(j,i).le.0.d0) then
         write(n6,*)
       write(n6,*) 'For ion heat flux qi//, metric check: hvsb(j,i)'
         write(n6,'(/2x,"  j, i, hvsb(j,i): ",2i4,f9.2)')
     >        j, i, hvsb(j,i)
         call wexit("pxviscN.f",
     >     "For ion heat flux qi//, hvsb is not correctly defined.")
      endif

!  Unit transfer to [ joule/(m2.s) ]
      heat_flux_density_para_temp_2
     &     = heat_flux_density_para_temp_1 * cev / hvsb(j,i)
!
      heat_flux_para_by_kappa_para(j,i) = heat_flux_density_para_temp_2
!
!=====================================================================
!
      if(set_hdsv .ne. 0) then
          Aie = hdsp(j,i,kce)*zxicl+hdsv(j,i,kce)*zxian
          Aiw = hdsp(j,i,kcw)*zxicl+hdsv(j,i,kcw)*zxian
          Ain = hdsp(j,i,kcn)*zxicl+hdsv(j,i,kcn)*zxian
          Ais = hdsp(j,i,kcs)*zxicl+hdsv(j,i,kcs)*zxian
      else ! as default
          Aie = hdsp(j,i,kce)*zxicl
          Aiw = hdsp(j,i,kcw)*zxicl
          Ain = hdsp(j,i,kcn)*zxicl
          Ais = hdsp(j,i,kcs)*zxicl
      endif
      Ainp = dmax1( Ain, 0.0d0 )
      Ainm = Ain - Ainp
      Aisp = dmax1( Ais, 0.0d0 )
      Aism = Ais - Aisp
      xicfE(jw) = Aie + 0.5d0*(Ainp-Aism)
      xicfW(jw) = Aiw + 0.5d0*(Aisp-Ainm)
!
!::d-coef Xe  at cell boundary (j+1/2,i)
!x    zxecl = 1.0d0/(hwtm(j,i)/xecl(j)+hwtp(j,i)/xecl(jp))
!x    zxean =  fcan/(hwtm(j,i)/xean(j)+hwtp(j,i)/xean(jp))
      zxecl =fundav(fccl,hwtm(j,i),hwtp(j,i),xecl(j),xecl(jp))
      zxean =fundav(fcan,hwtm(j,i),hwtp(j,i),xean(j),xean(jp))
!
!================================================================
!::flux limit of electron heat diffusivity
      Aee = hdsp(j,i,kce)*zxecl
      Aew = hdsp(j,i,kcw)*zxecl
      zqcl = - Aee*vte(jp,i) + Aew*vte(j,i)
      zqcl = zqcl*cev
!
!--new interpolation of flux at cell boundary
      zneo = hwtm(j,i)*vne(j,i)+hwtp(j,i)*vne(jp,i)
      zteo = hwtm(j,i)*vte(j,i)+hwtp(j,i)*vte(jp,i)
      zqfl = flime*zneo*zteo*cev*dsqrt(zteo*cev/cme)*hvsb(j,i)
!
      zfxe = 1.0d0
      if( flime.gt.0.0d0 ) zfxe = 1.0d0/(1.0d0+dabs(zqcl/zqfl))
      zxecl = zxecl*zfxe
      flmxe(j,i) = zfxe
!
!:: electron heat flux density qe// only //diffusive heat transfer
!   at cell boundary  j-jp, i.e. j+1/2. Unit [joule/(m2*s)]
!   Cf. Note page "sld mtc 21 & 24"             2018/10/08 Y.Homma
      Aee_limited = hdsp(j,i,kce)*zxecl
      Aew_limited = hdsp(j,i,kcw)*zxecl
      heat_flux_para_elec_temp_1
     &     = - Aee_limited*vte(jp,i) + Aew_limited*vte(j,i)

!  Non-zero check: hvsb
      if(hvsb(j,i).le.0.d0) then
         write(n6,*)
       write(n6,*) 'For ion heat flux qe//, metric check: hvsb(j,i)'
         write(n6,'(/2x,"  j, i, hvsb(j,i): ",2i4,f9.2)')
     >        j, i, hvsb(j,i)
         call wexit("pxviscN.f",
     >     "For ion heat flux qe//, hvsb is not correctly defined.")
      endif

!  Unit transfer to [ joule/(m2.s) ]
      heat_flux_para_elec_temp_2
     &     = heat_flux_para_elec_temp_1 * cev / hvsb(j,i)
!
      heat_flux_para_elec(j,i) = heat_flux_para_elec_temp_2

!=====================================================================
!
      if(set_hdsv .ne. 0) then
          Aee = hdsp(j,i,kce)*zxecl+hdsv(j,i,kce)*zxean
          Aew = hdsp(j,i,kcw)*zxecl+hdsv(j,i,kcw)*zxean
          Aen = hdsp(j,i,kcn)*zxecl+hdsv(j,i,kcn)*zxean
          Aes = hdsp(j,i,kcs)*zxecl+hdsv(j,i,kcs)*zxean
      else ! as default
          Aee = hdsp(j,i,kce)*zxecl
          Aew = hdsp(j,i,kcw)*zxecl
          Aen = hdsp(j,i,kcn)*zxecl
          Aes = hdsp(j,i,kcs)*zxecl
      endif
      Aenp = dmax1( Aen, 0.0d0 )
      Aenm = Aen - Aenp
      Aesp = dmax1( Aes, 0.0d0 )
      Aesm = Aes - Aesp
      xecfE(jw) = Aee + 0.5d0*(Aenp-Aesm)
      xecfW(jw) = Aew + 0.5d0*(Aesp-Aenm)
!
!::flux
      zf33 = - Aie*vti(jp,i) + Aiw*vti(j,i) - Ain*ztin + Ais*ztis
      zf43 = - Aee*vte(jp,i) + Aew*vte(j,i) - Aen*zten + Aes*ztes
      fx_cd(jw,m3) = zf33*cev
      fx_cd(jw,m4) = zf43*cev
 210  continue  ! loop(jw)
!
!
!-----------------------------------------------------------------------
!::jacobian
!-----------------------------------------------------------------------
!
      do jw = 2, jmax1
      j  = jcel(jw,it)
      jp = jcel(jw+1,it)
      jm = jcel(jw-1,it)
      jw1 = jw-1
!
!::clear
      nequ = 2*nion + 2
      do k2 = 1, nequ
      do k1 = 1, nequ
      zcc(k1,k2) = 0.0d0
      zdp(k1,k2) = 0.0d0
      zdm(k1,k2) = 0.0d0
      zee(k1,k2) = 0.0d0
      wcc(k1,k2) = 0.0d0
      wdp(k1,k2) = 0.0d0
      wdm(k1,k2) = 0.0d0
      wee(k1,k2) = 0.0d0
      enddo
      enddo

!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!
!::equation [q1a = ma*na]
!
!::G11 = -Da*grad(q1a)
!x    zc11 = -dncfW(jw1,ia)
!x    zd11 =  dncfE(jw1,ia) + dncfW(jw,ia)
!x    ze11 = -dncfE(jw, ia)
!
      wcc(m1a,m1a) = wcc(m1a,m1a) + dacfW(jw1,ia)
      wdm(m1a,m1a) = wdm(m1a,m1a) - dacfE(jw1,ia)
      wdp(m1a,m1a) = wdp(m1a,m1a) + dacfW(jw, ia)
      wee(m1a,m1a) = wee(m1a,m1a) - dacfE(jw, ia)
!
!::equation [q2a = ma*na*vpa]
!
!::G21 = -va*Da*grad(q1a) =>  va*dlt_[-Da*grad(q1a)] dlt_[va]*flx = 0.0
!x    zc21 = -xsva(jw1,ia)*dacfW(jw1,ia)
!x    zd21 =  xsva(jw1,ia)*dacfE(jw1,ia)+xsva(jw,ia)*dacfW(jw,ia)
!x    ze21 = -xsva(jw, ia)*dacfE(jw, ia)
!
      wcc(m2a,m1a) = wcc(m2a,m1a) + xsva(jw1,ia)*dacfW(jw1,ia)
      wdm(m2a,m1a) = wdm(m2a,m1a) - xsva(jw1,ia)*dacfE(jw1,ia)
      wdp(m2a,m1a) = wdp(m2a,m1a) + xsva(jw, ia)*dacfW(jw, ia)
      wee(m2a,m1a) = wee(m2a,m1a) - xsva(jw, ia)*dacfE(jw, ia)
!
!::G22 = -eta*grad(va)
!x    zc22 = -etcfW(jw1,ia)
!x    zd22 =  etcfE(jw1,ia)+etcfW(jw,ia)
!x    ze22 = -etcfE(jw, ia)
!
      zcc(m2a,m1a) = zcc(m2a,m1a) + etcfW(jw1,ia)*r1av(jm,ia)
      zdm(m2a,m1a) = zdm(m2a,m1a) - etcfE(jw1,ia)*r1av(j, ia)
      zdp(m2a,m1a) = zdp(m2a,m1a) + etcfW(jw, ia)*r1av(j, ia)
      zee(m2a,m1a) = zee(m2a,m1a) - etcfE(jw, ia)*r1av(jp,ia)
      zcc(m2a,m2a) = zcc(m2a,m2a) + etcfW(jw1,ia)*r2av(jm,ia)
      zdm(m2a,m2a) = zdm(m2a,m2a) - etcfE(jw1,ia)*r2av(j, ia)
      zdp(m2a,m2a) = zdp(m2a,m2a) + etcfW(jw, ia)*r2av(j, ia)
      zee(m2a,m2a) = zee(m2a,m2a) - etcfE(jw, ia)*r2av(jp,ia)
!
!::equation [q3 = energy_I ]
!
!::G31 = -sum[ (ea+pa)/roa*Da*grad(q1a) ]
!x    zc31 = -xsti(jw1,ia)*dacfW(jw1,ia)
!x    zd31 =  xsti(jw1,ia)*dacfE(jw1,ia)+xsti(jw,ia)*dacfW(jw,ia)
!x    ze31 = -xsti(jw, ia)*dacfE(jw, ia)
!
      wcc(m3,m1a) = wcc(m3,m1a) + xsti(jw1,ia)*dacfW(jw1,ia)
      wdm(m3,m1a) = wdm(m3,m1a) - xsti(jw1,ia)*dacfE(jw1,ia)
      wdp(m3,m1a) = wdp(m3,m1a) + xsti(jw, ia)*dacfW(jw, ia)
      wee(m3,m1a) = wee(m3,m1a) - xsti(jw, ia)*dacfE(jw, ia)
!
!::G32 = -sum[ eta*grad(1/2*va**2) ]
!x    zc32 = -ehcfW(jw1,ia)
!x    zd32 =  ehcfE(jw1,ia)+ehcfW(jw,ia)
!x    ze32 = -ehcfE(jw, ia)
!
      zcc(m3,m1a) = zcc(m3,m1a) + ehcfW(jw1,ia)*r1aw(jm,ia)
      zdm(m3,m1a) = zdm(m3,m1a) - ehcfE(jw1,ia)*r1aw(j, ia)
      zdp(m3,m1a) = zdp(m3,m1a) + ehcfW(jw, ia)*r1aw(j, ia)
      zee(m3,m1a) = zee(m3,m1a) - ehcfE(jw, ia)*r1aw(jp,ia)
      zcc(m3,m2a) = zcc(m3,m2a) + ehcfW(jw1,ia)*r2aw(jm,ia)
      zdm(m3,m2a) = zdm(m3,m2a) - ehcfE(jw1,ia)*r2aw(j, ia)
      zdp(m3,m2a) = zdp(m3,m2a) + ehcfW(jw, ia)*r2aw(j, ia)
      zee(m3,m2a) = zee(m3,m2a) - ehcfE(jw, ia)*r2aw(jp,ia)
!
!::G33 = -Ki*grad(Ti)
!x    zc33 = -xicfW(jw1)
!x    zd33 =  xicfE(jw1)+xicfW(jw)
!x    ze33 = -xicfE(jw )
!
      zcc(m3,m1a) = zcc(m3,m1a) + xicfW(jw1)*r1ai(jm,ia)
      zdm(m3,m1a) = zdm(m3,m1a) - xicfE(jw1)*r1ai(j, ia)
      zdp(m3,m1a) = zdp(m3,m1a) + xicfW(jw )*r1ai(j, ia)
      zee(m3,m1a) = zee(m3,m1a) - xicfE(jw )*r1ai(jp,ia)
      zcc(m3,m2a) = zcc(m3,m2a) + xicfW(jw1)*r2ai(jm,ia)
      zdm(m3,m2a) = zdm(m3,m2a) - xicfE(jw1)*r2ai(j, ia)
      zdp(m3,m2a) = zdp(m3,m2a) + xicfW(jw )*r2ai(j, ia)
      zee(m3,m2a) = zee(m3,m2a) - xicfE(jw )*r2ai(jp,ia)
!
!::equation [q4 = energy_e ]
!
!::G41 = -sum[ 5/2*Te*za/ma*Da*grad(q1a) ]
!x    zc41 = -xste(jw1,ia)*dacfW(jw1,ia)
!x    zd41 =  xste(jw1,ia)*dacfE(jw1,ia)+xste(jw,ia)*dacfW(jw,ia)
!x    ze41 = -xste(jw, ia)*dacfE(jw, ia)
!
      wcc(m4,m1a) = wcc(m4,m1a) + xste(jw1,ia)*dacfW(jw1,ia)
      wdm(m4,m1a) = wdm(m4,m1a) - xste(jw1,ia)*dacfE(jw1,ia)
      wdp(m4,m1a) = wdp(m4,m1a) + xste(jw, ia)*dacfW(jw, ia)
      wee(m4,m1a) = wee(m4,m1a) - xste(jw, ia)*dacfE(jw, ia)
!
!::G43 = -Ke*grad(Te)
!x    zc43 = -xecfW(jw1)
!x    zd43 =  xecfE(jw1)+xecfW(jw)
!x    ze43 = -xecfE(jw)
!
      zcc(m4,m1a) = zcc(m4,m1a) + xecfW(jw1)*r1ae(jm,ia)
      zdm(m4,m1a) = zdm(m4,m1a) - xecfE(jw1)*r1ae(j, ia)
      zdp(m4,m1a) = zdp(m4,m1a) + xecfW(jw )*r1ae(j, ia)
      zee(m4,m1a) = zee(m4,m1a) - xecfE(jw )*r1ae(jp,ia)
!
!::ff
!x      ff(m1a,jw) = ff(m1a,jw) - dflna(jw,ia) + dflna(jw1,ia)
!x      ff(m2a,jw) = ff(m2a,jw) - dfl2a(jw,ia) + dfl2a(jw1,ia)
!
      enddo  ! loop(ia)
!
!::G33
      zcc(m3,m3) = zcc(m3,m3) + xicfW(jw1)*r3i(jm)
      zdm(m3,m3) = zdm(m3,m3) - xicfE(jw1)*r3i(j )
      zdp(m3,m3) = zdp(m3,m3) + xicfW(jw )*r3i(j )
      zee(m3,m3) = zee(m3,m3) - xicfE(jw )*r3i(jp)
!
!::G43
      zcc(m4,m4) = zcc(m4,m4) + xecfW(jw1)*r4e(jm)
      zdm(m4,m4) = zdm(m4,m4) - xecfE(jw1)*r4e(j )
      zdp(m4,m4) = zdp(m4,m4) + xecfW(jw )*r4e(j )
      zee(m4,m4) = zee(m4,m4) - xecfE(jw )*r4e(jp)
!
!:cc,dd,ee
      do m = 1, m4
      do n = 1, m4
      cc(m,n,jw) = cc(m,n,jw) -wcc(m,n)-zcc(m,n)
      dd(m,n,jw) = dd(m,n,jw) -wdm(m,n)-zdm(m,n)+wdp(m,n)+zdp(m,n)
      ee(m,n,jw) = ee(m,n,jw) +wee(m,n)+zee(m,n)
      enddo
      enddo
!
!::ff
      do m = 1, m4
      ff(m,jw) = ff(m,jw) - fx_df(jw,m) + fx_df(jw1,m)
     >                    - fx_cd(jw,m) + fx_cd(jw1,m)
      enddo
      ff(m3,jw) = ff(m3,jw) - fx_vh(jw) + fx_vh(jw1)
!
!::bc-condition  ! %%% 2002/11/11
      if( jw.eq.2 ) then
      jwb = jw - 1
      do m = 1, m4
      do n = 1, m4
      dd(m,n,jwb) = dd(m,n,jwb) + zcc(m,n)
      ee(m,n,jwb) = ee(m,n,jwb) + zdm(m,n)
      enddo
      enddo
      endif
!
      if( jw.eq.jmax1 ) then
      jwb = jw + 1
      do m = 1, m4
      do n = 1, m4
      cc(m,n,jwb) = cc(m,n,jwb) - zdp(m,n)
      dd(m,n,jwb) = dd(m,n,jwb) - zee(m,n)
      enddo
      enddo
      endif
!
      enddo  ! loop(jw)
!
!::ff at boundary  ! %%% 2002/11/12
      jwb = 1
      jw  = 1
      do m = 1, m4
      ff(m,jwb) = ff(m,jwb) - fx_df(jw,m) - fx_cd(jw,m)
      enddo
      ff(m3,jwb) = ff(m3,jwb) - fx_vh(jw)
      jwb = jmax
      jw  = jmax1
      do m = 1, m4
      ff(m,jwb) = ff(m,jwb) + fx_df(jw,m) + fx_cd(jw,m)
      enddo
      ff(m3,jwb) = ff(m3,jwb) + fx_vh(jw)
!
!::qcon  ! %%%  2002/11/11
      do jw = 1, jmax1
      j = jcel(jw,it)
      do m = 1, m4
      qfx_df(j,i,m) = fx_df(jw,m)
      qfx_cd(j,i,m) = fx_cd(jw,m)
      enddo
      qfx_vh(j,i) = fx_vh(jw)
      enddo
!
      return
      end

