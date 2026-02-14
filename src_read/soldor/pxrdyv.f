!***********************************************************************
      subroutine pxrdyv(it)
!***********************************************************************
!
!      residuals of F^psi  +  flux
!
!       term eta_an*grad(vp) is dropped    ! %%% 2002/10/24
!       1st order up-wind scheme           ! %%% 2002/12/19  (pxrdyv)
!
!      hia = (vea+Na*Ti*cev)/q1a ==> (vea+Na*Ti*cevpr)/q1a
!      hea = c52*Za/Ma*Te*cev    ==> chfcv*Za/Ma*Te*cev
!-----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : ama, aza, cevpr, chfcv, gg, nion, q1a, vdet
     >    , vdxe, vdxi, vea, vna, vne, vni, vte, vteg, vti, vtig, vva
     >    , vvag, xfdfa, xwtam, xwtap
      use cplmet, only : gdsv, gwtm, gwtp, icel, jcel, jtmax, kce, kcn
     >    , kcs, kcw
      use cplqcn, only : qfy_cd, qfy_df, qfy_vh
      use cmeffz, only : xzflz
      use csonic, only : lfopt
      implicit none
!
      integer, intent(in) :: it
!
      integer  i, ip, im, ihp, ihm, jmax, jmax1, jw, j, jm
      integer  ia, m1a, m2a, m3, m4
      real*8   Atep, Atwp, Atnp, Atsp, Atem, Atwm, Atnm, Atsm
      real*8   Aiep, Aiwp, Ainp, Aisp, Aiem, Aiwm, Ainm, Aism
      real*8   Aeep, Aewp, Aenp, Aesp, Aeem, Aewm, Aenm, Aesm
      real*8   zf11p, zf21p, zf31p, zf41p, zf11m, zf21m, zf31m, zf41m
      real*8   zf22p, zf32p, zf33p, zf43p, zf22m, zf32m, zf33m, zf43m
      real*8   zf1p,  zf2p,  zf3p,  zf4p,  zf1m,  zf2m,  zf3m,  zf4m
      real*8   zetc,zetp,zetm,zethp,zethm
      real*8   zxic,zxip,zxim,zxihp,zxihm, zxec,zxep,zxem,zxehp,zxehm
      real*8   xvahp, xtihp, xtehp, xvahm, xtihm, xtehm
      real*8   zhic, zhip, zhim, zhec, zhep, zhem
      real*8   zvaep, zvawp, zvaem, zvawm
      real*8   ztiep, ztiwp, ztiem, ztiwm, zteep, ztewp, zteem, ztewm
      real(8) :: ztep, ztem, zf4zp, zf4zm
      real*8   zwa1, zwa2, zwb1, zwb2
!
      if( lfopt(2).eq.0 ) return
!
      i   = icel(1,it)
      ip  = i + 1
      im  = i - 1
      ihp = i
      ihm = i - 1
      jmax  = jtmax(it)
      jmax1 = jmax-1
!
      m3 = 2*nion + 1
      m4 = 2*nion + 2
!
!-----------------------------------------------------------------------
!::loop x-direction
!-----------------------------------------------------------------------
      do jw = 2, jmax1
        j  = jcel(jw,it)
        jm = jcel(jw-1,it)
!
        zf31p = 0.0d0
        zf32p = 0.0d0
        zf41p = 0.0d0
        zf31m = 0.0d0
        zf32m = 0.0d0
        zf41m = 0.0d0
!
        do ia = 1, nion
          m1a = 2*ia - 1
          m2a = 2*ia
!
!::d-coef Eta  at cell boundary (j,i+1/2) & (j,i-1/2)
          zetc = vdet(j,i, ia)*ama(ia)*vna(j,i, ia)
          zetp = vdet(j,ip,ia)*ama(ia)*vna(j,ip,ia)
          zetm = vdet(j,im,ia)*ama(ia)*vna(j,im,ia)
          zethp =1.0d0/(gwtm(j,ihp)/zetc+gwtp(j,ihp)/zetp)
          zethm =1.0d0/(gwtm(j,ihm)/zetm+gwtp(j,ihm)/zetc)
!
          Atep = gdsv(j,ihp,kce)*zethp
          Atwp = gdsv(j,ihp,kcw)*zethp
          Atnp = gdsv(j,ihp,kcn)*zethp
          Atsp = gdsv(j,ihp,kcs)*zethp
          Atem = gdsv(j,ihm,kce)*zethm
          Atwm = gdsv(j,ihm,kcw)*zethm
          Atnm = gdsv(j,ihm,kcn)*zethm
          Atsm = gdsv(j,ihm,kcs)*zethm
!
!::parameter (roa,va) at E(j+1/2,i+-1/2), W(j-1/2,i+-1/2)
!x    zraep = c14*(q1a(j,i,ia)+q1a(jp,i,ia)+q1a(j,ip,ia)+q1a(jp,ip,ia))
          zvaep = vvag(j, i,ia)
          zvawp = vvag(jm,i,ia)
          zvaem = vvag(j, im,ia)
          zvawm = vvag(jm,im,ia)
!
!::X-point   ( 00/10/30 )
!
!::X-interpolation (X) at O(j,i+-1/2)
          zhic  = 
     >     (vea(j,i, ia)+vna(j,i, ia)*vti(j,i )*cevpr)/q1a(j,i, ia)
          zhip  = 
     >     (vea(j,ip,ia)+vna(j,ip,ia)*vti(j,ip)*cevpr)/q1a(j,ip,ia)
          zhim  = 
     >     (vea(j,im,ia)+vna(j,im,ia)*vti(j,im)*cevpr)/q1a(j,im,ia)
          zhec  = chfcv*aza(ia)/ama(ia)*vte(j,i )*cev
          zhep  = chfcv*aza(ia)/ama(ia)*vte(j,ip)*cev
          zhem  = chfcv*aza(ia)/ama(ia)*vte(j,im)*cev
!
!::weit
          zwa1 = xwtam(j,i,ia)
          zwa2 = xwtap(j,i,ia)
          zwb1 = xwtam(j,im,ia)
          zwb2 = xwtap(j,im,ia)
!
!::X at i+1/2 and i-1/2
          xvahp = zwa1*vva(j,i, ia) + zwa2*vva(j,ip,ia)
          xtihp = zwa1*zhic         + zwa2*zhip
          xtehp = zwa1*zhec         + zwa2*zhep
          xvahm = zwb1*vva(j,im,ia) + zwb2*vva(j,i, ia)
          xtihm = zwb1*zhim         + zwb2*zhic
          xtehm = zwb1*zhem         + zwb2*zhec
!
!::flux  (j,i+1/2)
          zf11p = xfdfa(j,i,ia)
          zf21p = zf11p*xvahp
          zf31p = zf11p*xtihp + zf31p
          zf41p = zf11p*xtehp + zf41p
          zf22p = -Atep*zvaep+Atwp*zvawp
     >            -Atnp*vva(j,ip,ia)+Atsp*vva(j,i,ia)
          zf32p =(-Atep*zvaep**2+Atwp*zvawp**2-Atnp*vva(j,ip,ia)**2
     >        +Atsp*vva(j,i,ia)**2)*0.5d0 + zf32p
          zf1p  =        zf11p
          zf2p  =        zf21p + zf22p
!
!::flux  (j,i-1/2)
          zf11m = xfdfa(j,im,ia)
          zf21m = zf11m*xvahm
          zf31m = zf11m*xtihm + zf31m
          zf41m = zf11m*xtehm + zf41m
          zf22m = -Atem*zvaem+Atwm*zvawm
     >        -Atnm*vva(j,i,ia)+Atsm*vva(j,im,ia)
          zf32m =(-Atem*zvaem**2+Atwm*zvawm**2-Atnm*vva(j,i,ia)**2
     >        +Atsm*vva(j,im,ia)**2)*0.5d0 + zf32m
          zf1m  =        zf11m
          zf2m  =        zf21m + zf22m
!
!::residuals
          gg(m1a,jw) = gg(m1a,jw) - zf1p + zf1m
          gg(m2a,jw) = gg(m2a,jw) - zf2p + zf2m
!
!::qcon
          qfy_df(j,i,m1a) = zf11p
          qfy_df(j,i,m2a) = zf21p
          qfy_cd(j,i,m2a) = zf22p
          if( it.eq.2 ) then
            qfy_df(j,im,m1a) = zf11m
            qfy_df(j,im,m2a) = zf21m
            qfy_cd(j,im,m2a) = zf22m
          endif
!
        enddo  ! loop(ia)
!
!::d-coef Xi  at cell boundary (j,i+1/2) & (j,i-1/2)
        zxic  = vdxi(j,i )*vni(j,i )
        zxip  = vdxi(j,ip)*vni(j,ip)
        zxim  = vdxi(j,im)*vni(j,im)
        zxihp = 1.0d0/(gwtm(j,ihp)/zxic+gwtp(j,ihp)/zxip)
        zxihm = 1.0d0/(gwtm(j,ihm)/zxim+gwtp(j,ihm)/zxic)
!
        Aiep = gdsv(j,ihp,kce)*zxihp
        Aiwp = gdsv(j,ihp,kcw)*zxihp
        Ainp = gdsv(j,ihp,kcn)*zxihp
        Aisp = gdsv(j,ihp,kcs)*zxihp
        Aiem = gdsv(j,ihm,kce)*zxihm
        Aiwm = gdsv(j,ihm,kcw)*zxihm
        Ainm = gdsv(j,ihm,kcn)*zxihm
        Aism = gdsv(j,ihm,kcs)*zxihm
!
!::d-coef Xe  at cell boundary (j,i+1/2)
        zxec  = vdxe(j,i )*vne(j,i )
        zxep  = vdxe(j,ip)*vne(j,ip)
        zxem  = vdxe(j,im)*vne(j,im)
        zxehp = 1.0d0/(gwtm(j,ihp)/zxec+gwtp(j,ihp)/zxep)
        zxehm = 1.0d0/(gwtm(j,ihm)/zxem+gwtp(j,ihm)/zxec)
!
        Aeep = gdsv(j,ihp,kce)*zxehp
        Aewp = gdsv(j,ihp,kcw)*zxehp
        Aenp = gdsv(j,ihp,kcn)*zxehp
        Aesp = gdsv(j,ihp,kcs)*zxehp
        Aeem = gdsv(j,ihm,kce)*zxehm
        Aewm = gdsv(j,ihm,kcw)*zxehm
        Aenm = gdsv(j,ihm,kcn)*zxehm
        Aesm = gdsv(j,ihm,kcs)*zxehm
!
!::parameter (Ti,Te) at E(j+1/2,i+-1/2), W(j-1/2,i+-1/2)
!x    ztiep = c14*(vti(j,i)+vti(jp,i)+vti(j,ip)+vti(jp,ip))
!
        ztiep = vtig(j ,i)
        ztiwp = vtig(jm,i)
        zteep = vteg(j ,i)
        ztewp = vteg(jm,i)
        ztiem = vtig(j, im)
        ztiwm = vtig(jm,im)
        zteem = vteg(j, im)
        ztewm = vteg(jm,im)
!
!::5/2*Fz(i+1/2)*Te(i+1/2)
        ztep  = gwtp(j,i)*vte(j,i) + gwtm(j,i)*vte(j,ip)
        zf4zp = chfcv*xzflz(j,i)*ztep*cev
!
!::5/2*Fz(i-1/2)*Te(i-1/2)
        ztem  = gwtp(j,im)*vte(j,im) + gwtm(j,im)*vte(j,i)
        zf4zm = chfcv*xzflz(j,im)*ztem*cev
!
!::flux (j,i+1/2)
        zf33p =(-Aiep*ztiep+Aiwp*ztiwp
     >          -Ainp*vti(j,ip)+Aisp*vti(j,i))*cev
        zf43p =(-Aeep*zteep+Aewp*ztewp
     >          -Aenp*vte(j,ip)+Aesp*vte(j,i))*cev
        zf3p  = zf31p + zf32p + zf33p
        zf4p  = zf41p         + zf43p + zf4zp
!
!::flux (j,i-1/2)
        zf33m =(-Aiem*ztiem+Aiwm*ztiwm
     >          -Ainm*vti(j,i)+Aism*vti(j,im))*cev
        zf43m =(-Aeem*zteem+Aewm*ztewm
     >          -Aenm*vte(j,i)+Aesm*vte(j,im))*cev
        zf3m  = zf31m + zf32m + zf33m
        zf4m  = zf41m         + zf43m + zf4zm
!
!::residuals
        gg(m3,jw) = gg(m3,jw) - zf3p + zf3m
        gg(m4,jw) = gg(m4,jw) - zf4p + zf4m
!
!::qcon
        qfy_df(j,i,m3) = zf31p
        qfy_df(j,i,m4) = zf41p
        qfy_vh(j,i)    = zf32p
        qfy_cd(j,i,m3) = zf33p
        qfy_cd(j,i,m4) = zf43p
        if( it.eq.2 ) then
          qfy_df(j,im,m3) = zf31m
          qfy_df(j,im,m4) = zf41m
          qfy_vh(j,im)    = zf32m
          qfy_cd(j,im,m3) = zf33m
          qfy_cd(j,im,m4) = zf43m
        endif
!
      enddo  ! loop(jw)
      end