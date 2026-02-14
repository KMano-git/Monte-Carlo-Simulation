!***********************************************************************
      subroutine plmflx
!***********************************************************************
!
!     Fa^psi(j)  = - Da*grad(roa)
!     Qi^psi(j)  = - sum[(ea+pa)/roa*Da*grad(roa)]      - Ki*grad(Ti)
!     Qe^psi(j)  = - sum[(ee+pe)/ne*za/ma*Da*grad(roa)] - Ke*grad(Te)
!
!     Fa = - sum_(j,edge)[Fa^psi(j)] = 2.0d22
!     Qi = - sum_(j,edge)[Qi^psi(j)] = 4 MW
!     Qe = - sum_(j,edge)[Qe^psi(j)] = 8 MW
!
!      hia = (vea+Na*Ti*cev)/q1a ==> (vea+Na*Ti*cevpr)/q1a
!      hea = c52*Za/Ma*Te*cev    ==> chfcv*Za/Ma*Te*cev
!
!      hia(i+1/2) = gwtp(i+1/2)*X(i) + gwtm(i+1/2)*X(ip)  (i+1/2) = (i)
!
!-----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : ama, aza, cevpr, chfcv, nion, q1a, tflna, tflne
     >    , tflni, tflqe, tflqi, vdda, vdxe, vdxi, vea, vna, vnag, vne
     >    , vni, vte, vteg, vti, vtig
      use cplmet, only : gdsv, gwtm, gwtp, icmpe, itmpe, jcel, jtmax
     >    , kce, kcn, kcs, kcw
      use cmeffz, only : xzflz
      use csize,  only : ndsp
      use csonic, only : itim, time
      use cunit,  only : n6
      implicit none
!
      integer :: i, ip, it
      integer :: jmax, jmax1, jw, j, jm, ia
      real*8   Ade,Adw,Adn,Ads, Aie,Aiw,Ain,Ais, Aee,Aew,Aen,Aes
      real*8   zqiv, zqic, zqev, zqec, zda1, zda2, zda, zqevz
      real*8   zrae, zraw,  zhia1, zhia2, zhia, zhea1, zhea2, zhea, zhte
      real*8   zxi1, zxi2, zxi, ztie, ztiw, zxe1, zxe2, zxe, ztee, ztew
      real*8   zf11, zf33, zf43, tflnez
      real*8   zfa(ndsp)
!
!::index
      i  = icmpe - 1
      ip = i + 1
!
!::clear
      call setd( zfa, ndsp, 0.0d0 )
      zqiv = 0.0d0
      zqic = 0.0d0
      zqev = 0.0d0
      zqec = 0.0d0
      zqevz = 0.0d0
      tflnez = 0.0d0
!
!::integrate
      it = itmpe
      jmax  = jtmax(it)
      jmax1 = jmax-1
      do 210 jw = 2, jmax1
      j  = jcel(jw,it)
      jm = jcel(jw-1,it)
!
      do 220 ia = 1, nion
!
      zda1 = vdda(j,i, ia)
      zda2 = vdda(j,ip,ia)
      zda  = 1.0d0/(gwtm(j,i)/zda1+gwtp(j,i)/zda2)
      Ade  = gdsv(j,i,kce)*zda
      Adw  = gdsv(j,i,kcw)*zda
      Adn  = gdsv(j,i,kcn)*zda
      Ads  = gdsv(j,i,kcs)*zda
!x      zrae = c14*(q1a(j,i,ia)+q1a(jp,i,ia)+q1a(j,ip,ia)+q1a(jp,ip,ia))
!x      zraw = c14*(q1a(j,i,ia)+q1a(jm,i,ia)+q1a(j,ip,ia)+q1a(jm,ip,ia))
      zrae = ama(ia)*vnag(j, i,ia)
      zraw = ama(ia)*vnag(jm,i,ia)
      zf11 = -Ade*zrae + Adw*zraw - Adn*q1a(j,ip,ia) + Ads*q1a(j,i,ia)
      zfa(ia) = zfa(ia) + zf11/ama(ia)
!
      zhia1 = (vea(j,i ,ia)+vna(j,i ,ia)*vti(j, i)*cevpr)/q1a(j,i ,ia)
      zhia2 = (vea(j,ip,ia)+vna(j,ip,ia)*vti(j,ip)*cevpr)/q1a(j,ip,ia)
      zhia  = gwtp(j,i)*zhia1 + gwtm(j,i)*zhia2
      zqiv  = zqiv + zf11*zhia
!
      zhea1 = chfcv*aza(ia)/ama(ia)*vte(j,i )*cev
      zhea2 = chfcv*aza(ia)/ama(ia)*vte(j,ip)*cev
      zhea  = gwtp(j,i)*zhea1 + gwtm(j,i)*zhea2
      zqev  = zqev + zf11*zhea
!
 220  continue
!
!::impurity  xzflz = -Azn*dnz(j,ip) + Azs*dnz(j,i)
      zhte = gwtm(j,i)*vte(j,i) + gwtp(j,i)*vte(j,ip)
      zqevz = zqevz + chfcv*xzflz(j,i)*zhte*cev
      tflnez = tflnez + xzflz(j,i)
!
      zxi1 = vdxi(j,i )*vni(j,i )
      zxi2 = vdxi(j,ip)*vni(j,ip)
      zxi  = 1.0d0/(gwtm(j,i)/zxi1+gwtp(j,i)/zxi2)
      Aie  = gdsv(j,i,kce)*zxi
      Aiw  = gdsv(j,i,kcw)*zxi
      Ain  = gdsv(j,i,kcn)*zxi
      Ais  = gdsv(j,i,kcs)*zxi
      ztie = vtig(j,i)
      ztiw = vtig(jm,i)
      zf33 = -Aie*ztie + Aiw*ztiw - Ain*vti(j,ip) + Ais*vti(j,i)
      zqic = zqic + zf33*cev
!
      zxe1 = vdxe(j,i )*vne(j,i )
      zxe2 = vdxe(j,ip)*vne(j,ip)
      zxe  = 1.0d0/(gwtm(j,i)/zxe1+gwtp(j,i)/zxe2)
      Aee  = gdsv(j,i,kce)*zxe
      Aew  = gdsv(j,i,kcw)*zxe
      Aen  = gdsv(j,i,kcn)*zxe
      Aes  = gdsv(j,i,kcs)*zxe
      ztee = vteg(j,i)
      ztew = vteg(jm,i)
      zf43 = -Aee*ztee + Aew*ztew - Aen*vte(j,ip) + Aes*vte(j,i)
      zqec = zqec + zf43*cev
!
 210  continue
!
!
!::direction of flux
      zqev = zqev + zqevz
      zqiv = -zqiv
      zqic = -zqic
      zqev = -zqev
      zqec = -zqec
      tflnez = -tflnez
!
!::common variabes
      tflni = 0.0d0
      tflne = 0.0d0
      do 310 ia = 1, nion
      zfa(ia) = -zfa(ia)
      tflna(ia) = zfa(ia)
      tflni = tflni + tflna(ia)
      tflne = tflne + aza(ia)*tflna(ia)
 310  continue
      tflne = tflne + tflnez
      tflqi = zqiv + zqic
      tflqe = zqev + zqec
!
!::debug write
      j  = (jcel(2,it)+jcel(jmax1,it))/2
      ia = 1
!
      write(n6,'(2x,"*** plmflx",i7,1pe12.4," Fi,Qv,Qc =",1p3e10.3,
     >  "  Fe,Qv,Qc =",1p4e10.3,"  na,Ti,Te =",1p4e10.3)')
     >    itim, time, tflni, zqiv, zqic, tflne, tflnez, zqev,zqec,
     >    vna(j,ip,ia), vne(j,ip), vti(j,ip), vte(j,ip)
!
      return
      end
