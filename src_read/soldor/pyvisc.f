!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   subroutine pyvisc(nj)
      subroutine pyvisc(j)
!***********************************************************************
!
!    Ad_(E,W,N,S) = gdsv*D^ano
!
!    dacfN(j,i) = AdN + 1/2*(AdE(+)-AdW(-))
!    dacfS(j,i) = AdS + 1/2*(AdW(+)-AdE(-))
!
!-----------------------------------------------------------------------
!
!     G1a = - Da*grad(roa)
!     G2a = - va*Da*grad(roa) - eta*grad(va)
!     G3  = - sum[(ea+pa)/roa*Da*grad(roa)] - sum[eta*grad(1/2*va**2)]
!           - Ki*grad(Ti)
!     G4  = - sum[(ee+pe)/ne*za/ma*Da*grad(roa)] - Ke*grad(Te)
!
!         note.  roa = Q1a
!                hia = (ea+pa)/roa
!                hea = (ee+pe)/ne*za/ma = 5/2*Te*za/ma
!
!     new X-interpolation                                   00/12/14
!     Dif|i+1/2   interpolation method                      03/07/26
!
!     hia = (vea+Na*Ti*cev)/q1a ==> (vea+Na*Ti*cevpr)/q1a
!     hea = c52*Za/Ma*Te*cev    ==> chfcv*Za*Ma*Te*cev
!
!   impurity effect G4z = c52*Te*ceV*flez|(i+1/2)         2012/12/21
!     flez = -Dz*sum(z)Zz*grad(Nz)
!     Te(i+1/2) = gwtp*Te(i) + gwtm*Te(ip)   Not upwind  temporary
!
!     r1ae=dTe/dq1a,  r4e= dTe/dq4
!
!        You must modify pxrdyv.f and pyvisc.f simultaneously.
!
!-----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : ama, aza, c13, c23, cc, cevpr, chfcv, daan, dd
     >    , dfl1a, dfl2a, dfl3, dfl4, ee, etan, ff, hea, hia, nion, q1a
     >    , q2a, r1ae, r1ai, r1av, r1aw, r2ai, r2av, r2aw, r3i, r4e
     >    , vdda, vdet, vdxe, vdxi, vea, vna, vne, vni, vte, vteg, vti
     >    , vtig, vva, vvag, xean, xfdfa, xian, xste, xsti, xsva, xwtam
     >    , xwtap
      use cplmet, only : gdsv, gwtm, gwtp, icmax, jnxm, kce, kcn, kcs
     >    , kcw
      use cmeffz, only : xzflz
      use csize,  only : ndeq, ndsp, ndxy, ndy
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nj
      integer, intent(in) :: j
!
!::local variables
      integer k1, k2
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer jc, j, imax, imax1, iw, i, nd1, nd2, ia, ip, jp, jm
!ik   integer  im, iw1
      integer imax, imax1, i, ia, ip, jm
      integer  im
      real*8  Ade, Adw, Adn, Ads, Adep, Adem, Adwp, Adwm
      real*8  Ate, Atw, Atn, Ats, Atep, Atem, Atwp, Atwm
      real*8  Aie, Aiw, Ain, Ais, Aiep, Aiem, Aiwp, Aiwm
      real*8  Aee, Aew, Aen, Aes, Aeep, Aeem, Aewp, Aewm
      real*8  dacfN(ndy,ndsp),dacfS(ndy,ndsp)
      real*8  etcfN(ndy,ndeq),etcfS(ndy,ndeq)
      real*8  xicfN(ndy), xicfS(ndy), xecfN(ndy), xecfS(ndy)
      real(8), dimension(ndy) ::  xstz
!
      real*8  zdaan, zetan, zxian, zxean
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8  zrae, zraw, zvae, zvaw, ztie, ztiw, ztee, ztew
      real*8  zvae, zvaw, ztie, ztiw, ztee, ztew
      real*8  xc21, xd21, xe21
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8  zf11, zf21, zf1, zf41, zf22, zf31, zf32, zf33, zf43
      real*8  zf11, zf21, zf41, zf22, zf31, zf32, zf33, zf43
      real*8  zc11, zd11, ze11, zc21, zd21, ze21
      real*8  zc22, zd22, ze22, zc31, zd31, ze31
      real*8  zc32, zd32, ze32, zc33, zd33, ze33
      real*8  zc41, zd41, ze41, zc43, zd43, ze43
      real(8) :: xhte, zc4z, zd4z, ze4z, zf4z
      integer  m1a, m2a, m3, m4
!
      real*8  zwt1, zwt2
!
!::interpolation method for diffusion at j+1/2
      real*8  fundav, fdfc, fdwm, fdwp, fddm, fddp
! function
      fundav(fdfc,fdwm,fdwp,fddm,fddp) = fdfc/(fdwm/fddm+fdwp/fddp)
!
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   jc = nj
!ik   j  = nj
!ik   imax  = icmax(jc)
      imax  = icmax(j)
      imax1 = imax - 1
!
!x      write(n6,'(2x,"***  pyvisc  ***   jc =",i3,2i4)') jc,itim,nlp
!
!-----------------------------------------------------------------------
!::flux clear
!-----------------------------------------------------------------------
      do k2=1,nion
      do k1=1,ndxy
      dfl1a(k1,k2) = 0.0d0
      dfl2a(k1,k2) = 0.0d0
      enddo
      enddo
!
      do k1=1,ndxy
      dfl3(k1) = 0.0d0
      dfl4(k1) = 0.0d0
      enddo
!
!-----------------------------------------------------------------------
!::anomalous-diffusion coefficients
!::plasma parameter & trans-function at cell center (j,i)
!-----------------------------------------------------------------------
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   do 110 iw = 1, imax
!ik   i = iw
      do 110 i = 1, imax
      do 115 ia = 1, nion
      daan(i,ia) = vdda(j,i,ia)
      etan(i,ia) = vdet(j,i,ia)*ama(ia)*vna(j,i,ia)
      hia(i,ia)  = (vea(j,i,ia)+vna(j,i,ia)*vti(j,i)*cevpr)/q1a(j,i,ia)
      hea(i,ia)  =  chfcv*aza(ia)/ama(ia)*vte(j,i)*cev
!
      r1av(i,ia) = -q2a(j,i,ia)/q1a(j,i,ia)**2
      r2av(i,ia) = 1.0d0/q1a(j,i,ia)
      r1aw(i,ia) = vva(j,i,ia)*r1av(i,ia)
      r2aw(i,ia) = vva(j,i,ia)*r2av(i,ia)
      r1ai(i,ia) = (c13*vva(j,i,ia)**2-vti(j,i)*cev/ama(ia))/vni(j,i)
      r2ai(i,ia) = -c23*vva(j,i,ia)/vni(j,i)
      r1ae(i,ia) = -aza(ia)/ama(ia)*vte(j,i)*cev/vne(j,i)
 115  continue
      xian(i) = vdxi(j,i)*vni(j,i)
      xean(i) = vdxe(j,i)*vne(j,i)
!
      r3i(i) = c23/vni(j,i)
      r4e(i) = c23/vne(j,i)
 110  continue
!
!-----------------------------------------------------------------------
!::diff. coef. & flux at cell boundary (j,i+1/2)
!-----------------------------------------------------------------------
! modified 4/2 lines organize local variables and include files by kamata 2021/05/31
!ik   do 210 iw = 1, imax1
!ik   i   = iw
!ik   ip  = iw+1
!ik   jp  = jnxp(j,i)
      do 210 i = 1, imax1
      ip  = i+1
      jm  = jnxm(j,i)
!
      do 220 ia = 1, nion
!
!::d-coef (Da) at cell boundary (j,i+1/2)
!x    zdaan =1.0d0/(gwtm(j,i)/daan(i,ia)+gwtp(j,i)/daan(ip,ia))
      zdaan =fundav(1.0d0,gwtm(j,i),gwtp(j,i),daan(i,ia),daan(ip,ia))
!
      Ade = gdsv(j,i,kce)*zdaan
      Adw = gdsv(j,i,kcw)*zdaan
      Adn = gdsv(j,i,kcn)*zdaan
      Ads = gdsv(j,i,kcs)*zdaan
      Adep = dmax1( Ade, 0.0d0 )
      Adem = Ade - Adep
      Adwp = dmax1( Adw, 0.0d0 )
      Adwm = Adw - Adwp
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   dacfN(iw,ia) = Adn + 0.5d0*(Adep-Adwm)
!ik   dacfS(iw,ia) = Ads + 0.5d0*(Adwp-Adem)
      dacfN(i,ia) = Adn + 0.5d0*(Adep-Adwm)
      dacfS(i,ia) = Ads + 0.5d0*(Adwp-Adem)
!
!::d-coef Eta  at cell boundary (j,i+1/2)
!x    zetan =1.0d0/(gwtm(j,i)/etan(i,ia)+gwtp(j,i)/etan(ip,ia))
      zetan =fundav(1.0d0,gwtm(j,i),gwtp(j,i),etan(i,ia),etan(ip,ia))
!
      Ate = gdsv(j,i,kce)*zetan
      Atw = gdsv(j,i,kcw)*zetan
      Atn = gdsv(j,i,kcn)*zetan
      Ats = gdsv(j,i,kcs)*zetan
      Atep = dmax1( Ate, 0.0d0 )
      Atem = Ate - Atep
      Atwp = dmax1( Atw, 0.0d0 )
      Atwm = Atw - Atwp
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   etcfN(iw,ia) = Atn + 0.5d0*(Atep-Atwm)
!ik   etcfS(iw,ia) = Ats + 0.5d0*(Atwp-Atem)
      etcfN(i,ia) = Atn + 0.5d0*(Atep-Atwm)
      etcfS(i,ia) = Ats + 0.5d0*(Atwp-Atem)
!
!::parameter (roa,va) at E(j+1/2,i+1/2), W(j-1/2,i+1/2)
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   zrae= ama(ia)*vnag(j, i,ia)
!ik   zraw= ama(ia)*vnag(jm,i,ia)
      zvae= vvag(j, i,ia)
      zvaw= vvag(jm,i,ia)
!
!::X-interpolation (X) at O(j,i+1/2)
      zwt1 = xwtam(j,i,ia)
      zwt2 = xwtap(j,i,ia)
!
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   xsva(iw,ia) = zwt1*vva(j,i,ia)+zwt2*vva(j,ip,ia)
!ik   xsti(iw,ia) = zwt1*hia(i,ia)  +zwt2*hia(ip,ia)
!ik   xste(iw,ia) = zwt1*hea(i,ia)  +zwt2*hea(ip,ia)
      xsva(i,ia) = zwt1*vva(j,i,ia)+zwt2*vva(j,ip,ia)
      xsti(i,ia) = zwt1*hia(i,ia)  +zwt2*hia(ip,ia)
      xste(i,ia) = zwt1*hea(i,ia)  +zwt2*hea(ip,ia)
!
!::flux
      zf11 = xfdfa(j,i,ia)
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   zf21 = zf11*xsva(iw,ia)
!ik   zf31 = zf11*xsti(iw,ia)
!ik   zf41 = zf11*xste(iw,ia)
      zf21 = zf11*xsva(i,ia)
      zf31 = zf11*xsti(i,ia)
      zf41 = zf11*xste(i,ia)
      zf22 = - Ate*zvae + Atw*zvaw - Atn*vva(j,ip,ia) + Ats*vva(j,i,ia)
      zf32 = - Ate*zvae**2 + Atw*zvaw**2
     >       - Atn*vva(j,ip,ia)**2 + Ats*vva(j,i,ia)**2
! modified 4/4 lines organize local variables and include files by kamata 2021/05/31
!ik   dfl1a(iw,ia) =            zf11
!ik   dfl2a(iw,ia) =            zf21 + zf22
!ik   dfl3(iw)     = dfl3(iw) + zf31 + 0.5d0*zf32
!ik   dfl4(iw)     = dfl4(iw) + zf41
      dfl1a(i,ia) =           zf11
      dfl2a(i,ia) =           zf21 + zf22
      dfl3(i)     = dfl3(i) + zf31 + 0.5d0*zf32
      dfl4(i)     = dfl4(i) + zf41
!
 220  continue
!
!::d-coef Xi  at cell boundary (j,i+1/2)
!x    zxian = 1.0d0/(gwtm(j,i)/xian(i)+gwtp(j,i)/xian(ip))
      zxian =fundav(1.0d0,gwtm(j,i),gwtp(j,i),xian(i),xian(ip))
!
      Aie = gdsv(j,i,kce)*zxian
      Aiw = gdsv(j,i,kcw)*zxian
      Ain = gdsv(j,i,kcn)*zxian
      Ais = gdsv(j,i,kcs)*zxian
      Aiep = dmax1( Aie, 0.0d0 )
      Aiem = Aie - Aiep
      Aiwp = dmax1( Aiw, 0.0d0 )
      Aiwm = Aiw - Aiwp
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   xicfN(iw) = Ain + 0.5d0*(Aiep-Aiwm)
!ik   xicfS(iw) = Ais + 0.5d0*(Aiwp-Aiem)
      xicfN(i) = Ain + 0.5d0*(Aiep-Aiwm)
      xicfS(i) = Ais + 0.5d0*(Aiwp-Aiem)
!
!::d-coef Xe  at cell boundary (j,i+1/2)
!x    zxean = 1.0d0/(gwtm(j,i)/xean(i)+gwtp(j,i)/xean(ip))
      zxean =fundav(1.0d0,gwtm(j,i),gwtp(j,i),xean(i),xean(ip))
!
      Aee = gdsv(j,i,kce)*zxean
      Aew = gdsv(j,i,kcw)*zxean
      Aen = gdsv(j,i,kcn)*zxean
      Aes = gdsv(j,i,kcs)*zxean
      Aeep = dmax1( Aee, 0.0d0 )
      Aeem = Aee - Aeep
      Aewp = dmax1( Aew, 0.0d0 )
      Aewm = Aew - Aewp
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   xecfN(iw) = Aen + 0.5d0*(Aeep-Aewm)
!ik   xecfS(iw) = Aes + 0.5d0*(Aewp-Aeem)
      xecfN(i) = Aen + 0.5d0*(Aeep-Aewm)
      xecfS(i) = Aes + 0.5d0*(Aewp-Aeem)
!
!::parameter (Ti,Te) at E(j+1/2,i+1/2), W(j-1/2,i+1/2)
      ztie= vtig(j, i)
      ztiw= vtig(jm,i)
      ztee= vteg(j, i)
      ztew= vteg(jm,i)
!
!::flux
      zf33 = - Aie*ztie + Aiw*ztiw - Ain*vti(j,ip) + Ais*vti(j,i)
      zf43 = - Aee*ztee + Aew*ztew - Aen*vte(j,ip) + Aes*vte(j,i)
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   dfl3(iw)  = dfl3(iw) + zf33*cev
!ik   dfl4(iw)  = dfl4(iw) + zf43*cev
      dfl3(i)  = dfl3(i) + zf33*cev
      dfl4(i)  = dfl4(i) + zf43*cev
!
!::flux due to impurity
      xhte = gwtp(j,i)*vte(j,i)+gwtm(j,i)*vte(j,ip)
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   xstz(iw) = chfcv*xzflz(j,i)
!ik   zf4z = xstz(iw)*xhte
!ik   dfl4(iw) = dfl4(iw) + zf4z*cev
      xstz(i) = chfcv*xzflz(j,i)
      zf4z = xstz(i)*xhte
      dfl4(i) = dfl4(i) + zf4z*cev
 210  continue
!
!-----------------------------------------------------------------------
!::jacobian
!-----------------------------------------------------------------------
      m3 = 2*nion + 1
      m4 = 2*nion + 2
!
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   do 310 iw = 2, imax1
!ik   i   = iw
      do 310 i = 2, imax1
      ip  = i + 1
      im  = i - 1
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   iw1 = iw-1
!
      do 320 ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::G11 = -Da*grad(q1a)
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   zc11 = -dacfS(iw1,ia)
!ik   zd11 =  dacfN(iw1,ia) + dacfS(iw,ia)
!ik   ze11 = -dacfN(iw, ia)
      zc11 = -dacfS(im,ia)
      zd11 =  dacfN(im,ia) + dacfS(i,ia)
      ze11 = -dacfN(i, ia)
!
!::G21 = -va*Da*grad(q1a)
!  va*dlt_[-Da*grad(q1a)]
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   zc21 = -xsva(iw1,ia)*dacfS(iw1,ia)
!ik   zd21 =  xsva(iw1,ia)*dacfN(iw1,ia)+xsva(iw,ia)*dacfS(iw,ia)
!ik   ze21 = -xsva(iw, ia)*dacfN(iw, ia)
      zc21 = -xsva(im,ia)*dacfS(im,ia)
      zd21 =  xsva(im,ia)*dacfN(im,ia)+xsva(i,ia)*dacfS(i,ia)
      ze21 = -xsva(i, ia)*dacfN(i, ia)
!  dlt_[va]*flx = 0.0
      xc21 = 0.0d0; xd21 = 0.0d0; xe21 = 0.0d0   ! %%% 2002/11/11
!
!::G22 = -eta*grad(va)
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   zc22 = -etcfS(iw1,ia)
!ik   zd22 =  etcfN(iw1,ia)+etcfS(iw,ia)
!ik   ze22 = -etcfN(iw, ia)
      zc22 = -etcfS(im,ia)
      zd22 =  etcfN(im,ia)+etcfS(i,ia)
      ze22 = -etcfN(i, ia)
!
!::G31 = -sum[ (ea+pa)/roa*Da*grad(q1a) ]
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   zc31 = -xsti(iw1,ia)*dacfS(iw1,ia)
!ik   zd31 =  xsti(iw1,ia)*dacfN(iw1,ia)+xsti(iw,ia)*dacfS(iw,ia)
!ik   ze31 = -xsti(iw, ia)*dacfN(iw, ia)
      zc31 = -xsti(im,ia)*dacfS(im,ia)
      zd31 =  xsti(im,ia)*dacfN(im,ia)+xsti(i,ia)*dacfS(i,ia)
      ze31 = -xsti(i, ia)*dacfN(i, ia)
!
!::G32 = -sum[ eta*grad(1/2*va**2) ]
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   zc32 = -etcfS(iw1,ia)
!ik   zd32 =  etcfN(iw1,ia)+etcfS(iw,ia)
!ik   ze32 = -etcfN(iw, ia)
      zc32 = -etcfS(im,ia)
      zd32 =  etcfN(im,ia)+etcfS(i,ia)
      ze32 = -etcfN(i, ia)
!
!::G33 = -Ki*grad(Ti)
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   zc33 = -xicfS(iw1)
!ik   zd33 =  xicfN(iw1)+xicfS(iw)
!ik   ze33 = -xicfN(iw )
      zc33 = -xicfS(im)
      zd33 =  xicfN(im)+xicfS(i)
      ze33 = -xicfN(i )
!
!::G41 = -sum[ 5/2*Te*za/ma*Da*grad(q1a) ]
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   zc41 = -xste(iw1,ia)*dacfS(iw1,ia)
!ik   zd41 =  xste(iw1,ia)*dacfN(iw1,ia)+xste(iw,ia)*dacfS(iw,ia)
!ik   ze41 = -xste(iw, ia)*dacfS(iw, ia)
      zc41 = -xste(im,ia)*dacfS(im,ia)
      zd41 =  xste(im,ia)*dacfN(im,ia)+xste(i,ia)*dacfS(i,ia)
      ze41 = -xste(i, ia)*dacfS(i, ia)
!
!::G4z = 5/2*fleZ*(gwtp*Te(i)+gwtm*Te(ip))
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   zc4z = xstz(iw1)*gwtp(j,im)
!ik   zd4z = xstz(iw1)*gwtm(j,im) - xstz(iw)*gwtp(j,i)
!ik   ze4z = xstz(iw)*gwtm(j,i)
      zc4z = xstz(im)*gwtp(j,im)
      zd4z = xstz(im)*gwtm(j,im) - xstz(i)*gwtp(j,i)
      ze4z = xstz(i)*gwtm(j,i)
!
!::G43 = -Ke*grad(Te)
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   zc43 = -xecfS(iw1)
!ik   zd43 =  xecfN(iw1)+xecfS(iw)
!ik   ze43 = -xecfN(iw)
      zc43 = -xecfS(im)
      zd43 =  xecfN(im)+xecfS(i)
      ze43 = -xecfN(i)
!
!::equation [q1a = ma*na]
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m1a,m1a,iw) = cc(m1a,m1a,iw) + zc11
!ik   dd(m1a,m1a,iw) = dd(m1a,m1a,iw) + zd11
!ik   ee(m1a,m1a,iw) = ee(m1a,m1a,iw) + ze11
      cc(m1a,m1a,i) = cc(m1a,m1a,i) + zc11
      dd(m1a,m1a,i) = dd(m1a,m1a,i) + zd11
      ee(m1a,m1a,i) = ee(m1a,m1a,i) + ze11
!
!::equation [q2a = ma*na*vpa]
! modified 6/6 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m2a,m1a,iw) = cc(m2a,m1a,iw) + zc21 + (xc21+zc22)*r1av(im,ia)
!ik   dd(m2a,m1a,iw) = dd(m2a,m1a,iw) + zd21 + (xd21+zd22)*r1av(i, ia)
!ik   ee(m2a,m1a,iw) = ee(m2a,m1a,iw) + ze21 + (xe21+ze22)*r1av(ip,ia)
!ik   cc(m2a,m2a,iw) = cc(m2a,m2a,iw)        + (xc21+zc22)*r2av(im,ia)
!ik   dd(m2a,m2a,iw) = dd(m2a,m2a,iw)        + (xd21+zd22)*r2av(i, ia)
!ik   ee(m2a,m2a,iw) = ee(m2a,m2a,iw)        + (xe21+ze22)*r2av(ip,ia)
      cc(m2a,m1a,i) = cc(m2a,m1a,i) + zc21 + (xc21+zc22)*r1av(im,ia)
      dd(m2a,m1a,i) = dd(m2a,m1a,i) + zd21 + (xd21+zd22)*r1av(i, ia)
      ee(m2a,m1a,i) = ee(m2a,m1a,i) + ze21 + (xe21+ze22)*r1av(ip,ia)
      cc(m2a,m2a,i) = cc(m2a,m2a,i)        + (xc21+zc22)*r2av(im,ia)
      dd(m2a,m2a,i) = dd(m2a,m2a,i)        + (xd21+zd22)*r2av(i, ia)
      ee(m2a,m2a,i) = ee(m2a,m2a,i)        + (xe21+ze22)*r2av(ip,ia)
!
!::equation [q3 = energy_I ]
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m3,m1a,iw) = cc(m3,m1a,iw)
      cc(m3,m1a,i) = cc(m3,m1a,i)
     >                       + zc31 + zc32*r1aw(im,ia)+zc33*r1ai(im,ia)
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   dd(m3,m1a,iw) = dd(m3,m1a,iw)
      dd(m3,m1a,i) = dd(m3,m1a,i)
     >                       + zd31 + zd32*r1aw(i, ia)+zd33*r1ai(i, ia)
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   ee(m3,m1a,iw) = ee(m3,m1a,iw)
      ee(m3,m1a,i) = ee(m3,m1a,i)
     >                       + ze31 + ze32*r1aw(ip,ia)+ze33*r1ai(ip,ia)
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m3,m2a,iw) = cc(m3,m2a,iw) + zc32*r2aw(im,ia)+zc33*r2ai(im,ia)
!ik   dd(m3,m2a,iw) = dd(m3,m2a,iw) + zd32*r2aw(i, ia)+zd33*r2ai(i, ia)
!ik   ee(m3,m2a,iw) = ee(m3,m2a,iw) + ze32*r2aw(ip,ia)+ze33*r2ai(ip,ia)
      cc(m3,m2a,i) = cc(m3,m2a,i) + zc32*r2aw(im,ia)+zc33*r2ai(im,ia)
      dd(m3,m2a,i) = dd(m3,m2a,i) + zd32*r2aw(i, ia)+zd33*r2ai(i, ia)
      ee(m3,m2a,i) = ee(m3,m2a,i) + ze32*r2aw(ip,ia)+ze33*r2ai(ip,ia)
!
!::equation [q4 = energy_e ]
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m4,m1a,iw) = cc(m4,m1a,iw) + zc41 + zc43*r1ae(im,ia)
!ik   dd(m4,m1a,iw) = dd(m4,m1a,iw) + zd41 + zd43*r1ae(i, ia)
!ik   ee(m4,m1a,iw) = ee(m4,m1a,iw) + ze41 + ze43*r1ae(ip,ia)
      cc(m4,m1a,i) = cc(m4,m1a,i) + zc41 + zc43*r1ae(im,ia)
      dd(m4,m1a,i) = dd(m4,m1a,i) + zd41 + zd43*r1ae(i, ia)
      ee(m4,m1a,i) = ee(m4,m1a,i) + ze41 + ze43*r1ae(ip,ia)
!
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m4,m1a,iw) = cc(m4,m1a,iw) + zc4z*r1ae(im,ia)
!ik   dd(m4,m1a,iw) = dd(m4,m1a,iw) + zd4z*r1ae(i, ia)
!ik   ee(m4,m1a,iw) = ee(m4,m1a,iw) + ze4z*r1ae(ip,ia)
      cc(m4,m1a,i) = cc(m4,m1a,i) + zc4z*r1ae(im,ia)
      dd(m4,m1a,i) = dd(m4,m1a,i) + zd4z*r1ae(i, ia)
      ee(m4,m1a,i) = ee(m4,m1a,i) + ze4z*r1ae(ip,ia)
!
!::ff
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   ff(m1a,iw) = ff(m1a,iw) - dfl1a(iw,ia) + dfl1a(iw1,ia)
!ik   ff(m2a,iw) = ff(m2a,iw) - dfl2a(iw,ia) + dfl2a(iw1,ia)
      ff(m1a,i) = ff(m1a,i) - dfl1a(i,ia) + dfl1a(im,ia)
      ff(m2a,i) = ff(m2a,i) - dfl2a(i,ia) + dfl2a(im,ia)
 320  continue  ! loop(ia)
!
! modified 6/6 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m3,m3, iw) = cc(m3,m3, iw) + zc33*r3i(im)
!ik   dd(m3,m3, iw) = dd(m3,m3, iw) + zd33*r3i(i)
!ik   ee(m3,m3, iw) = ee(m3,m3, iw) + ze33*r3i(ip)
!ik   cc(m4,m4, iw) = cc(m4,m4, iw) + zc43*r4e(im)
!ik   dd(m4,m4, iw) = dd(m4,m4, iw) + zd43*r4e(i)
!ik   ee(m4,m4, iw) = ee(m4,m4, iw) + ze43*r4e(ip)
      cc(m3,m3, i) = cc(m3,m3, i) + zc33*r3i(im)
      dd(m3,m3, i) = dd(m3,m3, i) + zd33*r3i(i)
      ee(m3,m3, i) = ee(m3,m3, i) + ze33*r3i(ip)
      cc(m4,m4, i) = cc(m4,m4, i) + zc43*r4e(im)
      dd(m4,m4, i) = dd(m4,m4, i) + zd43*r4e(i)
      ee(m4,m4, i) = ee(m4,m4, i) + ze43*r4e(ip)
!
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m4,m4, iw) = cc(m4,m4, iw) + zc4z*r4e(im)
!ik   dd(m4,m4, iw) = dd(m4,m4, iw) + zd4z*r4e(i)
!ik   ee(m4,m4, iw) = ee(m4,m4, iw) + ze4z*r4e(ip)
      cc(m4,m4, i) = cc(m4,m4, i) + zc4z*r4e(im)
      dd(m4,m4, i) = dd(m4,m4, i) + zd4z*r4e(i)
      ee(m4,m4, i) = ee(m4,m4, i) + ze4z*r4e(ip)
!
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   ff(m3,iw) = ff(m3,iw) - dfl3(iw) + dfl3(iw1)
!ik   ff(m4,iw) = ff(m4,iw) - dfl4(iw) + dfl4(iw1)
      ff(m3,i) = ff(m3,i) - dfl3(i) + dfl3(im)
      ff(m4,i) = ff(m4,i) - dfl4(i) + dfl4(im)
!
 310  continue  ! loop(jw)
!
      return
      end
