!***********************************************************************
      subroutine pybcon(j)
!***********************************************************************
!
!   Yacobian at boundary
!
!   boundary condition at divertor plate
!      i=1  (sol wall) / i=imax (prv wall)
!
!  (1) na   e^psi*grad(na) = 1/alng*na(i=ib)
!  (2) va   va = 0.0                         when lbcsw/lbcpw = 1
!  (3) Ti   e^psi*grad(Ti) = 1/alng*Ti(i=ib)  or  Ti_pw = bwpti
!  (4) Te   e^psi*grad(Te) = 1/alng*Te(i=ib)  or  Te_sw = bwpte
!
!
!      alng : characteristics length
!               ITER  0.03;   SC  0.03;  C-mod  0.01
!
!-----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : ama, aza, bwpfn, bwpni, bwpte, bwpti, bwsfn
     >    , bwsni, bwste, bwsti, c13, c23, cc, dd, ee, ff, gwpni, gwpte
     >    , gwpti, gwsni, gwste, gwsti, lbcpw, lbcsw, lfbbc, nion, nlp
     >    , q1a, q2a, q3, q4, qb1a, qb2a, qb3, qb4, vna, vnag, vne, vni
     >    , vte, vteg, vti, vtig, vva
      use cplmet, only : gare, gdsv, icmax, jcxp1, jcxp2, jnxm, kce, kcn
     >    , kcs, kcw, icmin
      use csize,  only : ndy
      use csonic, only : itim
      use cunit,  only : n6
      implicit none

!::argument
      integer, intent(in) :: j
!
!::local variables
      integer  imax, iw, iwm, iwp, jm
      integer m1a, m2a, m3, m4, ia, kbc
      real*8  gnaw1, gtiw1, gtew1, Abe, Abw, Abn, Abs, Abep
      real*8  Ab1, Ab2, Abem, Abwp, Abwm, zgw, zd1, ze1, zf1
      real*8  znae, znaw, znan, znas
      real*8  ztie, ztiw, ztin, ztis, ztee, ztew, zten, ztes
      real*8  zd3, ze3, zf3, zd4, ze4, zf4
      real*8  gnaw2, gtiw2, gtew2
      real*8  vniw1, vtiw1, vtew1, vniw2, vtiw2, vtew2
      real*8  zc1, zc3, zc4
      integer lp; data lp/0/; save lp
!::KH variables for feedback boundary condition 20100426
      real*8  wrkb,wrkc(1:3),wrkd(1:3),wrke, wrkk(1:3)
      real*8,save :: lamda(1:3) = -0.01d0
      real*8  rst(ndy), ren(ndy), rhf(ndy)
      integer n,k,kmax,k2
!
      kmax=5
!
      if( lp.eq.0 ) then
        lp = 1
        write(n6,'(/2x,"***  pybcon  ***  gradient ")')
        write(n6,'(4x,"va = 0.0d0   10/14/2002")')
      endif
!
!::index
      imax  = icmax(j)
      m3   = 2*nion + 1
      m4   = 2*nion + 2
!
!-----------------------------------------------------------------------
!::boundary condition i = 1 (sol-wall)
!-----------------------------------------------------------------------
!::input
      gnaw1 = gwsni; gtiw1 = gwsti; gtew1 = gwste
      vniw1 = bwsni; vtiw1 = bwsti; vtew1 = bwste
!
      iw  = 1
      iwp = iw + 1
      jm  = jnxm(j,iw)
      if( j.eq.1 ) jm = j
!
      do ia = 1, nion
        m1a = 2*ia - 1
        m2a = 2*ia
!
!::metric
        Abe = 0.0d0 ! Abe = gdsv(j,iw,kce) ?
        Abw = 0.0d0 ! Abw = gdsv(j,iw,kcw) ?
!-----
        Abn = gdsv(j,iw,kcn)
        Abs = gdsv(j,iw,kcs)
        Abep = dmax1( Abe, 0.0d0 )
        Abem = Abe - Abep
        Abwp = dmax1( Abw, 0.0d0 )
        Abwm = Abw - Abwp
        Ab2  = Abn + 0.5d0*(Abep-Abwm)
        Ab1  = Abs + 0.5d0*(Abwp-Abem)
!
!::density (fixed value)
        if( lbcsw.eq.2 ) then
          zd1 = 1.0d0
          zf1 = (-vna(j,iw,ia) + vniw1*bwsfn(ia))*ama(ia)
          dd(m1a,m1a,iw) = zd1
          ff(m1a,iw)     = zf1
!::density (gradient)
        else
          zgw  =  gare(j,iw)/gnaw1
          zd1  =  Ab1 + zgw
          ze1  = -Ab2
          znae = vnag(j, iw,ia)
          znaw = vnag(jm,iw,ia)
          znan = vna(j,iwp,ia)
          znas = vna(j,iw,ia)
          zf1  =(Abe*znae - Abw*znaw + Abn*znan - Abs*znas
     >      - zgw*vna(j,iw,ia))*ama(ia)
          dd(m1a,m1a,iw) = zd1
          ee(m1a,m1a,iw) = ze1
          ff(m1a,iw)     = zf1
        endif
!
!::velociy (va=0.0)
        dd(m2a,m2a,iw) = 1.0d0
        ff(m2a,iw) = 0.0d0 - q2a(j,iw,ia)
      enddo
!
!::ion energy (fixed value)
      if( lbcsw.ge.1 ) then
        zd3  =  1.0d0
        zf3  = ( -vti(j,iw) + vtiw1 )*cev
        do ia = 1, nion
          m1a = 2*ia - 1
          m2a = 2*ia
          dd(m3,m1a,iw) = zd3*(c13*vva(j,iw,ia)**2
     >                      -vti(j,iw)*cev/ama(ia))/vni(j,iw)
          dd(m3,m2a,iw) =-zd3*c23*vva(j,iw,ia)/vni(j,iw)
        enddo
        dd(m3,m3, iw) = zd3*c23/vni(j,iw)
        ff(m3,iw)     = zf3
!::ion energy (gradient)
      else
        zgw =  gare(j,iw)/gtiw1
        zd3 =  Ab1 + zgw
        ze3 = -Ab2
        ztie = vtig(j, iw)
        ztiw = vtig(jm,iw)
        ztin = vti(j,iwp)
        ztis = vti(j,iw)
        zf3  = ( Abe*ztie - Abw*ztiw + Abn*ztin - Abs*ztis
     >         - zgw*vti(j,iw) )*cev
        do ia = 1, nion
          m1a = 2*ia - 1
          m2a = 2*ia
          dd(m3,m1a,iw) = zd3*(c13*vva(j,iw,ia)**2
     >                       -vti(j,iw)*cev/ama(ia))/vni(j,iw)
          dd(m3,m2a,iw) =-zd3*c23*vva(j,iw,ia)/vni(j,iw)
          ee(m3,m1a,iw) = ze3*(c13*vva(j,iwp,ia)**2
     >                       -vti(j,iwp)*cev/ama(ia))/vni(j,iwp)
          ee(m3,m2a,iw) =-ze3*c23*vva(j,iwp,ia)/vni(j,iwp)
        enddo
        dd(m3,m3,iw)  = zd3*c23/vni(j,iw)
        ee(m3,m3,iw)  = ze3*c23/vni(j,iwp)
        ff(m3,iw)     = zf3
      endif
!
!::electron energy (fixed value)
      if( lbcsw.ge.1 ) then
        zd4  =  1.0d0
        zf4  = ( -vte(j,iw) + vtew1 )*cev
        do ia = 1, nion
          m1a = 2*ia - 1
          m2a = 2*ia
          dd(m4,m1a,iw) =-zd4*aza(ia)/ama(ia)*vte(j,iw)*cev/vne(j,iw)
        enddo
        dd(m4,m4, iw) = zd4*c23/vne(j,iw)
        ff(m4,iw)     = zf4
!
!::electron energy (gradient)
      else
        zgw =  gare(j,iw)/gtew1
        zd4 =  Ab1 + zgw
        ze4 = -Ab2
        ztee = vteg(j, iw)
        ztew = vteg(jm,iw)
        zten = vte(j,iwp)
        ztes = vte(j,iw)
        zf4  = ( Abe*ztee - Abw*ztew + Abn*zten - Abs*ztes
     >        - zgw*vte(j,iw) )*cev
        do ia = 1, nion
          m1a = 2*ia - 1
          dd(m4,m1a,iw) =-zd4*aza(ia)/ama(ia)*vte(j,iw)*cev/vne(j,iw)
          ee(m4,m1a,iw) =-ze4*aza(ia)/ama(ia)*vte(j,iwp)*cev/vne(j,iwp)
        enddo
        dd(m4,m4,iw) = zd4*c23/vne(j,iw)
        ee(m4,m4,iw) = ze4*c23/vne(j,iwp)
        ff(m4,iw) = zf4
      endif
!
!-----------------------------------------------------------------------
!::boundary condition i = imax (prv-wall)
!-----------------------------------------------------------------------
      if( j.gt.jcxp1 .and. j.lt.jcxp2 ) goto 300
      gnaw2 = -gwpni; gtiw2 = -gwpti; gtew2 = -gwpte
      vniw2 =  bwpni; vtiw2 =  bwpti; vtew2 =  bwpte
!
      iw  = imax
      iwm = iw - 1
      jm  = jnxm(j,iw)
      if( j.eq.1 )     jm = j
!
      do ia = 1, nion
        m1a = 2*ia - 1
        m2a = 2*ia
!
!::metric
        Abe = 0.0d0 ! gdsv(j,iwm,kce)
        Abw = 0.0d0 ! gdsv(j,iwm,kcw)
!-----
        Abn = gdsv(j,iwm,kcn)
        Abs = gdsv(j,iwm,kcs)
        Abep = dmax1( Abe, 0.0d0 )
        Abem = Abe - Abep
        Abwp = dmax1( Abw, 0.0d0 )
        Abwm = Abw - Abwp
        Ab1  = Abn + 0.5d0*(Abep-Abwm)
        Ab2  = Abs + 0.5d0*(Abwp-Abem)
!
!KH 20100426 in sol2d.v2
!KH 201021111 add to sonicV3
!KH feedback boundary condition
        if( lfbbc .ne. 0)then
          if(mod(itim,10).eq.1 .and. nlp.eq.1) then
            call codrad(j,rst,ren,rhf,n,ndy,icmin(j),icmax(j))
            wrkb=0.0d0
            wrkc=0.0d0
            wrkd=0.0d0
            wrke=0.0d0
            wrkk=0.0d0
            do k=1,kmax
              wrkk(1)=log(vna(j,imax-k,ia))
              wrkk(2)=log(vti(j,imax-k))
              wrkk(3)=log(vte(j,imax-k))
              wrkb=wrkb+rhf(imax-k)**2
              do k2=1,3
                wrkc(k2)=wrkc(k2)+wrkk(k2)
                wrkd(k2)=wrkd(k2)+dabs(rhf(imax-k))*wrkk(k2)
              enddo
              wrke=wrke+dabs(rhf(imax-k))
            enddo
            do k2=1,3
              lamda(k2)=
     >          (5.0d0*wrkb-wrke**2)/(5.0d0*wrkd(k2)-wrkc(k2)*wrke)
              if(lamda(k2) .ge. 0.0d0) lamda(k2) = -1.0d0
              lamda(k2)=min(lamda(k2), -1.0d-3)
              lamda(k2)=max(lamda(k2), -0.5d0)
            enddo
            if(mod(itim,1000).eq.1)then
              if(j.eq.2) write(995,*)
              write(995,'(1h 4i6,3f12.6)')itim,nlp,j,imax,lamda
            endif
          endif
          gnaw2=lamda(1)
          if(lfbbc==1)then
            gtiw2=lamda(2)
            gtew2=lamda(3)
          endif
        endif
!
!::density (fixed value)
        if( lbcpw.eq.2 ) then
          zd1 = 1.0d0
          zf1 = (-vna(j,iw,ia) + vniw2*bwpfn(ia))*ama(ia)
          dd(m1a,m1a,iw) = zd1
          ff(m1a,iw)     = zf1
!
!::density (gradient)
        else
          zgw  =  gare(j,iwm)/gnaw2
          zd1  = -Ab1 + zgw
          zc1  =  Ab2
          znae = vnag(j, iwm,ia)
          znaw = vnag(jm,iwm,ia)
          znan = vna(j,iw,ia)
          znas = vna(j,iwm,ia)
          zf1   =(Abe*znae - Abw*znaw + Abn*znan - Abs*znas
     >      - zgw*vna(j,iw,ia))*ama(ia)
          dd(m1a,m1a,iw) = zd1
          cc(m1a,m1a,iw) = zc1
          ff(m1a,iw)     = zf1
        endif
!
!::velociy (va=0.0)
        dd(m2a,m2a,iw) = 1.0d0
        ff(m2a,iw) = 0.0d0 - q2a(j,iw,ia)
      enddo
!
!::ion energy (fixed value)
      if( lbcpw.ge.1 ) then
        zd3  =  1.0d0
        zf3  = ( -vti(j,iw) + vtiw2 )*cev
        do ia = 1, nion
          m1a = 2*ia - 1
          m2a = 2*ia
          dd(m3,m1a,iw) = zd3*(c13*vva(j,iw,ia)**2
     >                      -vti(j,iw)*cev/ama(ia))/vni(j,iw)
          dd(m3,m2a,iw) =-zd3*c23*vva(j,iw,ia)/vni(j,iw)
        enddo
        dd(m3,m3, iw) = zd3*c23/vni(j,iw)
        ff(m3,iw)     = zf3
!
!::ion energy (gradient)
      else
        zgw  =  gare(j,iwm)/gtiw2
        zd3  = -Ab1 + zgw
        zc3  =  Ab2
        ztie = vtig(j,iwm)
        ztiw = vtig(jm,iwm)
        ztin = vti(j,iw)
        ztis = vti(j,iwm)
        zf3   =(Abe*ztie - Abw*ztiw + Abn*ztin - Abs*ztis
     >      - zgw*vti(j,iw) )*cev
        do ia = 1, nion
          m1a = 2*ia - 1
          m2a = 2*ia
          dd(m3,m1a,iw) = zd3*(c13*vva(j,iw,ia)**2
     >                      -vti(j,iw)*cev/ama(ia))/vni(j,iw)
          dd(m3,m2a,iw) =-zd3*c23*vva(j,iw,ia)/vni(j,iw)
          cc(m3,m1a,iw) = zc3*(c13*vva(j,iwm,ia)**2
     >                      -vti(j,iwm)*cev/ama(ia))/vni(j,iwm)
          cc(m3,m2a,iw) =-zc3*c23*vva(j,iwm,ia)/vni(j,iwm)
        enddo
        dd(m3,m3, iw) = zd3*c23/vni(j,iw)
        cc(m3,m3, iw) = zc3*c23/vni(j,iwm)
        ff(m3,iw)     = zf3
      endif
!
!::electron energy (fixed value)
      if( lbcpw.ge.1 ) then
        zd4  =  1.0d0
        zf4  = ( -vte(j,iw) + vtew2 )*cev
        do ia = 1, nion
          m1a = 2*ia - 1
          m2a = 2*ia
          dd(m4,m1a,iw) =-zd4*aza(ia)/ama(ia)*vte(j,iw)*cev/vne(j,iw)
        enddo
        dd(m4,m4, iw) = zd4*c23/vne(j,iw)
        ff(m4,iw)     = zf4
!
!::electron energy (gradient)
      else
        zgw  =  gare(j,iwm)/gtew2
        zd4  = -Ab1 + zgw
        zc4  =  Ab2
        ztee = vteg(j, iwm)
        ztew = vteg(jm,iwm)
        zten = vte(j,iw)
        ztes = vte(j,iwm)
        zf4   =(Abe*ztee - Abw*ztew + Abn*zten - Abs*ztes
     >         - zgw*vte(j,iw) )*cev
        do ia = 1, nion
          m1a = 2*ia - 1
          m2a = 2*ia
          dd(m4,m1a,iw) =-zd4*aza(ia)/ama(ia)*vte(j,iw)*cev/vne(j,iw)
          cc(m4,m1a,iw) =-zc4*aza(ia)/ama(ia)*vte(j,iwm)*cev/vne(j,iwm)
        enddo
        dd(m4,m4, iw) = zd4*c23/vne(j,iw)
        cc(m4,m4, iw) = zc4*c23/vne(j,iwm)
        ff(m4,iw)     = zf4
      endif
!
      return
!
!-----------------------------------------------------------------------
!::boundary condition i = imax (edge plasma)
!-----------------------------------------------------------------------
 300  continue
      kbc = kcn
      iw  = imax
      do ia = 1, nion
        m1a = 2*ia - 1
        m2a = 2*ia
        dd(m1a,m1a,iw) = 1.0d0
        dd(m2a,m2a,iw) = 1.0d0
        ff(m1a,iw) = qb1a(j,ia,kbc)-q1a(j,iw,ia)
        ff(m2a,iw) = qb2a(j,ia,kbc)-q2a(j,iw,ia)
      enddo
      dd(m3,m3,  iw) = 1.0d0
      dd(m4,m4,  iw) = 1.0d0
      ff(m3,iw) = qb3(j,kbc)-q3(j,iw)
      ff(m4,iw) = qb4(j,kbc)-q4(j,iw)
!
      return
      end
