!**********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   subroutine pxpres(nt)
      subroutine pxpres(it)
!**********************************************************************
!
!     dq/dt + div(F) + Hp = S
!
!     dlt( dq/dt + divF + Hp - S ) = -( dq/dt + div(F) + Hp - S )^l
!
!      Cj*q(j-1) + Dj*q(j) + Ej*q(j+1) = fj
!
!      H1a = 0
!      H2a = - zdvb*pa                  [zf21]  qvl_pi(ia)
!            + (za*na/ne)*V*b*grad(pe)  [zf22]  qvl_pe(m2a)
!      H3  =  ve*V*b*grad(pe)           [zf3]   qvl_pe(m3)
!      H4  =  - H3                      [zf4]   qvl_pe(m4)
!
!   note.  b = unit vector along B
!          V = volume of cell (j,i)
!          zdvb = Integ[div(b)] = [S^sita*b]j+1/2-[S^sita*b]j-1/2
!
!        The term of va^(Da)*grad(pe)
!            will be included for comparison with B2-EIRENE code.
!
!     New pxpres   div ==> grd  2002/10/20
!         b*grad(pe)j =1/(dlp)j*(p(j+1/2)-p(j-1/2))
!         interpolation  pe(j+1/2) = Hm(j+1/2)*pe(j)+Hp(j+1/2)*pe(j+1)
!                                                        2006/02/02
!-----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : ama, aza, c13, c23, cc, crod, dd, ee, ff, nion
     >    , ss, vna, vne, vni, vte, vti, vva, vve
      use cplmet, only : hdxm, hdxp, hpit, hvol, hvsb, hwtm, hwtp, icel
     >    , jcel, jtmax
      use cplqcn, only : qvl_pe, qvl_pi
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nt
      integer, intent(in) :: it
!
!::local variables
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, i, jmax, jmax1, jw, j, jp, jm
!ik   integer  ia, ib, m1a, m2a, m3, m4, m1b, m2b
      integer  i, jmax, jmax1, jw, j, jp, jm
      integer  ia, ib, m2a, m3, m4, m1b, m2b
      real*8   zfna, zfpa, zdvb, zdlb, zpe, zpep, zpem
!-----
!xx   real*8   zgdm, zgdp
!-----02.02
      real*8   zhme, zhpe, zhmw, zhpw, zgdc, zgdd, zgde
!-----
      real*8   zr1b, zr2b, zr3, zf21
      real*8   zc22, zd22, ze22, zf22, zc3, zd3, ze3, zf3
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8   zgpe, zs3, zve
      real*8   zve
!
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   it = nt
      i  = icel(1,it)
      jmax  = jtmax(it)
      jmax1 = jmax - 1
!
      m3 = 2*nion + 1
      m4 = 2*nion + 2
!
      do ia = 1, nion
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   m1a = ia*2 - 1
      m2a = ia*2
!
      do jw = 2, jmax1
      j   = jcel(jw,it)
      jp  = jcel(jw+1,it)
      jm  = jcel(jw-1,it)
!
      zfna = vna(j,i,ia)/vni(j,i)
      zfpa = aza(ia)*vna(j,i,ia)/vne(j,i)
      zdvb = hvsb(j,i)-hvsb(jm,i)
      zdlb = hdxp(jm,i)/hpit(jm,i) + hdxm(j,i)/hpit(j,i)
!-----
!xx   zgdm = hwtp(jm,i)/zdlb
!xx   zgdp = hwtm(j,i)/zdlb
!-----02.02
      zhme = hwtp(j,i)        ! H-(j+1/2)
      zhpe = hwtm(j,i)        ! H+(j+1/2)
      zhmw = hwtp(jm,i)       ! H-(j-1/2)
      zhpw = hwtm(jm,i)       ! H+(j-1/2)
      zgdc = -zhmw/zdlb
      zgdd = (zhme-zhpw)/zdlb
      zgde = zhpe/zdlb
!-----
      zpe  = vne(j,i)*vte(j,i)*cev
      zpep = vne(jp,i)*vte(jp,i)*cev
      zpem = vne(jm,i)*vte(jm,i)*cev
!
      do ib = 1, nion
      m1b  = 2*ib - 1
      m2b  = 2*ib
!
!::H21 = -zdvb*pa   ! 2009/01/16
!xx   if( it.ge.itpvs .and. it.le.itpve ) zdvb = 0.0d0
      zr1b = zfna*(c13*vva(j,i,ib)**2-vti(j,i)*cev/ama(ib))
     >              +crod(ia,ib)*vti(j,i)*cev/ama(ia)
      zr2b = -c23*zfna*vva(j,i,ib)
      ss(m2a,m1b,jw) = ss(m2a,m1b,jw) - zdvb*zr1b
      ss(m2a,m2b,jw) = ss(m2a,m2b,jw) - zdvb*zr2b
      enddo  !  loop(ib)
      zr3 = c23*zfna
      ss(m2a,m3,jw) = ss(m2a,m3,jw) - zdvb*zr3
      zf21 = -zdvb*vna(j,i,ia)*vti(j,i)*cev
!
!::H22 = (za*na/ne)*V*b*grad(pe)
!-----
!xx      zc22 = -zfpa*hvol(j,i)*zgdm*c23
!xx      zd22 =  zfpa*hvol(j,i)*(zgdm-zgdp)*c23
!xx      ze22 =  zfpa*hvol(j,i)*zgdp*c23
!xx      zf22 = +zfpa*hvol(j,i)*(-zgdm*zpem+(zgdm-zgdp)*zpe+zgdp*zpep)
!-----02.02
      zc22 =  zfpa*hvol(j,i)*zgdc*c23
      zd22 =  zfpa*hvol(j,i)*zgdd*c23
      ze22 =  zfpa*hvol(j,i)*zgde*c23
      zf22 = +zfpa*hvol(j,i)*(zgdc*zpem+zgdd*zpe+zgde*zpep)
!-----
!
      cc(m2a,m4,jw) = cc(m2a,m4,jw) + zc22
      dd(m2a,m4,jw) = dd(m2a,m4,jw) + zd22
      ee(m2a,m4,jw) = ee(m2a,m4,jw) + ze22
!
      ff(m2a,jw) = ff(m2a,jw) - zf21 - zf22
!
!::qcon
      qvl_pi(j,i,ia)  = zf21
      qvl_pe(j,i,m2a) = zf22
!
      enddo    !  loop(jw)
      enddo    !  loop(ia)
!
!::H3 = ve*V*b*grad(pe)
      do jw = 2, jmax1
      j   = jcel(jw,it)
      jp  = jcel(jw+1,it)
      jm  = jcel(jw-1,it)
      zdlb = hdxp(jm,i)/hpit(jm,i) + hdxm(j,i)/hpit(j,i)
!-----
!xx   zgdm = hwtp(jm,i)/zdlb
!xx   zgdp = hwtm(j,i)/zdlb
!-----02.02
      zhme = hwtp(j,i)        ! H-(j+1/2)
      zhpe = hwtm(j,i)        ! H+(j+1/2)
      zhmw = hwtp(jm,i)       ! H-(j-1/2)
      zhpw = hwtm(jm,i)       ! H+(j-1/2)
      zgdc = -zhmw/zdlb
      zgdd = (zhme-zhpw)/zdlb
      zgde = zhpe/zdlb
!-----
      zpe  = vne(j,i)*vte(j,i)*cev
      zpep = vne(jp,i)*vte(jp,i)*cev
      zpem = vne(jm,i)*vte(jm,i)*cev
      zve  = vve(j,i)
!
!-----
!xx   zc3 = -zve*hvol(j,i)*zgdm*c23
!xx   zd3 =  zve*hvol(j,i)*(zgdm-zgdp)*c23
!xx   ze3 =  zve*hvol(j,i)*zgdp*c23
!xx   zf3 = +zve*hvol(j,i)*(-zgdm*zpem+(zgdm-zgdp)*zpe+zgdp*zpep)
!-----02.02
      zc3 =  zve*hvol(j,i)*zgdc*c23
      zd3 =  zve*hvol(j,i)*zgdd*c23
      ze3 =  zve*hvol(j,i)*zgde*c23
      zf3 = +zve*hvol(j,i)*(zgdc*zpem+zgdd*zpe+zgde*zpep)
!-----
!
      cc(m3,m4,jw) = cc(m3,m4,jw) + zc3
      dd(m3,m4,jw) = dd(m3,m4,jw) + zd3
      ee(m3,m4,jw) = ee(m3,m4,jw) + ze3
      ff(m3,jw)    = ff(m3,jw)    - zf3
!
!::H4 = -H3
      cc(m4,m4,jw) = cc(m4,m4,jw) - zc3
      dd(m4,m4,jw) = dd(m4,m4,jw) - zd3
      ee(m4,m4,jw) = ee(m4,m4,jw) - ze3
      ff(m4,jw)    = ff(m4,jw)    + zf3
!
!::qcon
      qvl_pe(j,i,m3) =  zf3
      qvl_pe(j,i,m4) = -zf3
!
      enddo
!
!-----------------------------------------------------------------------
!::H3 = -H4 = sum[za/ma*(Te-<pe>/ne)*(-Da)*grad(roa)]  new-term
!-----------------------------------------------------------------------
      return
      end
