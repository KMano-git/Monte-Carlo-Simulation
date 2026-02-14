!**********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   subroutine pyvpin(nj)
      subroutine pyvpin(j)
!**********************************************************************
!
!        G1a = roa*vp
!        G2a = v//a*roa*vp
!        G3  = sum[(ea+Pa)/roa*roa*vp] = (Q3+Pi)*vp
!        G4  = (ee+Pe)/roe*roe*vp      = 5/3*Q4*vp
!
!          Note   vp = vpinch for all sepcies  vp > 0 inward velocity
!
!                                                        2007/03/24
!----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : c13, c23, c53, cc, dd, dfl1a, dfl2a, dfl3, dfl4
     >    , ee, ff, nion, q1a, q2a, q3, q4, r1ai, r2ai, r3i, vni, vti
     >    , vva
      use cplmet, only : gare, gwtm, gwtp, icmax
      use cplvpn, only : vdvp
      use csize,  only : ndxy
      implicit none
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nj
      integer, intent(in) :: j
!
!::local variables
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  jc, j, imax, imax1, k1, k2
!ik   integer  iw, i, ip, im, ihp, ihm, iw1
      integer  imax, imax1, k1, k2
      integer  i, ip, im
      integer  m1a, m2a, m3, m4, ia
      real*8   zwt1, zwt2, zf1, zf2, zf3, zf4
      real*8   zq1, zq2, zq3, zq4, zpi, zc, zd, ze
!
!xx   write(n6,'(/2x,"*** pyvpin ***  jc =",i3)') jc
!
! modified 3/1 lines organize local variables and include files by kamata 2021/05/31
!ik   jc = nj
!ik   j  = nj
!ik   imax  = icmax(jc)
      imax  = icmax(j)
      imax1 = imax - 1
!
!----------------------------------------------------------------------
!::flux clear
!----------------------------------------------------------------------
      do k2 = 1, nion
      do k1 = 1, ndxy
      dfl1a(k1,k2) = 0.0d0
      dfl2a(k1,k2) = 0.0d0
      enddo
      enddo
!
      do k1 = 1, ndxy
      dfl3(k1) = 0.0d0
      dfl4(k1) = 0.0d0
      enddo
!
!----------------------------------------------------------------------
!::trans-function of Pi at cell center (j,i)
!----------------------------------------------------------------------
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   do iw = 1, imax
!ik   i = iw
      do i = 1, imax
      do ia = 1, nion
      r1ai(i,ia) = c13*vva(j,i,ia)**2
      r2ai(i,ia) = -c23*vva(j,i,ia)
      enddo
      r3i(i) = c23
      enddo
!
!----------------------------------------------------------------------
!::flux at cell boundary (j,i+1/2)
!----------------------------------------------------------------------
! modified 6/4 lines organize local variables and include files by kamata 2021/05/31
!ik   do iw = 1, imax1
!ik   i   = iw
!ik   ip  = iw + 1
!ik   ihp = iw
!ik   zwt1 = gwtp(j,ihp)
!ik   zwt2 = gwtm(j,ihp)
      do i = 1, imax1
      ip  = i + 1
      zwt1 = gwtp(j,i)
      zwt2 = gwtm(j,i)
!
      do ia = 1, nion
      zq1 = zwt1*vdvp(j,i)*q1a(j,i,ia)+zwt2*vdvp(j,ip)*q1a(j,ip,ia)
      zq2 = zwt1*vdvp(j,i)*q2a(j,i,ia)+zwt2*vdvp(j,ip)*q2a(j,ip,ia)
! modified 4/4 lines organize local variables and include files by kamata 2021/05/31
!ik   zf1 = gare(j,ihp)*zq1
!ik   zf2 = gare(j,ihp)*zq2
!ik   dfl1a(iw,ia) = zf1
!ik   dfl2a(iw,ia) = zf2
      zf1 = gare(j,i)*zq1
      zf2 = gare(j,i)*zq2
      dfl1a(i,ia) = zf1
      dfl2a(i,ia) = zf2
      enddo
      zq3 = zwt1*vdvp(j,i)*q3(j,i)+zwt2*vdvp(j,ip)*q3(j,ip)
      zpi = zwt1*vdvp(j,i)*vni(j,i)*vti(j,i)*cev
     >     +zwt2*vdvp(j,ip)*vni(j,ip)*vti(j,ip)*cev
      zq4 = zwt1*vdvp(j,i)*q4(j,i)+zwt2*vdvp(j,ip)*q4(j,ip)
! modified 4/4 lines organize local variables and include files by kamata 2021/05/31
!ik   zf3 = gare(j,ihp)*(zq3+zpi)
!ik   zf4 = gare(j,ihp)*c53*zq4
!ik   dfl3(iw) = zf3
!ik   dfl4(iw) = zf4
      zf3 = gare(j,i)*(zq3+zpi)
      zf4 = gare(j,i)*c53*zq4
      dfl3(i) = zf3
      dfl4(i) = zf4
      enddo
!
!----------------------------------------------------------------------
!::jacobian
!----------------------------------------------------------------------
      m3 = 2*nion + 1
      m4 = 2*nion + 2
!
! modified 7/3 lines organize local variables and include files by kamata 2021/05/31
!ik   do iw = 2, imax1
!ik   i   = iw
!ik   ip  = i + 1
!ik   im  = i - 1
!ik   ihp = i
!ik   ihm = i - 1
!ik   iw1 = iw - 1
      do i = 2, imax1
      ip  = i + 1
      im  = i - 1
!
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   zc = -gare(j,ihm)*gwtp(j,ihm)*vdvp(j,im)
!ik   zd = (gare(j,ihp)*gwtp(j,ihp)-gare(j,ihm)*gwtm(j,ihm))*vdvp(j,i)
!ik   ze =  gare(j,ihp)*gwtm(j,ihp)*vdvp(j,ip)
      zc = -gare(j,im)*gwtp(j,im)*vdvp(j,im)
      zd = (gare(j,i)*gwtp(j,i)-gare(j,im)*gwtm(j,im))*vdvp(j,i)
      ze =  gare(j,i)*gwtm(j,i)*vdvp(j,ip)
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::equation [q1a=ma*na]
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m1a,m1a,iw) = cc(m1a,m1a,iw) + zc
!ik   dd(m1a,m1a,iw) = dd(m1a,m1a,iw) + zd
!ik   ee(m1a,m1a,iw) = ee(m1a,m1a,iw) + ze
      cc(m1a,m1a,i) = cc(m1a,m1a,i) + zc
      dd(m1a,m1a,i) = dd(m1a,m1a,i) + zd
      ee(m1a,m1a,i) = ee(m1a,m1a,i) + ze
!
!::equation [q2a=ma*na*vpa]
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m2a,m2a,iw) = cc(m2a,m2a,iw) + zc
!ik   dd(m2a,m2a,iw) = dd(m2a,m2a,iw) + zd
!ik   ee(m2a,m2a,iw) = ee(m2a,m2a,iw) + ze
      cc(m2a,m2a,i) = cc(m2a,m2a,i) + zc
      dd(m2a,m2a,i) = dd(m2a,m2a,i) + zd
      ee(m2a,m2a,i) = ee(m2a,m2a,i) + ze
!
!::equation [q3=energy_I]
! modified 6/6 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m3,m1a,iw) = cc(m3,m1a,iw) + zc*r1ai(im,ia)
!ik   dd(m3,m1a,iw) = dd(m3,m1a,iw) + zd*r1ai(i, ia)
!ik   ee(m3,m1a,iw) = ee(m3,m1a,iw) + ze*r1ai(ip,ia)
!ik   cc(m3,m2a,iw) = cc(m3,m2a,iw) + zc*r2ai(im,ia)
!ik   dd(m3,m2a,iw) = dd(m3,m2a,iw) + zd*r2ai(i, ia)
!ik   ee(m3,m2a,iw) = ee(m3,m2a,iw) + ze*r2ai(ip,ia)
      cc(m3,m1a,i) = cc(m3,m1a,i) + zc*r1ai(im,ia)
      dd(m3,m1a,i) = dd(m3,m1a,i) + zd*r1ai(i, ia)
      ee(m3,m1a,i) = ee(m3,m1a,i) + ze*r1ai(ip,ia)
      cc(m3,m2a,i) = cc(m3,m2a,i) + zc*r2ai(im,ia)
      dd(m3,m2a,i) = dd(m3,m2a,i) + zd*r2ai(i, ia)
      ee(m3,m2a,i) = ee(m3,m2a,i) + ze*r2ai(ip,ia)
!
!::ff
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   ff(m1a,iw) = ff(m1a,iw) - dfl1a(iw,ia) + dfl1a(iw1,ia)
!ik   ff(m2a,iw) = ff(m2a,iw) - dfl2a(iw,ia) + dfl2a(iw1,ia)
      ff(m1a,i) = ff(m1a,i) - dfl1a(i,ia) + dfl1a(im,ia)
      ff(m2a,i) = ff(m2a,i) - dfl2a(i,ia) + dfl2a(im,ia)
      enddo !  loop(ia)
!
! modified 6/6 lines organize local variables and include files by kamata 2021/05/31
!ik   cc(m3,m3,iw) = cc(m3,m3,iw) + zc*c53
!ik   dd(m3,m3,iw) = dd(m3,m3,iw) + zd*c53
!ik   ee(m3,m3,iw) = ee(m3,m3,iw) + ze*c53
!ik   cc(m4,m4,iw) = cc(m4,m4,iw) + zc*c53
!ik   dd(m4,m4,iw) = dd(m4,m4,iw) + zd*c53
!ik   ee(m4,m4,iw) = ee(m4,m4,iw) + ze*c53
      cc(m3,m3,i) = cc(m3,m3,i) + zc*c53
      dd(m3,m3,i) = dd(m3,m3,i) + zd*c53
      ee(m3,m3,i) = ee(m3,m3,i) + ze*c53
      cc(m4,m4,i) = cc(m4,m4,i) + zc*c53
      dd(m4,m4,i) = dd(m4,m4,i) + zd*c53
      ee(m4,m4,i) = ee(m4,m4,i) + ze*c53
!
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   ff(m3,iw) = ff(m3,iw) - dfl3(iw) + dfl3(iw1)
!ik   ff(m4,iw) = ff(m4,iw) - dfl4(iw) + dfl4(iw1)
      ff(m3,i) = ff(m3,i) - dfl3(i) + dfl3(im)
      ff(m4,i) = ff(m4,i) - dfl4(i) + dfl4(im)
!
      enddo  ! loop(i)
!
      return
      end
