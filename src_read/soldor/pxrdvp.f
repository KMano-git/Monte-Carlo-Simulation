!**********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   subroutine pxrdvp(nt)
      subroutine pxrdvp(it)
!**********************************************************************
!
!       pinch term    vp > 0 inward velocity
!       flux > 0      inward  (Note i-direction)
!
!        roa*Vdf => roa*Vdf + roa*Vpinch
!
!----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : c53, gg, nion, q1a, q2a, q3, q4, vni, vti 
      use cplmet, only : gare, gwtm, gwtp, icel, itmps, itsle, jcel
     >    , jtmax
      use cplqcn, only : qfy_cd, qfy_df
      use cplvpn, only : qfy_vp, vdvp
      implicit none
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  nt
      integer, intent(in) :: it
!
!::local variables
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  it, i, ip, im, ihp, ihm, jmax, jmax1, i6
!ik   integer  jw, j, jp, jm
      integer  i, ip, im, ihp, ihm, jmax, jmax1, i6
      integer  jw, j
      integer  m1a, m2a, m3, m4, ia
      real*8   zwt1, zwt2, zq1a, zq2a, zq3, zq4, zpi
      real*8   zf1p, zf2p, zf3p, zf4p, zf1m, zf2m, zf3m, zf4m
!
!xx   write(n6,'(/2x,"*** pxrdvp ***  nt =",i3)') nt
!
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   it = nt
      i  = icel(1,it)
      ip = i + 1
      im = i - 1
      ihp = i
      ihm = i - 1
      jmax  = jtmax(it)
      jmax1 = jmax - 1
!
      m3 = 2*nion + 1
      m4 = 2*nion + 2
!
!::loop x-direction
      do jw = 2, jmax1
      j  = jcel(jw,it)
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   jp = jcel(jw+1,it)
!ik   jm = jcel(jw-1,it)
!
!::loop ia
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
!
!::flux-(q1a,q2a) at (j,i+1/2)
      zwt1 = gwtp(j,ihp)
      zwt2 = gwtm(j,ihp)
      zq1a = zwt1*vdvp(j,i)*q1a(j,i,ia) + zwt2*vdvp(j,ip)*q1a(j,ip,ia)
      zq2a = zwt1*vdvp(j,i)*q2a(j,i,ia) + zwt2*vdvp(j,ip)*q2a(j,ip,ia)
      zf1p = gare(j,ihp)*zq1a
      zf2p = gare(j,ihp)*zq2a
!
!::flux-(q1a,q2a) at (j,i-1/2)
      zwt1 = gwtp(j,ihm)
      zwt2 = gwtm(j,ihm)
      zq1a = zwt1*vdvp(j,im)*q1a(j,im,ia) + zwt2*vdvp(j,i)*q1a(j,i,ia)
      zq2a = zwt1*vdvp(j,im)*q2a(j,im,ia) + zwt2*vdvp(j,i)*q2a(j,i,ia)
      zf1m = gare(j,ihm)*zq1a
      zf2m = gare(j,ihm)*zq2a
!
!::residuals
      gg(m1a,jw) = gg(m1a,jw) - zf1p + zf1m
      gg(m2a,jw) = gg(m2a,jw) - zf2p + zf2m
!
!::qcon
      qfy_vp(j,i,m1a) = zf1p
      qfy_vp(j,i,m2a) = zf2p
      if( it.eq.2 ) then
      qfy_vp(j,im,m1a) = zf1m
      qfy_vp(j,im,m2a) = zf2m
      endif
!
      enddo   !  loop(ia)
!
!::flux-(q3,q4) at (j,i+1/2)
      zwt1 = gwtp(j,ihp)
      zwt2 = gwtm(j,ihp)
      zq3  = zwt1*vdvp(j,i)*q3(j,i) + zwt2*vdvp(j,ip)*q3(j,ip)
      zpi  = zwt1*vdvp(j,i)*vni(j,i)*vti(j,i)*cev +
     >       zwt2*vdvp(j,ip)*vni(j,ip)*vti(j,ip)*cev
      zq4  = zwt1*vdvp(j,i)*q4(j,i) + zwt2*vdvp(j,ip)*q4(j,ip)
      zf3p = gare(j,ihp)*(zq3+zpi)
      zf4p = gare(j,ihp)*c53*zq4
!
!::flux-(q3,q4) at (j,i-1/2)
      zwt1 = gwtp(j,ihm)
      zwt2 = gwtm(j,ihm)
      zq3  = zwt1*vdvp(j,im)*q3(j,im)  + zwt2*vdvp(j,i)*q3(j,i)
      zpi  = zwt1*vdvp(j,im)*vni(j,im)*vti(j,im)*cev +
     >       zwt2*vdvp(j,i)*vni(j,i)*vti(j,i)*cev
      zq4  = zwt1*vdvp(j,im)*q4(j,im) + zwt2*vdvp(j,i)*q4(j,i)
      zf3m = gare(j,ihm)*(zq3+zpi)
      zf4m = gare(j,ihm)*c53*zq4
!
!::residuals
      gg(m3,jw) = gg(m3,jw) - zf3p + zf3m
      gg(m4,jw) = gg(m4,jw) - zf4p + zf4m
!
!::qcon
      qfy_vp(j,i,m3) = zf3p
      qfy_vp(j,i,m4) = zf4p
      if( it.eq.2 ) then
      qfy_vp(j,im,m3) = zf3m
      qfy_vp(j,im,m4) = zf4m
      endif
!
      enddo  ! loop(jw)
!
      return
!
!::debug write
      if( it.eq.itsle .or. it.eq.itmps ) then
      i6 = 112
      ia = 1
      m1a = 2*ia-1
      write(i6,'(/2x,"*** pxrdvp ***  it =",i3)') it
      write(i6,'(5x,"jw",2x,"j",3x,"i",4x,"df(m1)",5x,"vp(m1)",5x,
     > "df(m3)",5x,"cd(m3)",5x,"vp(m3)",5x,"df(m4)",5x,"cd(m4)",5x,
     > "vp(m4)")')
!
      do jw = 1, jmax
      j = jcel(jw,it)
      write(i6,'(2x,3i4,1p8e11.2)')
     >   jw, j, i, qfy_df(j,i,m1a),  qfy_vp(j,i,m1a),
     >  qfy_df(j,i,m3), qfy_cd(j,i,m3), qfy_vp(j,i,m3),
     >  qfy_df(j,i,m4), qfy_cd(j,i,m4), qfy_vp(j,i,m4)
      enddo
!
      if( it.eq.itmps ) call wexit("pxrdvp","test")
      endif
      return
      end
