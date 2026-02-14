!**********************************************************************
! modified 1/1 lines integrated calculation with TOPICS by kamata 2020/11/16
!ik   subroutine plinit
      subroutine plinit_t
!**********************************************************************
!
!::input (cpmpls)
!       real*8  anamp(ndy,ndsp), atimp(ndy), atemp(ndy)
!       real*8  anasl(ndy,ndsp), atisl(ndy), atesl(ndy)
!       real*8  anapv(ndy,ndsp), atipv(ndy), atepv(ndy)
!
!::output(cplcom)
!       vna(j,i,ia), vva(j,i,ia), vti(j,i), vte(j,i)
!       q1a(j,i,ia), q2a(j,i,ia), q3(j,i),  q4(j,i)
!
!
!  Mach  -1.0        0.0                 0.0         +1.0
!         |-----------x-------------------x-----------|
!        zld1        zlx1                zlx2        zld2
!
!                  profile of vp                 !  %%% 2002/10/07
!-----------------------------------------------------------------------
      use cphcns, only : cev
      use cplcom, only : ama, anamp, anasl, anapv, atemp, atepv, atesl
     >    , atimp, atipv, atisl, aza, nad0, nas0, nion, q1a, q2a, q3, q4
     >    , qw1a, qw2a, qw3, qw4, ted0, tes0, tid0, tis0, vna, vne, vte
     >    , vti, vva
      use cplmet, only : icel, icspx, itmax, itmpe, itmps, itpve, itpvs
     >    , itsle, itsls, jcel, jtmax, jtmin, kreg
      use csize,  only : ndx
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  iy, ia, i, it, jw, jtst, jten, jt, jc, ic, icx, icy
      integer  ia, i, it, jw, jtst, jten, jt, jc, ic
      integer  k, j, mpr, kbc, kty
      integer  jwd1, jwd2, jwhf, jwx1, jwx2
      real*8   fwna(ndx), fwti(ndx), fwte(ndx), alnx(ndx)
      real*8   zne, zpr, x, zro, zcs, zcs1, zcs2, frg, zld1, zld2
      real*8   zlx1, zlx2, zl, zs, zv, zq3
!
!::input  (outer and inner divertot)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   iy = icspx
      ia = 1
! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   nas0 = anasl(iy,ia)
!ik   tis0 = atisl(iy)
!ik   tes0 = atesl(iy)
      nas0 = anasl(icspx,ia)
      tis0 = atisl(icspx)
      tes0 = atesl(icspx)
!
      nad0 = 5.0d19
      tid0 = 10.0d0
      ted0 = 10.0d0
!
! modified 1/1 lines integrated calculation with TOPICS by kamata 2020/11/16
!ik   write(n6,'(/2x,"*** plinit ***")')
      write(n6,'(/2x,"*** plinit_t ***")')
      write(n6,'(2x,"nas0, tis0, tes0 =",1p3e12.3)') nas0,tis0,tes0
      write(n6,'(2x,"nad0, tid0, ted0 =",1p3e12.3)') nad0,tid0,ted0
!
!----------------------------------------------------------------------
!::poloidal dependence at separatrix   fwna, fwti, fwte
!----------------------------------------------------------------------
      it = itsle
      jwd1 = jtmin(it)
      jwd2 = jtmax(it)
      jwhf = (jwd1+jwd2)/2
      do jw = jwd1, jwhf
      j = jcel(jw,it)
      i = icel(jw,it)
      if( kreg(j,i).eq.2 ) exit
      jwx1 = jw
      enddo
      do jw = jwd2, jwhf, -1
      j = jcel(jw,it)
      i = icel(jw,it)
      if( kreg(j,i).eq.2 ) exit
      jwx2 = jw
      enddo
!
      write(n6,'(2x,"jwd1 =",i4,"  jwx1 =",i4,"  jwx2 =",i4,"  jwd2 =",
     >   i4)') jwd1, jwx1, jwx2, jwd2
!
!::default
      fwna(jwd1:jwd2) = 1.0d0
      fwti(jwd1:jwd2) = 1.0d0
      fwte(jwd1:jwd2) = 1.0d0
!
!::outer divertor
      do jw = jwd1, jwx1
        x = dfloat(jw-jwd1)/dfloat(jwx1-jwd1)
        fwna(jw) = nad0/nas0*(1.0d0-x) + x
        fwti(jw) = tid0/tis0*(1.0d0-x) + x
        fwte(jw) = ted0/tes0*(1.0d0-x) + x
      enddo
!
!::inner divertor
      do jw = jwx2, jwd2
        x = dfloat(jwd2-jw)/dfloat(jwd2-jwx2)
        fwna(jw) = nad0/nas0*(1.0d0-x) + x
        fwti(jw) = tid0/tis0*(1.0d0-x) + x
        fwte(jw) = ted0/tes0*(1.0d0-x) + x
      enddo
!
      write(n6,'(2x,"jw   =",10i12)') (jw,jw=jwd1,jwd2)
      write(n6,'(2x,"fwna =",1p10e12.3)') (fwna(jw),jw=jwd1,jwd2)
      write(n6,'(2x,"fwti =",1p10e12.3)') (fwti(jw),jw=jwd1,jwd2)
      write(n6,'(2x,"fwte =",1p10e12.3)') (fwte(jw),jw=jwd1,jwd2)
!
!----------------------------------------------------------------------
!::na,ti,te in sol region
!----------------------------------------------------------------------
      do it = itsls,itsle
      jtst = jtmin(it)
      jten = jtmax(it)
!
      do jt = jtst, jten
      jc = jcel(jt,it)
      ic = icel(jt,it)
      zne = 0.0
      do ia = 1, nion
      vna(jc,ic,ia) = anasl(ic,ia)*fwna(jc)
      zne = zne + aza(ia)*vna(jc,ic,ia)
      enddo
      vne(jc,ic) = zne
      vti(jc,ic) = atisl(ic)*fwti(jc)
      vte(jc,ic) = atesl(ic)*fwte(jc)
      enddo
!
      enddo
!
!----------------------------------------------------------------------
!::na,ti,te in private region
!----------------------------------------------------------------------
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   icx = jcel(1,itsle)
      do 320 it = itpvs, itpve
      jtst = jtmin(it)
      jten = jtmax(it)
      do 330 jt = jtst, jten
      jc = jcel(jt,it)
      ic = icel(jt,it)
      zne = 0.0d0
!
      do ia = 1, nion
      vna(jc,ic,ia) = anapv(ic,ia)*fwna(jc)
      zne = zne + aza(ia)*vna(jc,ic,ia)
      enddo
      vne(jc,ic) = zne
      vti(jc,ic) = atipv(ic)*fwti(jc)
      vte(jc,ic) = atepv(ic)*fwte(jc)
 330  continue
!
 320  continue
!
!----------------------------------------------------------------------
!::pararell velocity
!----------------------------------------------------------------------
      do it = itsls, itpve
      jtst = jtmin(it)
      jten = jtmax(it)
!
!::sound velocity at divertor plate
      do k = 1, 2
      if( k.eq.1 ) jt = jtst
      if( k.eq.2 ) jt = jten
      jc = jcel(jt,it)
      ic = icel(jt,it)
      zpr = vne(jc,ic)*vte(jc,ic)
      zro = 0.0d0
      do ia = 1, nion
      zpr = zpr + vna(jc,ic,ia)*vti(jc,ic)
      zro = zro + ama(ia)*vna(jc,ic,ia)
      enddo
      zcs = sqrt(zpr*cev/zro)
      if( k.eq.1 ) zcs1 = zcs
      if( k.eq.2 ) zcs2 = zcs
      enddo
!
      if( it.le.itsle ) then
        frg = 5.0d0; mpr = 20
      else
        frg = 3.0d0; mpr = 5
      endif
      call plenc(it,alnx)
      zld1 = alnx(jtst)
      zld2 = alnx(jten)
      zlx1 = zld1 + 1.0d0/frg*(zld2-zld1)
      zlx2 = zld2 - 1.0d0/frg*(zld2-zld1)
!
      do jt = jtst, jten
      jc = jcel(jt,it)
      ic = icel(jt,it)
      zl = alnx(jt)
      if( zl.le.zlx1 ) then
      zs = dabs((zl-zlx1)/(zld1-zlx1))
      zv = -zs**mpr*zcs1
      elseif( zl.le.zlx2 ) then
      zs = 0.0d0
      zv = 0.0d0
      else
      zs = dabs((zl-zlx2)/(zld2-zlx2))
      zv = +zs**mpr*zcs2
      endif
      if( it.eq.itsls ) zv = 0.0d0  ! sol-side wall
      if( it.eq.itpve ) zv = 0.0d0  ! prv-side wall
      do ia = 1, nion
      vva(jc,ic,ia) = zv
      enddo   ! loop (ia)
      enddo   ! loop (jt)
!
      enddo   ! loop (it)
!
!----------------------------------------------------------------------
!::na,va,ti,te in main plasma
!----------------------------------------------------------------------
      do 410 it = itmps, itmpe
      jtst = jtmin(it)
      jten = jtmax(it)
      do 420 jt = jtst, jten
      jc = jcel(jt,it)
      ic = icel(jt,it)
      zne = 0.0d0
      do ia = 1, nion
      vna(jc,ic,ia) = anamp(ic,ia)
      vva(jc,ic,ia) = 0.0d0
      zne = zne + aza(ia)*vna(jc,ic,ia)
      enddo
      vti(jc,ic) = atimp(ic)
      vte(jc,ic) = atemp(ic)
      vne(jc,ic) = zne
 420  continue
 410  continue
!
!-----------------------------------------------------------------------
!::conservative variables & vne
!-----------------------------------------------------------------------
      do 510 it = 1, itmax
      do 520 jt = jtmin(it), jtmax(it)
      jc  = jcel(jt,it)
      ic  = icel(jt,it)
      zne = 0.0d0
      zq3 = 0.0
      do 530 ia = 1,nion
      q1a(jc,ic,ia) = ama(ia)*vna(jc,ic,ia)
      q2a(jc,ic,ia) = ama(ia)*vna(jc,ic,ia)*vva(jc,ic,ia)
      zne = zne + aza(ia)*vna(jc,ic,ia)
      zq3 = zq3 + 1.5d0*vna(jc,ic,ia)*vti(jc,ic)*cev
     >          + 0.5d0*ama(ia)*vna(jc,ic,ia)*vva(jc,ic,ia)**2
 530  continue
      vne(jc,ic) = zne
      q3(jc,ic)  = zq3
      q4(jc,ic)  = 1.5d0*vne(jc,ic)*vte(jc,ic)*cev
 520  continue
 510  continue
!
!-----------------------------------------------------------------------
!::boundary condition at sol-wall and prv-wall
!-----------------------------------------------------------------------
      kbc = 1
      it  = itsls
      jt  = (jtmin(it)+jtmax(it))/2
      jc = jcel(jt,it)
      ic = icel(jt,it)
      do ia = 1, nion
      qw1a(ia,kbc) = q1a(jc,ic,ia)
      qw2a(ia,kbc) = q2a(jc,ic,ia)
      enddo
      qw3(kbc) = q3(jc,ic)
      qw4(kbc) = q4(jc,ic)
!
      kbc = 2
      it = itpve
      jt = (jtmin(it)+jtmax(it))/2
      jc = jcel(jt,it)
      ic = icel(jt,it)
      do ia = 1, nion
      qw1a(ia,kbc) = q1a(jc,ic,ia)
      qw2a(ia,kbc) = q2a(jc,ic,ia)
      enddo
      qw3(kbc) = q3(jc,ic)
      qw4(kbc) = q4(jc,ic)
!
!-----------------------------------------------------------------------
!::non conservative variables
!-----------------------------------------------------------------------
!      write(n6,'(2x,"CHECK plinit_t  izmax =",i5)') izmax
      call plauxv
!
!----------------------------------------------------------------------
!::diffusion coef.
!----------------------------------------------------------------------
      call pldcof
!
!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
      write(n6,'(5x,"plinit_t normal end")')
!
      return
!
!  kty = 2 (Na,Va,Ti,Te), kty = 4 (Na,N1,N2,Ve,V1,V2)
      kty = 4
      call ploutp(kty,itsle)
      call ploutp(kty,itpvs)
      call ploutp(kty,itmps)
!
      return
!
      ia = 1
      it = itsle
!
      write(n6,'(2x,"ia =",i2,"  it =",i3)') ia,it
      call plprnt(n6,"plc",ia,it,1.0d0, 1)
      call plprnt(n6,"vna",ia,it,1.0d19,1)
      call plprnt(n6,"vva",ia,it,0.0d0, 1)
      call plprnt(n6,"vti",ia,it,0.0d0, 1)
      call plprnt(n6,"vte",ia,it,0.0d0, 1)
      call plprnt(n6,"vcs",ia,it,0.0d0, 1)
!
      return
      end
