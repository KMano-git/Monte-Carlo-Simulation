!**********************************************************************
      subroutine plinit
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
      use cplcom, only : ama, anamp, anasl, atemp, atesl, atimp, atisl
     >    , aza, nad0, nas0, nion, q1a, q2a, q3, q4, qw1a, qw2a, qw3
     >    , qw4, tes0, tis0, ted0, tid0, vna, vne, vni, vte, vti, vva
      use cplmet, only : icel, itmax, itmpe, itmps, itpve, itpvs, itsle
     >    , itsls, jcel, jtmax, jtmin, kreg
      use csize,  only : ndsp, ndx
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  it, jt, jtst, jten, jc, ic, ia, kbc, i, k, mpr
      integer  kty
      real*8   zl, zne, zpr, zro, zcs
      real*8   zni, zq3, zcs1, zcs2, zld1, zld2, zlx1, zlx2, zs, zv
      real*8   alnx(ndx), frg
!
      integer  j, jwd1, jwd2, jwhf, jw, jwx1, jwx2, icx
      real*8   fratio(ndsp), x
      real*8   fwna(ndx), fwti(ndx), fwte(ndx), fcprv
!
      write(n6,'(/2x,"*** plinit ***")')
!
!----------------------------------------------------------------------
!::constant type
!----------------------------------------------------------------------
      it = itsle
      jt = 1
      ia = 1
      j  = jcel(jt,it)
      i  = icel(jt,it)
      if( nas0.eq.0.0d0 ) nas0 = anasl(i,ia)
      if( tis0.eq.0.0d0 ) tis0 = atisl(i)
      if( tes0.eq.0.0d0 ) tes0 = atesl(i)
!
      if( nad0.eq.0.0d0 ) nad0 = nas0
      if( tid0.eq.0.0d0 ) tid0 = tis0
      if( ted0.eq.0.0d0 ) ted0 = tes0
!
      do ia = 1, nion
      fratio(ia) = anasl(i,ia)/anasl(i,1)
      enddo
!
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
      write(n6,'(2x,"  nas0 =",1pe12.4,"  tis0 =",1pe12.4,"  tes0 =",
     >  1pe12.4)') nas0, tis0, tes0
      write(n6,'(2x,"  nad0 =",1pe12.4,"  tid0 =",1pe12.4,"  ted0 =",
     >  1pe12.4)') nad0, tid0, ted0
      write(n6,'(2x,"jwd1 =",i4,"  jwx1 =",i4,"  jwx2 =",i4,"  jwd2 =",
     >   i4)') jwd1, jwx1, jwx2, jwd2
!
      fwna(jwd1:jwd2) = 1.0d0
      fwti(jwd1:jwd2) = 1.0d0
      fwte(jwd1:jwd2) = 1.0d0
!
      do jw = jwd1, jwx1
        x = dfloat(jw-jwd1)/dfloat(jwx1-jwd1)
        fwna(jw) = nad0/nas0*(1.0d0-x) + x
        fwti(jw) = tid0/tis0*(1.0d0-x) + x
        fwte(jw) = ted0/tes0*(1.0d0-x) + x
      enddo
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
      vna(jc,ic,ia) = nas0*fwna(jt)*fratio(ia)
      zne = zne + aza(ia)*vna(jc,ic,ia)
      enddo
      vne(jc,ic) = zne
      vti(jc,ic) = tis0*fwti(jt)
      vte(jc,ic) = tes0*fwte(jt)
      enddo
!
      enddo
!
!----------------------------------------------------------------------
!::na,ti,te in private region
!----------------------------------------------------------------------
      write(n6,'(2x,"*** plinitL ***  prv region")')
      icx = jcel(1,itsle)
      do 320 it = itpvs, itpve
      jtst = jtmin(it)
      jten = jtmax(it)
      fcprv = 0.9d0-0.3d0*(it-itpvs)/(itpve-itpvs)  ! KH111108
      do 330 jt = jtst, jten
      jc = jcel(jt,it)
      ic = icel(jt,it)
      zne = 0.0d0
      do 335 ia = 1, nion
      vna(jc,ic,ia) = vna(jc,icx,ia)*fcprv
      vna(jc,ic,ia) = dmax1(vna(jc,ic,ia),0.3d0*nad0*fratio(ia)*fcprv)
      zne = zne + aza(ia)*vna(jc,ic,ia)
 335  continue
      vne(jc,ic) = zne
      vti(jc,ic) = vti(jc,icx)*fcprv
      vte(jc,ic) = vte(jc,icx)*fcprv
!::debug
!----------
      if( jt.eq.jtst ) then
      write(n6,'(2x,"jc,ic =",2i5,"  icx =",i5,"  vti =",1pe12.3,
     >   "  fcprv =",f8.3)') jc, ic, icx, vti(jc,ic), fcprv
      endif
!----------
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
      zni = 0.0d0
      zne = 0.0d0
      do 430 ia = 1, nion
      vna(jc,ic,ia) = anamp(ic,ia)
      vva(jc,ic,ia) = 0.0d0
      zni = zni + vna(jc,ic,ia)
      zne = zne + aza(ia)*vna(jc,ic,ia)
 430  continue
      vti(jc,ic) = atimp(ic)
      vte(jc,ic) = atemp(ic)
      vni(jc,ic) = zni
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
      call plauxv
!
!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
      write(n6,'(5x,"plinit normal end")')
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
