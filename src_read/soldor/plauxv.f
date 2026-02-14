!***********************************************************************
      subroutine plauxv
!***********************************************************************
!
!       zne ==> zne0 = sum(Za*Na) for nion
!       zne ==>  sum(Za*Na) + sum(Zz*Nz)      !  KS 20121030
!
!CHK_ia  plauxv   [do ia]   2012/12/19
!-----------------------------------------------------------------------
      use cmeffz, only : vsez, vsezz
      use cphcns, only : cev
      use cplcom, only : ama, aza, gcse, gcsi, nimin_aux, nion, q1a, q2a
     >    , q3, q4, temin_aux, timin_aux, vcs, vea, vna, vne, vnezef
     >    , vni, vva, vte, vve, vti, vzf 
      use cplmet, only : icel, itmax, jcel, jtmax, jtmin
      use csize,  only : ndx
      use csonic, only : itim, lstop, time
      use cunit,  only : n6
      implicit none
!
      integer  nter, it, jtst, jten, jerr, jt, j, i, kerr(ndx), ia
      real*8   zni, zne, zzf, zve, zki, zro, zpe, zpi, zne0
      real*8   zpemin, zpimin, zge, zgi
      real*8   zte, zti, tmp, zdnz
!
!-----------------------------------------------------------------------
!::plasma parameter at cell center
!-----------------------------------------------------------------------
      nter = 0
      do 100 it = 1, itmax
!
      jtst = jtmin(it)
      jten = jtmax(it)
      jerr = 0
      do 110 jt = jtst, jten
      j = jcel(jt,it)
      i = icel(jt,it)
      kerr(j) = 0
      zni = 0.0d0
      zne0 = 0.0d0
      zzf = 0.0d0
      zve = 0.0d0
      zki = 0.0d0
      zro = 0.0d0
!
      do ia = 1,nion
      vna(j,i,ia) = q1a(j,i,ia)/ama(ia)
      if(vna(j,i,ia).lt.nimin_aux) then
        if( vna(j,i,ia).le.0.0d0) then
          jerr = jerr + 1
          kerr(j) = kerr(j) + 1
        else
          write(n6,'(1h ,x, "*** plauxv * "," itim = ", i6,
     >      "   na(",2i4,i2," ) =",1pe11.3, " -> nimin_aux =",1pe11.3)')
     >      itim,j,i,ia, vna(j,i,ia),nimin_aux
        endif
        vna(j,i,ia)=nimin_aux
        q1a(j,i,ia)=nimin_aux*ama(ia)
      endif
      vva(j,i,ia) = q2a(j,i,ia)/q1a(j,i,ia)
      zni = zni + vna(j,i,ia)
      zne0 = zne0 + aza(ia)*vna(j,i,ia)
      zzf = zzf + aza(ia)**2*vna(j,i,ia)
      zve = zve + aza(ia)*q2a(j,i,ia)/ama(ia)
      zki = zki + 0.5d0*q2a(j,i,ia)*vva(j,i,ia)
      zro = zro + q1a(j,i,ia)
      enddo
!
!::effect of impurity   see plfimp.f
!::KSFUJI  define mcel for dummy cell  see.
      zdnz = 0.0d0
      zzf = zzf + vsezz(j,i)
!
!::limit of impurity density  (see plfimp)
      zdnz   = vsez(j,i)
      zne = zne0 + zdnz
      vnezef(j,i) = zne
!
      zpe = q4(j,i)/1.5d0
      zpi = (q3(j,i)-zki)/1.5d0
      zpemin = temin_aux*zne*cev
      zpimin = timin_aux*zni*cev
!
!::Temin
      if( zpe.lt.zpemin ) then
          write(n6,'(1h ,x, "*** plauxv * "," itim = ", i6,
     >      "   Te(",2i4" ) =",f9.4,"   -> temin_aux =",f9.4)')
     >      itim,j,i, zpe/zne/cev, temin_aux
          zpe = zpemin
          q4(j,i) = 1.5d0*zpe
          zpe = q4(j,i)/1.5d0
      endif
!::Timin
      if( zpi.lt.zpimin ) then
        if(zpi.le.0.0d0) then
          jerr = jerr + 1
          kerr(j) = kerr(j) + 1000
        else
          write(n6,'(1h ,x, "*** plauxv * "," itim = ", i6,
     >      "   Ti(",2i4" ) =",f9.4,"   -> timin_aux =",f9.4)')
     >      itim,j,i, zpi/zni/cev, timin_aux
          zpi = zpimin
          q3(j,i) = 1.5d0*zpi + zki
          zpi = (q3(j,i)-zki)/1.5d0
        endif
      endif
!
!::check of sound velocity  Note. 2009/09/08
      zge = gcse(it)
      zgi = gcsi(it)
      zti = zpi/(zni*cev)
      zte = zpe/(zne*cev)   ! <== impurity effect
!
      if( zro.le.0.0d0) then
        jerr = jerr + 1
        kerr(j) = kerr(j) + 1
        vcs(j,i) = sqrt((zti+zte)*cev/ama(1))
      else
        tmp = (zge*zpe+zgi*zpi)/zro
        if( tmp.le.0.0d0 ) then
          jerr = jerr + 1
          kerr(j) = kerr(j) + 1
          tmp = dabs(tmp)
        endif
        vcs(j,i) = sqrt(tmp)
      endif
!
      vti(j,i) = zpi/(zni*cev)
      vve(j,i) = zve/zne
      vzf(j,i) = zzf/zne
      vte(j,i) = zpe/(zne*cev)
      vne(j,i) = zne
      vni(j,i) = zni
!
      do ia = 1, nion
      vea(j,i,ia) = 1.5d0*vna(j,i,ia)*vti(j,i)*cev
     >             +0.5d0*ama(ia)*vna(j,i,ia)*vva(j,i,ia)**2
      enddo
!
!::check vva < Mach
!x      if( dabs(vva(j,i,ia)/vcs(j,i)).gt.1.08d0 ) then
!::check Te,Ti
      if( vte(j,i).lt.0.9*temin_aux ) then
        jerr = jerr + 1
        kerr(j) = kerr(j) + 100
      endif
      if( vti(j,i).lt.0.9*timin_aux ) then
        jerr = jerr + 1
        kerr(j) = kerr(j) + 1000
      endif
!
 110  continue
!
!::output
      if( jerr.gt.0 ) then
      nter = nter + 1
      if( nter.eq.1 ) then
      write(n6,'(2x,"serious error  *** plauxv ***   time =",1pe11.3,
     >  "  Ti,Te,Cs**2 < 0 .or. Cs>1")') time
      write(n6,'( 2x,1x,"it",2x,"jt",2x,"ic",3x,"jc",3x,"eimc",2x,"ni",
     >  9x,"vp",9x,"ti",9x,"te",9x,"cs",9x,"mch")')
      endif
      ia = 1
      do 220 jt = jtst, jten
      j = jcel(jt,it)
      i = icel(jt,it)
      if( kerr(j).eq.0 ) goto 220
      write(n6,'(1x,4i4,i7,1p6e11.3)') it,jt,i,j,kerr(j)
     >  ,vni(j,i),vva(j,i,ia),vti(j,i),vte(j,i),vcs(j,i)
     >  ,vva(j,i,ia)/vcs(j,i)
 220  continue
      endif
!
 100  continue
!
!-----------------------------------------------------------------------
!::plasma parameter at grid point
!-----------------------------------------------------------------------
      call plauxg
!
!-----------------------------------------------------------------------
!::flag
!-----------------------------------------------------------------------
      if(nter.gt.0) then
        lstop = 1
      endif
!
      return
      end
