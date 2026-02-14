!***********************************************************************
      subroutine plbman
!***********************************************************************
!
!  plasma parameter at edge
!      flxna, flxqi, flxq2 are specified.
!      vp = 0.0
!
!             see sub. plrprf
!
!       c52 => chfcv
!       include impurity effect   hne = hne0 + dmin1(hdnz, 0.2*hne0)
!         hdnz = sum(z) z*Nz   hne0 = sum(a) a*Nz   2012/12/11
!
!CHK_ia  plbman   [do ia]   2012/12/19
!-----------------------------------------------------------------------
      use cmeffz, only : vsez, xzflz
      use cphcns, only : cev
      use cplcom, only : ama, anamp, aza, chfcv, daan, nion, qb1a, qb2a
     >    , qb3, qb4, vdda, vdxe, vdxi, vna, vne, vni, vte, vti, xean
     >    , xian
      use cplmet, only : gdsv, gwtm, gwtp, icmpe, jcxp1, jcxp2, kcn, kcs
      use cpmpls, only : fflna, flxni, flxqe, flxqi, lnedg
      use csize,  only : ndsp
      use csonic, only : itim
      use cunit,  only : n6
      implicit none
!
!::local common
      real*8  bna(ndsp), bva(ndsp), bti, bte
! deleted 1 line replace all include files with module files by kamata 2021/08/18
!ik   common /com_plbman/ bna, bva, bti, bte
!
      real*8   zfla(ndsp)
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  ic, i1, i2, isf, jcst, jcen, ia, j, kbc
!ik   real*8   zni1, zne1, zflxa, sum1, sum2, zdf1, zdf2, zfli, zfle
      integer  ic, i1, i2, isf, jcst, jcen, ia, j
      real*8   zni1, zne1, zflxa, sum1, sum2, zdf1, zdf2, zfli
      real*8   sum, sum3, sum4, zfi, zfe, zfa, zne, zq3, smfe, smfz
      real(8) :: zne1a, zfea, zflea, qec, qev

! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   integer  jh, i
!
!x      write(n6,'(/2x,"***  plbman  ***   time =",1pe12.3)') time
!
!::index
      ic  = icmpe
      i1  = ic       ! bc at core edge
      i2  = ic - 1
      isf = i2
!
      jcst = jcxp1 + 1
      jcen = jcxp2 - 1
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   jh = (jcst+jcen)/2
!
!----------------------------------------------------------------------
!::density
!----------------------------------------------------------------------
      zni1 = 0.0d0
      zne1a = 0.0d0
      do 110 ia = 1, nion
!
!::flux
      if( lnedg.eq.0 ) then
      zflxa = flxni*fflna(ia)
      sum1  = 0.0d0
      sum2  = 0.0d0
      do 120 j = jcst, jcen
      zdf1 = vdda(j,i1,ia)
      zdf2 = vdda(j,i2,ia)
      daan(j,ia)  = 1.0d0/(gwtm(j,isf)/zdf2+gwtp(j,isf)/zdf1)
      sum1 = sum1 + gdsv(j,isf,kcn)*daan(j,ia)
      sum2 = sum2 + gdsv(j,isf,kcs)*daan(j,ia)*vna(j,i2,ia)
 120  continue
      bna(ia) = (zflxa+sum2)/sum1
!
      else
!::density
      do j = jcst, jcen
      zdf1 = vdda(j,i1,ia)
      zdf2 = vdda(j,i2,ia)
      daan(j,ia)  = 1.0d0/(gwtm(j,isf)/zdf2+gwtp(j,isf)/zdf1)
      enddo
!x    bna(ia) = bnamp(ia)
      ic = icmpe
      bna(ia) = anamp(ic,ia)
      endif
!
      bva(ia) = 0.0d0
      zni1 = zni1 + bna(ia)
      zne1a = zne1a + aza(ia)*bna(ia)
 110  continue
!
!::flux   (debug No meaning)
      zfli = 0.0d0
      zflea = 0.0d0
      do 130 ia = 1, nion
      sum  = 0.0d0
      do 140 j = jcst, jcen
      sum  = sum + gdsv(j,isf,kcn)*daan(j,ia)*bna(ia)
     >           - gdsv(j,isf,kcs)*daan(j,ia)*vna(j,i2,ia)
 140  continue
      zfla(ia) = sum
      zfli = zfli + zfla(ia)
      zflea = zflea + aza(ia)*zfla(ia)
 130  continue
!
!----------------------------------------------------------------------
!::ion temperature   NEW TYPE
!----------------------------------------------------------------------
      sum1  = 0.0d0
      sum2  = 0.0d0
      sum3  = 0.0d0
      sum4  = 0.0d0
      do 220 j = jcst, jcen
      zdf1 = vdxi(j,i1)*zni1
      zdf2 = vdxi(j,i2)*vni(j,i2)
      xian(j)  = 1.0d0/(gwtm(j,isf)/zdf2+gwtp(j,isf)/zdf1)
      sum1 = sum1 + gdsv(j,isf,kcn)*xian(j)
      sum2 = sum2 + gdsv(j,isf,kcs)*xian(j)*vti(j,i2)*cev
      zfi  = 0.0d0
      do 225 ia = 1, nion
      zfi  = zfi + gdsv(j,isf,kcn)*daan(j,ia)*bna(ia)
     >           - gdsv(j,isf,kcs)*daan(j,ia)*vna(j,i2,ia)
 225  continue
      sum3 = sum3 + gwtm(j,isf)*zfi
      sum4 = sum4 + gwtp(j,isf)*zfi*vti(j,i2)*cev
 220  continue
      bti  = (flxqi+sum2-chfcv*sum4)/(chfcv*sum3+sum1)/cev
!
!x      zcnv = chfcv*(sum3*bti*cev+sum4)
!x      write(n6,'(2x,"< plbman >  ",1pe12.3,"  bti =",1p2e12.3,
!x      "  qcnv =",1p2e12.3)') time, bti, bti2, zcnv, chfcv*zfli*bti2*cev
!
!----------------------------------------------------------------------
!::electron temperature   NEW TYPE
!----------------------------------------------------------------------
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      sum4 = 0.0d0
      smfe = 0.0d0
      smfz = 0.0d0
!
      do j = jcst, jcen
      zne1 = zne1a + vsez(j,i1)   ! impurity effect Ne*Xe
      zdf1 = vdxe(j,i1)*zne1
      zdf2 = vdxe(j,i2)*vne(j,i2)
      xean(j)  = 1.0d0/(gwtm(j,isf)/zdf2+gwtp(j,isf)/zdf1)
      sum1 = sum1 + gdsv(j,isf,kcn)*xean(j)
      sum2 = sum2 + gdsv(j,isf,kcs)*xean(j)*vte(j,i2)*cev
      zfea  = 0.0d0
      do ia = 1, nion
      zfa  = gdsv(j,isf,kcn)*daan(j,ia)*bna(ia)
     >     - gdsv(j,isf,kcs)*daan(j,ia)*vna(j,i2,ia)
      zfea  = zfea + aza(ia)*zfa
      enddo  ! loop(ia)
!
!::impuirty effect Fz = -D*grad(Nz) = -Azn*dnez(j,ip) + Azs*dnez(jj,i)
      zfe  = zfea - xzflz(j,isf)
      sum3 = sum3 + gwtm(j,isf)*zfe
      sum4 = sum4 + gwtp(j,isf)*zfe*vte(j,i2)*cev
      smfe = smfe + zfe
      smfz = smfz - xzflz(j,isf)
      enddo  ! loop(j)
      bte  = (flxqe+sum2-chfcv*sum4)/(chfcv*sum3+sum1)/cev
!
      j = jcst+10
      qec = sum1*bte*cev - sum2
      qev = chfcv*(sum3*bte*cev+sum4)
!
      if( itim.le.5 .or. mod(itim,100).eq.0 ) then
      write(n6,'(2x,"plbman  Te =",1p3e12.3,"  smfe,smfz =",1p2e12.3,
     >  "  qec,qev =",1p2e12.3)') bte, vte(j,i1), vte(j,i2),
     >  smfe, smfz, qec,qev
      endif
!
!----------------------------------------------------------------------
!::plasma parameter at edge (north side)
!----------------------------------------------------------------------
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   kbc = kcn
      do 310 j = jcst, jcen
      zne = 0.0d0
      zq3 = 0.0d0
      do 320 ia = 1, nion
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   qb1a(j,ia,kbc) = ama(ia)*bna(ia)
!ik   qb2a(j,ia,kbc) = ama(ia)*bna(ia)*bva(ia)
      qb1a(j,ia,kcn) = ama(ia)*bna(ia)
      qb2a(j,ia,kcn) = ama(ia)*bna(ia)*bva(ia)
      zne = zne + aza(ia)*bna(ia)
      zq3 = zq3 + 1.5d0*bna(ia)*bti*cev
     >          + 0.5d0*ama(ia)*bna(ia)*bva(ia)**2
 320  continue
      zne = zne + vsez(j,i1)
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   qb3(j,kbc)  = zq3
!ik   qb4(j,kbc)  = 1.5d0*zne*bte*cev
      qb3(j,kcn)  = zq3
      qb4(j,kcn)  = 1.5d0*zne*bte*cev
 310  continue
!
!::debug write
!x      if( itim.le.5 .or. mod(itim,100).eq.0 ) then
!x      j  = jcst + 10
!x      ic = icmpe
!x      i  = ic
!x      ia = 1
!x      write(n6,'(2x,"*** plbman-1 ",i8,2i5,i3,2x,a,1p4e12.3,2x,a,2x,
!x     > ,1p5e12.3)') itim, j, i, kbc,  "bq1a,bq2a,bq3,bq4 =",
!x     >   qb1a(j,ia,kbc), qb2a(j,ia,kbc), qb3(j,kbc), qb4(j,kbc),
!x     >   " <== Na,Va,Ti,Te =", bna(ia), bva(ia), bti, bte, zne
!x      call dbg_bcval
!x      endif
!
!
!::debug write
!x      write(n6,'(2x,"=== plbman ===   itim =",i3)') itim
!x      call plmflx
!x      ia = 1
!x      write(n6,'(2x,"na,Ti,Te =",1p3e12.3,"  Fa,Fa,Qi,Qe =",1p4e12.3)')
!x     >  bna(ia), bti, bte, zfla(ia), tflna(ia), tflqi, tflqe
!x      write(n6,'(2x,"jcst,jcen =",2i4,"  qb1a,qb2a,qb3,qb4 =",
!x     >  1p4e12.3)') jcst, jcen, qb1a(jcst,ia,kbc), qb2a(jcst,ia,kbc),
!x     >  qb3(jcst,kbc), qb4(jcst,kbc)
!
      return
      end
!
!***********************************************************************
      subroutine dbg_bcval
!***********************************************************************
      use cmeffz, only : vdnz
      use cplcom, only : nion, qb1a, qb2a, qb3, qb4
      use cplimp, only : ismaxl, wmc_nty
      use cplmet, only : icmpe, jcxp1, kcn
      use csize,  only : ndmis, ndsp, nzsmx
      use csonic, only : itim
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: jcst, jcen, jh, j, ic, i, kbc, ia, iz, nty
      integer :: jcst, j, ia, iz, nty
!
      real(8), dimension(ndsp)    :: bq1a, bq2a
      real(8) :: bq3, bq4
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/05/08
!ik   real(8), dimension(0:ndis2L, nzsmx) :: bnz
      real(8), dimension(0:ndmis, nzsmx) :: bnz
      real(8), dimension(ndsp)    :: bvna, bvva
      real(8) :: bvti, bvte, bvni, bvne, bvve, bvzf
!
      jcst = jcxp1 + 1
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   jcen = jcxp2 - 1
!ik   jh = (jcst+jcen)/2
      j  = jcst + 10
! modified 7/4 lines organize local variables and include files by kamata 2021/05/31
!ik   ic = icmpe
!ik   i  = ic
!ik   kbc = kcn
!ik   bq1a(1:nion) = qb1a(j,1:nion,kbc)
!ik   bq2a(1:nion) = qb2a(j,1:nion,kbc)
!ik   bq3 = qb3(j,kbc)
!ik   bq4 = qb4(j,kbc)
      bq1a(1:nion) = qb1a(j,1:nion,kcn)
      bq2a(1:nion) = qb2a(j,1:nion,kcn)
      bq3 = qb3(j,kcn)
      bq4 = qb4(j,kcn)
      do nty = 1, wmc_nty
      do iz = 0, ismaxL(nty)
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik     bnz(iz,nty) = vdnz(j,i,iz,nty)
        bnz(iz,nty) = vdnz(j,icmpe,iz,nty)
      enddo
      enddo
!
      call dbg_auxv(bq1a,bq2a,bq3,bq4,bnz,
     >              bvna, bvva, bvti, bvte, bvni, bvne, bvve, bvzf)
!
      ia = 1
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   write(n6,'(2x,"*** plbman-2 ",i8,2i5,i3,2x,a,1p4e12.3,2x,a,2x,
!ik  >  1p5e12.3)') itim, j, i, kbc,
      write(n6,'(2x,"*** plbman-2 ",i8,2i5,i3,2x,a,1p4e12.3,2x,a,2x,
     >  1p5e12.3)') itim, j, icmpe, kcn,
     >  "bq1a,bq2a,bq3,bq4 =",  bq1a(ia), bq2a(ia), bq3, bq4,
     >   " ==> Na,Va,Ti,Te =", bvna(ia), bvva(ia), bvti, bvte, bvne
!
      return
      end
!
!***********************************************************************
      subroutine dbg_auxv(bq1a, bq2a, bq3, bq4, bnz,
     >              bvna, bvva, bvti, bvte, bvni, bvne, bvve, bvzf)
!***********************************************************************
      use cphcns, only : cev
      use cplcom, only : ama, aza, nion
      use cplimp, only : ismaxl, wmc_nty
      use csize,  only : ndmis, ndsp, nzsmx
      implicit none
!
!::argument
      real(8), dimension(ndsp)    :: bq1a, bq2a
      real(8) :: bq3, bq4
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/05/08
!ik   real(8), dimension(0:ndis2L,nzsmx) :: bnz
      real(8), dimension(0:ndmis,nzsmx) :: bnz
      real(8), dimension(ndsp)    :: bvna, bvva
      real(8) :: bvti, bvte, bvni, bvne, bvve, bvzf
!
!::local variables
      real(8) :: zni, zne0, zne, zzf, zve, zki, zro
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real(8) :: zdnz, za, anz, znezef, zpe, zpi, zti, zte
      real(8) :: zdnz, za, anz, zpe, zpi
      integer :: ia, iz, nty
!
      zni = 0.0
      zne0 = 0.0
      zzf = 0.0
      zve = 0.0
      zki = 0.0
      zro = 0.0
!
      do ia = 1,nion
        bvna(ia) = bq1a(ia)/ama(ia)
        bvva(ia) = bq2a(ia)/bq1a(ia)
        zni = zni + bvna(ia)
        zne0 = zne0 + aza(ia)*bvna(ia)
        zzf = zzf + aza(ia)**2*bvna(ia)
        zve = zve + aza(ia)*bq2a(ia)/ama(ia)
        zki = zki + 0.5d0*bq2a(ia)*bvva(ia)
        zro = zro + bq1a(ia)
      enddo
!
!::effect of impurity   see plfimp.f
      zdnz = 0.0d0
      do nty = 1, wmc_nty
      do iz = 1, ismaxL(nty)
        za   = iz
        anz  = bnz(iz,nty)       ! see plfimp
        zdnz = zdnz + za*anz
        zzf  = zzf +  za**2*anz
      enddo
      enddo
!
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   znezef = zne0 + zdnz
!ik   zne = znezef
      zne = zne0 + zdnz
!
      zpe = bq4/1.5d0
      zpi = (bq3-zki)/1.5d0
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   zti = zpi/(zni*cev)
!ik   zte = zpe/(zne*cev)
!
      bvti = zpi/(zni*cev)
      bvve = zve/zne
      bvzf = zzf/zne
      bvte = zpe/(zne*cev)
      bvne = zne
      bvni = zni
!
      return
      end
