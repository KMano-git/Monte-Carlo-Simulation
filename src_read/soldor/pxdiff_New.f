!***********************************************************************
      subroutine pxdiff(it)
!***********************************************************************
!
!       classical & anomalous diffusion coefficients
!          Xi enhance at initial state
!          include impurity effect
!
!  eta calculated by wtip = dmax1( wti, 2.5d0 )   14/2/18  K.Shimizu
!  minimun Ti for eta timneta  by hoshino  140514
!  2019.10.04 YAMOTO
!  Balescu model implemented (activate when mdl_bale.eq.1 @ inppls)
!-----------------------------------------------------------------------
      use cphcns, only : cev, cme, cmp
      use cplcom, only : aion, alnr, c43, daan, dacl, etan, etcl, feqp
     >    , kappa_helander_test, ldluc, mdl_bale, mdl_ini, nion, nlp
     >    , temndf, timndf, timneta, vdda, vdet, vdxe, vdxi, vna, vne
     >    , vni, vte, vti, vzf, xean, xecl, xian, xicl, xicl_test
      use cplimp, only : ismaxl, wmc_nty
      use cplmet, only : icel, itmps, itsle, jcel, jtmax, kreg
      use cmeffz, only : mion, vdnz, wai, wfab, wfeb, wkab, wma 
      use csize,  only : mdsp, ndx
      use csonic, only : itim
      use cunit,  only : n6
      implicit none
!
      integer, intent(in) :: it
!
      real*8  wna(ndx,mdsp), wni(ndx), wne(ndx)
      real*8  wti(ndx), wte(ndx), wzf(ndx)
      real(8) :: wtip(ndx)
!
      integer  i, jmax, jw, j, ia, ib
      integer  ii, iz, irg, nty
      real*8   zte32, zti32,  fke, fki, feq, feta
      real(8) :: zti32p
      real*8   zeta, zke, zki, fka
      real*8   zfce, zfci, flce(ndx), flci(ndx)
      real*8   hzeff, hAi, hZi, hne, hte, hti, zlne, coulog
      real*8   prefac_ke, prefac_ki, prefac_eta
!
      i  = icel(1,it)
      jmax = jtmax(it)
!
!-----------------------------------------------------------------------
!::minimun Te and Ti for diffusion coefficients  temndf, timndf
!-----------------------------------------------------------------------
!::define at plinpt
!
!-----------------------------------------------------------------------
!::impurity effect
!-----------------------------------------------------------------------
!::clear
      do j  = 1, ndx
      do ia = 1, mdsp
      wna(j,ia) = 0.0d0
      enddo
      enddo
!
      do jw = 1, jmax
      j = jcel(jw,it)
!
      ii = 0
      do ia = 1, nion
      ii = ii + 1
      wna(j,ii) = vna(j,i,ia)
      enddo
!
      do nty = 1, wmc_nty
      do iz = 1, ismaxL(nty)
      ii = ii + 1
      wna(j,ii) = vdnz(j,i,iz,nty)
      enddo
      enddo
!
      wni(j) = vni(j,i)
      wne(j) = vne(j,i)
      wzf(j) = vzf(j,i)
!
      wti(j) = vti(j,i)
      wte(j) = vte(j,i)
      wtip(j) = vti(j,i)
      wte(j) = dmax1( wte(j), temndf )   !  <== Optcheck
      wti(j) = dmax1( wti(j), timndf )
      wtip(j) = dmax1( wtip(j), timneta )  !  use for eta
      enddo !jw
!
!::reduction in Xe, Xi in Div region for initial calculation
      if( mdl_ini.eq.1 ) then
      do jw = 1, jmax
      j = jcel(jw,it)
      irg = kreg(j,i)
      if( irg.eq.1 .or. irg.eq.3 ) then
      wte(j) = dmin1( wte(j), 100.0d0 )
      wti(j) = dmin1( wti(j), 150.0d0 )
      endif
      enddo
      endif
!
!-----------------------------------------------------------------------
!::classical-diffusion coefficients at cell center   B2 version
!-----------------------------------------------------------------------
      do jw = 1, jmax
      j  = jcel(jw,it)
!
      zte32 = wte(j)*sqrt(wte(j))
      zti32 = wti(j)*sqrt(wti(j))
      zti32p = wtip(j)*sqrt(wtip(j))
!
!::coulomb logarithm
!   set lne = 12  2002/12/3
!   subroutine made by M. Honda  2011/06/17
!
!xx   zlne = 12.0d0     !  <=== Optcheck
      hzeff = wzf(j)
      hAi = aion(1)
      hZi = 1.0d0
      hne = wne(j)
      hte = wte(j)
      hti = wti(j)
      zlne = coulog( hzeff, hne, hte, hti, hAi, hZi )
      alnr(j) = zlne
!
!::frequency, coefficient
      fke = 0.0d0
      fki = 0.0d0
      feq = 0.0d0
      do ia = 1, mion    ! <== Note 2012/12/17
!
      feta = 0.0d0
      fka = 0.0d0
      do ib = 1, mion
      feta = feta + wfab(ia,ib)*wna(j,ib)
      fka  = fka  + wkab(ia,ib)*wna(j,ib)
      enddo
      feta = feta*zlne/zti32p
      if(mdl_bale.eq.1)then
         prefac_eta = 1.81d0
      else
         prefac_eta = 0.96d0
      end if
!      zeta = 0.96*wna(j,ia)*wtip(j)*cev/feta
      zeta = prefac_eta*wna(j,ia)*wtip(j)*cev/feta
      if( ia.le.nion ) etcl(j,ia) = c43*zeta
!
      fke = fke + wfeb(ia)*wna(j,ia)
      fka = fka*zlne/zti32
      fki = fki + wna(j,ia)/fka
      feq = feq + wfeb(ia)/wai(ia)*wna(j,ia)
      enddo
!
      fke = fke*zlne/zte32
!      zke = 3.16*wne(j)*wte(j)*cev/cme/fke
!      zki = 3.90*wti(j)*cev/cmp*fki
      if(mdl_bale.eq.1)then
         prefac_ke = 2.5d0*2.16d0*wzf(j)/(1.0d0+0.27d0*wzf(j))
         prefac_ki = 3.98d0
      else
         prefac_ke = 3.16d0
         prefac_ki = 3.90d0
      end if
      zke = prefac_ke*wne(j)*wte(j)*cev/cme/fke
      zki = prefac_ki*wti(j)*cev/cmp*fki
      feq = feq*zlne/zte32*3.00*cme/cmp
!
      xecl(j) = zke
      xicl(j) = zki
      feqp(j) = feq

!     kappa // test 2018/05/29 Y.Homma
      xicl_test(j,i) = xicl(j)
! HERE CALCULATE KAPPA HELANDER, you can use zlne above.
      call kappa_para_helander( wni(j), wti(j),
     >     alnr(j), kappa_helander_test(j,i) )
      enddo
!
!-----------------------------------------------------------------------
!::Luciani factor in the main plasma
!-----------------------------------------------------------------------
      if( ldluc.eq.1 .and. it.ge.itmps ) then
      do jw = 2, jmax-1
      j = jcel(jw,it)
      zfce = 0.3d16*wte(j)**2/(wzf(j)*wne(j))
      zfci = 0.3d16*wti(j)**2/(wzf(j)*wne(j))
      flce(j) = 1.0d0/(1.0d0+zfce)
      flci(j) = 1.0d0/(1.0d0+zfci)
      do ia = 1, nion
      etcl(j,ia) = etcl(j,ia)*flci(j)
      enddo
      xecl(j) = xecl(j)*flce(j)
      xicl(j) = xicl(j)*flci(j)

!     kappa // test 2018/05/29 Y.Homma
      xicl_test(j,i) = xicl(j)
      enddo
      endif
!
!-----------------------------------------------------------------------
!::anomalous-diffusion coefficients at cell center
!-----------------------------------------------------------------------
      do jw = 1, jmax
      j = jcel(jw,it)
      do ia = 1, nion
      dacl(j,ia) = 0.0d0
      daan(j,ia) = vdda(j,i,ia)
      etan(j,ia) = vdet(j,i,ia)*wma(ia)*wna(j,ia)
      enddo
      xian(j) = vdxi(j,i)*wni(j)
      xean(j) = vdxe(j,i)*wne(j)
      enddo
!
 100  continue
      return
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      if( it.eq.itsle ) then
      if( itim.eq.0 .or. nlp.eq.5 ) then
      write(n6,'(/2x,"*** pxdiff_New ***  itim =",i5,"  nlp =",i2)')
     >    itim, nlp
      do jw = 1, jmax
      j = jcel(jw,it)
      write(n6,'(3i5,1p6e12.3)')
     >  jw, j, i, vne(j,i), wne(j), wzf(j), xecl(j), xicl(j), feqp(j)
      enddo
      endif
      endif
!
      return
      end
!
!fj del >
!      subroutine test_pxdiff
!      real*8  dacl(2,2)
!      integer jno(3)
!!
!!-----------------------------------------------------------------------
!!::dacl in private region
!!-----------------------------------------------------------------------
!      ii = 0
!      do ix = 1, ndxy
!      jno(ix) = 0
!      do ia = 1, ndsp
!        dacl(ix,ia) = 0.0d0
!      enddo
!      enddo
!      goto 100
!!
!      if( it.ge.itpvs .and. it.le.itpve ) then
!      ii = 0
!      do jw = 5, jmax-5
!      j = jcel(jw,it)
!      if( vne(j,i).ge.1.2d20 ) cycle
!      jno(ii) = j
!      ii = ii + 1
!      do ia = 1, nion
!        dacl(j,ia) = 3.0d6
!      enddo
!      enddo
!      endif
!!
! 100  continue
!      if(  mod(itim,50).eq.0 .and. nlp.eq.5  .and. ii.gt.0 ) then
!
!      i6 = 992
!      write(i6,'(2x,"pxdiff ",i7,i3,i4,1pe11.2,100i4)')
!     >   itim, nlp, it, dacl(jno(1),1), (jno(j),j=1,ii)
!      endif
!!
!      return
!      end
!fj del <

