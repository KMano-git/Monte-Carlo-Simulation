!**********************************************************************
      subroutine ntpump(ktyp)
!**********************************************************************
!
!    Conductance   Q = C * P   [Pa*m3/s] = [Watt] = [m3/s] * [Pa]
!
!         Q = 1/4*v*N * S * T = F * S * T
!
!        [Pa*m3/s] = [1/s] * [m2] * [Pa/m2*m^3]
!
!         Q(mole) = 1/2 * (flux of atom) * 300/11600 * 1.6021e-19
!                 = 1/2 * (flux of atom) * 4.143e-21
!
!          flux of D  1 Pam3/s = 1/2.07e-21 = 0.5e21 1/s
!
!         zcnd(i,1) : cond. for atom
!         zcnd(i,2) : cond. for mole
!         zcnd(i,3) : cond. for atom & mole
!
!        lout = 1     cond
!        lout = 2     cond & flux
!        lout = 3     cond & flux & density
!        lout = 4     pressure at wall
!
!----------------------------------------------------------------------
      use clocal, only : vden0, vdeng, vpre0, vpreg
      use cntcom, only : ccnd, cflx, chwl, cpre, cqfl, cstyp, cvol, den0
     >    , deng, eng0, engg, icwl, igtw, ihwl, ipmp, migx, migy, ngas
     >    , nhwl, nopm, nowl2, npew, npmp, npsw,  npwl2, nwpm, rcywl
     >    , temwl, volm
      use cntmnt, only : dotn
      use cntwcn, only : swnrm, swsum, whta, whtm, wnrm, wtot
      use cphcns, only : cev
      use csize,  only : ndgs, ndpm, ndwp
      use csonic, only : itim, kpcn, time
      use cunit,  only : lmspe, lmype, n6
      implicit none
!
!::argument
      integer, intent(in) :: ktyp
!
!::local variables
      real*8  fc
      integer ih, iw, ig, ic, isiz, j, nw, iws, iwe, i, k
!
      real*8  vhta(30,3), vhtm(30,3)
      real(8)    ze0, zeg
      save       vhta, vhtm    ! SGI-5/1
! vden0, vdeng, vpre0, vpreg moved to clocal.f

      real*8  pmvl(ndpm), pmn0(ndpm,ndgs), pmp0(ndpm,ndgs)
      real*8  pmng(ndpm,2), pmpg(ndpm,2)
      save    pmvl, pmn0, pmp0, pmng, pmpg
      real*8  ztmdg, ztmwl
      real*8  zf0, zn0, zt0, zp0, zc0, zfg, zng, ztg, zpg, zcg, zct
!
      integer lout; data lout/3/; save lout
!
      if( npmp.le.0 ) return

!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!::clear
!-----------------------------------------------------------------------
      if( ktyp.eq.1 ) then
      if( nhwl+1.gt.30 ) call wexit("ntpump","nhwl+1.gt.30")
      isiz = 30*3
      call setd( vhta, isiz, 0.0d0 )
      call setd( vhtm, isiz, 0.0d0 )
      isiz = ndwp*ndgs
      call setd( vden0, isiz, 0.0d0 )
      call setd( vpre0, isiz, 0.0d0 )
      isiz = ndwp*2
      call setd( vdeng, isiz, 0.0d0 )
      call setd( vpreg, isiz, 0.0d0 )
      isiz = ndpm*ndgs
      call setd( pmn0, isiz, 0.0d0 )
      call setd( pmp0, isiz, 0.0d0 )
      isiz = ndpm*2
      call setd( pmng, isiz, 0.0d0 )
      call setd( pmpg, isiz, 0.0d0 )
!
      if( lmype.eq.lmspe .and. itim.le.1 ) then
      write(68,'(3x,"itim",2x,"time",7x,"wnam",6x,"ccnd",8x,"cpre",
     >  8x,"cflx",6x,"wnam",6x,"ccnd",8x,"cpre",8x,"cflx")')
      endif
!
!-----------------------------------------------------------------------
!::dplate, wall, recomb
!-----------------------------------------------------------------------
      elseif( ktyp.eq.2 ) then
!
      fc = dotn/swnrm
      do ih = 1, nhwl+1
      do j  = 1, 3
      vhta(ih,j) = vhta(ih,j) + fc*whta(ih,j)
      vhtm(ih,j) = vhtm(ih,j) + fc*whtm(ih,j)
      enddo
      enddo
!
      call ntdens( dotn )
!
!::density & pressure at albdeo
      call setd( pmvl, ndpm, 0.0d0 )
      do i = 1, npmp
      iw = ipmp(i)
      ic = icwl(iw)
      k = -igtw(iw)
      if( k.le.0 ) cycle
      pmvl(k) = pmvl(k) + volm(ic)
      do ig = 1, ngas
      pmn0(k,ig) = pmn0(k,ig) + den0(ic,ig)*volm(ic)
      pmp0(k,ig) = pmp0(k,ig) + den0(ic,ig)*eng0(ic,ig)*volm(ic)
      enddo
      do ig = 1, 2
      pmng(k,ig) = pmng(k,ig) + deng(ic,ig)*volm(ic)
      pmpg(k,ig) = pmpg(k,ig) + deng(ic,ig)*engg(ic,ig)*volm(ic)
      enddo
      enddo
!
!::density & pressure at wall
      do iw = 1, npwl2
      ih = ihwl(iw)
      ic = icwl(iw)
      if( ic.le.0 ) cycle
      do ig = 1, ngas
      vden0(iw,ig) = vden0(iw,ig) + den0(ic,ig)
      vpre0(iw,ig) = vpre0(iw,ig) + den0(ic,ig)*eng0(ic,ig)
      enddo
      do ig = 1, 2
      vdeng(iw,ig) = vdeng(iw,ig) + deng(ic,ig)
      vpreg(iw,ig) = vpreg(iw,ig) + deng(ic,ig)*engg(ic,ig)
      enddo
      enddo
!
!xx   if( kpcn.eq.0 ) return
      if( mod(itim,1000).ne.0 ) return
!
!::debug write
!::conduc
      if( lout.ge.1 ) then
      write(n6,'(/2x,"*** ntpump(",i2," ) ***  itim =",i6,"  time =",
     >  1pe12.3,2x,i2,2x,a)') ktyp,itim,time,ktyp,cstyp
      write(n6,'(3x,"ih",3x,"cwal",4x,"vol",9x,"Ctot",8x
     >  "Fa",10x,"Ca",10x,"N0",10x,"T0",10x,"P0",10x,
     >  "Fm",10x,"Cm",10x,"Ng",10x,"Tg",10x,"Pg")')
      do i = 1, nopm
      ih = nwpm(i)
      zn0 = 0.0d0; zt0 = 0.0d0; zp0 = 0.0d0; zc0 = 0.0d0
      zng = 0.0d0; ztg = 0.0d0; zpg = 0.0d0; zcg = 0.0d0
      zf0 = vhta(ih,3)
      zn0 = pmn0(i,1)/pmvl(i)
      zp0 = pmp0(i,1)/pmvl(i)*cev
      if( pmn0(i,1).gt.0.0d0 ) zt0 = pmp0(i,1)/pmn0(i,1)
      if( zp0.gt.0.0d0 ) zc0 = zf0*zt0*cev/zp0
      zfg = vhtm(ih,3)
      zng = pmng(i,1)/pmvl(i)
      zpg = pmpg(i,1)/pmvl(i)*cev
      if( pmng(i,1).gt.0.0d0 ) ztg = pmpg(i,1)/pmng(i,1)
      if( zpg.gt.0.0d0 ) zcg = 0.5d0*zfg*ztg*cev/zpg
      zct = 0.0d0
      if( zp0+zpg.gt.0.0d0 )
     >zct = (zf0*zt0*cev+0.5d0*zfg*ztg*cev)/(zp0+zpg)
      write(n6,'(2x,i3,2x,a,2x,1p12e12.3)')
     >  ih, chwl(ih), pmvl(i), zct,
     >  zf0, zc0, zn0, zt0, zp0, zfg, zcg, zng, ztg, zpg
      enddo
      endif
!
!::flux
      if( kpcn.eq.0 ) return
!
      if( lout.ge.2 ) then
      write(n6,'(2x)')
      write(n6,'(4x,"flux to wall   dotn =",1pe12.3,"  wtot =",1p2e12.3,
     >  "  wnrm =",1p2e12.3,3x,a)') dotn,wtot,swsum,wnrm,swnrm
     >  ," 1:incident, 2:reflect,  3:absorp/pump"
      write(n6,'(3x,"ih",3x,"chwl",4x,"atm1",8x,"atm2",8x,"atm3",
     >  8x,"atmc",8x,"mol1",8x,"mol2",8x,"mol3",8x,"molc",8x,"fin ",
     >  8x,"frfl",8x,"fabs")')
      do ih = 1, nhwl+1
      if( ih.eq.nhwl+1 .or.
     >    whta(ih,3)+whtm(ih,3).gt.0.0 ) then
      write(n6,'(2x,i3,2x,a,2x,1p12e12.3)')
     >  ih, chwl(ih)
     > ,fc*whta(ih,1), fc*whta(ih,2), fc*whta(ih,3)
     > ,fc*(whta(ih,1)-whta(ih,2)-whta(ih,3))
     > ,fc*whtm(ih,1), fc*whtm(ih,2), fc*whtm(ih,3)
     > ,fc*(whtm(ih,1)-whtm(ih,2)-whtm(ih,3))
     > ,fc*(whta(ih,1)+whtm(ih,1))
     > ,fc*(whta(ih,2)+whtm(ih,2))
     > ,fc*(whta(ih,3)+whtm(ih,3))
      endif
      enddo
      endif
!
!::density at albedo
      if( lout.ge.3 ) then
      ig = 1
      write(n6,'(2x)')
      write(n6,'(2x,"density and pressure just before albedo")')
      write(n6,'(4x,"iw",3x,"ih",3x,"chwl",3x"ic",4x,"icx",2x,"icy",4x
     > ,"d0",10x,"e0",10x,"dg",10x,"eg",10x,"pd2",9x,"ptot")')
      do i = 1, npmp
      iw = ipmp(i)
      ic = icwl(iw)
      ih = ihwl(iw)
      if( ic.le.0 ) cycle
      write(n6,'(2x,2i5,2x,a,3i6,1p6e12.3)')
     >  iw,ih,chwl(ih),ic,migx(ic),migy(ic),den0(ic,ig),eng0(ic,ig)
     > ,deng(ic,1),engg(ic,1),deng(ic,1)*engg(ic,1)*cev
     > ,(den0(ic,ig)*eng0(ic,ig)+deng(ic,1)*engg(ic,1))*cev
      enddo
      endif
!
!-----------------------------------------------------------------------
!::total
!-----------------------------------------------------------------------
      elseif( ktyp.eq.3 ) then
!
!::conductance  Q = C * P
      call setd( cflx, ndpm, 0.0d0 )
      call setd( cqfl, ndpm, 0.0d0 )
      call setd( cpre, ndpm, 0.0d0 )
      call setd( ccnd, ndpm, 0.0d0 )
      call setd( cvol, ndpm, 0.0d0 )
!
      do iw = 1, npwl2
      k = -igtw(iw)
      if( k.le.0 ) cycle
      ic = icwl(iw)
      if( ic.eq.0 ) cycle    !  Fuji  2010/03/06
      cpre(k) = cpre(k) + vpreg(iw,1)*cev*volm(ic)
      cvol(k) = cvol(k) + volm(ic)
      enddo
!
      do i = 1, nopm
      ih = nwpm(i)
      ztmdg = temwl(ih)
      ztmwl = (ztmdg+273.0d0)/11600.0d0
      cflx(i) = vhtm(ih,3)
      cqfl(i) = 0.5d0*cflx(i)*ztmwl*cev
      ccnd(i) = 0.0d0
      if( cvol(i).gt.0.0d0 ) cpre(i) = cpre(i)/cvol(i)
      if( cpre(i).gt.0.0d0 ) ccnd(i) = cqfl(i)/cpre(i)
      enddo
!
!::debug write
      if( lmype.eq.lmspe ) then
      write(68,'(2x,i5,1pe12.4,5(2x,a,2x,1p3e12.3))')
     >   itim, time, (chwl(nwpm(i)),ccnd(i),cpre(i),cflx(i),i=1,nopm)
      endif
!
      if( kpcn.eq.0 ) return
!
!::cond
      if( lout.ge.1 ) then
      write(n6,'(/2x,"*** ntpump(",i2," ) ***  itim =",i6,"  time =",
     >  1pe12.3,2x,i2,2x,a)') ktyp,itim,time,ktyp,"Total"
      write(n6,'(2x,"conductance",3x,"rcyw",6x,"temw",6x,"pflx",8x,
     >  "Qflx",8x,"Pre",9x,"Cond")')
      do i = 1, nopm
      ih = nwpm(i)
      write(n6,'(2x,i3,2x,a,2x,f8.3,2x,f8.2,2x,1p4e12.3)')
     >  ih,chwl(ih),rcywl(ih),temwl(ih),cflx(i),cqfl(i),cpre(i),ccnd(i)
      enddo
      endif
!
!::flux
      if( lout.ge.2 ) then
      write(n6,'(3x,"ih",3x,"chwl",4x,"atm1",8x,"atm2",8x,"atm3",
     >  8x,"atmc",8x,"mol1",8x,"mol2",8x,"mol3",8x,"molc",8x,"fin ",
     >  8x,"frfl",8x,"fabs")')
      do ih = 1, nhwl+1
      write(n6,'(2x,i3,2x,a,2x,1p12e12.3)')
     >  ih, chwl(ih)
     > ,vhta(ih,1),vhta(ih,2),vhta(ih,3)
     > ,vhta(ih,1)-vhta(ih,2)-vhta(ih,3)
     > ,vhtm(ih,1),vhtm(ih,2),vhtm(ih,3)
     > ,vhtm(ih,1)-vhtm(ih,2)-vhtm(ih,3)
     > ,vhta(ih,1)+vhtm(ih,1)
     > ,vhta(ih,2)+vhtm(ih,2)
     > ,vhta(ih,3)+vhtm(ih,3)
      enddo
      endif
!
!::denisty
      if( lout.ge.3 ) then
      ig = 1
      write(n6,'(/2x,"density and pressure just before albedo")')
      write(n6,'(4x,"iw",3x,"ih",3x,"chwl",3x"ic",4x,"icx",2x,"icy",4x
     > ,"d0",10x,"e0",10x,"dg",10x,"eg",10x,"pd2",9x,"ptot")')
      do i = 1, npmp
      iw = ipmp(i)
      ic = icwl(iw)
      ih = ihwl(iw)
      if( ic.le.0 ) cycle
      ze0 = 0.0d0
      zeg = 0.0d0
      if( vden0(iw,ig).gt.0.0d0 ) ze0 = vpre0(iw,ig)/vden0(iw,ig)
      if( vdeng(iw,1).gt.0.0d0  ) zeg = vpreg(iw,1)/vdeng(iw,1)
      write(n6,'(2x,2i5,2x,a,3i6,1p6e12.3)')
     >  iw,ih,chwl(ih),ic,migx(ic),migy(ic)
     > ,vden0(iw,ig),ze0
     > ,vdeng(iw,1), zeg, vpreg(iw,1)*cev
     > ,(vpre0(iw,ig)+vpreg(iw,1))*cev
      enddo
      endif
!
!-----------------------------------------------------------------------
!::pressure at wall
!-----------------------------------------------------------------------
      elseif( ktyp.eq.4 ) then
      ig = 1
      write(n6,'(2x)')
      write(n6,'(2x,"*** ntpump(4) ***")')
      write(n6,'(2x,"density and pressure at wall")')
      write(n6,'(4x,"iw",3x,"ih",3x,"chwl",3x"ic",4x,"icx",2x,"icy",4x
     > ,"d0",10x,"e0",10x,"dg",10x,"eg",10x,"pd2",9x,"ptot")')
      do nw = 1, nowl2
      iws = npsw(nw)
      iwe = npew(nw)
      do iw = iws, iwe
      ic = icwl(iw)
      ih = ihwl(iw)
      ze0 = 0.0d0
      zeg = 0.0d0
      if( vden0(iw,ig).gt.0.0d0 ) ze0 = vpre0(iw,ig)/vden0(iw,ig)
      if( vdeng(iw,1).gt.0.0d0  ) zeg = vpreg(iw,1)/vdeng(iw,1)
      write(n6,'(2x,2i5,2x,a,3i6,1p6e12.3)')
     >  iw,ih,chwl(ih),ic,migx(ic),migy(ic)
     > ,vden0(iw,ig),ze0
     > ,vdeng(iw,1), zeg, vpreg(iw,1)*cev
     > ,(vpre0(iw,ig)+vpreg(iw,1))*cev
      enddo
      enddo
      endif
!
      return
      end
