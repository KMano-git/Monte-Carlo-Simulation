!***********************************************************************
      subroutine plfimp
!***********************************************************************
!
!    plfimp : soldor(pl) flux(f) impurity(imp)
!
!     vsez0:  dNe = sum(z)[z*Nz] at cell center
!     vsez:   dNe = dmin1( dNe, clmdnz*Ne0 )  Ne0 = sum(a)[za*Na]
!     xzflz:  Fyz = -Dz*grad(dNe) = -Azn*dNe(ip) + Azs*dNe(i)
!     vdnz:   Nz
!
!   dzan(j,i)  : diffusion coefficient of impurity ions  (0.3 m2/s)
!   vsez(j,i)  : v(..) s(source) e(electron density)  due to z(impurity)
!   vsezg(j,i) : dNe at grid points
!   xzwtm(j,i), xzwtp(j,i) : weit of Te at (j,i+1/2)
!               Ti(j,i+1/2) = xzwtm(j,i)*vti(j,i) + xzwtp(j,i)*vti(j,ip)
!   xzflz(j,i) : partcile flux due to impurity   Fz(j,i+1/2)
!   vdnz(j,i,iz) : v(soldor non-cons. variables) dn(density) z(impurity)
!
!-----------------------------------------------------------------------
      use cmeffz, only :  dzan, vdnz, vdnz0, vsez, vsez0, vsezg, vsezz
     >    , xzflz, xzwtm, xzwtp
      use cntcom, only : mcel
      use cphcns, only : cev
      use cplcom, only : ama, aza, c12, c14, chfcv, clmdnz, mdl_fimp
     >    , nion, nlp, nlpmn, vna, vte, xfdfa
      use cplimp, only : cdifl, ismaxl, tdnzl, wmc_nty
      use cplmet, only : gdsv, gwtm, gwtp, icel, itmax, itmps, itpve
     >    , jcdp1, jcdp2, jcel, jtmax, jtmin, kce, kcn, kcs, kcw
      use cplwrd, only : wfac
      use csize,  only : ndx, ndy
      use csonic, only : itim, limp, lmstd, time
      use cunit,  only : mype
      implicit none
!
!xx   real(8), dimension(ndx,ndy,0:ndis2) :: vdnz
!xx   real(8), dimension(ndx,ndy) :: dzan, vsez, vsezg
!xx   real(8), dimension(ndx,ndy) :: xzwtm, xzwtp, xzflz
!
!::local variables
      integer :: it, jwst, jwen, jw, j, i, ic, iz, jmax
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: jp, jm, ip, j1, j2, nty
      integer :: jp, jm, ip, nty
      real(8) :: sez, sezz,za, anz, zdzan, zf11, zwtm, zwtp
      real(8) :: Aze, Azw, Azn, Azs
!
!::debug
!     itsls=1 itsle=25 itpvs=26 itpve=33 itmps=34 itmpe=43
      integer,dimension(10) ::
     >   tbit = (/25, 20, 15, 27, 34, 35, 36, 41, 42, 43/)
      integer :: ia, ift, m
      real(8) :: hfcv, hqcv, hte, hne0, hfac, fac
      character(20) :: dsn
!
      if( lmstd.eq.1 ) return
      if( limp.eq.0  ) return
!
!
!::impurity density in soldor
!::increase in electron density due to impurity
!
      do it = 2, itmax
        if( it.eq.itpve ) cycle
        jwst = jtmin(it)
        jwen = jtmax(it)
        if( it.ge.itmps ) then
        jwst = jwst + 1
        jwen = jwen - 1
        endif
!
        do jw = jwst, jwen
          j = jcel(jw,it)
          i = icel(jw,it)
          ic = mcel(j,i)
          if( ic.eq.0 ) cycle
!
          hne0 = 0.0d0
          do ia = 1, nion
            hne0 = hne0 + aza(ia)*vna(j,i,ia)
          enddo
!
          sez = 0.0d0
          sezz = 0.0d0
          do nty = 1, wmc_nty
            fac=1.0d0
            if(nty==1)fac=wfac
            do iz = 0, ismaxL(nty)
              za = dble(iz)
              anz  = tdnzL(iz,ic,nty)*fac
              sez  = sez  + za*anz
              sezz = sezz + za**2*anz
! modified 1/1 lines density of SONIC impurities to TOPICS by kamata 2022/03/12
!ik           vdnz(j,i,iz,nty) = anz
              vdnz0(j,i,iz,nty) = anz
            enddo
          enddo
!
          vsez0(j,i) = sez
          vsez(j,i)  = dmin1( sez, clmdnz*hne0 )
          hfac = 1.0d0
          if( sez.gt.0.0d0 ) hfac = vsez(j,i)/vsez0(j,i)
!
          vsezz(j,i) = sezz*hfac
          do nty = 1, wmc_nty
            vdnz(j,i,0:ismaxL(nty),nty)
! modified 1/1 lines density of SONIC impurities to TOPICS by kamata 2022/03/12
!ik  >            = vdnz(j,i,0:ismaxL(nty),nty)*hfac
     >            = vdnz0(j,i,0:ismaxL(nty),nty)*hfac
          enddo
!
        enddo  ! loop(jw)
      enddo    ! loop(it)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!:: impurity effect on convection flux
      if( mdl_fimp > 1) then  ! default is 0
!::cdif
        dzan (1:ndx,1:ndy) = cdifL(1) ! tentative set, to be developed
!::vsez at grid points (see plauxg)
        do it = 1, itmax
          jmax  = jtmax(it)
          jwst  = 1
          jwen  = jmax-1
          if( it.ge.itmps ) jwst = 2
!
          do jw = jwst, jwen
            j  = jcel(jw,it)
            jp = jcel(jw+1,it)
            i  = icel(jw,it)
            ip = i + 1
            vsezg(j,i) = c14*(vsez(j,i)+vsez(j,ip)
     >       +vsez(jp,i)+vsez(jp,ip))
          enddo   ! loop(jw)
        enddo       ! loop(it)
!
!::d-plate (j=3/2, j=jmax-1/2)   (see sub. plauxg)
! modified 4/1 lines organize local variables and include files by kamata 2021/05/31
!ik   j1 = jcdp1
!ik   j2 = jcdp2
!ik   do it = 1, itpve-1
!ik     i  = it
        do i = 1, itpve-1
          ip = i + 1
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik     vsezg(j1,  i) = c12*(vsez(j1,i)+vsez(j1,ip))
!ik     vsezg(j2-1,i) = c12*(vsez(j2,i)+vsez(j2,ip))
          vsezg(jcdp1,  i) = c12*(vsez(jcdp1,i)+vsez(jcdp1,ip))
          vsezg(jcdp2-1,i) = c12*(vsez(jcdp2,i)+vsez(jcdp2,ip))
        enddo
!
!::y-flux across cell boundary(j,i+1/2) due to impurity
        do it = 1, itmax-1
          if( it.eq.itpve ) cycle
          jmax = jtmax(it)
!
          do jw = 2, jmax-1
            j  = jcel(jw,it)
            i  = icel(jw,it)
            ip = i + 1
            jp = jcel(jw+1,it)
            jm = jcel(jw-1,it)
!
!::d-coef (Dz) at cell boundary (j,i+1/2) & (j,i-1/2)
            zdzan =1.0d0/(gwtm(j,i)/dzan(j,i)+gwtp(j,i)/dzan(j,ip))
            Aze = gdsv(j,i,kce)*zdzan
            Azw = gdsv(j,i,kcw)*zdzan
            Azn = gdsv(j,i,kcn)*zdzan
            Azs = gdsv(j,i,kcs)*zdzan
!
!::velocity across cell boundary
            zf11 = -Aze*vsezg(j,i)+Azw*vsezg(jm,i)
     >             -Azn*vsez(j,ip)+Azs*vsez(j,i)
!
!::up-wind
            zwtm = 1.0d0
            zwtp = 0.0d0
            if( zf11.lt.0.0d0 ) then
              zwtm = 0.0d0
              zwtp = 1.0d0
            endif
!
            xzwtm(j,i) = zwtm
            xzwtp(j,i) = zwtp
            xzflz(j,i) = zf11
          enddo  ! loop(jw)
        enddo    ! loop(it)
!
      else  ! mdl_fimp = 0 -> zf11=0
! modified 3/3 lines i and j are undefined (bug) by kamata 2021/08/18
!ik       xzwtm(j,i) = 1.0d0
!ik       xzwtp(j,i) = 0.0d0
!ik       xzflz(j,i) = 0.0d0
        xzwtm(1:ndx,1:ndy) = 1.0d0
        xzwtp(1:ndx,1:ndy) = 0.0d0
        xzflz(1:ndx,1:ndy) = 0.0d0
      endif
!
      return
!
!::debug write
      if( mod(itim,100).eq.0 .and. nlp.eq.nlpmn+1 ) then
      ift = 180000 + 500 + mype
!
      ia = 1
      do m = 1, 10
        it = tbit(m)
        write(dsn,'(a,i2.2,a)') "plfimp_",it,".txt"
        open(unit=ift,file=dsn)
        write(ift,'(2x,"time =",1pe12.5,"  itim =",i7,"  nlp =",i2)')
     >   time,  itim, nlp
        write(ift,'(4x,"it",5x,"j",3x,"i",3x,"vna",9x,"vsez",8x,"flz",
     >   9x,"fla",9x,"wtm",9x,"wtp",9x,"Te(i)",7x,"Te(ip)",6x,"Te(ih)",
     >   6x,"hqcv",8x,"hfcv")')
!
        hfcv = 0.0d0
        hqcv = 0.0d0
        jwst = 2
        jwen = jtmax(it)-1
!
        do jw = jwst, jwen
          j = jcel(jw,it)
          i = icel(jw,it)
          ip = i + 1
!xx       hte = xzwtm(j,i)*vte(j,i) + xzwtp(j,i)*vte(j,ip)
          hte = gwtp(j,i)*vte(j,i) + gwtm(j,i)*vte(j,ip)
          hfcv = hfcv + xzflz(j,i)
          hqcv = hqcv + chfcv*xzflz(j,i)*hte*cev
          write(ift,'(2x,i4,i6,i4,1p11e12.3)') it, j, i,
     >     vna(j,i,ia), vsez(j,i), xzflz(j,i), xfdfa(j,i,ia)/ama(ia),
     >     gwtp(j,i), gwtm(j,i), vte(j,i), vte(j,ip), hte,
     >     hqcv, hfcv
        enddo   ! loop(jw)
!
      close(ift)
      enddo     ! loop(m)
!
      endif
!
      return
      end
