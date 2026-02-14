!*********************************************************************
      subroutine imp_phyv
!*********************************************************************
!)
!)   phyv : physical variables
!)   MC scoreing (denz, temz) => (Nz, Tz)  time averaged quantities
!)   Nz = vznz,  Tz = vztz,  Wrad = vzwrd,  Wci = vzwci  c: collision
!)   Nz = tdnz,  Tz = n/a, Wrad = twrd SY
!)
!)   time averaged density and flux at rho = 0.98
!)    mypcn at tst + 0.5*dtm
!)
!) Old type  t: total c: colision  twne: Ne
!)     twrd(ndmc), tdnz(0:ndis2,ndmc), twci(ndmc), twne(ndmc)
!) New type  vz: header  z: impurity
!)     vznz(0:ndis2,ndmc), vztz(0:ndis2,ndmc),
!)     vzwrd(ndmc), vzwci(ndmc)
!)
!)  denZ(iz,ic) = denZ(iz,ic) + pemt(ip)*wght(ip)*dt0
!)  temZ(iz,ic) = temZ(iz,ic) + pemt(ip)*wght(ip)*dt0*vv(ip)
!)  Tz = temz/denz*0.5d0*amz/cev*2.0d0/3.0d0
!)  wsct(ic) = wsct(ic) - pemt(ip)*wght(ip)*dww
!)
!)--------------------------------------------------------------------
      use cimcom, only : amz, cdi, cdlz_i, cdlz_r, denz, dtscr, dtsiz
     >    , friz, hregw, ionz, ismax, phdnz, phfrz, phionz, phrecz
     >    , phthz, phvlz, phwci, phwrd, recz, sflux, sptyc, thfz, timez
     >    , vzpz, wexh, whit, wsct
      use cimden, only : eipot
      use cntcom, only : chwl, mrgn, ncmax, ncmax2, nhwl, volm
      use cntpls, only : dene, deni
      use cphcns, only : cev
      use csize,  only : ndmc, ndmis
      use cunit,  only : n6
      implicit none

! modified 3/3 lines organize local variables and include files by kamata 2021/07/04
!ik   real(8) :: fav, fwc, ftz, zsum
!ik   real(8) :: zsm_ion, zsm_lin, zsm_rec, zsm_cx, znz, zne, ptz
!ik   integer :: it, itm, irg, ia, nv, ic, iz, jsf, nvp, i, ih, ik
      real(8) :: fav, fwc, zsum
      real(8) :: zsm_ion, zsm_lin, zsm_rec, znz, zne
      integer :: irg, ia, nv, ic, iz, i, ih, ik
      real(8) :: fv_in, fv_out, fv_abs, fv_pmp, fv_vac, fv_err
      real(8), dimension(0:10) :: dzrg, dirg, dvrg, dnrg
      real(8), dimension(30,4) :: fzwln, fzwli
      real(8), dimension(30)   :: fzexh
      real(8) :: fc, dtsim, tnptl

!::check
      dtsiz = dtscr
      if( dtsiz /= dtscr ) then
        write(n6,'(/2x,"*** stop at sub. imp_phyv  dtsiz,dtscr =",
     >   1p2e12.4)') dtsiz,dtscr
        call wexit("imp_phyv","dtsiz /= dtscr")
      endif

!::time step in imp_step
      fav = 1.0d0/dtsiz

!::density [1/m3], temperature [eV]
      zsum = 0.0d0
      phdnz(0:ndmis,1:ndmc) = 0.0d0
      phfrz(0:ndmis,1:ndmc) = 0.0d0
      phthz(0:ndmis,1:ndmc) = 0.0d0
      phionz(0:ndmis,1:ndmc) = 0.0d0
      phrecz(0:ndmis,1:ndmc) = 0.0d0
      phvlz(0:ndmis,1:ndmc) = 0.0d0
      fwc = 0.5d0*amz
! deleted 1 line organize local variables and include files by kamata 2021/07/04
!ik   ftz = 0.5d0*amz/cev*2.0d0/3.0d0
      nv  = ncmax2   ! incl. vac-region
      do ic = 1, nv
        phwci(ic) = fav*fwc*wsct(ic)/volm(ic)
        do iz = 0, ismax
          phdnz(iz,ic) = fav*denz(iz,ic)/volm(ic)
          phionz(iz,ic) = fav*ionz(iz,ic)/volm(ic)
          phrecz(iz,ic) = fav*recz(iz,ic)/volm(ic)
          if(denz(iz,ic).gt.0.0d0)then
             phfrz(iz,ic) = friZ(iz,ic)/denz(iz,ic)
             phthz(iz,ic) = thfZ(iz,ic)/denz(iz,ic)
             phvlz(iz,ic) = vzpz(iz,ic)/denz(iz,ic)
          else
             phfrz(iz,ic) = 0.0d0
             phthz(iz,ic) = 0.0d0
             phvlz(iz,ic) = 0.0d0
          endif
!$$$          vztz(iz,ic) = 0.0d0
!$$$          if( vztz(iz,ic) > 0.0d0 ) then
!$$$            vztz(iz,ic) = ftz*temz(iz,ic)/denz(iz,ic)
!$$$          endif
          zsum = zsum + phdnz(iz,ic)
        enddo
      enddo
      !DBG
      write(n6,*) 'dbg: phdnz sum = ', zsum
!::flux
      write(n6,'(2x,"imp_phyv   dtsiz =",1pe12.3,"  fav =",1pe12.3)')
     >    dtsiz, fav

!$$$      do jsf = 1, 5
!$$$      do iz  = 0, ismax
!$$$        fbc_flxI(iz,jsf) = fav*fbc_flxI(iz,jsf)
!$$$        fbc_flxO(iz,jsf) = fav*fbc_flxO(iz,jsf)
!$$$      enddo
!$$$      enddo

!::wrad
      phwrd(1:ndmc) = 0.0d0
! modified 2/1 lines organize local variables and include files by kamata 2021/07/04
!ik   nvp = ncmax
!ik   do ic = 1, nvp
      do ic = 1, ncmax
        zne = dene(ic)
        zsm_ion = 0.0d0
        zsm_lin = 0.0d0
        zsm_rec = 0.0d0
! deleted 1 line organize local variables and include files by kamata 2021/07/04
!ik     zsm_cx  = 0.0d0

!::ion/lin/rec
        do iz = 0, ismax
          znz = phdnz(iz,ic)
          zsm_ion = zsm_ion + zne*znz*cdLz_i(iz,ic)
          zsm_lin = zsm_lin + znz*cdi(iz,ic)*eipot(iz)*cev
          zsm_rec = zsm_rec + zne*znz*cdLz_r(iz,ic)
!::cx
!xx       zn0 = tden0(ic,1)
!xx       zsm_cx  = zsm_cx  + zn0*znz*cdLz_cx(iz,ic)
        enddo

        phwrd(ic) = zsm_ion + zsm_lin + zsm_rec
      enddo

!::particle number
      do irg = 0, 10
        dzrg(irg) = 0.0d0
        dirg(irg) = 0.0d0
        dvrg(irg) = 0.0d0
        dnrg(irg) = 0.0d0
      enddo
! Yamoto
      ia = 1
      do ic = 1, ncmax2
        irg = mrgn(ic)
        dirg(irg) = dirg(irg) + deni(ic,ia)*volm(ic)
        dvrg(irg) = dvrg(irg) + volm(ic)
        do iz = 0, ismax
!        do iz = 1, ismax  !   2012/09/10  by Asakura
          dzrg(irg) = dzrg(irg) + phdnz(iz,ic)*volm(ic)
        enddo
      enddo
      do irg = 1, 8
        if(dvrg(irg).ne.0.0d0)then
        dnrg(irg) = dzrg(irg)/dvrg(irg)
        else
        dnrg(irg) = 0.0d0
        end if
      enddo
      dtsim = timeZ
      tnptl = 0.0d0
      do i = 1, 8
      tnptl = tnptl + dzrg(i)
      enddo
      write(n6,'(2x)')
      write(n6,'(2x,"sput = ",a,2x,"DotN =",1pe12.4,"  dtsim =",
     >  1pe12.4,"  Nzin =",1pe12.4,"  Nptl =",1pe12.4,2x,a)')
     >  sptyc, sflux, dtsim, sflux*dtsim, tnptl, "Nzpt excep is=0"
      write(n6,'(2x,"Nzpt : odv =",1pe12.4,"  sol =",1pe12.4,"  idv =",
     >  1pe12.4,"  opv =",1pe12.4,"  ipv =",1pe12.4,"  edg =",1pe12.4,
     >  "  cor =",1pe12.4,"  vac =",1pe12.4 )') (dzrg(i),i=1,8)
      write(n6,'(2x,"nzpt : odv =",1pe12.4,"  sol =",1pe12.4,"  idv =",
     >  1pe12.4,"  opv =",1pe12.4,"  ipv =",1pe12.4,"  edg =",1pe12.4,
     >  "  cor =",1pe12.4,"  vac =",1pe12.4 )') (dnrg(i),i=1,8)
      write(n6,'(2x,"dvrg : odv =",1pe12.4,"  sol =",1pe12.4,"  idv =",
     >  1pe12.4,"  opv =",1pe12.4,"  ipv =",1pe12.4,"  edg =",1pe12.4,
     >  "  cor =",1pe12.4,"  vac =",1pe12.4 )') (dvrg(i),i=1,8)
      write(n6,'(2x,"Nipt : odv =",1pe12.4,"  sol =",1pe12.4,"  idv =",
     >  1pe12.4,"  opv =",1pe12.4,"  ipv =",1pe12.4,"  edg =",1pe12.4,
     >  "  cor =",1pe12.4,"  vac =",1pe12.4 )') (dirg(i),i=1,8)
! Compression factor of impurity Yamoto
!$$$      write(n6,'(2x,"Compression: =", 1pe12.4)')
!$$$     > dnrg(2)/((dzrg(1)+dzrg(3))/(dvrg(1)+dvrg(3)))
!::flux to hit wall
      fzwln(1:30,1:4) = 0.0d0
      fzwli(1:30,1:4) = 0.0d0
      do ih = 1, nhwl+1
      do ik = 1, 4
        iz = 0
        fzwln(ih,ik) = whit(iz,ih,ik)
        do iz = 1, ismax
          fzwli(ih,ik) = fzwli(ih,ik) + whit(iz,ih,ik)
        enddo
      enddo
      enddo
!::absorption
      fzexh(1:30) = 0.0d0
      do ih = 1, nhwl+1
      do iz = 0, ismax
        fzexh(ih) = fzexh(ih) + wexh(iz,ih)
      enddo
      enddo
!
      write(n6,'(2x)')
      write(n6,'(3x,"ih",2x,"chwl",5x,
     >  "fin-n",7x,"frfl_n",6x,"fpth_n",6x,"fabs_n",6x,"ferr_n",6x,
     >  "fin-i",7x,"frfl_i",6x,"fpth_i",6x,"fabs_i",6x,"ferr_i",6x,
     >  "fexh")')
      fc = 1.0d0
      do ih = 1, nhwl+1
      write(n6,'(2x,i3,2x,a,2x,1p11e12.4)') ih, chwl(ih)
     > ,fc*fzwln(ih,1),fc*fzwln(ih,2),fc*fzwln(ih,3),fc*fzwln(ih,4)
     > ,fc*(fzwln(ih,1)-fzwln(ih,2)-fzwln(ih,3)-fzwln(ih,4))
     > ,fc*fzwli(ih,1),fc*fzwli(ih,2),fc*fzwli(ih,3),fc*fzwli(ih,4)
     > ,fc*(fzwli(ih,1)-fzwli(ih,2)-fzwli(ih,3)-fzwli(ih,4))
     > ,fc*fzexh(ih)
      enddo
!::conservation in vac-region
      fv_in  = 0.0d0
      fv_out = 0.0d0
      fv_abs = 0.0d0
      fv_pmp = 0.0d0
      fv_vac = 0.0d0
!
      do ih = 1, nhwl+1
      if( chwl(ih)(1:1).eq."W" .and. chwl(ih)(3:3).eq."g" ) then
        fv_in = fv_in + fzwln(ih,1) + fzwli(ih,1)
      elseif( chwl(ih)(1:1).eq."g" ) then
        fv_out = fv_out + fzwln(ih,1) + fzwli(ih,1)
      elseif( chwl(ih)(1:1).eq."V" .and. chwl(ih)(3:3).eq." " ) then
        fv_abs = fv_abs + fzwln(ih,3) + fzwli(ih,3)
      elseif( chwl(ih)(1:1).eq."a" ) then
        fv_pmp = fv_pmp + fzwln(ih,4) + fzwli(ih,4)  ! 3=>4
      endif
      enddo
!
      irg = 8
      fv_vac = hregw(irg)
      fv_err = fv_in - (fv_out + fv_abs + fv_pmp + fv_vac)
!
      write(n6,'(2x,"conservation in Vac  in =",1pe12.4,
     >  "  out =",1pe12.4,"  abs =",1pe12.4,"  pmp =",1pe12.4,
     >  "  vac =",1pe12.4,"  err =",1pe12.4)')
     >   fv_in, fv_out, fv_abs, fv_pmp, fv_vac, fv_err
      return
      end
