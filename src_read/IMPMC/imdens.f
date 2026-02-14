!**********************************************************************
      subroutine imdens
!**********************************************************************
!
!  fin(W3g1/2_in), fout(g1/2), fabs(V1_abs), fpmp(a1_abs), fsty(rgn8)
!      fin = fout + fabs + fpmp + fsty
!
!   path for Shevron and Partition
!
!   whit(ih,k)
!    hit(1), rfl(2), path/abs(3), err
!    hit(1), rfl(2), path(3), abs(4), err  ! 2012/05/28  imwrfn
!
!----------------------------------------------------------------------
      use cimcom, only : dtimz, hrgpc, ismax, ndis, sdtmz, sflux, sitmp
     >    , sptcl, sptyc, stimp, stimz, swtot, timez, wexh, whit, wpemt
     >    , wtemt
      use cntcom, only : chwl, mrgn, ncmax2, nhwl, volm
      use cntpls, only : deni
      use csize,  only : ndmc
      use csonic, only : itim, time
      use cunit,  only : n6
      use mod_shexe, only : impmc_model
      implicit none
!
!::local variables
      integer :: nv_dummy, irg, ia, ic, iz, i, ih, ik, lpr
      character  zwl*4
      real(8) :: fv_in, fv_out, fv_abs, fv_pmp, fv_vac, fv_err
      real(8), dimension(0:10) :: dzrg, dirg, dvrg, dnrg
      real(8), dimension(30,4) :: fzwln, fzwli
      real(8), dimension(30)   :: fzexh
      real(8) :: fc, dtsim, tnptl
      real(8) :: wrd_dummy(ndmc), wci_dummy(ndmc), dnz(0:ndis,ndmc)
     > ,engz_dummy(0:ndis,ndmc)

      write(n6,'(/2x,"*** imdens ***")')

!::total weight of sputtered impurities  wtemt in imwcon
      if( impmc_model == 0 ) then
      swtot = wtemt    ! set at imwcon(1) total PE  at imwcon(0) 1PE
      sptcl = wpemt
      endif
      stimp = time
      sitmp = itim
      stimz = timeZ
      sdtmz = dtimZ
!
!::impurity density
      call calwrd("ionrec",nv_dummy,wrd_dummy,wci_dummy,dnz,engz_dummy)
!
!::particle number
      do irg = 0, 10
        dzrg(irg) = 0.0d0
        dirg(irg) = 0.0d0
        dvrg(irg) = 0.0d0
        dnrg(irg) = 0.0d0
      enddo
!
! Yamoto
      ia = 1
      do ic = 1, ncmax2
        irg = mrgn(ic)
        dirg(irg) = dirg(irg) + deni(ic,ia)*volm(ic)
        dvrg(irg) = dvrg(irg) + volm(ic)
        do iz = 0, ismax
          dzrg(irg) = dzrg(irg) + dnz(iz,ic)*volm(ic)
        enddo
      enddo
      do irg = 1, 8
        if(dvrg(irg).ne.0.0d0)then
        dnrg(irg) = dzrg(irg)/dvrg(irg)
        else
        dnrg(irg) = 0.0d0
        end if
      enddo
!
      dtsim = timeZ
      tnptl = 0.0d0
      do i = 1, 8
      tnptl = tnptl + dzrg(i)
      enddo
!
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
!
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
!
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
      if( impmc_model == 1 ) then
        lpr = 0
        zwl = chwl(ih)
        if( zwl(1:1) == "W" .and. zwl(3:3) == "g" ) lpr = 1
        if( index("gpsa",zwl(1:1)) > 0 ) lpr = 1
        if( lpr == 0 ) cycle
      endif
      write(n6,'(2x,i3,2x,a,2x,1p11e12.4)') ih, chwl(ih)
     > ,fc*fzwln(ih,1),fc*fzwln(ih,2),fc*fzwln(ih,3),fc*fzwln(ih,4)
     > ,fc*(fzwln(ih,1)-fzwln(ih,2)-fzwln(ih,3)-fzwln(ih,4))
     > ,fc*fzwli(ih,1),fc*fzwli(ih,2),fc*fzwli(ih,3),fc*fzwli(ih,4)
     > ,fc*(fzwli(ih,1)-fzwli(ih,2)-fzwli(ih,3)-fzwli(ih,4))
     > ,fc*fzexh(ih)
      enddo
      if( impmc_model == 1 ) return
!
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
      fv_vac = hrgpc(irg)
      fv_err = fv_in - (fv_out + fv_abs + fv_pmp + fv_vac)
!
      write(n6,'(2x,"conservation in Vac  in =",1pe12.4,
     >  "  out =",1pe12.4,"  abs =",1pe12.4,"  pmp =",1pe12.4,
     >  "  vac =",1pe12.4,"  err =",1pe12.4)')
     >   fv_in, fv_out, fv_abs, fv_pmp, fv_vac, fv_err
!
      return
      end

! dbg_denz and dbg_vzpz are defind here as IMPMC_TD. (moved from immont.f)
!***********************************************************************
      subroutine dbg_denz(cmsg)
!***********************************************************************
      use cimcom, only : denz, ismax, sflux, sptyc, swtot
      use cntcom, only : ncmax2
      use cunit,  only : n6
      implicit none

      character(*), intent(in) :: cmsg

!::local variables
      integer :: ic, iz
      real(8) :: zsum

      zsum = 0.0d0
      do ic = 1, ncmax2
      do iz = 0, ismax
      zsum = zsum + denz(iz,ic)
      enddo
      enddo

      write(n6,'(2x,"### dbg_denz",2x,a,2x,a,2x,
     >  "sflux, swtot, tot_denz =",1p3e12.3)')
     >  trim(sptyc), trim(cmsg), sflux, swtot, zsum

      return
      end

!***********************************************************************
      subroutine dbg_vzpz(cmsg)
!***********************************************************************
      use cimcom,    only : ismax, recz, sflux, sptyc, swtot, vzpz
      use cntcom,    only : ncmax2
      use cunit,     only : n6
      use mod_shexe, only : impmc_model
      implicit none

      character(*), intent(in) :: cmsg

!::local variables
      integer :: ic, iz
      real(8) :: zsum
      select case( impmc_model )
      case( 0 )
        zsum = 0.0d0
        do ic = 1, ncmax2
          do iz = 0, ismax
            zsum = zsum + recz(iz,ic)
          enddo
        enddo

        write(n6,'(2x,"### dbg_vzpz",2x,a,2x,a,2x,
     >    "sflux, swtot, tot_recz =",1p3e12.3)')
     >    trim(sptyc), trim(cmsg), sflux, swtot, zsum
      case( 1 )
        zsum = 0.0d0
        do ic = 1, ncmax2
          do iz = 0, ismax
            zsum = zsum + vzpZ(iz,ic)
          enddo
        enddo

        write(n6,'(2x,"### dbg_vzpz",2x,a,2x,a,2x,
     >    "sflux, swtot, tot_vzpz =",1p3e12.3)')
     >    trim(sptyc), trim(cmsg), sflux, swtot, zsum
      end select

      return
      end
