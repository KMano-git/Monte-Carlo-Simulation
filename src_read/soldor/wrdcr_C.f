!--------------------------------------------------------------------
!
!xxx    subroutine plwrdcrG
!
!       swrd = 0.0
!       do i = 1, wcr_nty
!         iz = i
!         cimz(1:10) = wcr_cnc(1:10,iz)
!         kz = wcr_ity(i)
!         if( trim(cnam) == "C"  ) call wrdcr_C(cimz,wrd)
!         if( trim(cnam) == "W"  ) call wrdcr_W(cimz,wrd)
!         if( trim(cnam) == "W2" ) call wrdcr_W2(cimz,wrd)
!         if( trim(cnam) == "Ar" ) call wrdcr_Ar(cimz,wrd)
!         else
!xxx        call wexit("plwrdcrG","no found imp-name "//trim(cnam))
!        endif
!         swrd = swrd + wrd
!       enddo
!
!**********************************************************************
      subroutine wrdcr_C(cimz,wrd)
!**********************************************************************
!
!   omit the model (mdl_cgen = 0)
!
!::impurity speces to be programed
!      prg_typ = (/"C", "W", "W2", "Ar", "Ne"/)
!
!----------------------------------------------------------------------
      use cntcom, only : mcel
      use cplcom, only : nlp, nlpmn, vne, vte
      use cplmet, only : icel, itmax, itpve, itpvs, jcel, jtmax, kreg
      use csize,  only : ndmc
      use csonic, only : itim, lpost
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/05/31
!ik   real(8) :: cimz(10), wrd(ndmc)
      real(8), intent(in)  :: cimz(10)
      real(8), intent(out) :: wrd(ndmc)
!
!::local common
! modified 2/12 lines replace all include files with module files by kamata 2021/08/18
!ik   real*8  xc(0:5),a0(5),a1(5),a2(5),a3(5),a4(5)
!ik   common /blk_Lz_carbon/ xc, a0, a1, a2, a3, a4
      real(8) :: xc(0:5) = (/ 0.588d0
     >    , 1.122d0, 4.786d0, 15.49d0, 144.5d0, 20.0d3 /)
      real(8) :: a0(5) =
     >   (/ -4.4938d-34,  3.2630d-33, -4.5050d-32,  1.3574d-31, 0.0d0 /)
      real(8) :: a1(5) =
     >   (/  2.2734d-33, -6.0158d-33,  1.6000d-32, -4.0782d-33, 0.0d0 /)
      real(8) :: a2(5) =
     >   (/ -4.1241d-33,  3.2423d-33, -3.0887d-34,  5.6257d-35, 0.0d0 /)
      real(8) :: a3(5) =
     >   (/  3.0825d-33, -3.8419d-34, -2.0692d-35, -3.4589d-37, 0.0d0 /)
      real(8) :: a4(5) =
     >   (/ -7.3691d-34,  3.0320d-35,  5.7628d-37,  7.8564d-40, 0.0d0 /)
!
!::local variables
      character(4), save :: atmz = "C"

      integer :: it, jmax, lprv, jw, icx, icy, ir, ic, ldbg
      real(8) :: zne_min(10)
      real(8) :: zte, zne, zne2, zx, zlz, zwr, wrsm, wrgn(10)

!::clear
      wrd(1:ndmc) = 0.0d0

!::flux tube
      do it = 2, itmax-1
      if( it.eq.itpve ) cycle
      jmax = jtmax(it)
      lprv = 0
      if( it.ge.itpvs .and. it.le.itpve ) lprv = 1
!
!::ne_min
      if( lprv.eq.1 ) then
      call setd( zne_min, 10, 1.0d50 )
      do jw = 2, jmax-1
      icx = jcel(jw,it)
      icy = icel(jw,it)
      ir  = kreg(icx,icy)
      zne_min(ir) = dmin1(zne_min(ir),vne(icx,icy))
      enddo
      endif
!
!::omit the model
!::impurity generation proportional to Ti for carbon
!::local impurity contamination
!::new sputtering model
!
!::cal. corona model
      do jw = 2, jmax-1
      icx = jcel(jw,it)
      icy = icel(jw,it)
!
!::ir   div(1,3) sol(2), prv(4,5), main(6), vacm(7)
      ir  = kreg(icx,icy)
      zte = vte(icx,icy)
      zne = vne(icx,icy)
!
!
!-----------------------------------
!::Lz of Carbon
      zx = zte
      if( zx.le.xc(0) ) then
        zlz = 1.0d-60
      elseif( zx.le.xc(1) ) then
        zlz = (((a4(1)*zx+a3(1))*zx+a2(1))*zx+a1(1))*zx+a0(1)
      elseif( zx.le.xc(2) ) then
        zlz = (((a4(2)*zx+a3(2))*zx+a2(2))*zx+a1(2))*zx+a0(2)
      elseif( zx.le.xc(3) ) then
        zlz = (((a4(3)*zx+a3(3))*zx+a2(3))*zx+a1(3))*zx+a0(3)
      elseif( zx.le.xc(4) ) then
        zlz = (((a4(4)*zx+a3(4))*zx+a2(4))*zx+a1(4))*zx+a0(4)
      else
        zlz = 2.0d-32
      endif

!::electron density :  Wr = Cimp*zne_nz * zne_wr * Lz(Te)
      zne2 = zne
      if( lprv.eq.1 ) zne2 = zne_min(ir)
      zwr = -cimz(ir)*zne2*zne*zlz

      ic = mcel(icx,icy)
      wrd(ic) = zwr
      enddo    ! loop (jw)
      enddo    ! loop (it)
      return

!::debug write
      ldbg = 0
      if( mod(itim,50).eq.0  .and. nlp.eq.nlpmn ) ldbg = 1
      if( lpost.eq.1 ) ldbg = 1

      if( ldbg == 1 ) then
      write(n6,'(2x,a,2x,1p10e11.3)')
     >   atmz, (cimz(ir),ir=1,6)
      call wrintg(wrd,wrsm,wrgn)
      write(n6,'(2x,"wrd",a2,2x,1pe14.6,1x,1p12e11.3)')
     >   atmz(1:2),wrsm,(wrgn(ir),ir=1,7)
      endif

      return
      end

!
! deleted 29 lines replace all include files with module files by kamata 2021/08/18
!**********************************************************************
!ik   block data bkdat_Lz_Carbon
!**********************************************************************
!
!       simplified non-cornal equilibrium
!              Carbon  ne*tau = 4.0e15
!
!----------------------------------------------------------------------
!ik   implicit none
!
!::local common  (blkdat_Lz_Carbon)
!ik   real*8  xc(0:5),a0(5),a1(5),a2(5),a3(5),a4(5)
!ik   common /blk_lz_carbon/ xc, a0, a1, a2, a3, a4
!
!ik   data xc/0.588d0, 1.122d0, 4.786d0, 15.49d0, 144.5d0, 20.0d3/
!
!ik   data a0(1),a1(1),a2(1),a3(1),a4(1)
!ik  >  /-4.4938d-34, 2.2734d-33, -4.1241d-33, 3.0825d-33, -7.3691d-34/

!ik   data a0(2),a1(2),a2(2),a3(2),a4(2)
!ik  >  /3.263d-33, -6.0158d-33, 3.2423d-33, -3.8419d-34, 3.032d-35/
!
!ik   data a0(3),a1(3),a2(3),a3(3),a4(3)
!ik  >  /-4.505d-32, 1.60d-32, -3.0887d-34, -2.0692d-35, 5.7628d-37/
!
!ik   data a0(4),a1(4),a2(4),a3(4),a4(4)
!ik  >  /1.3574d-31, -4.0782d-33, 5.6257d-35, -3.4589d-37, 7.8564d-40/
!
!ik   end
