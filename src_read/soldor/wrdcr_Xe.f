!**********************************************************************
      subroutine wrdcr_Xe(cimz,wrd)
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

!::local variables
      character(4), save ::  atmz = "Xe"
      integer :: it, jmax, lprv, jw, icx, icy, ir, ic, ldbg
      real(8) :: zne_min(10)
      real(8) :: zte, zne, zne2, zx, zy, zlz, zwr, wrsm, wrgn(10)

!
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
!-----------------------------------
!::Lz(Te) of Xe 190416 KH rlos=1e-16, ne=1e20, n0/ne=0
      zx = zte
      zx = dmin1(zx,1.0d4)
      zx = dmax1(zx,1.0d0)
      zx = dlog10(zx)
      if(zx<1.6055d0)then
        zy = -35.122d0 + 8.813d0  *zx    - 0.92995d0*zx**2
     >                 - 13.948d0 *zx**3 + 18.751d0 *zx**4
     >                 - 10.152d0 *zx**5 + 2.0015d0 *zx**6
      else
        zy = -20.601d0 - 23.689d0 *zx    + 22.407d0  *zx**2
     >                 - 10.807d0 *zx**3 + 2.8366d0  *zx**4
     >                 - 0.38908d0*zx**5 + 0.021955d0*zx**6
      endif
      zlz = 10.0d0**zy

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
     >    atmz, (cimz(ir),ir=1,6)
      call wrintg(wrd,wrsm,wrgn)
      write(n6,'(2x,"wrd",a2,2x,1pe14.6,1x,1p12e11.3)')
     >   atmz(1:2),wrsm,(wrgn(ir),ir=1,7)
      endif

      return
      end
