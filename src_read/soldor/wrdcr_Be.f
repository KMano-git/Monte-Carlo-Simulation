!**********************************************************************
      subroutine wrdcr_Be(cimz,wrd)
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
      real(8), intent(in)  :: cimz(10)
      real(8), intent(out) :: wrd(ndmc)


!::local variables
      character(4), save ::  atmz = "Be"
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
!KH140930  non-coronal Be rlos = 1e-18, N=1e20
      zx = zte
      zx = dmin1(zx,1.0d4)
      zx = dmax1(zx,1.0d0)
      zx = dlog10(zx)
      if(zte.le.10.0d0)then
      zy = -32.308d0 + 5.7225d0*zx -      10.763*zx**2
     >               - 73.776d0*zx**3 + 217.34d0*zx**4
     >               - 211.96d0*zx**5 + 71.284d0*zx**6
      else
      zy =  127.26d0 - 602.67d0*zx    + 926.94d0*zx**2
     >               - 772.93d0*zx**3 + 385.65d0*zx**4
     >               - 118.74d0*zx**5 + 22.151d0*zx**6
     >               - 2.2982d0*zx**7 + 0.10186d0*zx**8
      endif

      zlz = 10.0d0**zy
!
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
      if( mod(itim,100).eq.0  .and. nlp.eq.nlpmn ) ldbg = 1
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
