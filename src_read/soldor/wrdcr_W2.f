!**********************************************************************
      subroutine wrdcr_W2(cimz,wrd)
!**********************************************************************
!
!        ADPAK (Te<10 eV)  ADAS (Te>30 eV)
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
      character(4), save :: atmz = "W2"

      integer :: it, jmax, lprv, jw, icx, icy, ir, ic, ldbg
      real(8) :: zne_min(10)
      real(8) :: zte, zne, zt, zne2, zx, zy, zlz, zwr, wrsm, wrgn(10)

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
!::Lz of W2
      zt  = zte
      zt  = dmax1( zt, 1.0d0 )
      zt  = dmin1( zt, 1.0d4 )
      zx  = dlog10(zt)

      zy = -31.9910d0 +zx*(0.422792d0  +zx*(7.11035d0 + zx*( -9.65662d0
     >     +zx*(4.86245d0 +zx*(-1.06964d0 + zx*0.0864512d0)))))
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
     >   atmz, (cimz(ir),ir=1,6)
      call wrintg(wrd,wrsm,wrgn)
      write(n6,'(2x,"wrd",a2,2x,1pe14.6,1x,1p12e11.3)')
     >   atmz(1:2),wrsm,(wrgn(ir),ir=1,7)
      endif

      return
      end

!
!      call fit_W2
!      stop
!      end
!
!**********************************************************************
      Subroutine fit_W2
!**********************************************************************
      implicit none

!::local variables
      integer, parameter :: ndim = 301
      real(8) :: temn, temx, dx, zx, zte, zy, zlz, zt
      integer :: np, i

      temn = 0.10d0
      temx = 1.0d5
      np = 61
      dx = dlog10(temx/temn)/dfloat(np-1)

      write(6,'(5x,"i",3x,"Te",10x,"Lz")')

      do i = 1, np
      zx  = dx*dfloat(i-1)
      zte = temn*10.0d0**zx

      zt  = zte
      zt  = dmax1( zt, 1.0d0 )
      zt  = dmin1( zt, 1.0d4 )
      zx  = dlog10(zt)

      zy = -31.9910d0 +zx*(0.422792d0  +zx*(7.11035d0 + zx*( -9.65662d0
     >     +zx*(4.86245d0 +zx*(-1.06964d0 + zx*0.0864512d0)))))
      zlz = 10.0d0**zy

      write(6,'(2x,i4,1p2e12.3)') i, zte, zlz
      enddo

      return
      end
