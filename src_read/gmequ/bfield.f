!***********************************************************************
      subroutine bfield
!***********************************************************************
      use com_eqdat, only : dr, dz, nr, nz, psi, rbt0, rg, ubx, uby, ubz
     >    , zg
      use cunit,     only : n6
      implicit none
!-----------------------------------------------------------------------
!     common /com_eqdat/ rg, zg, psi, nr, nz
!     common /com_eqdt2/ dr, dz, rbt0,
!    >    raxs, zaxs, paxs, rsep, zsep, psep
!
!     common /com_unitb/ ubx, uby, ubz
!
!     br = -1.0d0/r*sy
!     bz =  1.0d0/r*sx
!     bt =  rbt0/r         (see sub. mfangl)
!
!-----------------------------------------------------------------------
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer  i, j, ir, iz, ii, ii0
      integer  i, j, ir, iz, ii
      real*8   br, bz, bt, bb
! function
      integer iino
      iino(i,j) = nr*(j-1) + i
!
      write(n6,'(/2x,"*** bfield ***")')
!
      do ir = 2, nr-1
      do iz = 2, nz-1
      ii = iino(ir,iz)
      br = -(psi(iino(ir,iz+1))-psi(iino(ir,iz-1)))/(2.0d0*dz*rg(ir))
      bz =  (psi(iino(ir+1,iz))-psi(iino(ir-1,iz)))/(2.0d0*dr*rg(ir))
      bt =  rbt0/rg(ir)
      bb = dsqrt(br**2+bz**2+bt**2)
      ubx(ii) = br/bb
      uby(ii) = bz/bb
      ubz(ii) = bt/bb
      enddo
      enddo
!
!::BC-value
      call bfbcon
!
!::debug write
      iz = nz/2
      write(n6,'(2x,"zz  =",10f11.4)') zg(iz)
      write(n6,'(2x,"rr  =",10f11.4/(2x,5x,10f11.4))') (rg(ir),ir=1,nr)
      write(n6,'(2x,"ubx =",10f11.4/(2x,5x,10f11.4))')
     >    (ubx(iino(ir,iz)),ir=1,nr)
      write(n6,'(2x,"uby =",10f11.4/(2x,5x,10f11.4))')
     >    (uby(iino(ir,iz)),ir=1,nr)
      write(n6,'(2x,"ubz =",10f11.4/(2x,5x,10f11.4))')
     >    (ubz(iino(ir,iz)),ir=1,nr)
!
      return
      end
!
!xx      call test_bfbcon
!xx      stop
!xx      end
!
!***********************************************************************
      subroutine test_bfbcon
!***********************************************************************
      use com_eqdat, only : ndr, ndz, nr, nz, ubx
      implicit none
!
! function
      integer :: iino
      integer :: i, j
      iino(i,j) = nr*(j-1) + i
!
!::local variables
      integer :: nmax, ir, iz, ii
!
      nr = 129
      nz = 257
!
      nmax = ndr*ndz
      ubx(1:nmax) = -3.1415
!
      do ir = 2, nr-1
      do iz = 2, nz-1
      ii = iino(ir,iz)
      ubx(ii) = 2.718
      enddo
      enddo
!
      call bfbcon
!
      stop
      end
!
!***********************************************************************
      subroutine bfbcon
!***********************************************************************
      use com_eqdat, only : ndr, ndz, nr, nz, ubx, uby, ubz
      use cunit,     only : n6
      implicit none
!
!::local variables
      integer :: nmax
      integer :: ir, iz, ir0, iz0, ii, ii0, nbc, ner, k
! function
      integer :: iino
      integer :: i, j
      iino(i,j) = nr*(j-1) + i
!
!
      nmax = ndr*ndz
!
      write(n6,'(2x,"*** bfbcon ***  define BC-value")')
      write(n6,'(2x,"nr,ndr =",2i4,"  nz,ndz =",2i4,"  ndr*ndz =",i6)')
     >  nr, ndr, nz, ndz, nmax
!
      nbc = 0
      do k  = 1, 2
      if( k.eq.1 ) ir =  1
      if( k.eq.2 ) ir = nr
      do iz = 1, nz
      if( ir.eq.1  ) ir0 = 2
      if( ir.eq.nr ) ir0 = nr-1
      if( iz.eq.1  ) iz0 = 2
      if( iz.eq.nz ) iz0 = nz-1
      ii  = iino(ir,iz)
      ii0 = iino(ir0,iz0)
      ubx(ii) = ubx(ii0)
      uby(ii) = uby(ii0)
      ubz(ii) = ubz(ii0)
      nbc = nbc + 1
!xx      write(n6,'(2x,"def ubx_bc ",
!xx     >    2i5,2x,1pe12.3,2x,2i5,2x,i6)') ir, iz, ubx(ii), ir0, iz0, nbc
      enddo
      enddo
!
      do k  = 1, 2
      if( k.eq.1 ) iz =  1
      if( k.eq.2 ) iz = nz
      do ir = 1, nr
      if( ir.eq.1  ) ir0 = 2
      if( ir.eq.nr ) ir0 = nr-1
      if( iz.eq.1  ) iz0 = 2
      if( iz.eq.nz ) iz0 = nz-1
      ii  = iino(ir,iz)
      ii0 = iino(ir0,iz0)
      ubx(ii) = ubx(ii0)
      uby(ii) = uby(ii0)
      ubz(ii) = ubz(ii0)
      nbc = nbc + 1
!xx      write(n6,'(2x,"def ubx_bc ",
!xx     >    2i5,2x,1pe12.3,2x,2i5,2x,i6)') ir, iz, ubx(ii), ir0, iz0, nbc
      enddo
      enddo
      write(n6,'(2x,"nbc =",i6)') nbc
!
      return
!-----------------------------------------------------------------------
!
!::check  NO pass
      ner = 0
      do ir = 1, ndr
      do iz = 1, ndz
      ii = iino(ir,iz)
      if( ubx(ii).lt.0.0d0 .or. uby(ii).lt.0.0d0 .or.
     >    ubz(ii).lt.0.0d0 ) then
        ner = ner + 1
        if( mod(ner,50).eq.0 ) write(n6,'(2x,"ubx < 0 ",
     >    i8,2i5,2x,1pe12.3,2x,i6)') ii, ir, iz, ubx(ii), ner
      endif
      enddo
      enddo
      write(n6,'(2x,"ner =",i6)') ner
!
      return
      end
