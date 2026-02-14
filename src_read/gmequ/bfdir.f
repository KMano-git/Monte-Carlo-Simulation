!***********************************************************************
      subroutine bfdir(bdir)
!***********************************************************************
      use cgdcom, only : grdx, grdy
      use cntcom, only : mcel
      use cplmet, only : icel, itsle, jcel, jtmax
      use cunit,  only : n6
      implicit none
!
!::argument
      real(8), intent(out) :: bdir

      integer  jt, jte, j, i, ic
      integer  mx1, mx2, mx3, mx4, my1, my2, my3, my4
      real*8   x1, y1, x2, y2, zln, esx, esy, zlb, ebx, eby, sgn
      real*8   bx, by, bz

!::local variables
!
      write(n6,'(/2x,"*** bfdir ***")')
      jte = jtmax(itsle)
      do jt = 2, jte-1
      j  = jcel(jt,itsle)
      i  = icel(jt,itsle)
      call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
!
      x1 = grdx(mx1,my1)
      y1 = grdy(mx1,my1)
      x2 = grdx(mx2,my2)
      y2 = grdy(mx2,my2)
      zln = dsqrt( (x2-x1)**2+(y2-y1)**2 )
      esx = (x2-x1)/zln
      esy = (y2-y1)/zln
!
      call mfdir(x1,y1,bx,by,bz)
      zlb = dsqrt( bx**2+by**2 )
      ebx = bx/zlb
      eby = by/zlb
!
      sgn = esx*ebx+esy*eby
      if( jt.eq.2 ) bdir = sgn
!
      ic = mcel(j,i)
      write(n6,'(2x,2i5,i7,6f10.5)')
     >  j, i, ic, x1, y1, bx, by, bz, sgn
      enddo
!
      return
      end
!
!***********************************************************************
      subroutine mfdir(xp,yp,bx,by,bz)
!***********************************************************************
!
!       xp : R-poistion
!       yp : Z-position
!       bx : Br/B
!       by : Bz/B
!       bz : Bf/B  (toroidal direction)
!
!----------------------------------------------------------------------
      use com_eqdat, only : dr, dz, nr, rg, ubx, uby, ubz, zg
      implicit none
!
!::argument
      real(8), intent(in)  :: xp, yp
      real(8), intent(out) :: bx, by, bz
!
!::local variables
      real*8   tx, ty, txy, tc0, tc1, tc2
      real*8   tbx, tby, tbz
      integer  itx, ity, itx1, ity1
!
      integer iino
      integer ir, iz
      iino(ir,iz) = nr*(iz-1) + ir
!
      tx   = (xp-rg(1))/dr + 1.0d0
      ty   = (yp-zg(1))/dz + 1.0d0
      itx  = int(tx)
      ity  = int(ty)
      itx1 = itx + 1
      ity1 = ity + 1
      tx   = tx - dble(itx)
      ty   = ty - dble(ity)
      txy  = tx * ty
      tc0  = 1.0d0 - tx - ty + txy
      tc1  = tx - txy
      tc2  = ty - txy
      tbx  = tc0*ubx(iino(itx,ity))  + tc1*ubx(iino(itx1,ity))
     >     + tc2*ubx(iino(itx,ity1)) + txy*ubx(iino(itx1,ity1))
      tby  = tc0*uby(iino(itx,ity))  + tc1*uby(iino(itx1,ity))
     >     + tc2*uby(iino(itx,ity1)) + txy*uby(iino(itx1,ity1))
      tbz  = tc0*ubz(iino(itx,ity))  + tc1*ubz(iino(itx1,ity))
     >     + tc2*ubz(iino(itx,ity1)) + txy*ubz(iino(itx1,ity1))
!
      bx = tbx
      by = tby
      bz = tbz
!
      return
      end
