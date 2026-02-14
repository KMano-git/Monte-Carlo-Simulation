!***********************************************************************
      subroutine velion(xp,yp,zp,vlx,vly,vlz,vlb,ierr)
!***********************************************************************
!
!         xp  (i)  x-position
!         yp  (i)  y-position
!         zp  (i)  z-position
!         vlx (i)  x-velocity
!         vly (i)  y-velocity
!         vlz (i)  z-velocity
!         vlb (o)  parallel velocity   b*vel
!
!-----------------------------------------------------------------------
      use cimcom,    only : ipcm, lpcm
      use com_eqdat, only : dr, dz, nr, nz, rg, ubx, uby, ubz, zg
      use cunit,     only : lmype
      implicit none
!
!::argument
      real(8), intent(in)  :: xp, yp, zp, vlx, vly, vlz
      real(8), intent(out) :: vlb
      integer, intent(out) :: ierr
!
!::local common
      real*8  tbr, tbz, tbt, tvr, tvz, tvt
!
!::local variables
      real*8  r1, z1
      real*8  tx, ty, tc0, tc1, tc2, txy
      integer itx, ity, itx1, ity1
      real*8  csf, snf
      integer  lper
      data  lper/0/; save
! function
      integer  iino, i, j
      iino(i,j) = nr*(j-1) + i
!
!::(r,z)
      r1 = sqrt(xp**2+yp**2)
      z1 = zp
!
!::(br,bz,bt)
      tx   = (r1-rg(1))/dr + 1.0d0
      ty   = (z1-zg(1))/dz + 1.0d0
      itx  = int(tx)
      ity  = int(ty)
      if( itx.le.1 .or. itx.ge.nr .or.
     >    ity.le.1 .or. ity.ge.nz ) goto 910
      itx1 = itx + 1
      ity1 = ity + 1
      tx   = tx - itx
      ty   = ty - ity
      txy  = tx * ty
      tc0  = 1.0d0 - tx - ty + txy
      tc1  = tx - txy
      tc2  = ty - txy
!
      tbr  = tc0*ubx(iino(itx,ity))  + tc1*ubx(iino(itx1,ity))
     >     + tc2*ubx(iino(itx,ity1)) + txy*ubx(iino(itx1,ity1))
      tbz  = tc0*uby(iino(itx,ity))  + tc1*uby(iino(itx1,ity))
     >     + tc2*uby(iino(itx,ity1)) + txy*uby(iino(itx1,ity1))
      tbt  = tc0*ubz(iino(itx,ity))  + tc1*ubz(iino(itx1,ity))
     >     + tc2*ubz(iino(itx,ity1)) + txy*ubz(iino(itx1,ity1))
!
!::(vr,vz,vt)
      csf = xp/r1
      snf = yp/r1
      tvr =  vlx*csf + vly*snf
      tvt = -vlx*snf + vly*csf
      tvz =  vlz
!
!::vlb
      vlb = tvr*tbr + tvz*tbz + tvt*tbt
      ierr = 0
!
      return
!
!::error
 910  continue
      lper = lper + 1
      write(0,'(2x,"error velion lper =",i3,"  lpcm,ipcm =",i10,i8,
     > "  itx,ity =",2i6,"  r1,z1=",1p2e12.3,"  eq_reg =",1p4e11.3
     > ," at lmype =", i8)')
     >  lper,lpcm,ipcm, itx,ity, r1,z1, rg(1),rg(nr),zg(1),zg(nz)
     > ,lmype
      vlb = 0.0_8
      ierr = 1
      if( lper.gt.50 ) call wexit("velion","lper.gt.50")
!
      end
