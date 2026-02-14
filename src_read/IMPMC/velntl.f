!***********************************************************************
      subroutine velntl(ri,zi,fi,viv,vip,xp,yp,zp,vlx,vly,vlz)
!***********************************************************************
!
!         ri  (i)  r-position at neutralized point
!         zi  (i)  z-position at neutralized point
!         fi  (i)  fai (toroidal angle)  (0.0 < fi < 1.0)
!         viv (i)  vertical velocity   sqrt(vel**2-(b*vel)**2)
!         vip (i)  parallel velocity   b*vel
!         xp  (0)  x-position of neutral
!         yp  (0)  y-position of neutral
!         zp  (0)  z-position of neutral
!         vlx (0)  x-velocity of neutral
!         vly (0)  y-velocity of neutral
!         vlz (0)  z-velocity of neutral
!
!-----------------------------------------------------------------------
      use cntcom,    only : fnfi, tcfi, tsfi
      use com_eqdat, only : dr, dz, nr, rg, ubx, uby, ubz, zg
      implicit none
!
!::argument
      real(8), intent(in)  :: ri, zi, fi, viv, vip
      real(8), intent(out) :: xp, yp, zp, vlx, vly, vlz
!
!::local variables
      integer  ifi, ifa
      real*8   zvr2, zvz2, zvf2, zvr1, zvz1, zvf1, zvr, zvz, zvf
      real*8   csth, snth, csfi, snfi, csa, sna
!
      real*8   tx, ty, txy, tc0, tc1, tc2, tbr, tbz, tbt
      integer  itx, ity, itx1, ity1
! function
      real(8)    random
      integer    iino, i, j
      iino(i,j) = nr*(j-1) + i

!
!::(r2,z2,f2)
      ifi = int(fnfi*random(0)+1.0d0)
      zvr2 =  viv*tcfi(ifi)
      zvz2 = -viv*tsfi(ifi)
      zvf2 =  vip
!
!::(br,bz,bt)
      tx   = (ri-rg(1))/dr + 1.0d0
      ty   = (zi-zg(1))/dz + 1.0d0
      itx  = int(tx)
      ity  = int(ty)
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
!::(r1,z1,f1)
      csth = tbt
      snth = sqrt(tbr**2+tbz**2)
      zvr1 = zvr2*csth + zvf2*snth
      zvz1 = zvz2
      zvf1 = -zvr2*snth + zvf2*csth
!
!::(r,z,f)
      csfi = tbr/snth
      snfi = tbz/snth
      zvr  = zvr1*csfi - zvz1*snfi
      zvz  = zvr1*snfi + zvz1*csfi
      zvf  = zvf1
!
!::(x,y,z)
      ifa  = int(fnfi*fi+1.0d0)
      csa  = tcfi(ifa)
      sna  = tsfi(ifa)
      xp   = ri*csa
      yp   = ri*sna
      zp   = zi
      vlx  = zvr*csa - zvf*sna
      vly  = zvr*sna + zvf*csa
      vlz  = zvz
      return
      end
