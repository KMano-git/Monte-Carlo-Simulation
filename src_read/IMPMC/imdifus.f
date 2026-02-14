!***********************************************************************
      subroutine imdifus(ip,dtunt)
!***********************************************************************
      use cimcns,    only : fgau, gaus0
      use cimcom,    only : amz, dfcf, ien, ir, is, rr, vv, vz, zz
     > , dtimZ, v_imp_conv
      use cimptl,    only : wrr, wzz
      use cntcom,    only : mgrd, mseg, xpnt, ypnt
      use com_eqdat, only : dr, dz, nr, nz, rg, ubx, uby, zg
      use cphcns,    only : cev
      use cunit,     only : mype, n6
      implicit none
!
!::argument
      integer, intent(in) :: ip
      real(8), intent(in) :: dtunt
!
!::local variables
      real*8   tx, ty, tc0, tc1, tc2, txy, tbx, tby
      integer  itx, ity, itx1, ity1
      real*8   tbb, bnx, bny, dln, zdf, conv_diff, bdir
      integer  i, j, ix, ic, k
!
      integer iino
! function
      integer    imox, imoy
      real(8)    random
!
      iino(i,j) = nr*(j-1) + i
!
      ic  = ir(ip)   !  imdifus  2008/10/01
      wrr = rr(ip)
      wzz = zz(ip)
!
      tx  = (wrr-rg(1))/dr + 1.0d0
      ty  = (wzz-zg(1))/dz + 1.0d0
!-----
      if( wrr.le.rg(1) .or. wrr.ge.rg(nr) .or.
     >    wzz.le.zg(1) .or. wzz.ge.zg(nz) ) goto 900
!-----
      itx  = int(tx)
      ity  = int(ty)
      itx1 = itx + 1
      ity1 = ity + 1
      tx   = tx - dfloat(itx)
      ty   = ty - dfloat(ity)
      txy  = tx * ty
      tc0  = 1.0d0 - tx - ty + txy
      tc1  = tx - txy
      tc2  = ty - txy
      tbx  = tc0*ubx(iino(itx,ity))  + tc1*ubx(iino(itx1,ity))
     >     + tc2*ubx(iino(itx,ity1)) + txy*ubx(iino(itx1,ity1))
      tby  = tc0*uby(iino(itx,ity))  + tc1*uby(iino(itx1,ity))
     >     + tc2*uby(iino(itx,ity1)) + txy*uby(iino(itx1,ity1))
!xx   tbz  = tc0*ubz(iino(itx,ity))  + tc1*ubz(iino(itx1,ity))
!xx  >     + tc2*ubz(iino(itx,ity1)) + txy*ubz(iino(itx1,ity1))
!
      tbb  = dsqrt(tbx*tbx+tby*tby)
      bnx  = tbx/tbb
      bny  = tby/tbb
!
      zdf  = dsqrt(2.0d0*dtunt*dfcf(is(ip),ic))   !  imdifus_new
      ix   = int(fgau*random(0) + 1.0d0)
      call bfdir_light(bdir)
      conv_diff = -1.0d0*dsign(1.0d0,bdir)*v_imp_conv(is(ip),ic)*dtimZ ! term of convection diffusion. v_imp_conv>0:move to core from SOL. (default:v_imp_conv = 0)
      dln  = zdf*gaus0(ix)+conv_diff
!
      wrr = wrr - dln*bny
      wzz = wzz + dln*bnx
!
      rr(ip) = wrr
      zz(ip) = wzz
      return
!
 900  continue
      ic = ir(ip)
      write(n6,'(2x,"error imdifus  mype =",i5,"  ip =",i6,"  is =",
     >  i3,"  ic =",i7,"  ix,iy =",2i5,"  rr,zz =",1p2e10.2,
     >  "  vv,vz =", 1p2e10.2,"  Tz =",1pe10.2)')
     >   mype, ip, is(ip), ic, imox(ic), imoy(ic),
     >   rr(ip), zz(ip), vv(ip), vz(ip), 0.5d0*amz*vv(ip)/cev
      write(n6,'(2x,"posit ",2x,2f14.8)') rr(ip), zz(ip)
      do k = 1, mseg(ic)+1
      write(n6,'(2x,"gpoint",2x,2f14.8)')
     >   xpnt(mgrd(ic,k)), ypnt(mgrd(ic,k))
      enddo
      ien(ip) = 9
      return
      end
