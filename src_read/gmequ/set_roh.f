!***********************************************************************
      subroutine set_roh
!***********************************************************************
      use clocal,    only : tpsi, troh, tvol
      use cntcom,    only : mcel, mgrd, mseg, xpnt, ypnt
      use com_eqdat, only : nr, nz, rg, zg
      use com_trclm, only : trxlm1, trxlm2, trylm1, trylm2
      use cphcns,    only : cpi
      use cplmet,    only : icaxs, icspx, jcxp1, jcxp2
      use cpmpls,    only : vlmp
      use csize,     only : ndy
      use cunit,     only : n6
      implicit none
!
!::local common
! deleted 2 lines replace all include files with module files by kamata 2021/08/18
!ik   real*8  tvol(ndy), tpsi(ndy), troh(ndy)
!ik   common /com_set_roh/ tvol, tpsi, troh
!
!::local variables
      integer  jst, jen, i, j, ic, ip1, ip2, ii
      real*8   zvl,  zp1, zp2, zpv
      real*8   xp, yp, sp, sx, sy, sxy, sxx, syy
      real*8   zpmn(ndy), zpmx(ndy)
!
      write(n6,'(/2x,"*** set_roh ***")')
!
!-----------------------------------------------------------------------
!::eqinit
!-----------------------------------------------------------------------
!
!::spline
      call splset
      write(n6,'(2x)')
!
!::trac-limit
      trxlm1 = rg(2)
      trxlm2 = rg(nr-1)
      trylm1 = zg(2)
      trylm2 = zg(nz-1)
!
!-----------------------------------------------------------------------
!::volume
!-----------------------------------------------------------------------
      jst = jcxp1 + 1
      jen = jcxp2 - 1
!
      do i = icspx, icaxs
      zvl = 0.0d0
      zp1 = -1.0d20
      zp2 = +1.0d20
      zpv = 0.0d0
      ii  = 0
!
      do j = jst, jen
      ic = mcel(j,i)
      ip1 = mgrd(ic,4)
      ip2 = mgrd(ic,3)
      if( mseg(ic).eq.3 ) ip1 = mgrd(ic,3)
      zvl = zvl+0.5d0*(xpnt(ip1)**2+xpnt(ip2)**2)*(ypnt(ip2)-ypnt(ip1))
!
      xp = xpnt(ip1)
      yp = ypnt(ip1)
      call spln2( xp, yp, sp, sx, sy, sxy, sxx, syy )
      zp1 = dmax1(zp1,sp)
      zp2 = dmin1(zp2,sp)
      zpv = zpv + sp
      ii = ii + 1
!
      enddo
      zvl = cpi*zvl
      zpv = zpv/dfloat(ii)
!
      tvol(i) = zvl
      tpsi(i) = zpv
      zpmx(i) = zp1
      zpmn(i) = zp2
      enddo
!
!::roh
      zvl = tvol(icspx)
      do i = icspx, icaxs
      troh(i) = sqrt(tvol(i)/zvl)
      enddo
!
!::debug write
      write(n6,'(2x,2x,"i",3x,"vlmp",8x,"vol",9x,"roh",9x,"psi",9x,
     > "pmx",9x,"pmn")')
      do i = icspx, icaxs
      write(n6,'(2x,i3,1p6e12.3)')
     >  i, vlmp(i), tvol(i), troh(i), tpsi(i), zpmx(i), zpmn(i)
      enddo
!
      return
      end
!
!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   function ptoroh(xp,yp)
      real(8) function ptoroh(xp,yp)
!***********************************************************************
      use clocal, only : troh, tpsi
      use cplmet, only : icaxs, icspx
      use csize,  only : ndy
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   real*8  xp, yp, ptoroh
      real(8), intent(in)  :: xp, yp
!
!::local common
! deleted 2 lines replace all include files with module files by kamata 2021/08/18
!ik   real*8  tvol(ndy), tpsi(ndy), troh(ndy)
!ik   common /com_set_roh/ tvol, tpsi, troh
!
!::local variavles
      integer  i, iy
      real*8  sp, sx, sy, sxy, sxx, syy
!
      call spln2( xp, yp, sp, sx, sy, sxy, sxx, syy )
!
      do i = icspx, icaxs-1
      iy = i
      if( (sp-tpsi(i))*(sp-tpsi(i+1)).le.0.0d0 ) goto 120
      enddo
!
      ptoroh = 1.5d0
      return
!
 120  continue
      ptoroh = troh(iy) + (troh(iy+1)-troh(iy))/(tpsi(iy+1)-tpsi(iy))
     >      *(sp-tpsi(iy))
      return
      end
