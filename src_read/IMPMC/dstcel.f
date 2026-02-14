!***********************************************************************
      subroutine dstset
!***********************************************************************
      use cimp_loc, only : dlbp, dldf, dlpt
      use cntcom,   only : mcel, migx, migy, ncmax
      use cplmet,   only : icel, itmax, itmps, itpve, itpvs, itsle, jcel
     >    , jtmax
      use csize,    only : ndx
      use cunit,    only : n6
      implicit none
!
      integer  ic, it, jt, jmax, icx, icy, md
      real*8   zdlbp, zdldf, zdlpt, dlnx(ndx), slbp, slbb
!
      write(n6,'(/2x,"*** dstset ***")')
!
      do ic = 1, ncmax
      call dstcel(ic,zdlbp,zdldf,zdlpt)
      dlbp(ic) = zdlbp
      dldf(ic) = zdldf
      dlpt(ic) = zdlpt
      enddo
!
!::debug write
      do it = 2, itmax-1
      if( it.eq.itpve ) cycle
      if( it.eq.2 .or. it.eq.itsle .or. it.eq.itpvs .or.
     >    it.eq.itmps ) then
      write(n6,'(2x,"it =",i3)') it
      write(n6,'(5x,"jt",3x,"it",3x,"ic",5x,"ix",3x,"iy",4x,"pitch",
     >  6x,"dlbp",7x,"slbp",7x,"dlbb",7x,"slbb",7x,"dlnx",7x,"dldf")')
!
      call gdlen(it,dlnx)
      jmax = jtmax(it)
      slbp = 0.0d0
      slbb = 0.0d0
!
      md = jmax/5
      if( md.le.0 ) md = 1
      do jt = 2, jmax-1
      icx = jcel(jt,it)
      icy = icel(jt,it)
      ic  = mcel(icx,icy)
      slbp = slbp + dlbp(ic)
      slbb = slbb + dlbp(ic)/dlpt(ic)
      if( mod(jt-2,md).eq.0 ) then
      write(n6,'(2x,2i5,i7,2i5,7f11.6)')
     >  jt,it,ic,migx(ic),migy(ic)
     > ,dlpt(ic),dlbp(ic),slbp,dlbp(ic)/dlpt(ic),slbb,dlnx(jt),dldf(ic)
      endif
      enddo
      endif
      enddo
!
      return
      end
!
!***********************************************************************
      subroutine dstcel(ic,dlbp,dldf,dlpt)
!***********************************************************************
!
!       side        : mknd(ic,k)     see monte/ntgrid.f
!
!    psi   mknd = 3
!          D---r----C
!     A    |        |              xpnt(ipmg),ypnt(ipmg)
!     |   4s        q2         A :  ipmg = mgrd(ic,1)
!     |    |        |  dldf    B :  ipmg = mgrd(ic,2)
!          A---p----B          C :  ipmg = mgrd(ic,3)
!              1  dlbp         D :  ipmg = mgrd(ic,4)
!                 ==> Bp       A :  ipmg = mgrd(ic,5)
!
!          pit = sqrt(bx^2+by^2)/sqrt(bx^2+by^2+bz^2)
!
!-----------------------------------------------------------------------
      use cntcom, only : bvx, bvy, bvz, iplx, iply, mgrd, migx, migy
     >    , mseg, next, xpnt, ypnt
      use cplmet, only : hpit, icmpe
      use cunit,  only : n6
      implicit none
!
      integer, intent(in)  :: ic
      real(8), intent(out) :: dlbp, dldf, dlpt
!
!::local variables
      integer  ja, jb, jc, jd
      real*8   xp, xq, xr, xs, yp, yq, yr, ys
      real*8   vx, vy, wx, wy, vv, ww, sn
      real*8   bx, by, bz, pt
!
      integer  ix,iy
      integer  ic2, ix2, iy2
      real*8   pt1, pt2, ptav
!
!::unclear dummy cell ?
      dlbp = 0.1d-5
      dldf = 0.1d-5
      dlpt = 1.0d0
      sn   = 1.0d0
!
!::unclear
!
      if( mseg(ic).ne.3 .and. mseg(ic).ne.4 ) then
        write(n6,'(2x,"Warning  dstcel mseg.ne.4  ",4i5)')
     >   ic, migx(ic), migy(ic), mseg(ic)
        call wexit("dstcel","mseg(ic).ne.3/4")
      endif
!
      ja = mgrd(ic,1)
      jb = mgrd(ic,2)
      jc = mgrd(ic,3)
      jd = mgrd(ic,4)
      if( mseg(ic).eq.3 ) jd = jc
!
      xp = 0.5d0*(xpnt(ja)+xpnt(jb))
      xq = 0.5d0*(xpnt(jb)+xpnt(jc))
      xr = 0.5d0*(xpnt(jc)+xpnt(jd))
      xs = 0.5d0*(xpnt(jd)+xpnt(ja))
!
      yp = 0.5d0*(ypnt(ja)+ypnt(jb))
      yq = 0.5d0*(ypnt(jb)+ypnt(jc))
      yr = 0.5d0*(ypnt(jc)+ypnt(jd))
      ys = 0.5d0*(ypnt(jd)+ypnt(ja))
!
      vx = xq - xs
      vy = yq - ys
      wx = xr - xp
      wy = yr - yp
!
      vv = dsqrt(vx*vx+vy*vy)
      ww = dsqrt(wx*wx+wy*wy)
      sn = (vx*wy-vy*wx)/(vv*ww)
!
      dlbp = vv
      dldf = ww*sn
!
      bx = bvx(ic)
      by = bvy(ic)
      bz = bvz(ic)
      pt = dsqrt(bx*bx+by*by)/dsqrt(bx*bx+by*by+bz*bz)
      dlpt = pt
!
      return
      end
!
!**********************************************************************
      subroutine gdlen(it,dlnx)
!**********************************************************************
!
!       |.*.|..x..|..x..|..x..|..x..|..x..|
!   J     1 |  2     3    j-1    j    j+1         jhm = j-1/2
!           |                 |-><--|             dxp = hdxm(jhm)
!           |                 dxp dxm             dxm = hdxp(jhp)
!                             <----->
!
!          dlnx(j) = dxp(j-1)/hpit(j-1) + dxm(j)/hpit(j)
!          dlnx(j) = (dxp(j-1)+dxm(j))/(hpit(j-1)+hpit(j))/2
!                                                         02/06/12
!
!                copy from ~/sol2d/soc/soldor/plenc.f     00/05/03
!----------------------------------------------------------------------
      use cgdcom, only : hbr, hbt, hbz
      use cplmet, only : hdxm, hdxp, icel, jcel, jtmax, jtmin
      use csize,  only : ndx
      implicit none
!
!::argument
      integer, intent(in)  :: it
      real(8), intent(out) :: dlnx(ndx)
!
!::local
      real*8    tln, tlp, zhm, zhp, zln, zlp, zh
      integer   ic, j, jsm, jsp, jt, jte, jts
!
      do 110 j = 1, ndx
      dlnx(j) = 0.0
 110  continue
!
      tlp = 0.0d0
      tln = 0.0d0
!
      jts = jtmin(it) + 1
      jte = jtmax(it)
      do jt = jts, jte
        jsp = jcel(jt,it)
        ic  = icel(jt,it)
        jsm = jcel(jt-1,it)
!-----
        zhm = sqrt(hbr(jsm,ic)**2+hbz(jsm,ic)**2)
     >     /sqrt(hbr(jsm,ic)**2+hbz(jsm,ic)**2+hbt(jsm,ic)**2)
        zhp = sqrt(hbr(jsp,ic)**2+hbz(jsp,ic)**2)
     >     /sqrt(hbr(jsp,ic)**2+hbz(jsp,ic)**2+hbt(jsp,ic)**2)
        zh  = 0.5d0*(zhm+zhp)
!-----
        zlp = hdxp(jsm,ic) + hdxm(jsp,ic)
        zln = zlp/zh
!xx     zln = hdxp(jsm,ic)/zhm + hdxm(jsp,ic)/zhp
        dlnx(jt) = zln
        tlp = tlp + zlp
        tln = tln + zln
      enddo
!
      return
      end
