!***********************************************************************
      subroutine test_getwal
!***********************************************************************
!
!      write(clin(mj+1:),'(3x,a1,1x,a1,5x)') char(9), char(9)
!
!-----------------------------------------------------------------------
      use csize
      use cntcom
      use cplmet
      use cunit
      implicit none

!
      integer  nd
      parameter  (nd=1000)
      real*4  wx(nd,5),wy(nd,5)   ! wall data
!
      character clin*256
      integer  nsw, nvw, nms, nsl, npv
      integer  np, i, mj, lenx, i6
!
      call getwal("swal",nsw,wx(1,1),wy(1,1),nd)
      call getwal("vwal",nvw,wx(1,2),wy(1,2),nd)
      call getwal("mspx",nms,wx(1,3),wy(1,3),nd)
      call getwal("psol",nsl,wx(1,4),wy(1,4),nd)
      call getwal("pprv",npv,wx(1,5),wy(1,5),nd)
!
      i6 = 21
      open(unit=i6,file="@wall")
!
      write(i6,'(4x,"xsw,",6x,"ysw,",6x,"xvw,",6x,"yvw,",6x,
     >   "xsp,",6x,"ysp,",6x,"xsl,",6x,"ysl,",6x,"xpv,",6x,"ypv,")')
      np = max0(nsw,nvw,nms,nsl,npv)
!
      do i = 1, np
      clin = "  "
      mj = 2
      if( i.le.nsw ) then
        write(clin(mj+1:),'(f9.4,a1,f9.4,a1)') wx(i,1),",",wy(i,1),","
      else
        write(clin(mj+1:),'(9x,a1,9x,a1)') ",",","
      endif
      mj = mj+20
      if( i.le.nvw ) then
        write(clin(mj+1:),'(f9.4,a1,f9.4,a1)') wx(i,2),",",wy(i,2),","
      else
        write(clin(mj+1:),'(9x,a1,9x,a1)') ",",","
      endif
      mj = mj+20
      if( i.le.nms ) then
        write(clin(mj+1:),'(f9.4,a1,f9.4,a1)') wx(i,3),",",wy(i,3),","
      else
        write(clin(mj+1:),'(9x,a1,9x,a1)') ",",","
      endif
      mj = mj+20
      if( i.le.nsl ) then
        write(clin(mj+1:),'(f9.4,a1,f9.4,a1)') wx(i,4),",",wy(i,4),","
      else
        write(clin(mj+1:),'(9x,a1,9x,a1)') ",",","
      endif
      mj = mj+20
      if( i.le.npv ) then
        write(clin(mj+1:),'(f9.4,a1,f9.4,a1)') wx(i,5),",",wy(i,5),","
      else
        write(clin(mj+1:),'(9x,a1,9x,a1)') ",",","
      endif
      write(i6,'(a)') clin(1:lenx(clin))
      enddo
!
      close(i6)
      return
      end 

!***********************************************************************
      subroutine getwal(ctyp,np,px,py,nd)
!***********************************************************************
!
!       ctyp (i)  "swal", "vwal", "mspx", "psol", "pprv"
!       np   (o)  number of data
!       px   (o)  x-point
!       py   (o)  y-point
!       nd   (i)  dimension size
!
!-----------------------------------------------------------------------
      use csize
      use cntcom
      use cplmet
      use cunit
      implicit none

!
!::argument
      character ctyp*(*) 
      integer   np, nd
      real*4    px(nd), py(nd)   ! wall data
!
!::local variables
      integer  ii, k, i, j, ic, ip, n, kmax, it, jt
      integer  jxa, jxb, jxd, iya, iyb, iyd
      integer  imox, imoy, ip1
!
!-----------------------------------------------------------------------
!::sol wall
!-----------------------------------------------------------------------
      if( ctyp.eq."swal" ) then
!
      ii = 0
      do k = 1, 4
      if( k.eq.1 ) then
      jxa = jdp1;  jxb = jdp2-1;  jxd = +1
      iya = ivs1;  iyb = ivs1;    iyd = +1
      elseif( k.eq.2 ) then
      jxa = jdp2-1;  jxb = jdp2-1;  jxd = +1
      iya = ivs1;    iyb = ivs2-1;  iyd = +1
      elseif( k.eq.3 ) then
      jxa = jdp2-1;  jxb = 1;       jxd = -1
      iya = ivs2-1;  iyb = ivs2-1;  iyd = +1
      elseif( k.eq.4 ) then
      jxa = 1;       jxb = 1;       jxd = +1
      iya = ivs2-1;  iyb = 1;       iyd = -1
      endif
!
      do i = iya, iyb, iyd
      do j = jxa, jxb, jxd
      if( k.eq.3 ) then
        if( j.ge.jxp1 .and. j.le.jsl2 ) cycle
      endif
!
      ic = nocl(j,i)
      if( ic.le.0 ) cycle
!
      ip = mgrd(ic,k)
      ii = ii + 1
      if( ii.gt.nd ) goto 910
      px(ii) = xpnt(ip)
      py(ii) = ypnt(ip)
      enddo
      enddo
      enddo
      ii = ii + 1
      if( ii.gt.nd ) goto 910
      px(ii) = px(1)
      py(ii) = py(1)
      np = ii
      endif
!
!-----------------------------------------------------------------------
!::vacume vessel
!-----------------------------------------------------------------------
      if( ctyp.eq."vwal" ) then
      ii = 0
      i = ivb1
      do j = jvb1, jvb2
      ip = nogdv(j,i)
      ii = ii + 1
      if( ii.gt.nd ) goto 910
      px(ii) = xpnt(ip)
      py(ii) = ypnt(ip)
      enddo
      j = jvb2
      do i = ivb1, ivb2
      ip = nogdv(j,i)
      ii = ii + 1
      if( ii.gt.nd ) goto 910
      px(ii) = xpnt(ip)
      py(ii) = ypnt(ip)
      enddo
      i = ivb2
      do j = jvb2, jvb1, -1
      ip = nogdv(j,i)
      ii = ii + 1
      if( ii.gt.nd ) goto 910
      px(ii) = xpnt(ip)
      py(ii) = ypnt(ip)
      enddo
      j = jvb1
      do i = ivb2, ivb1, -1
      ip = nogdv(j,i)
      ii = ii + 1
      if( ii.gt.nd ) goto 910
      px(ii) = xpnt(ip)
      py(ii) = ypnt(ip)
      enddo
      np = ii
      endif
!
!-----------------------------------------------------------------------
!::separatrix
!-----------------------------------------------------------------------
      if( ctyp.eq."mspx" ) then
      ii = 0
      it = itsle
      do jt = jtmin(it), jtmax(it)
      j = jcel(jt,it)
      i = icel(jt,it)
      ic = mcel(j,i)
      if( ic.le.0 ) cycle
      ip = mgrd(ic,4)
      ii = ii + 1
      if( ii.gt.nd ) goto 910
      px(ii) = xpnt(ip)
      py(ii) = ypnt(ip)
      ip1 = mgrd(ic,3)
      enddo
      ii = ii + 1
      px(ii) = xpnt(ip1)
      py(ii) = ypnt(ip1)
      np = ii
      endif
!
!-----------------------------------------------------------------------
!::sol
!-----------------------------------------------------------------------
      if( ctyp.eq."psol" ) then
      ii = 0
      it = itsls+1
      do jt = jtmin(it), jtmax(it)
      j = jcel(jt,it)
      i = icel(jt,it)
      ic = mcel(j,i)
      if( ic.le.0 ) cycle
      ip = mgrd(ic,1)
      ii = ii + 1
      if( ii.gt.nd ) goto 910
      px(ii) = xpnt(ip)
      py(ii) = ypnt(ip)
      ip1 = mgrd(ic,2)
      enddo
      ii = ii + 1
      px(ii) = xpnt(ip1)
      py(ii) = ypnt(ip1)
      np = ii
      endif
!
!-----------------------------------------------------------------------
!::prv
!-----------------------------------------------------------------------
      if( ctyp.eq."pprv" ) then
      ii = 0
      it = itpve-1
      do jt = jtmin(it), jtmax(it)
      j = jcel(jt,it)
      i = icel(jt,it)
      ic = mcel(j,i)
      if( ic.le.0 ) cycle
      ip = mgrd(ic,4)
      ii = ii + 1
      if( ii.gt.nd ) goto 910
      px(ii) = xpnt(ip)
      py(ii) = ypnt(ip)
      ip1 = mgrd(ic,3)
      enddo
      ii = ii + 1
      px(ii) = xpnt(ip1)
      py(ii) = ypnt(ip1)
      np = ii
      endif
!
!::debug write
!xx      write(n6,'(2x,"*** getwal ***  ",a,2x,i5)') ctyp, np
!xx      write(n6,'(2x,"wx =",15f9.4)') (px(i),i=1,np)
!xx      write(n6,'(2x,"wy =",15f9.4)') (py(i),i=1,np)
!
      return
!
!::dimension error
 910  continue
      write(n6,'(/2x,"*** pltwal ***  ii.gt.nd ",2i7)') ii,nd
      call wexit("pltwal","ii.gt.nd")
      end
