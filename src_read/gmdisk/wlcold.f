!***********************************************************************
      subroutine wlcold(xst,yst,xen,yen,ncol,xcl,ycl,tcl,scl,
     >                  nsu,xsu,ysu)
!***********************************************************************
!
!    O  i4  ncol    : number of collision point
!    O  r8  xcl(2)  : x position
!    O  r8  ycl(2)  : y position
!    O  r8  tcl(2)  : time    if t.lt.1.0,  on line
!    O  r8  scl(2)  : position 
!    I  i4  nsu     : polinominal line
!    I  r8  xsu(nsu): x position of p-line
!    I  r8  ysu(nsy): y position of p-line
!
!-----------------------------------------------------------------------
      implicit none
!
!::argument
      integer ncol, nsu
      real*8  xst, yst, xen, yen, xcl(2),ycl(2),tcl(2), scl(2)
      real*8  xsu(nsu),ysu(nsu)
!
!::local variables
      real*8  vex, vey, dx, dy, px, py, aa, ss, tt, xsv, ysv, tsv, ssv
      integer is, i
!
      vex=xen-xst
      vey=yen-yst
!
      is=0
      do 100 i=1,nsu-1
      dx=xsu(i+1)-xsu(i)
      dy=ysu(i+1)-ysu(i)
      px=xst-xsu(i)
      py=yst-ysu(i)
      aa=dx*vey-dy*vex
      if(aa.eq.0.0) go to 100
      ss=(px*vey-py*vex)/aa
      tt=(px*dy -py*dx )/aa
      if(ss.lt.0.0 .or. ss.gt.1.0) go to 100
      if(ss.eq.1.0) go to 100
      if(tt.le.0.0) go to 100
      is=is+1
      xcl(is)=xsu(i)+ss*dx
      ycl(is)=ysu(i)+ss*dy
      tcl(is)=tt
      scl(is)= float(i) + ss
      if(is.eq.2) go to 120
  100 continue
      if(is.eq.0) then
      ncol=0
      return
      endif
!
  120 continue
      ncol=is
      if(ncol.eq.2 .and. tcl(2).lt.tcl(1)) then
      xsv = xcl(1)
      ysv = ycl(1)
      tsv = tcl(1)
      ssv = scl(1)
      xcl(1) = xcl(2)
      ycl(1) = ycl(2)
      tcl(1) = tcl(2)
      scl(1) = scl(2)
      xcl(2) = xsv
      ycl(2) = ysv
      tcl(2) = tsv
      scl(2) = ssv
      endif
!
      return
      end
