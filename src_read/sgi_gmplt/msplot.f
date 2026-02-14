c**********************************************************************
      subroutine msplot(kk,kw)
c**********************************************************************
      use com_gmsiz
      use cplmet
      use cgdcom
      use cntcom
      use cntmnt
      use csonic
      use cunit
      implicit none
c
c::argument
      integer kk, kw
c
c::local variables
      integer  ista(3), iend(3), jsta(3), jend(3)
      integer  l, i, j, n, ist, ien, jst, jen, ii, np, lftmp
      real*4   wpen, gx(200), gy(200)
      character cno*2, cinp*80
      integer  mji, no
c
      call gmsiz(rmi,rmx,zmi,zmx,xln,yln)
      call pslfct(1.5)
c
c::ista,iend
      ista(1) = mpw1; iend(1) = mpw2; jsta(1) = mqd1; jend(1) = mqx1
      ista(2) = mpw1; iend(2) = mpax; jsta(2) = mqs1; jend(2) = mqs2
      ista(3) = mpw1; iend(3) = mpw2; jsta(3) = mqx2; jend(3) = mqd2
c
      if( kw.eq.0 ) then
      ista(1) = mpw1+1;  iend(1) = mpw2-1
      ista(2) = mpw1+1 
      ista(3) = mpw1+1;  iend(3) = mpw2-1
      endif
!
      write(n6,'(2x,"msplot  mpw1,mpw2,mpax =",3i5)') mpw1,mpw2,mpax
      write(n6,'(2x,"msplot  mqx1,mqx2      =",2i5)') mqx1,mqx2
c
c::region
      do l = 1, 3
      ist = ista(l); ien = iend(l)
      jst = jsta(l); jen = jend(l)
      write(n6,'(2x,"msplot  reg =",i2,2x,4i5)') l,jst,jen,ist,ien
c
c::p-line
      do i = ist, ien
      wpen = 0.2
      call pslfct(wpen)
c
      ii = 0
      do j = jst, jen
      ii = ii + 1
      gx(ii) = grdx(j,i)
      gy(ii) = grdy(j,i)
      enddo
      np = ii
      call gdplt(1,np,gx,gy,1,' ')
      enddo  
c
c::q-line
      do j = jst, jen
      wpen = 0.2
      if( l.eq.1 .and. j.eq.jst ) wpen = 1.5
      if( l.eq.3 .and. j.eq.jen ) wpen = 1.5
      call pslfct(wpen)
c
      ii = 0
      do i = ist, ien
      ii = ii + 1
      gx(ii) = grdx(j,i)
      gy(ii) = grdy(j,i)
      enddo
      np = ii
      call gdplt( 1,np,gx,gy,1,'  ' )
      enddo
c
      enddo
      call pslfct(1.0)
c
      if( kk.eq.1 ) then
      call gdpag(1)
      endif
c
      return
      end
