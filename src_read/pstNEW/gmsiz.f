c***********************************************************************
      subroutine gmsiz(rmi,rmx,zmi,zmx,xln,yln)
c***********************************************************************
      implicit none
!xx      include 'clattr'
      integer  ldev, ltdev, lsymfg, lpaper
      common /cmsym/ ldev, ltdev, lsymfg, lpaper

! lpapaer : papaer = 1 : Landscape,  = 2 : Portrait
! lsymfg No use
c
c::argument
      real*4  rmi, rmx, zmi, zmx, xln, yln
c
c::local variables
      real*4  xmi, xmx, ymi, ymx
      real*4  xst, yst
!
      xmi = rmi
      xmx = rmx
      ymi = zmi
      ymx = zmx
c
cx      write(6,'(2x,"gmsiz   rmi,rmx,zmi,zmx =",4f9.4)')
cx     >   rmi,rmx,zmi,zmx
c
      xst = 2.0
      yst = 3.0
      xln = 15.0
      yln = 20.0

c::gdsiz
      xst = 2.2
      yst = 1.5
      if( lpaper.eq.1 ) yln = 16.0    ! Landscape
      if( lpaper.eq.2 ) yln = 21.0    ! Portlate
!
      if( lpaper.eq.1 ) yln = 15.0
c
      xln = yln/(ymx-ymi)*(xmx-xmi)
      if( lpaper.eq.2 .and. xln.gt.16.0 ) then
      xln = 16.0
      yln = xln/(xmx-xmi)*(ymx-ymi)
      endif
c
      call gdsiz(xst,yst,xln,yln)
      call gdaxs(0,xmi,xmx,"r (m)")
      call gdays(1,ymi,ymx,"z (m)")
c
cx      call gdgtfr(xst,yst,xln,yln)
cx      write(6,'(2x,"gdgtfr  ",4f10.4)') xst,yst,xln,yln
c
      return
      end
