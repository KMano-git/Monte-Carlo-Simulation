!***********************************************************************
      subroutine pst_wimp
!***********************************************************************
      use catcom, only : catmz
      use cimcom, only : denflxztowall, eneflxztowall, ismax, ndis
     >    , weflx, wpflx
      use cimctl, only : nprg
      use cimden, only : nsput, nzmx
      use cntcom, only : icwl, ipwl, npsw, npew, xpnt, ypnt
     > , ncmax
      use cntwfl, only : jxwl, iywl, xare, xdeg
      use cplwrd, only : wfac
      use csize,  only : ndwp
      use csonic, only : itim, time
      use cunit,  only : lmspe, lmype, n6
      use mod_externalgrid
      implicit none
!
!::local variables
      integer  nw, iw, iws, iwe, ixy, iz
      integer  i, i6
      integer  ic, ip, ip1
      real*8   x1, y1, x2, y2, x0, y0, th
      real*8   lwfac
      real(8) :: zTotDenFlxZ(ndwp) ,zTotEneFlxZ(ndwp)
      real(8) :: wTotDenFlxZ(0:ndis) ,wTotEneFlxZ(0:ndis)
      real(8) :: wholeDenFlxZ, wholeEneFlxZ,
     &           totDenTmp, totEneTmp
      character(10) :: fnam
!
!::   assume ismax <= 74, write format set "200" real numbers
!

      zTotEneFlxZ = 0.0_8 ! for avoid NaN error
      zTotDenFlxZ = 0.0_8 ! for avoid NaN error

      if( lmype.ne.lmspe ) return

      write(n6,'(/2x,"*** pst_wimp ***"  )')
      write(n6,*) ' === denFlxZtoWall'

      denFlxZtoWall = 0.0d0
      eneFlxZtoWall = 0.0d0
      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe
          if(.not.use_exdata .and. iw==iwe) cycle
          do i = 1, nsput
            denFlxZtoWall(0:nzmx,iw)
     &          = denFlxZtoWall(0:nzmx,iw) + wpflx(0:nzmx,iw,i)
            eneFlxZtoWall(0:nzmx,iw)
     &          = eneFlxZtoWall(0:nzmx,iw) + weflx(0:nzmx,iw,i)
          enddo
          denFlxZtoWall(0:nzmx,iw) = denFlxZtoWall(0:nzmx,iw) * xare(iw)
          eneFlxZtoWall(0:nzmx,iw) = eneFlxZtoWall(0:nzmx,iw) * xare(iw)
        enddo
      enddo

      if(nprg.eq.1)then
        write(n6,'(4x,a,f7.3)')"wfac is multiplied: wfac =",wfac
        denFlxZtoWall = denFlxZtoWall * wfac
        eneFlxZtoWall = eneFlxZtoWall * wfac
        lwfac=wfac
      else
        lwfac=1.0d0
      endif

      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe-1, 10
          if( nw.eq.1 .or. nw.eq.3 ) ixy = jxwl(iw) ! jxwl is not defined for exgrid
          if( nw.eq.2 .or. nw.eq.4 ) ixy = iywl(iw) ! iywl is not defined for exgrid
          write(n6,"(4x,3i5,1p2e12.4,9e12.4)")
     &      nw, iw, ixy, xare(iw), lwfac,
     &      (denFlxZtoWall(iz,iw),iz=0,3),
     &      (eneFlxZtoWall(iz,iw),iz=0,3)
        enddo
      enddo

! Sum up over charges
      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe
          if(.not.use_exdata .and. iw==iwe) cycle
          totDenTmp=0.d0
          totEneTmp=0.d0
          do iz=0, ismax
            totDenTmp = totDenTmp + denFlxZtoWall(iz,iw)
            totEneTmp = totEneTmp + eneFlxZtoWall(iz,iw)
          enddo
          zTotDenFlxZ(iw) = totDenTmp
          zTotEneFlxZ(iw) = totEneTmp
        enddo
      enddo

      wholeDenFlxZ = SUM(zTotDenFlxZ)
      wholeEneFlxZ = SUM(zTotEneFlxZ)

! Sum up over wall
      do iz=0, ismax
        totDenTmp=0.d0
        totEneTmp=0.d0
        do nw = 1, 4
          iws = npsw(nw)
          iwe = npew(nw)
          do iw = iws, iwe
            if(.not.use_exdata .and. iw==iwe) cycle
            totDenTmp = totDenTmp + denFlxZtoWall(iz,iw)
            totEneTmp = totEneTmp + eneFlxZtoWall(iz,iw)
          enddo
        enddo
        wTotDenFlxZ(iz) = totDenTmp
        wTotEneFlxZ(iz) = totEneTmp
      enddo
      write(n6,'(a,2e12.4)') '  [pst_wimp] TOTAL(wall): den/ene =',
     &     SUM(wTotDenFlxZ),SUM(wTotEneFlxZ)

!::power density

!::write
      i6 = 921
      write(fnam,'(i5)') nprg
      fnam = adjustl( fnam )
      fnam = 'Wimp_' // trim( fnam ) // '.txt'
      open(unit=i6,file=fnam)
      write(i6,'(2x,"impurity flux to wall: ",a2,"  wfac =", f7.3,
     >   "  itim =",i8,"  time =",1pe14.6)') catmz(1),lwfac,itim, time

      do nw = 1, 4
        write(i6,"(/A2,4A6,3(A8,2x),A8,4x,200(3x,i2.2,A3,i2.2,2x))")
     &      '# ','1:nw','2:iw','3:ixy','4:ic',
     &      '5:x0','6:y0','7:th','8:xare',
     &      9,       ':Ft',0,(iz+10,        ':Fz',iz,iz=0,ismax),
     &      11+ismax,':Qt',0,(iz+10+ismax+2,':Qz',iz,iz=0,ismax)
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe
          if(.not.use_exdata .and. iw==iwe) cycle
          if( nw.eq.1 .or. nw.eq.3 ) ixy = jxwl(iw) ! jxwl is not defined for exgrid
          if( nw.eq.2 .or. nw.eq.4 ) ixy = iywl(iw) ! iywl is not defined for exgrid
          ip  = ipwl(iw)
          if(use_exdata) then
            ip1 = ipwl2(iw)
          else
            ip1 = ipwl(iw+1)
          endif
          ic  = icwl(iw)
          if(ic.le.ncmax .or. .not.use_exdata) then
            x1 = xpnt(ip)
            y1 = ypnt(ip)
            x2 = xpnt(ip1)
            y2 = ypnt(ip1)
          elseif(ic .le. ncmax+vac_ele_size)then
            x1 = vac_grid_x(ip)
            x2 = vac_grid_x(ip1)
            y1 = vac_grid_y(ip)
            y2 = vac_grid_y(ip1)
          elseif(ic .le. ncmax+vac_ele_size+pri_ele_size)then
            x1 = pri_grid_x(ip)
            x2 = pri_grid_x(ip1)
            y1 = pri_grid_y(ip)
            y2 = pri_grid_y(ip1)
          else
            x1 = vgx_EX(ip)
            x2 = vgx_EX(ip1)
            y1 = vgy_EX(ip)
            y2 = vgy_EX(ip1)
          endif
          x0  = 0.5d0*(x1+x2)
          y0  = 0.5d0*(y1+y2)
          th  = xdeg(iw)
!
          write(i6,"(2x,4i6, 2f10.5,f10.3,1pe12.4,200e12.4)")
     &           nw, iw, ixy,
     &           ic, x0, y0, th,
     &           xare(iw),
     &           zTotDenFlxZ(iw)/xare(iw),
     &           (denFlxZtoWall(iz,iw)/xare(iw),iz=0,ismax),
     &           zTotEneFlxZ(iw)/xare(iw),
     &           (eneFlxZtoWall(iz,iw)/xare(iw),iz=0,ismax)
        enddo !iw
      enddo !nw
      write(i6,"(/61x,'#total ',200e12.4)")
     &           wholeDenFlxZ,
     &           (wTotDenFlxZ(iz),iz=0,ismax),
     &           wholeEneFlxZ,
     &           (wTotEneFlxZ(iz),iz=0,ismax)
      close(i6)
!
      return
      end
