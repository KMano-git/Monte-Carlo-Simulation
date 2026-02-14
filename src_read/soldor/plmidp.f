!**********************************************************************
      subroutine plmidp
!**********************************************************************
      use cgdcom, only : grdx, grdy, mqs1, npmx
      use cplcom, only : dmid, nmid, xmid, ymid
      use cplmet, only : gwtm, hgdx, hgdy, icspx, icwl1, jcxp1, jcxp2
      use csize,  only : ndy
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  ix, iy, ic, jcs, jce, jc, jcm, np, n
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real*8   xaxs, yaxs, eps, zwt, zl
      real*8   yaxs, eps, zwt, zl
      real*8   zgx(ndy),zgy(ndy)
! deleted 2 lines organize local variables and include files by kamata 2021/05/31
!ik   character creg*80
!ik   real*4 rmi,rmx,zmi,zmx   ! plot routine
      integer ldbg; parameter ( ldbg = 0 )
!
      write(n6,'(/2x,"*** plmidp ***")')
!
!::axis
      ix = mqs1
      iy = npmx
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   xaxs = grdx(ix,iy)
      yaxs = grdy(ix,iy)
      write(n6,'(2x,"axis point  xaxs =",f8.4,"  yaxs =",f8.4)')
     >  grdx(ix,iy),grdy(ix,iy)
!
!::inner-mid-plane
      ic  = icspx
      eps = 1.0e20
      jce = jcxp2
      jcs = (jcxp1+jcxp2)/2
      do jc = jcs, jce
      if( dabs(hgdy(jc,ic)-yaxs).lt.eps ) then
      eps  = dabs(hgdy(jc,ic)-yaxs)
      jcm = jc
      endif
      enddo
      write(n6,'(2x,"mid-plane (in)   xmid =",f8.4,"  ymid =",f8.4,
     >  "  jc,ic =",2i5)')  hgdx(jcm,ic),hgdy(jcm,ic),jcm,ic
!
!::outer-mid-plane
      ic  = icspx
      eps = 1.0e20
      jcs = jcxp1
      jce = (jcxp1+jcxp2)/2
      do jc = jcs, jce
      if( dabs(hgdy(jc,ic)-yaxs).lt.eps ) then
      eps  = dabs(hgdy(jc,ic)-yaxs)
      jcm = jc
      endif
      enddo
      write(n6,'(2x,"mid-plane (out)  xmid =",f8.4,"  ymid =",f8.4,
     >  "  jc,ic =",2i5)')  hgdx(jcm,ic),hgdy(jcm,ic),jcm,ic
!
!::point
      np  = 1
      zwt = gwtm(jcm,ic)
      zgx(np) = hgdx(jcm,ic) + zwt*(hgdx(jcm,ic+1)-hgdx(jcm,ic))
      zgy(np) = hgdy(jcm,ic) + zwt*(hgdy(jcm,ic+1)-hgdy(jcm,ic))
      do 120 ic = icspx, icwl1, -1
      np = np + 1
      zgx(np) = hgdx(jcm,ic)
      zgy(np) = hgdy(jcm,ic)
 120  continue
!
!::distance from separatrix
      zl = 0.0
      ic = icspx + 1
      xmid(ic) = zgx(1)
      ymid(ic) = zgy(1)
      dmid(ic) = zl
      do 210 n = 2, np
      ic = ic - 1
      zl = zl + sqrt( (zgx(n)-zgx(n-1))**2+(zgy(n)-zgy(n-1))**2 )
      xmid(ic) = zgx(n)
      ymid(ic) = zgy(n)
      dmid(ic) = zl
 210  continue
      nmid = np
!
!::debug
!x      write(n6,'(2x,"icspx =",i3,"  nmid =",i3)') icspx,nmid
!x      write(n6,'(2x,"xmid =",8f9.4)') (xmid(ic),ic=1,icspx+1)
!x      write(n6,'(2x,"ymid =",8f9.4)') (ymid(ic),ic=1,icspx+1)
      write(n6,'(2x,"dmid =",8f9.4)') (dmid(ic),ic=1,icspx+1)
!
!x      if( ldbg.gt.0 ) then
!x      call gdmod(1)
!x      creg = '3.5,4.5,-0.2,0.2'
!x      call plsiz(creg,rmi,rmx,zmi,zmx)
!x      call mggrid(2)
!x      call gdpltd(1,np,xmid,ymid,-7,' ')
!x      call gdpltd(1,1,xmid(np),ymid(np),-1,' ')
!x      call gdpagh(1)
!x      call gdmod(0)
!x      endif
!
      return
      end
