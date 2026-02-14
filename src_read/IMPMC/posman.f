!***********************************************************************
      subroutine posman(aro,np,xp,yp,icp)
!***********************************************************************
      use cntcom, only : mcel, mgrd, xpnt, ypnt
      use cplmet, only : icel, icspx, icaxs, itmps, jcel, jcxp1, jcxp2
     >    , jtmin
      use cpmpls, only : romp
      use csize,  only : ndx
      use cunit,  only : n6
      implicit none
!
!::argument
      integer, intent(in)  :: np
      integer, intent(out) :: icp(np)
      real(8), intent(in)  :: aro
      real(8), intent(out) :: xp(np), yp(np)
!
!::local variables
      integer  i, iro, ii, jst, jen, j, ic, ip1, ip2
      integer  nt, jt, ic0, ier, ko, ic1
      integer  kou, ierr, md
      real*8   tx(ndx), ty(ndx), ts(ndx)
      real*8   zro1, zro2, sum, ds, s
      real*8   ss, tt, x0, y0, x1, y1
      real(8)  zro
! function
      real(8)    dintp, ptoroh
      integer    imox, imoy
!
      integer :: ndmv
      integer    icmv(ndx*4), nmv
      real(8)    dtmv(ndx*4)
!
      write(n6,'(/2x,"*** posman ***")')
!
      do i = icspx, icaxs
        iro = i
        if( aro.gt.romp(i) ) goto 110
      enddo
      write(n6,'(2x,"no found iro  ",f9.5)') aro
      call wexit( 'posman', 'no found iro' )
!
 110  continue
      write(n6,'(2x,"aro =",f9.6,"  romp =",i4,f9.6,"  romp =",
     >   i4,f9.6)') aro, iro-1, romp(iro-1), iro, romp(iro)
!
      zro1 = romp(iro-1)
      zro2 = romp(iro)
!
      ii = 0
      jst = jcxp1 + 1
      jen = jcxp2 - 1
      do j = jst, jen
        ic  = mcel(j,iro)
        ip2 = mgrd(ic,4)
        ip1 = mgrd(ic,1)
        if( j.eq.jen ) then
          ip2 = mgrd(ic,3)
          ip1 = mgrd(ic,2)
        endif
        ii = ii + 1
        tx(ii) = xpnt(ip1) + 
     >    (xpnt(ip2)-xpnt(ip1))/(zro2-zro1)*(aro-zro1)
        ty(ii) = ypnt(ip1) + 
     >    (ypnt(ip2)-ypnt(ip1))/(zro2-zro1)*(aro-zro1)
      enddo
      nt = ii
!
      sum = 0.0d0
      ts(1) = sum
      do i = 2, nt
        sum = sum + sqrt( (tx(i)-tx(i-1))**2 + (ty(i)-ty(i-1))**2 )
        ts(i) = sum
      enddo
      do i = 1, nt
        ts(i) = ts(i)/sum
      enddo
!
!::position
      ds = 1.0d0/dfloat(np)
      ierr = 0
      do i = 1, np
        s = ds*(dfloat(i)-0.5d0)
        xp(i) = dintp(s,nt,ts,tx)
        yp(i) = dintp(s,nt,ts,ty)
!-----
        ic = 0
        do j = jst, jen
          ic1 = mcel(j,iro)
          call mchkin(xp(i),yp(i),ic1,kou)
          if( kou.eq.0 ) then
            ic = ic1
            cycle
          endif
        enddo
        if( ic.eq.0 ) ierr = ierr + 1
        icp(i) = ic
      enddo
!
!::debug write
      write(n6,'(2x,2x,"i",3x,"xp",8x,"yp",8x,"ro",8x,"ic",5x,"ix",
     <   3x,"iy",2x,"kou")')
      md = np/5
      if( md.le.0 ) md = 1
      do i = 1, np
        kou = 0
        ic = icp(i)
        if( ic.gt.0 ) call mchkin(xp(i),yp(i),ic,kou)
        if( mod(i,md).eq.0 .or. kou.ne.0 ) then
          zro = ptoroh(xp(i),yp(i))
          write(n6,'(2x,i4,3f10.5,i7,i5,i4,i3)')
     >     i, xp(i), yp(i), zro, ic, imox(ic), imoy(ic), kou
        endif
      enddo
      write(n6,'(2x,"ierr =",i3)') ierr
!
      if( ierr.gt.0 ) call wexit("posman","ierr > 0")
      return
      end
