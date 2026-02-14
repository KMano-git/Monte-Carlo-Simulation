!***********************************************************************
      subroutine test_crspnt
!***********************************************************************
      implicit none
      real*8  x1, y1, x2, y2, x3, y3, x4, y4, xp, yp
      data  x1,y1/1.0d0,5.0d0/, x2,y2/3.0d0,9.0d0/
      data  x3,y3/1.0d0,-2.0d0/, x4,y4/3.0d0,-8.0d0/
!
      call crspnt(x1,y1,x2,y2,x3,y3,x4,y4,xp,yp)
      write(6,'(2x,"xp,yp =",2f10.5)') xp,yp
!
      stop
      end
!
!***********************************************************************
      subroutine crspnt(x1,y1,x2,y2,x3,y3,x4,y4,xp,yp)
!***********************************************************************
!
!                       ---3
!                   ----   |
!               4---       |
!            -- |          |
!          --   |          |
!        --     |          |
!      p------  1----------2      cross point of line-12 and line-34
!
!-----------------------------------------------------------------------
      implicit none
!
!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/28
!ik   real*8  x1, y1, x2, y2, x3, y3, x4, y4, xp, yp
      real(8), intent(in)  :: x1, y1, x2, y2, x3, y3, x4, y4
      real(8), intent(out) :: xp, yp
!
!::local variables
      real*8  a, b
!
      a = (y2-y1)/(x2-x1)
      b = (y3-y4)/(x3-x4)
!
      xp = (a*x1-b*x4-y1+y4)/(a-b)
      yp = (a*b*(x1-x4)-b*y1+a*y4)/(a-b)
!
      return
      end
!
!***********************************************************************
      real(8) function funroh(ic,xp,yp)
!***********************************************************************
      use cntcom, only : mgrd, mrgn, mseg,  xpnt, ypnt, ncmax
      use cpmpls, only : romp
      use mod_externalgrid
      implicit none
!
!::argument
      integer, intent(in) :: ic
      real(8), intent(in) :: xp, yp
!
!::local variables
      integer  iy, kmax, k, ic_tmp
      real*8  x(0:4), y(0:4)
      real*8  xt, yt, tp, a, b, zro1, zro2
! function
      integer    imoy
!
      if( mrgn(ic).eq.6 .or. mrgn(ic).eq.7 ) goto 100
!
!::sol/div
      funroh = 1.5d0
      return
!
!::main region
 100  continue
!
      if(ic.le.ncmax .or. .not.use_exdata) then
        kmax = mseg(ic)
        do k = 1, kmax
         x(k) = xpnt(mgrd(ic,k))
         y(k) = ypnt(mgrd(ic,k))
        enddo
      elseif(ic .le. ncmax+vac_ele_size)then
        ic_tmp = ic-ncmax
        kmax = mseg_vacume(ic_tmp)
        do k = 1, kmax
          x(k) = vac_grid_x(vac_element(ic_tmp,k))
          y(k) = vac_grid_y(vac_element(ic_tmp,k))
        enddo
      elseif(ic .le. ncmax+vac_ele_size+pri_ele_size)then
        ic_tmp = ic-(ncmax+vac_ele_size)
        kmax = mseg_pri(ic_tmp)
        do k = 1, kmax
          x(k) = pri_grid_x(pri_element(ic_tmp,k))
          y(k) = pri_grid_y(pri_element(ic_tmp,k))
        enddo
      else
        ic_tmp = ic-(ncmax+vac_ele_size+pri_ele_size)
        kmax = mseg_subdiv(ic_tmp)
        do k = 1, kmax
          x(k) = vgx_EX(subdiv_cell(ic_tmp,k))
          y(k) = vgy_EX(subdiv_cell(ic_tmp,k))
        enddo
      endif
!
!::segment 4
      if( kmax.eq.4 ) then
        a = (y(2)-y(1))/(x(2)-x(1))
        b = (y(3)-y(4))/(x(3)-x(4))
        x(0) = (a*x(1)-b*x(4)-y(1)+y(4))/(a-b)
        y(0) = (a*b*(x(1)-x(4))-b*y(1)+a*y(4))/(a-b)
!
!::xt,yt
        a = (y(4)-y(1))/(x(4)-x(1))
        b = (y(0)-yp)/(x(0)-xp)
        xt = (a*x(1)-b*xp-y(1)+yp)/(a-b)
        yt = (a*b*(x(1)-xp)-b*y(1)+a*yp)/(a-b)
!
!::tp
        if( dabs(x(1)-x(4)).ge.dabs(y(1)-y(4)) ) then
          tp = (xt-x(1))/(x(4)-x(1))
        else
          tp = (yt-y(1))/(y(4)-y(1))
        endif
!
!::segment 3
      else
        x(0) = xp + x(2)-x(1)
        y(0) = yp + y(2)-y(1)
!
        a = (y(3)-y(1))/(x(3)-x(1))
        b = (y(0)-yp)/(x(0)-xp)
        xt = (a*x(1)-b*xp-y(1)+yp)/(a-b)
        yt = (a*b*(x(1)-xp)-b*y(1)+a*yp)/(a-b)
!
        if( dabs(x(1)-x(3)).ge.dabs(y(1)-y(3)) ) then
          tp = (xt-x(1))/(x(3)-x(1))
        else
          tp = (yt-y(1))/(y(3)-y(1))
        endif
      endif
!
!::roh
      iy = imoy(ic)
      zro1 = romp(iy-1)
      zro2 = romp(iy)
      funroh = zro1 + (zro2-zro1)*tp
!
      return
      end
