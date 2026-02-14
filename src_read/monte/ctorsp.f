!**********************************************************************
      subroutine ctorsp(ic,ss,tt,xx,yy,ier)
!**********************************************************************
!
!     purpose : change position from com. space to real space
!                                      (ss,tt) ==> (xx,yy)
!
!       copy ~/pjsonic/soc/gmapp/stsolv.f       2001/10/15
!
!        ic (i) :  cell number
!        ss (i) :  s-position in computational space
!        tt (i) :  t-position in computational space
!        xx (o) :  x-position in real space
!        yy (o) :  y-position in real space
!        ier(o) :  error flag
!
!----------------------------------------------------------------------
      use cntcom, only : mgrd, mseg, xpnt, ypnt, ncmax
      use mod_externalgrid
      implicit none
!
!::argument
      integer, intent(in)  :: ic
      integer, intent(out) :: ier
      real(8), intent(in)  :: ss, tt
      real(8), intent(out) :: xx, yy
!
!::local variables
      integer np, km, ic_tmp, k
      real*8  p0, q0, px, qx, py, qy, pxy, qxy
      real*8  xc(5),yc(5)
!
!::segment of cell
      if(ic.le.ncmax .or. .not.use_exdata) then
        np = mseg(ic)
        do k = 1, np
          xc(k) = xpnt(mgrd(ic,k))
          yc(k) = ypnt(mgrd(ic,k))
        enddo
      elseif(ic .le. ncmax+vac_ele_size)then
        ic_tmp = ic-ncmax
        np = mseg_vacume(ic_tmp)
        do k = 1, np
          xc(k) = vac_grid_x(vac_element(ic_tmp,k))
          yc(k) = vac_grid_y(vac_element(ic_tmp,k))
        enddo
      elseif(ic .le. ncmax+vac_ele_size+pri_ele_size)then
        ic_tmp = ic-(ncmax+vac_ele_size)
        np = mseg_pri(ic_tmp)
        do k = 1, np
          xc(k) = pri_grid_x(pri_element(ic_tmp,k))
          yc(k) = pri_grid_y(pri_element(ic_tmp,k))
        enddo
      else
        ic_tmp = ic-(ncmax+vac_ele_size+pri_ele_size)
        np = mseg_subdiv(ic_tmp)
        do k = 1, np
          xc(k) = vgx_EX(subdiv_cell(ic_tmp,k))
          yc(k) = vgy_EX(subdiv_cell(ic_tmp,k))
        enddo
      endif
      if( np.ne.3 .and. np.ne.4 ) then
        call wexit("ctorsp","mseg .ne. 3,4")
      endif
!
!::triangle
      if( np.eq.3 ) then
        p0  = xc(1)
        q0  = yc(1)
        px  = xc(2) - xc(1)
        qx  = yc(2) - yc(1)
        py  = xc(3) - xc(1) - px*ss
        qy  = yc(3) - yc(1) - qx*ss
!
        xx  = p0 + px*ss + py*tt
        yy  = q0 + qx*ss + qy*tt
        ier = 0
!::rectangle
      else
        p0  = xc(1)
        q0  = yc(1)
        px  = xc(2) - xc(1)
        qx  = yc(2) - yc(1)
        py  = xc(4) - xc(1)
        qy  = yc(4) - yc(1)
        pxy = xc(1) + xc(3) - xc(2) - xc(4)
        qxy = yc(1) + yc(3) - yc(2) - yc(4)
!
        xx  = p0 + px*ss + py*tt + pxy*ss*tt
        yy  = q0 + qx*ss + qy*tt + qxy*ss*tt
        ier = 0
      endif
!
      return
      end
