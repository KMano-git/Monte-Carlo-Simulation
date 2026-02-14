!***********************************************************************
      subroutine mchkin(xp,yp,ic,ko)
!***********************************************************************
!
!   ic : cell number
!   ko : flag of in(0)/out(1)    internal/external point in cell
!
!-----------------------------------------------------------------------
      use cntcom, only : mgrd, mseg, xpnt, ypnt, ncmax
      use mod_externalgrid
      implicit none
!
!::argument
      real(8), intent(in)  :: xp, yp
      integer, intent(in)  :: ic
      integer, intent(out) :: ko
!
!::local variables
      integer kmax, k, ic_tmp, k1
      real*8  x1, y1, x2, y2, ds, ds_save
      integer, parameter :: grid_max = 4
      real*8  x1_list(grid_max),x2_list(grid_max)
     > ,y1_list(grid_max),y2_list(grid_max)
!
      if( ic.le.0 ) goto 120
!
      if(ic.le.ncmax .or. (.not. use_exdata))then
          kmax = mseg(ic)
          if(kmax > grid_max) call wexit("mchkin","kmax>4")
          do k = 1, kmax
            x1_list(k) = xpnt(mgrd(ic,k))
            x2_list(k) = xpnt(mgrd(ic,k+1))
            y1_list(k) = ypnt(mgrd(ic,k))
            y2_list(k) = ypnt(mgrd(ic,k+1))
          enddo
        elseif(ic.le. ncmax+vac_ele_size)then
          ic_tmp = ic-ncmax
          kmax = mseg_vacume(ic_tmp)
          if(kmax > grid_max) call wexit("mchkin","kmax>4")
          do k = 1,kmax
            k1 = k+1
            if( k.eq.kmax ) k1 = 1
            x1_list(k) = vac_grid_x(vac_element(ic_tmp,k))
            y1_list(k) = vac_grid_y(vac_element(ic_tmp,k))
            x2_list(k) = vac_grid_x(vac_element(ic_tmp,k1))
            y2_list(k) = vac_grid_y(vac_element(ic_tmp,k1))
          enddo
        elseif(ic .le. ncmax+vac_ele_size+pri_ele_size) then ! private
          ic_tmp = ic-(ncmax+vac_ele_size)
          kmax = mseg_pri(ic_tmp)
          if(kmax > grid_max) call wexit("fndway","kmax>4")
          do k = 1,kmax
            k1 = k+1
            if( k.eq.kmax ) k1 = 1
            x1_list(k) = pri_grid_x(pri_element(ic_tmp,k))
            y1_list(k) = pri_grid_y(pri_element(ic_tmp,k))
            x2_list(k) = pri_grid_x(pri_element(ic_tmp,k1))
            y2_list(k) = pri_grid_y(pri_element(ic_tmp,k1))
          enddo
        else
          ic_tmp = ic-(ncmax+vac_ele_size+pri_ele_size)
          kmax = mseg_subdiv(ic_tmp)
          if(kmax > grid_max) call wexit("fndway","kmax>4")
          do k = 1,kmax
            k1 = k+1
            if( k.eq.kmax ) k1 = 1
            x1_list(k) = vgx_EX(subdiv_cell(ic_tmp,k))
            y1_list(k) = vgy_EX(subdiv_cell(ic_tmp,k))
            x2_list(k) = vgx_EX(subdiv_cell(ic_tmp,k1))
            y2_list(k) = vgy_EX(subdiv_cell(ic_tmp,k1))
          enddo
        endif

      do k = 1, kmax
        x1 = x1_list(k)
        y1 = y1_list(k)
        x2 = x2_list(k)
        y2 = y2_list(k)
        ds = (x1-xp)*(y2-yp)-(x2-xp)*(y1-yp) ! same as ds = (x2-x1)*(yp-y1)-(y2-y1)*(xp-x1)
        if(k.eq.1) then
          ds_save = ds
        else
          if( ds_save*ds.lt.0.0d0 ) goto 120
        endif
      enddo
!
      ko = 0
      return
!
 120  continue
      ko = 1
      return
      end