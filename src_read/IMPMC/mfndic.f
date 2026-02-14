!***********************************************************************
      subroutine mfndic(xp,yp,ic,ko)
!***********************************************************************
!
!   ic : cell number
!   ko : flag of in(0)/out(1)    internal/external point in cell
!
!-----------------------------------------------------------------------
      use cimcom, only : i6_trac, ipcm, isw2tr, iswtr, lpcm
      use cntcom, only : mgrd, mseg, xpnt, ypnt, ncmax
      use cunit,  only : n6
      use mod_externalgrid
      implicit none
!
!::argument
      real(8), intent(in)    :: xp, yp
      integer, intent(inout) :: ic
      integer, intent(out)   :: ko
!
!::local variables
      integer  kmax, k, icen, ic_tmp
      real*8   xc, yc
!
      integer, parameter :: ndmv = 200
      integer :: nmv
      integer, dimension(ndmv) :: icmv
      real*8,  dimension(ndmv) :: dtmv
!
      integer, parameter :: grid_max = 4
      real*8  x1_list(grid_max),y1_list(grid_max)
!
! function
      integer    imox, imoy
!
      if( ic.eq.0 ) then
        call wexit("mfndic","ic.eq.0")
      endif
!
      call mchkin(xp,yp,ic,ko)
      if( ko.eq.0 ) return
!
!::center of cell
      xc = 0.0d0
      yc = 0.0d0
      if(ic.le.ncmax .or. (.not. use_exdata))then
          kmax = mseg(ic)
          if(kmax > grid_max) call wexit("mfndic","kmax>4")
          do k = 1, kmax
            x1_list(k) = xpnt(mgrd(ic,k))
            y1_list(k) = ypnt(mgrd(ic,k))
          enddo
      elseif(ic.le. ncmax+vac_ele_size)then
          ic_tmp = ic-ncmax
          kmax = mseg_vacume(ic_tmp)
          if(kmax > grid_max) call wexit("mfndic","kmax>4")
          do k = 1,kmax
            x1_list(k) = vac_grid_x(vac_element(ic_tmp,k))
            y1_list(k) = vac_grid_y(vac_element(ic_tmp,k))
          enddo
      elseif(ic .le. ncmax+vac_ele_size+pri_ele_size) then ! private
          ic_tmp = ic-(ncmax+vac_ele_size)
          kmax = mseg_pri(ic_tmp)
          if(kmax > grid_max) call wexit("mfndic","kmax>4")
          do k = 1,kmax
            x1_list(k) = pri_grid_x(pri_element(ic_tmp,k))
            y1_list(k) = pri_grid_y(pri_element(ic_tmp,k))
          enddo
      else
          ic_tmp = ic-(ncmax+vac_ele_size+pri_ele_size)
          kmax = mseg_subdiv(ic_tmp)
          if(kmax > grid_max) call wexit("mfndic","kmax>4")
          do k = 1,kmax
            x1_list(k) = vgx_EX(subdiv_cell(ic_tmp,k))
            y1_list(k) = vgy_EX(subdiv_cell(ic_tmp,k))
          enddo
      endif
!
      do k = 1, kmax
        xc = xc + x1_list(k)
        yc = yc + y1_list(k)
      enddo
      xc = xc/dfloat(kmax)
      yc = yc/dfloat(kmax)
!-----
      call fndway(xc,yc,xp,yp,ic,icen,ndmv,nmv,icmv,dtmv)
!-----
!::hit wall
      if( icen.lt.0 ) then
        ko = icen
        return
      endif
!
      if( iswtr.gt.0 .and. isw2tr.gt.0 ) then
        write(i6_trac,'(2x,"fndw3",1p2e14.6,i6,i7,i6,i5,i3,i7)')
     >    xp, yp, icen, imox(icen), imoy(icen), nmv, icmv(nmv)
      endif
!
      if( nmv.eq.0 ) then
        call wexit("mfndic","nmv = 0")
      endif
!
      ic = icmv(nmv)
      call mchkin(xp,yp,ic,ko)

      if( ko.ne.0 )then
        write(n6,'(/2x,"Warning  mfndic after edit  lpcm =",i10,
     >    "  ipcm =",i7,"  ic =",i7,"  ix,iy =",i5,i4)')
     >    lpcm, ipcm, ic, imox(ic), imoy(ic)
        write(n6,'(2x,"point",2x,2f15.8)') xp, yp
!
        x1_list = 0.0d0
        y1_list = 0.0d0
        if(ic.le.ncmax .or. (.not. use_exdata))then
          kmax = mseg(ic)
          do k = 1, kmax
            x1_list(k) = xpnt(mgrd(ic,k))
            y1_list(k) = ypnt(mgrd(ic,k))
          enddo
        elseif(ic.le. ncmax+vac_ele_size)then
          ic_tmp = ic-ncmax
          kmax = mseg_vacume(ic_tmp)
          do k = 1,kmax
            x1_list(k) = vac_grid_x(vac_element(ic_tmp,k))
            y1_list(k) = vac_grid_y(vac_element(ic_tmp,k))
          enddo
        elseif(ic .le. ncmax+vac_ele_size+pri_ele_size) then ! private
          ic_tmp = ic-(ncmax+vac_ele_size)
          kmax = mseg_pri(ic_tmp)
          do k = 1,kmax
            x1_list(k) = pri_grid_x(pri_element(ic_tmp,k))
            y1_list(k) = pri_grid_y(pri_element(ic_tmp,k))
          enddo
        else
          ic_tmp = ic-(ncmax+vac_ele_size+pri_ele_size)
          kmax = mseg_subdiv(ic_tmp)
          do k = 1,kmax
            x1_list(k) = vgx_EX(subdiv_cell(ic_tmp,k))
            y1_list(k) = vgy_EX(subdiv_cell(ic_tmp,k))
          enddo
        endif
!
        do k = 1, kmax
          write(n6,'(2x,"grid ",2x,2f15.8)')
     >     x1_list(k), y1_list(k)
        enddo
      endif
!
      return
      end
