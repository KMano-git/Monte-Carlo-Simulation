!***********************************************************************
      subroutine fndway(xst,yst,xen,yen,icst,icen,ndmv,nmv,icmv,dtmv)
!***********************************************************************
!
!         xst   (i)   x-coordinate of start point
!         yst   (i)   y-coordinate of start point
!         xen   (i/o) x-coordinate of end point
!         yen   (i/o) y-coordinate of end point
!         icst  (i)   cell number of start point
!         icen  (o)   cell number of end point
!         ndmv  (i)   dimension size of movement
!         nmv   (o)   number of movement
!         icmv  (o)   cell number
!         dtmv  (o)   passing time in the cell
!
!    Note.  icmv(1) = icst
!           icen = 0    out of system
!           icen = lon (when hit wall)  iwp = -lon
!
!           sum_[i=1,nmv]dtmv(i) = 1.0d0  (icen > 0)
!                                < 1.0d0  (icen =< 0)
!
!-----------------------------------------------------------------------
      use mod_externalgrid
      use cunit, only : mype
      implicit none
!
!::argument
      integer, intent(in)    :: ndmv
      real(8), intent(in)    :: xst, yst
      real(8), intent(inout) :: xen, yen
      real(8), intent(inout) :: dtmv(ndmv)
      integer, intent(in)    :: icst
      integer, intent(out)   :: icen, nmv
      integer, intent(inout) :: icmv(ndmv)
!
!::local variables
      integer  lon, ii, ken, kou
      real*8   zpx, zpy
!
!::same cell
      call mchkin(xen,yen,icst,kou)
      if( kou.eq.0 ) then
        nmv  = 1
        icen = icst
        icmv(1) = icen
        dtmv(1) = 1.0d0
        return
      endif
!
      if(.not. use_exdata) then
        call trac_way(xst,yst,xen,yen,icst
     >       ,ndmv,icmv,dtmv,ii,lon,ken,zpx,zpy)
      else
        call trac_way_ex(xst,yst,xen,yen,icst
     >       ,ndmv,icmv,dtmv,ii,lon,ken,zpx,zpy)
      endif
!
      if(ii.gt.0)then
        nmv  = ii
        icen = icmv(ii)
        if( ken.eq.2 ) then
          icen = lon   ! <== wall number
        endif
        xen = zpx
        yen = zpy
      else
        nmv  = 0
        icen = 0
        icmv = 0
        dtmv = 0.0d0
      endif
!
      end subroutine fndway

!***********************************************************************
      subroutine trac_way(xst,yst,xen,yen,icst
     >       ,ndmv,icmv,dtmv,ii,lon,ken,zpx,zpy)
!***********************************************************************
      use cimcom, only : ipcm, lpcm
      use cntcom, only : mgrd, mknd, mseg, next, xpnt, ypnt
      use cunit,  only : n6
      implicit none
!
!::argument
      integer, intent(in)  :: ndmv
      real(8), intent(in)  :: xst, yst
      real(8), intent(in)  :: xen, yen
      real(8), intent(out) :: dtmv(ndmv)
      integer, intent(in)  :: icst
      integer, intent(out) :: icmv(ndmv)
      integer, intent(out) :: ii, lon, ken
      real(8), intent(out) :: zpx, zpy
!
!::local variables
      integer  ic, k
      integer  ipmg1, ipmg2, ipmg
      real*8   vex, vey, zdx, zdy, zex, zey, zaa, zss, ztt, ztt0
      integer  nxon(4); data nxon/3,4,1,2/
      save     nxon
      real*8   eps, epsh, ztmin
      integer  jj, j, jno
      real*8   svtm(10), svss(10), svpx(10), svpy(10)
      integer  svkl(10)
! function
      integer  imox, imoy
!
!::start
      vex = xen - xst
      vey = yen - yst
      ztt0= 0.0d0
      ic  = icst
      lon = 0
      ken = 0
      ii  = 0
      eps = 1.0d-10
      epsh = 0.5d0*eps
!
      do ! main loop
        ii = ii + 1
        if( ii.gt.ndmv ) call wexit("fndway","ii.gt.ndmv")
!
        ztmin = 1.0d20
        jno = 0
        jj = 0
!
        do k = 1, mseg(ic)
          ipmg1 = mgrd(ic,k)
          ipmg2 = mgrd(ic,k+1)
          zdx = xpnt(ipmg2)-xpnt(ipmg1)
          zdy = ypnt(ipmg2)-ypnt(ipmg1)
          zex = xst - xpnt(ipmg1)
          zey = yst - ypnt(ipmg1)
          zaa = zdx*vey - zdy*vex
          if( zaa.eq.0.0d0 ) cycle
          zss = (zex*vey-zey*vex)/zaa
          ztt = (zex*zdy-zey*zdx)/zaa
          jj = jj + 1
          svtm(jj) = ztt
          svss(jj) = zss
          svpx(jj) = xpnt(ipmg1) + zdx*zss
          svpy(jj) = ypnt(ipmg1) + zdy*zss
          svkl(jj) = k
          if( ztt.gt.ztt0+eps .and. ztt.lt.ztmin ) then
            ztmin = ztt
            jno = jj
          endif
        enddo

!::ztt = ztt0   ! %KS  2009/11/16
        if( jno.le.0 ) then
          do j = 1, jj
            if( (svss(j).ge.0.0d0 .and. svss(j).le.1.0d0) 
     >          .and. (dabs(svtm(j)-ztt0).le.5.0d-12 ) ) then
              jno = j
              ztt = svtm(jno)
              write(n6,'(2x,"ztt ~ ztt0  jno =",i3," ztt0,ztt =",
     >            1p2e16.8)') jno, ztt0, ztt
              ztt0 = dmin1( ztt0, ztt-epsh )
              ztmin = ztt
              exit
            endif
          enddo
!
!::error
          if( jno.eq.0 ) then  ! lon = 0 inside a cell.
            write(n6,'(2x,"*** fndway ***  no found col-point lpcm ="
     >       ,i10,"  ipcm =",i5)') lpcm, ipcm
            write(n6,'(2x,"ic =",i7,"  ix,iy =",2i5,"  ln0 =",i6)')
     >       ic, imox(ic), imoy(ic), lon
            write(n6,'(2x,"start",1p3e16.8)') xst, yst
            write(n6,'(2x,"end  ",1p3e16.8)') xen, yen
            do k = 1, mseg(ic)+1
              ipmg = mgrd(ic,k)
              write(n6,'(2x,"grid ",1p2e16.8)') xpnt(ipmg),ypnt(ipmg)
            enddo
            do j = 1, jj
              write(n6,'(2x,"collp",1p5e16.8,i4)')
     >        svpx(j),svpy(j),svtm(j),svss(j),svtm(j)-ztt0,svkl(j) ! KSFUJI
            enddo
            ii = 0
            return
          endif ! error
        endif ! not found
!
        ztt = svtm(jno)
        k   = svkl(jno)
        lon = mknd(ic,k)   ! <== Note
        if( ztt.lt.1.0d0 ) then
          if( lon.le.0 ) ken = 2
          zpx = svpx(jno)
          zpy = svpy(jno)
        else
          ken = 1
          zpx = xen   !  xst + vex*1.0d0
          zpy = yen   !  yst + vey*1.0d0
          ztt = 1.0d0
        endif
!
        icmv(ii) = ic
        dtmv(ii) = ztt-ztt0
        ztt0 = ztt
        if( ken.ne.0 ) exit
        ic  = next(ic,k)
        lon = nxon(lon)
      enddo ! main loop
      end subroutine trac_way

!***********************************************************************
      subroutine trac_way_ex(xst,yst,xen,yen,icst
     >       ,ndmv,icmv,dtmv,ii,lon,ken,zpx,zpy)
!***********************************************************************
      use cimcom, only : ipcm, lpcm
      use cntcom, only : mgrd, mknd, mseg, next, xpnt, ypnt, ncmax
      use cunit,  only : n6
      use mod_externalgrid
      implicit none
!
!::argument
      integer, intent(in)  :: ndmv
      real(8), intent(in)  :: xst, yst
      real(8), intent(in)  :: xen, yen
      real(8), intent(out) :: dtmv(ndmv)
      integer, intent(in)  :: icst
      integer, intent(out) :: icmv(ndmv)
      integer, intent(out) :: ii, lon, ken
      real(8), intent(out) :: zpx, zpy
!
!::local variables
      integer  ic, k, ipmg1, ipmg2, kmax, ic_tmp, k1, ic_new
      real*8   vex, vey, zdx, zdy, zex, zey, zaa, zss, ztt, ztt0
      real*8   eps, epsh, ztmin, x1, x2, y1, y2
      integer  jj, j, jno, ker
      real*8   svtm(10), svss(10), svpx(10), svpy(10)
      integer  svkl(10)
      integer, save ::  nxon(4); data nxon/3,4,1,2/
      integer, parameter :: grid_max = 4
      real*8  x1_list(grid_max),x2_list(grid_max)
     > ,y1_list(grid_max),y2_list(grid_max)
      real*8 distance,distance_max,x_sect,y_sect
!
!::start
      vex = xen - xst
      vey = yen - yst
      ztt0= 0.0d0
      ic  = icst
      lon = 0
      ken = 0
      ii  = 0
      eps = 1.0d-10
      epsh = 0.5d0*eps
!------------------------------------------------------------------
      do ! main loop
        ii = ii + 1
        if( ii.gt.ndmv ) call wexit("trac_way_ex","ii.gt.ndmv")
!
        x1_list = 0.0d0
        x2_list = 0.0d0
        y1_list = 0.0d0
        y2_list = 0.0d0
!
        ztmin = 1.0d20
        distance_max = -1.0d0
        jno = 0
        jj = 0
!------------------------------------------------------------------
        if(ic.le.ncmax)then
          kmax = mseg(ic)
          if(kmax > grid_max) call wexit("fndway","kmax>4")
          do k = 1, kmax
            ipmg1 = mgrd(ic,k)
            ipmg2 = mgrd(ic,k+1)
            x1_list(k) = xpnt(ipmg1)
            x2_list(k) = xpnt(ipmg2)
            y1_list(k) = ypnt(ipmg1)
            y2_list(k) = ypnt(ipmg2)
          enddo
        elseif(ic.le. ncmax+vac_ele_size)then
          ic_tmp = ic-ncmax
          kmax = mseg_vacume(ic_tmp)
          if(kmax > grid_max) call wexit("fndway","kmax>4")
          do k = 1, kmax
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
          do k = 1, kmax
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
          do k = 1, kmax
            k1 = k+1
            if( k.eq.kmax ) k1 = 1
            x1_list(k) = vgx_EX(subdiv_cell(ic_tmp,k))
            y1_list(k) = vgy_EX(subdiv_cell(ic_tmp,k))
            x2_list(k) = vgx_EX(subdiv_cell(ic_tmp,k1))
            y2_list(k) = vgy_EX(subdiv_cell(ic_tmp,k1))
          enddo
        endif
!------------------------------------------------------------------
!::trac point
        do k = 1, kmax
          x1 = x1_list(k)
          x2 = x2_list(k)
          y1 = y1_list(k)
          y2 = y2_list(k)
          zdx = x2 - x1
          zdy = y2 - y1
          zex = xst - x1
          zey = yst - y1
          zaa = zdx*vey - zdy*vex
          if( zaa.eq.0.0d0 ) cycle
          zss = (zex*vey-zey*vex)/zaa
          ztt = (zex*zdy-zey*zdx)/zaa
          jj = jj + 1
          svtm(jj) = ztt
          svss(jj) = zss
          svpx(jj) = x1 + zdx*zss
          svpy(jj) = y1 + zdy*zss
          svkl(jj) = k
          if(zss.gt.0.0d0 .and. zss.le.1.0d0) then
            if(ztt.gt.ztt0)then
              ! calculate the distance between intersection and start point
              x_sect = -zex+zss*zdx
              y_sect = -zey+zss*zdy
              distance = sqrt(x_sect*x_sect+y_sect*y_sect)
              if(distance_max < distance)then
                distance_max = distance
                ztmin = ztt
                jno = jj
              endif
            endif
          endif
        enddo !k
!------------------------------------------------------------------
!::ztt = ztt0   ! %KS  2009/11/16
        if( jno.le.0 ) then
          do j = 1, jj
            if( (svss(j).ge.0.0d0 .and. svss(j).le.1.0d0) 
     >          .and. (dabs(svtm(j)-ztt0).le.5.0d-12 ) ) then
              jno = j
              ztt = svtm(jno)
              write(n6,'(2x,"ztt ~ ztt0  jno =",i3," ztt0,ztt =",
     >            1p2e16.8)') jno, ztt0, ztt
              ztt0 = dmin1( ztt0, ztt-epsh )
              ztmin = ztt
              exit
            endif
          enddo
!
!::error
          if( jno.eq.0 ) then  ! lon = 0 inside a cell.
            write(n6,'(2x,"*** fndway ***  no found col-point lpcm ="
     >       ,i10,"  ipcm =",i5)') lpcm, ipcm
            write(n6,'(2x,"ic =",i7,"  ln0 =",i6)') ic, lon
            write(n6,'(2x,"start",1p3e16.8)') xst, yst
            write(n6,'(2x,"end  ",1p3e16.8)') xen, yen
            ii = 0
            return
          endif ! error
        endif ! not found
!
        ztt = svtm(jno)
        k   = svkl(jno)
        lon = mknd(ic,k)   ! <== Note
        if( ztt.lt.1.0d0 ) then
          if( lon.le.0 ) ken = 2
          zpx = svpx(jno)
          zpy = svpy(jno)
        else
          ken = 1
          zpx = xen   !  xst + vex*1.0d0
          zpy = yen   !  yst + vey*1.0d0
          ztt = 1.0d0
        endif
!
        icmv(ii) = ic
        dtmv(ii) = ztt-ztt0
        ztt0 = ztt
        if( ken.ne.0 ) return
        ic_new = next(ic,k)
        if(ic_new.eq.0 .and. lon.ge.0) then
          ker = 0
          call crossBorder(zpx,zpy,ic,ic_new,ker)
        endif
        ic = ic_new
        lon = nxon(lon)
      enddo ! main loop
      end subroutine trac_way_ex