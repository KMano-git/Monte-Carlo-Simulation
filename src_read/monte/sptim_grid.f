!**********************************************************************
      subroutine sptim_ex(xp,yp,vx,vy,ic,ln,wt,wx,wy,ker,kmax)
!**********************************************************************
      use cntcom, only:mgrd,xpnt,ypnt
      implicit none
!
!::argument
      real*8, intent(in) :: xp,yp
      real*8, intent(inout) :: wt, wx, wy,vx,vy
      integer, intent(inout) :: ic,ln
      integer, intent(in) :: kmax
      integer, intent(inout) :: ker
!
!::local variables
      integer  k, ip1, ip2
      real*8 x1_list(kmax),y1_list(kmax)
     > ,x2_list(kmax),y2_list(kmax)

      do k = 1,kmax
          ip1 = mgrd(ic,k)
          ip2 = mgrd(ic,k+1)
          x1_list(k) = xpnt(ip1)
          y1_list(k) = ypnt(ip1)
          x2_list(k) = xpnt(ip2)
          y2_list(k) = ypnt(ip2)
      enddo
      call sptim_main(xp,yp,vx,vy,ic,ln
     > ,wt,wx,wy,ker,kmax
     > ,x1_list,y1_list,x2_list,y2_list)
      end subroutine sptim_ex

!**********************************************************************
      subroutine sptim_grid(xp,yp,vx,vy,ic,ln
     > ,wt,wx,wy,ker
     > ,element,ele_size,kmax
     > ,grid_size,grid_x,grid_y,ic_in)
!**********************************************************************
      use csize, only:ndms
      implicit none
!
!::argument
      real*8, intent(in) :: xp, yp
      real*8, intent(inout) :: wt, wx, wy,vx, vy
      integer, intent(in) :: ele_size,grid_size,kmax,ic_in
      integer, intent(inout) :: ker
      integer, intent(inout) :: ln,ic
      integer, intent(in) :: element(ele_size,ndms)
      real*8, intent(in) :: grid_x(grid_size),grid_y(grid_size)
!
!::local variables
      integer  k, k1
      real*8 x1_list(kmax),y1_list(kmax)
     >   ,x2_list(kmax),y2_list(kmax)

      do k = 1,kmax
          k1 = k+1
          if( k.eq.kmax ) k1 = 1
          x1_list(k) = grid_x(element(ic_in,k))
          y1_list(k) = grid_y(element(ic_in,k))
          x2_list(k) = grid_x(element(ic_in,k1))
          y2_list(k) = grid_y(element(ic_in,k1))
      enddo
      call sptim_main(xp,yp,vx,vy,ic,ln
     > ,wt,wx,wy,ker,kmax
     > ,x1_list,y1_list,x2_list,y2_list)
      end subroutine sptim_grid

!**********************************************************************
      subroutine sptim_main(xp,yp,vx,vy,ic,ln
     > ,wt,wx,wy,ker,kmax
     > ,x1_list,y1_list,x2_list,y2_list)
!**********************************************************************
      use cntcom, only: mknd, next, weit, iptl
      use cunit,  only: mype, n6
      use mod_externalgrid ! for debug
      implicit none
!
!::argument
      real*8, intent(in) :: xp, yp  ! verocity vector and its start point
      real*8, intent(out) :: wt
      real*8, intent(inout) :: vx,vy,wx,wy ! wx,wy:intersection point
      integer, intent(in) :: kmax
      integer, intent(out) :: ker
      integer, intent(inout) :: ln,ic
      real*8, intent(in) :: x1_list(kmax),y1_list(kmax)
     > ,x2_list(kmax),y2_list(kmax)
!
!::local variables
      integer  k, ip1, ip2, kn, k1, ic_new
      real*8   zdx, zdy, zpx, zpy, zaa, zss, ztt
     > ,x1,y1,x2,y2,x_sect,y_sect,distance_max
      real*8 distances(kmax),wt_save(kmax)
     > ,wx_save(kmax),wy_save(kmax)

      ! initialization
      distances(1:kmax) = -10.0 ! Distance from the starting point() to the intersection
      wt_save(1:kmax) = 0.0_8
      wx_save(1:kmax) = 0.0_8
      wy_save(1:kmax) = 0.0_8

      ! determine where to move
      do k = 1, kmax
        x1 = x1_list(k)
        y1 = y1_list(k)
        x2 = x2_list(k)
        y2 = y2_list(k)

        zdx = x2 - x1
        zdy = y2 - y1
        zpx = xp - x1
        zpy = yp - y1
        zaa = zdx*vy - zdy*vx
        if( zaa.eq.0.0d0 ) goto 10
        zss = (zpx*vy-zpy*vx)/zaa
        ! Velocity vector intersects line segment (when 0<zss<=1) 
        if(zss.gt.0.0d0 .and. zss.le.1.0d0) then
          ztt = (zpx*zdy-zpy*zdx)/zaa
          if(ztt.ge.0.0d0) then
            ! calculate the distance between intersection and start point
            x_sect = -zpx+zss*zdx
            y_sect = -zpy+zss*zdy
            distances(k) = sqrt(x_sect*x_sect+y_sect*y_sect)
            ! save outputs
            wt_save(k) = ztt
            wx_save(k) = x1 + zdx*zss
            wy_save(k) = y1 + zdy*zss
          endif !ztt
        endif !zss
      enddo ! k
!::
      distance_max = -1.0
      do k = 1, kmax
        if(distance_max < distances(k)) then
          distance_max = distances(k)
          kn = k
        endif
      enddo
!::
      if(distance_max > 0)then
      !::next cell
        wt = wt_save(kn)
        wx = wx_save(kn)
        wy = wy_save(kn)
        ln = mknd(ic,kn)
        ic_new = next(ic,kn)
        ker = 0
        if(ic_new.eq.0 .and. ln.ge.0) then
          call crossBorder(wx,wy,ic,ic_new,ker)
        endif
        ic = ic_new
        return
      endif
      ker = 1
      return
!::
 10   continue
      write(n6,'(2x,"*** sptim  zaa.eq.0.0    ",1pe12.3)') zaa
      call wexit("sptim_grid","zaa.eq.0.0")
      end subroutine sptim_main

!**********************************************************************
      subroutine crossBorder(wx,wy,ic_from,ic_to,ker)
!**********************************************************************
      use mod_externalgrid
      use cntcom, only:ncmax,migx,migy,mgrd
      use cgdcom, only:mpw2,mqs1,mqs2,mqd2
      implicit none
!
!::argument
      real*8, intent(in) :: wx, wy
      integer, intent(in) :: ic_from
      integer, intent(inout) :: ic_to,ker
!
!::local variables
      integer k,jset,i,ic_offset
      logical ispri

      !***************
      ! Vacume to SOL
      !***************
      if(ic_from.gt.ncmax .and. ic_from.le.ncmax+vac_ele_size) then
          jset = 2
          ic_to = -1
          ispri = .false.
          call whereToGoToSOL(wx,wy,ic_to,jset,ispri,ic_from)
          if(ic_to.gt.0) then
              ker = 0
              return
          endif
      !***************
      ! private to SOL
      !***************
      elseif(ic_from.gt.ncmax+vac_ele_size) then
          jset = mpw2-1
          ic_to = -1
          ispri = .true.
          call whereToGoToSOL(wx,wy,ic_to,jset,ispri,ic_from)
          if(ic_to.gt.0) then
              ker = 0
              return
          endif
      !***************
      ! SOL to Vacume
      !***************
      elseif(ic_from.le.ncmax .and. 
     >       i_save(mgrd(ic_from,1)).eq.2) then
          ic_to = -1
          ic_offset = ncmax
          call whereToGoFromSOL(wx,wy
     >        ,vac_ele_size,vac_grid_size
     >        ,vac_grid_x,vac_grid_y
     >        ,mseg_vacume,vac_element
     >        ,boundary_vac,ic_offset,ic_from,ic_to)
          if(ic_to.gt.0) then
              ker = 0
              return
          endif
      !***************
      ! SOL to private
      !***************
      elseif(ic_from.le.ncmax 
     >    .and. i_save(mgrd(ic_from,3)).eq.mpw2-1
     >    .and. (j_save(mgrd(ic_from,3)).le.mqs1
     >           .or. j_save(mgrd(ic_from,3)).ge.mqs2)) then
          ic_to = -1
          ic_offset = ncmax+vac_ele_size
          call whereToGoFromSOL(wx,wy
     >        ,pri_ele_size,pri_grid_size
     >        ,pri_grid_x,pri_grid_y
     >        ,mseg_pri,pri_element
     >        ,boundary_pri,ic_offset,ic_from,ic_to)
          if(ic_to.gt.0) then
              ker = 0
              return
          endif
      endif

      ! not found where to go
      ker = 1
      call wexit("sptim_grid","not found where to go")
      return
      end subroutine crossBorder
!
!**********************************************************************
      subroutine whereToGoToSOL(wx,wy
     > ,ic_out,jset,ispri
     > ,ic_from)
!**********************************************************************
      use cgdcom, only:nqmx,grdx,grdy,mqs1,mqs2
      use cntcom, only:mcel,mclx,mcly,mcly2,nocl
      use cunit, only:n6, mype ! for debug
      use mod_externalgrid
      implicit none
!
!::argument
      real*8,  intent(in) :: wx,wy
      integer, intent(in) :: jset
      logical, intent(in) :: ispri
      integer, intent(out) :: ic_out
      integer, intent(in) :: ic_from
!
!::local variables
      integer i
      real*8 x1,x2,y1,y2, distance, distance_min
      real*8, parameter :: zero = 1e-10
!debug
      real*8, allocatable :: distance_allay(:)
!debug
      allocate(distance_allay(nqmx))
      distance_allay = 0.0d0

      distance_min = 1.0d10
      ic_out = -1
      do i = 1, nqmx-1
          if(ispri .and. (i+1>mqs1 .and. i<mqs2))cycle
          x1 = grdx(i,jset)
          y1 = grdy(i,jset)
          x2 = grdx(i+1,jset)
          y2 = grdy(i+1,jset)
          if(abs(x1-x2)<zero .and. abs(y1-y2)<zero) cycle
          call calc_dist(x1,y1,x2,y2,wx,wy,distance)
          ! debug
          distance_allay(i) = distance
          !
          if(distance.ge.0.0d0.and.distance<distance_min)then
            distance_min = distance
            if(ispri) then
              ic_out = nocl(i,jset-1)
            else
              ic_out = mcel(mclx(i),mcly(jset))
            endif
          endif
      enddo
!::
      if(ic_out.gt.0) then
        return
      endif

!:: not find where to go
      write(n6,'(2x,"whereToGoToSOL, not found ",
     >         " ispri =",l3)') ispri
      write(20000+mype,'("whereToGoToSOL, not found,is_pri"
     >      ,l3)') ispri
      write(20000+mype,'(2x,"wx:",1pe20.10,
     >     " wy:",1pe20.10," ic_from",i6)') wx,wy,ic_from
      do i = 1, nqmx
        write(20000+mype,'(2x,1pe10.2)') distance_allay(i)
      enddo
      call wexit("whereToGoToSOL","not found")
      end subroutine whereToGoToSOL

!**********************************************************************
      subroutine whereToGoFromSOL(wx,wy,ele_size,mygrid_size
     >    ,grid_x,grid_y,my_mseg,elements
     >    ,boundary,ic_offset,ic_from,ic_out)
!**********************************************************************
      use csize, only:ndms
      use mod_externalgrid!, only:label_SOLbound
      use cunit, only:n6,mype ! for debug
      implicit none
!
!::argument
      real*8, intent(in) :: wx,wy
      integer, intent(in) :: ele_size,mygrid_size,ic_offset,ic_from
      integer, intent(in) :: my_mseg(ele_size),elements(ele_size,ndms)
      real*8, intent(in) :: grid_x(mygrid_size),grid_y(mygrid_size)
      character(len=2), intent(in) :: boundary(ele_size,ndms)
      integer, intent(out) :: ic_out
!
!::local variables
      integer ic,k,kmax,ip1,ip2,i
      real*8 x1,x2,y1,y2, distance, distance_min
!debug
      real*8, allocatable :: distance_allay(:)
!debug
      allocate(distance_allay(ele_size))
      distance_allay = 0.0d0
!::
      distance_min = 1.0d10
      ic_out = -1
      do ic = 1,ele_size
          kmax = my_mseg(ic)
          do k = 1, kmax
              if(boundary(ic,k).eq.label_SOLbound) then
                  ip1 = elements(ic,k)
                  if(k==kmax) then
                      ip2 = elements(ic,1)
                  else
                      ip2 = elements(ic,k+1)
                  endif
                  x1 = grid_x(ip1)
                  y1 = grid_y(ip1)
                  x2 = grid_x(ip2)
                  y2 = grid_y(ip2)
                  call calc_dist(x1,y1,x2,y2,wx,wy,distance)
                  ! debug
                  distance_allay(ic) = distance
                  !
                  if(distance.ge.0.0d0.and.distance<distance_min)then
                     distance_min = distance
                     ic_out = ic+ic_offset
                  endif
              endif
          enddo
      enddo
!::
      if(ic_out.gt.0) then
        return
      endif
!:: not find where to go
      write(n6,'(2x,"whereToGoFromSOL, not found ")') 
      write(20000+mype,'("whereToGoFromSOL, not found")')
      write(20000+mype,'(2x,"wx:",1pe20.10,
     >     " wy:",1pe20.10)') wx,wy
      do i = 1, ele_size
        write(20000+mype,'(2x,1pe10.2)') distance_allay(i)
      enddo
      call wexit("whereToGoFromSOL","not found ic")
      end subroutine whereToGoFromSOL

!**********************************************************************
      subroutine calc_dist(x1,y1,x2,y2,wx,wy,distance)
!**********************************************************************
      implicit none
!
!::argument
      real*8, intent(in) :: x1,y1,x2,y2,wx,wy
      real*8, intent(out) :: distance
!
!::local variables
      real*8 a,b,c
      real*8, parameter :: zero = 1.0d-10
!
      distance = -1.0d0
!::   out of range
      if(min(x1,x2)-zero>wx .or. max(x1,x2)+zero<wx) return
      if(min(y1,y2)-zero>wy .or. max(y1,y2)+zero<wy) return
!
!::   ax+by+c=0:line equation
      a = y2-y1
      b = -x2+x1
      c = -x1*y2+x2*y1
!
!::   difference between point(wx,wy) to line
      distance = abs(a*wx+b*wy+c)/sqrt(a**2+b**2)
!
      end subroutine calc_dist