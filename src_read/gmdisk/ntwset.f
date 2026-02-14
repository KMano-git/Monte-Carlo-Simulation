!**********************************************************************
      subroutine ntwset
!**********************************************************************
!
!   set recycling and temperature at each wall segment
!   set exit number in vac-region  ihgw(ih)  = 0/ih < ndwh = 40
!       ihgw = (-2:W3g2, -1:W3g1, 0:Wwal, 1:g1, 2:g2, 3:a1-a2, 4:Vwal)
!             no cal. flux  when ihgw(ih) <=0  0:total flux into vac
!   set gate number igtw(iw)  ndwp = 600
!            = (-1:pump)/(0:wall)/(1,2:gate)/(3,4:shevron,partition)
!
!----------------------------------------------------------------------
      use cntcom, only : chwl, e0wl, ihgw, ipmp, irfw, nhwl, nogt, 
     >   nowl2, npew, npmp, npsw, rcwl, rcywl, temwl, tywl
      use csize,  only : ndwh, ndwp
      use cunit,  only : n6
      implicit none
!
!:: local variables
      integer  i, ih, nw, iw, ii, mj, n, iws, iwe
      integer  tnogt
      character ctyp*4
      real*8   zrcy, ztem, ze0
!
      write(n6,'(/2x,"*** ntwset ***  ",a)')
     > "set recy and temp at each wall segment"
!
!::clear
      do i = 1, ndwp
        rcwl(i) = 0.0d0
        e0wl(i) = 0.0d0
      enddo
!
!::recycling, wall temperature
      write(n6,'(2x,1x,"ih",2x,"chwl",5x,"rcy",5x,"twl",5x,"irfw")')
      do ih = 1, nhwl
        ctyp = chwl(ih)
        zrcy = rcywl(ih)
        ztem = temwl(ih)
        ze0  = 1.5d0*(ztem+273.0d0)/11600.0d0
        write(n6,'(2x,i3,2x,a,2f10.3,i4)') 
     >      ih, ctyp, zrcy, ztem, irfw(ih)
        do nw = 1, nowl2
          iws = npsw(nw)
          iwe = npew(nw)
          do iw = iws, iwe
            if( tywl(iw).eq.ctyp ) then
              rcwl(iw) = zrcy
              e0wl(iw) = ze0
            endif
          enddo
        enddo !nw
      enddo ! ih
!
!::pump region
      ii = 0
      do nw = 1, nowl2
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe
          if( index("ac",tywl(iw)(1:1)).gt.0 ) then
            ii = ii + 1
            ipmp(ii) = iw
          endif
        enddo
      enddo
      npmp = ii
      write(n6,'(2x,"npmp =",i3)') npmp
!
!::ihgw (new index for conductance)
      ihgw(1:ndwh) = 0
!
      tnogt = 0
      do ih = 1, nhwl
        ctyp = chwl(ih)
        if( ctyp(1:1).ne."W" ) cycle
        mj = index(ctyp,"g")
        if( mj.gt.0 ) then
          read(ctyp(mj+1:mj+1),*) n
          tnogt = max0( tnogt, n )
          ihgw(ih) = -n
        endif
      enddo
!
      do ih = 1, nhwl
        ctyp = chwl(ih)
        if( ctyp(1:1).eq."W" ) cycle
        mj = index(ctyp,"g")
        if( mj.gt.0 .and. ctyp(1:1).ne."V" ) then ! No use V1g1
          read(ctyp(mj+1:mj+1),*) n
          ihgw(ih) = n
        endif
!
        mj = index(ctyp,"a")
        if( mj.gt.0 .and. ctyp(1:1).ne."V" ) then
          ihgw(ih) = nogt + 1
        endif
!
        if( ctyp.eq."V1  " ) then
          ihgw(ih) = nogt + 2
        endif
      enddo
!
      write(n6,'(/2x,"nhwl =",i3,"  nogt =",2i3)') nhwl, tnogt, nogt
      write(n6,'(3x,"ih",2x,"chwl",4x,"rcywl",6x,"temwl",5x,
     >   "irfw",1x,"ihgw")')
      do ih = 1, nhwl
        write(n6,'(2x,i3,2x,a,2x,1p2e11.3,2i5)')
     >  ih, chwl(ih), rcywl(ih), temwl(ih), irfw(ih), ihgw(ih)
      enddo
!
      return
      end
