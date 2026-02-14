!***********************************************************************
      subroutine edtgate
!***********************************************************************
!
!      kg = 1, 2  number of gate
!      ks = 1 (plasma)  2(vacuum region)
!
!         0   dsgt(2)  (3)        (4)        5 cm   distance
!         |.......|.....|.....x....|..........|
!         *-------*-----*----------*----------*
!
!         *---------------*-------------------*
!        5 cm  vacuum region                  0.0
!
!-----------------------------------------------------------------------
      use cntcom, only : dsgt, ipgt, iwgt, negt, nogt, nsgt, xpnt, ypnt
      use csize,  only : ndgtp
      use cunit,  only : n6
      use mod_externalgrid, only:use_exdata
      implicit none
!
!::local variables
      integer :: kg, ks, ig, ngs, nge, ips, ipe
      integer :: ip1, ip2, ii, i
      integer :: iw
      real(8) :: sx, xp, yp, sum
      real(8), dimension(200) :: xtmp, ytmp
      integer, dimension(200) :: itmp
!
      if(use_exdata) then
        call edtgate_ex ! fix me ! may not be needed
        return
      endif
!
      write(n6,'(/2x,"***  edtgate ***")')
!
      ii = 0
      do kg = 1, nogt
      do ks = 1, 2
        ngs = nsgt(kg,ks)
        nge = negt(kg,ks)
        ips = ipgt(ngs)
        ipe = ipgt(nge)
        do ig = ngs, nge
          ip1 = ipgt(ig)
          iw  = iwgt(ig)
          if(ip1==ips) then
            xp  =  xpnt(ips)
            yp  =  ypnt(ips)
          else
            sx = (xpnt(ip1)-xpnt(ips))/(xpnt(ipe)-xpnt(ips))
            xp  =  xpnt(ip1)
            yp  =  ypnt(ips) + sx*(ypnt(ipe)-ypnt(ips))
          endif
          ii = ii + 1
          if( ii.gt.200 ) goto 910
          itmp(ii) = ip1
          xtmp(ii) = xp
          ytmp(ii) = yp
          write(n6,'(1x,i3,i3,i5,i5,i7,2f14.9,"  ==>",2f14.9,1p2e12.3)')
     >      kg, ks, ig, iw, ip1, xpnt(ip1), ypnt(ip1), xp, yp,
     >      xpnt(ip1)-xp, ypnt(ip1)-yp
        enddo   ! loop(ig)  gate point
      enddo     ! loop(ks)  side (pla,vac)
      enddo     ! loop(kg)  gate
!
!::new data
      do i = 1, ii
        ip1 = itmp(i)
        xpnt(ip1) = xtmp(i)
        ypnt(ip1) = ytmp(i)
      enddo
!
!::define dsgt
      dsgt(1:ndgtp) = 0.0d0
      do kg = 1, nogt
      do ks = 1, 2
        ngs = nsgt(kg,ks)
        nge = negt(kg,ks)
        sum = 0.0d0
        do ig = ngs, nge-1
          dsgt(ig) = sum
          ip1 = ipgt(ig)
          ip2 = ipgt(ig+1)
          sum = sum + sqrt((xpnt(ip2)-xpnt(ip1))**2
     >         +(ypnt(ip2)-ypnt(ip1))**2)
        enddo
        ig = nge
        dsgt(ig) = sum
      enddo
      enddo
!
      return
!
!::error
 910  continue
      call wexit("edtgate","ii.gt.200")
      end

!***********************************************************************
      subroutine edtgate_ex
!***********************************************************************
      use cntcom, only : dsgt, ipgt, negt, nogt, nsgt
      use csize,  only : ndgtp
      use cunit,  only : n6
      use mod_externalgrid
      implicit none
!
!::local variables
      integer kg, ks, ii, ip, ip2, ngs, nge
      real(8) :: sum, x1, x2, y1, y2
!::
      write(n6,'(/2x,"***  edtgate_ex ***")')
!
!::define dsgt
!::(for external grid, dsgt is not available because ipgt is not sort as spatially order)
      dsgt(1:ndgtp) = 0.0d0
      do kg = 1, nogt
        do ks = 1, 2
          ngs = nsgt(kg,ks)
          nge = negt(kg,ks)
          sum = 0.0d0
          do ii = ngs, nge-1
            dsgt(ii) = sum
            ip  = ipgt(ii)
            ip2 = ipgt2(ii)
            if(ks == 1) then ! plasma to vacume
              x1 = pri_grid_x(ip)
              y1 = pri_grid_y(ip)
              x2 = pri_grid_x(ip2)
              y2 = pri_grid_y(ip2)
            else ! vacume to plasma
              x1 = vgx_EX(ip)
              y1 = vgy_EX(ip)
              x2 = vgx_EX(ip2)
              y2 = vgy_EX(ip2)
            endif
            sum = sum + sqrt((x2-x1)**2+(y2-y1)**2)
          enddo ! ii
          dsgt(nge) = sum
        enddo
      enddo        
      end subroutine edtgate_ex
!
!***********************************************************************
      subroutine lstgate
!***********************************************************************
      use cntcom, only : dsgt, icgt, ipgt, iwgt, mgrd, mknd, mseg, negt
     >    , nogt, nsgt, xpnt, ypnt
      use cunit,  only : n6
      use mod_externalgrid, only:use_exdata
      implicit none
!
!::local variables
      integer :: kg, ks, ig, ngs, nge, ips, ipe
      integer :: ig1, ip1, ip2, ip3, ip4
      integer :: ic, iw, kln, k, ln
      real(8) :: sx, sy
! function
      integer    imox, imoy
!:
      ! ignore when exdata is used because this subroutine only write log
      if(use_exdata) return
!
      write(n6,'(/2x,"***  lstgate ***")')
      write(n6,'(2x,"kg",1x,"ks",3x,"ig",3x,"iw",5x,"ic",4x,"ix",3x,
     >   "iy",1x,"kln",3x,"ip1",3x,"ip2",3x,"ip3",3x,"ip4",3x,
     >   "xpnt",10x,"ypnt",10x,"dsgt",10x,"sx",12x,"sy",12x,"sx-sy")')
!
      do kg = 1, nogt
      if( kg.ne.1 ) write(n6,'(2x)')
      do ks = 1, 2
        ngs = nsgt(kg,ks)
        nge = negt(kg,ks)
        ips = ipgt(ngs)
        ipe = ipgt(nge)
        do ig = ngs, nge
          ig1 = ig + 1
          ig1 = min0( ig1, nge )
          ip1 = ipgt(ig)
          ip2 = ipgt(ig1)
          if(xpnt(ipe) .eq. xpnt(ips)) then
            if(xpnt(ip1) .eq. xpnt(ips)) then
              sx = 1.0d0
            else
              call wexit("lstgate","xpnt(ip1)-xpnt(ips)=0")
            endif
          else
            sx  = (xpnt(ip1)-xpnt(ips))/(xpnt(ipe)-xpnt(ips))
          endif
          if(ypnt(ipe) .eq. ypnt(ips)) then
            if(ypnt(ip1) .eq. ypnt(ips)) then
              sy = 1.0d0
            else
              call wexit("lstgate","ypnt(ipe)-ypnt(ips)=0")
            endif
          else
            sy  = (ypnt(ip1)-ypnt(ips))/(ypnt(ipe)-ypnt(ips))
          endif
          ic = icgt(ig)
!---------
          if( ic.eq.0 ) then
            iw = iwgt(ig)
            kln = 0
            ip3 = 0
            ip4 = 0
          else
            iw = iwgt(ig)
            kln = 0
            do k = 1, mseg(ic)
              ln  = mknd(ic,k)
              if( ln.eq.(-iw) ) then
                kln = k
                exit
              endif
            enddo
            ip3 = mgrd(ic,kln)
            ip4 = mgrd(ic,kln+1)
          endif
!
          write(n6,'(1x,i3,i3,i5,i5,i7,i6,i5,i4,4i6,
     >      2f14.9,f14.9,2f14.9,1pe12.3)')
     >      kg, ks, ig, iw, ic, imox(ic), imoy(ic),
     >      kln, ip1, ip2, ip3, ip4,
     >      xpnt(ip1), ypnt(ip1), dsgt(ig), sx, sy, sx-sy
!
      enddo   ! loop(ig)  gate point
      enddo     ! loop(ks)  side (pla,vac)
      enddo     ! loop(kg)  gate
!
      return
      end
