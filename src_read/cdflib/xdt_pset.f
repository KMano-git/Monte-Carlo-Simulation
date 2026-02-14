!**********************************************************************
      subroutine xdt_pset
!**********************************************************************
!
!        set varaiables for intepolation
!
!----------------------------------------------------------------------
      use cxdcom, only : ndxkn, nxs, pin, plg, prk, pxa, pxd, pxi, pxs
     >    , xdnam, xdrnk, xjsta, xwmax, xwmin, xwnum, xwspc, xwvar
      use cunit,  only : n6
      implicit none

!::local varaiables
      integer i, j, n, ner
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer jst, nw, iw, lenx, ixs, rnk
      integer nw, iw, rnk
      character  zspc*40
! added 2 lines organize local variables and include files by kamata 2021/06/16
! function
      integer    lenx
!
!-----------------------------------------------------------------------
!
!::clear
      do j = 1, ndxkn
        prk(j) = -1
      do i = 1, 3
        pin(i,j) = 0; plg(i,j) = 0
      enddo
      do i = 1, 2
        pxi(i,j)=0.0d0; pxa(i,j)=0.0d0; pxs(i,j)=0.0d0; pxd(i,j)=0.0d0
      enddo
      enddo
!
! deleted 2 lines organize local variables and include files by kamata 2021/06/16
!ik   do j = 1, nxs
!ik   ixs = j
!
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!k    n   = ixs
      do n = 1, nxs
      rnk = xdrnk(n)
      nw  = xdrnk(n) + 1
!
!::rank, pointer of y-data, number of x-data
!
      prk(n)   = xdrnk(n)
      pin(1,n) = xjsta(n) - 1
      pin(2,n) = xwnum(2,n)
      pin(3,n) = xwnum(3,n)
!
!::log flag   Note  plg(1,2,3) = y, x1, x2
      do i = 1, nw
      zspc = xwspc(i,n)
      plg(i,n) = 0
      if( zspc(1:lenx(zspc)).eq."log" ) plg(i,n) = 1
      enddo
!
!::min,max,start,dlt
      do i = 1, rnk
      iw = i + 1
      pxi(i,n) = xwmin(iw,n)
      pxa(i,n) = xwmax(iw,n)
      pxs(i,n) = xwmin(iw,n)
      pxd(i,n) =(xwmax(iw,n)-xwmin(iw,n))/dfloat(xwnum(iw,n)-1)
      if( plg(iw,n).eq.1 ) then
        pxi(i,n) = exp(pxi(i,n)); pxa(i,n) = exp(pxa(i,n))
      endif
      enddo
!
      enddo
!
!::debug write
      ner = 0
      write(n6,'(/2x,"*** xdt_pset ***")')
! modified 2/1 lines organize local variables and include files by kamata 2021/06/16
!ik   do j = 1, nxs
!ik   n = j
      do n = 1, nxs
      write(n6,'(i3,2x,a16,2x,a20,2x,i2,i8,i3,2(i6,i3,1p4e11.3))')
     >  n, xdnam(n)(1:lenx(xdnam(n))),xwvar(1,n)(1:lenx(xwvar(1,n)))
     > ,prk(n),pin(1,n),plg(1,n), (pin(i+1,n),plg(i+1,n)
     > ,pxi(i,n),pxa(i,n),pxs(i,n),pxd(i,n),i=1,prk(n))
      if( prk(n).lt.0 .or. prk(n).gt.2 ) ner = ner + 1
      enddo
      if( ner.gt.0 ) then
        call wexit("xdt_pset","incorrect rank")
      endif
!
      return
      end
