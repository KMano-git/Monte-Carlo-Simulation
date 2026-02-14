!**********************************************************************
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   subroutine atm_pset
      subroutine atm_pset( kis )
!**********************************************************************
!
!        set varaiables for intepolation
!
!----------------------------------------------------------------------
      use catcom, only : ndxkn, ndxrw, nxs, pin, plg, prk, pxa, pxd, pxi
     >    , pxs, xdnam, xdrnk, xjsta, xwmax, xwmin, xwnum, xwspc, xwvar
      use cunit,  only : n6
      implicit none
! added 3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
! arguments
      integer, intent(in) :: kis
! kis : impurity number

!::local varaiables
! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   integer i, j, n, ner, k
!ik   integer jst, nw, iw, lenx, ixs, rnk
      integer  i, j, n, ner
      integer  nw, iw, rnk
      character  zspc*12
! added 2 lines organize local variables and include files by kamata 2021/05/31
! function
! deleted 1 line treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   integer  lenx
!
!-----------------------------------------------------------------------
!
!::clear
! modified 9/7 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   do j = 1, ndxkn
!ik     prk(j) = -1
!ik   do i = 1, 3
!ik     pin(i,j) = 0; plg(i,j) = 0
!ik   enddo
!ik   do i = 1, 2
!ik     pxi(i,j)=0.0d0; pxa(i,j)=0.0d0; pxs(i,j)=0.0d0; pxd(i,j)=0.0d0
!ik   enddo
!ik   enddo
      pin(1:ndxrw,1:ndxkn,kis) = 0 
      plg(1:ndxrw,1:ndxkn,kis) = 0 
      prk(1:ndxkn,kis) = -1 
      pxa(1:2,1:ndxkn,kis) = 0.0_8
      pxd(1:2,1:ndxkn,kis) = 0.0_8
      pxi(1:2,1:ndxkn,kis) = 0.0_8
      pxs(1:2,1:ndxkn,kis) = 0.0_8
!
! modified 4/1 lines organize local variables and include files by kamata 2021/05/31
!ik   do j = 1, nxs
!ik   ixs = j
!
!ik   n   = ixs
! modified 3/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   do n = 1, nxs
!ik   rnk = xdrnk(n) - 1        ! exclude Iz
!ik   nw  = xdrnk(n) + 1
      do n = 1, nxs(kis)
      rnk = xdrnk(n,kis) - 1        ! exclude Iz
      nw  = xdrnk(n,kis) + 1
!
!::rank, pointer of y-data, number of x-data
! modified 5/5 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   prk(n)   = rnk
!ik   pin(1,n) = xjsta(n) - 1   ! <sig*v>
!ik   pin(2,n) = xwnum(2,n)     ! Ne
!ik   pin(3,n) = xwnum(3,n)     ! Te
!ik   pin(4,n) = xwnum(4,n)     ! Iz
      prk(n,kis)   = rnk
      pin(1,n,kis) = xjsta(n,kis) - 1   ! <sig*v>
      pin(2,n,kis) = xwnum(2,n,kis)     ! Ne
      pin(3,n,kis) = xwnum(3,n,kis)     ! Te
      pin(4,n,kis) = xwnum(4,n,kis)     ! Iz
!
!::log flag   Note  plg(1,2,3) = y, x1, x2
      do i = 1, nw
! modified 4/4 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   zspc = xwspc(i,n)
!ik   plg(i,n) = 0
!ik   if( zspc(1:5).eq."log  " ) plg(i,n) = 1
!ik   if( zspc(1:5).eq."log10" ) plg(i,n) = 2
      zspc = xwspc(i,n,kis)
      plg(i,n,kis) = 0
      if( zspc(1:5).eq."log  " ) plg(i,n,kis) = 1
      if( zspc(1:5).eq."log10" ) plg(i,n,kis) = 2
      enddo
!
!::min,max,start,dlt
      do i = 1, rnk
      iw = i + 1
! modified 6/8 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   pxi(i,n) = xwmin(iw,n)
!ik   pxa(i,n) = xwmax(iw,n)
!ik   pxs(i,n) = xwmin(iw,n)
!ik   pxd(i,n) =(xwmax(iw,n)-xwmin(iw,n))/dfloat(xwnum(iw,n)-1)
!ik   if( plg(iw,n).eq.1 ) then
!ik     pxi(i,n) = exp(pxi(i,n)); pxa(i,n) = exp(pxa(i,n))
      pxi(i,n,kis) = xwmin(iw,n,kis)
      pxa(i,n,kis) = xwmax(iw,n,kis)
      pxs(i,n,kis) = xwmin(iw,n,kis)
      pxd(i,n,kis) =(xwmax(iw,n,kis)-xwmin(iw,n,kis))
     >              /dfloat(xwnum(iw,n,kis)-1)
      if( plg(iw,n,kis) == 1 ) then
        pxi(i,n,kis) = exp(pxi(i,n,kis))
        pxa(i,n,kis) = exp(pxa(i,n,kis))
      endif
! modified 2/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   if( plg(iw,n).eq.2 ) then
!ik     pxi(i,n) = 10.0d0**pxi(i,n); pxa(i,n) = 10.0d0**pxa(i,n)
      if( plg(iw,n,kis) == 2 ) then
        pxi(i,n,kis) = 10.0d0**pxi(i,n,kis)
        pxa(i,n,kis) = 10.0d0**pxa(i,n,kis)
      endif
      enddo
!
      enddo
!
!-----------------------------------------------------------------------
!::check data
!      rank = 2,  spacing   vx1 & vx2 & vy   log10,  is < 100
!-----------------------------------------------------------------------
      ner = 0
! modified 1/2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   write(n6,'(/2x,"*** atm_pset ***  nxs =",i3)') nxs
      write(n6,'(/2x,"*** atm_pset ***  nxs =",i3, " kis =",i3 )')
     >     nxs(kis), kis
! modified 2/1 lines organize local variables and include files by kamata 2021/05/31
!ik   do j = 1, nxs
!ik   n = j
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   do n = 1, nxs
      do n = 1, nxs(kis)
      write(n6,'(i3,2x,a,2x,a,2x,i2,i8,i3,2(i6,i3,1p4e11.3))')
! modified 3/4 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik  >  n, xdnam(n)(1:lenx(xdnam(n))),xwvar(1,n)(1:lenx(xwvar(1,n)))
!ik  > ,prk(n),pin(1,n),plg(1,n), (pin(i+1,n),plg(i+1,n)
!ik  > ,pxi(i,n),pxa(i,n),pxs(i,n),pxd(i,n),i=1,prk(n))
     >   n, trim( xdnam(n,kis) ), trim( xwvar(1,n,kis) )
     > , prk(n,kis), pin(1,n,kis), plg(1,n,kis)
     > , ( pin(i+1,n,kis), plg(i+1,n,kis), pxi(i,n,kis), pxa(i,n,kis)
     >   , pxs(i,n,kis), pxd(i,n,kis), i=1,prk(n,kis) )
!-----
! modified 5/5 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   if( prk(n).ne.2 ) ner = ner + 1
!ik   if( pin(2,n).le.0 .or. pin(3,n).le.0 )   ner = ner + 1
!ik   if( pin(4,n).le.0 .or. pin(4,n).gt.100 ) ner = ner + 1
!ik   if( plg(1,n).ne.2 .or. plg(2,n).ne.2 .or.
!ik  >    plg(3,n).ne.2 ) ner = ner + 1
      if( prk(n,kis) /= 2 ) ner = ner + 1
      if( pin(2,n,kis) <= 0 .or. pin(3,n,kis) <= 0 )   ner = ner + 1
      if( pin(4,n,kis) <= 0 .or. pin(4,n,kis) >  100 ) ner = ner + 1
      if( plg(1,n,kis) /= 2 .or. plg(2,n,kis) /= 2 .or.
     >    plg(3,n,kis) /= 2 ) ner = ner + 1
!-----
      enddo
      if( ner.gt.0 ) then
        call wexit("xdt_pset","incorrect rank or spacing or data size")
      endif
!
      return
      end
