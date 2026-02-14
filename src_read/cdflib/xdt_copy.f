!**********************************************************************
      subroutine xdt_copy(cxs)
!**********************************************************************
      use cdfcom, only : ndat_max
      use cxdcom, only : ndxdt, ndxkn, ndxrw, nxs, nyp, xdaty, xdevl
     >    , xdnam, xdrnk, xjend, xjnum, xjsta, xwmax, xwmin, xwmlt
     >    , xwnum, xwspc, xwunt, xwvar
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   character cxs*(*)
      character, intent(in) :: cxs*(*)

      character cmsg*80
      integer  nymx, n, i, j, ndrw, ndty, nw
! modified 1/3 lines organize local variables and include files by kamata 2021/06/16
!ik   integer  ynum, ynum2, lenx
      integer  ynum, ynum2
! function
      integer    lenx
!
!-----------------------------------------------------------------------
!
!::pointer
      nxs = nxs + 1
      nyp = nyp + 1
      nymx= nyp + ndat_max
!
!::dimension error
      if( nxs.gt.ndxkn ) then
      write(cmsg,'("too many variables   nxs > ndxkn ",2i5)') nxs,ndxkn
      call wexit("xdt_copy",cmsg)
      endif
      if( nymx.gt.ndxdt ) then
      write(6,'(/2x,"Warning ! nyp+ndat_max > ndxdt ",2i8)') nymx,ndxdt
      endif
!
      n = nxs
      j = nyp
      ndrw = ndxrw
      ndty = ndxdt - nyp
!
      call cdf_dat(cxs,xdnam(n),xdevl(n),xdrnk(n),xdaty(j)
     >       ,xwvar(1,n),xwnum(1,n),xwmin(1,n),xwmax(1,n)
     >       ,xwmlt(1,n),xwunt(1,n),xwspc(1,n),ndrw,ndty)
!
!::pointer
      nw = xdrnk(n) + 1
      ynum  = xwnum(1,n)
      ynum2 = 1
      do i = 2, nw
      ynum2 = ynum2*xwnum(i,n)
      enddo
      xjsta(n) = j
      xjend(n) = j + ynum - 1
      xjnum(n) = ynum
!
!::pointer
      nyp = xjend(n)
!
!::debug write
      write(6,'(/2x,"*** xdt_copy ***  [",a,"]  nxs =",i3,"  cvar =",a
     > ,2x,"devl =",a12,"  drnk =",i2)')  xdnam(n)(1:lenx(xdnam(n)))
     > ,nxs,cxs(1:lenx(cxs)),xdevl(n),xdrnk(n)
      do i = 1, nw
      write(6,'(4x,"wvar =",a18,"  wnum =",i6,"  wmin =",1pe10.3,
     >  "  wmax =",1pe10.3,"  wmlt =",1pe10.3,"  wunt =",a8,
     >  "  wspc =",a8)') xwvar(i,n),xwnum(i,n),xwmin(i,n),xwmax(i,n)
     > ,xwmlt(i,n),xwunt(i,n),xwspc(i,n)
      enddo
      if( ynum.gt.1 ) then
      write(6,'(4x,"ynum =",2i6,
     >     "  daty =",i6,1x,1pe10.3,"  daty =",i6,1x,1pe10.3)')
     >   ynum,ynum2,xjsta(n),xdaty(xjsta(n)),xjend(n),xdaty(xjend(n))
      else
      write(6,'(4x,"ynum =",2i6,
     >     "  daty =",i6,1x,1pe10.3)')
     >   ynum,ynum2,xjsta(n),xdaty(xjsta(n))
      endif
!
!::linear/log
      do i = 2, nw
      if( xwspc(i,n)(1:lenx(xwspc(i,n))).eq."log" ) then
      xwmin(i,n) = log(xwmin(i,n))
      xwmax(i,n) = log(xwmax(i,n))
      endif
      enddo
      if( xwspc(1,n)(1:lenx(xwspc(1,n))).eq."log" ) then
      do j = xjsta(n), xjend(n)
      xdaty(j) = log(xdaty(j))
      enddo
      endif
!
      return
      end
