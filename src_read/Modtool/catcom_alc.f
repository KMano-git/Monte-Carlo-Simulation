! added treat 4 or more impurities with IMPMC by kamata 2022/04/15
! allocate variables in catcom
      subroutine catcom_alc( kk )
      use csize,  only : ndmis, nzsmx
      use catcom, only : aionz, catmz, dratm, ip_cxr, ip_ion, ip_plt
     >    , ip_prb, ip_prc, ip_rad, ip_rec, natmz, ndxdt, ndxkn, ndxrw
     >    , ndxx1, ndxx2, nxs, nyp, pin, plg, prk, pxa, pxd, pxi, pxs
     >    , xdaty, xdatz, xdevl, xdnam, xdrnk, xdsn, xdtx1, xdtx2, xjend
     >    , xjnum, xjsta, xwmax, xwmin, xwmlt, xwnum, xwspc, xwunt
     >    , xwvar
      implicit none
! arguments
      integer, intent(in) :: kk
! kk : flag = 1: allocate,    = 2: deallocate

! local variables
      integer    istat
      character  cmsg*80

      select case( kk )
      case( 1 )
! allocate
        if( .not. allocated( aionz ) ) then
          allocate( aionz(nzsmx), catmz(nzsmx), dratm(nzsmx)
     >        , ip_cxr(nzsmx), ip_ion(nzsmx), ip_plt(nzsmx)
     >        , ip_prb(nzsmx), ip_prc(nzsmx), ip_rad(nzsmx)
     >        , ip_rec(nzsmx), natmz(nzsmx), nxs(nzsmx), nyp(nzsmx)
     >        , pin(ndxrw,ndxkn,nzsmx), plg(ndxrw,ndxkn,nzsmx)
     >        , prk(ndxkn,nzsmx), pxa(2,ndxkn,nzsmx)
     >        , pxd(2,ndxkn,nzsmx), pxi(2,ndxkn,nzsmx)
     >        , pxs(2,ndxkn,nzsmx)
     >        , xdaty(ndxdt,nzsmx), xdatz(ndmis,ndxkn,nzsmx)
     >        , xdevl(ndxkn,nzsmx), xdnam(ndxkn,nzsmx)
     >        , xdrnk(ndxkn,nzsmx), xdsn(ndxkn,nzsmx)
     >        , xdtx1(ndxx1,ndxkn,nzsmx), xdtx2(ndxx2,ndxkn,nzsmx)
     >        , xjend(ndxkn,nzsmx), xjnum(ndxkn,nzsmx)
     >        , xjsta(ndxkn,nzsmx), xwmax(ndxrw,ndxkn,nzsmx)
     >        , xwmin(ndxrw,ndxkn,nzsmx), xwmlt(ndxrw,ndxkn,nzsmx)
     >        , xwnum(ndxrw,ndxkn,nzsmx), xwspc(ndxrw,ndxkn,nzsmx)
     >        , xwunt(ndxrw,ndxkn,nzsmx), xwvar(ndxrw,ndxkn,nzsmx)
     >        , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'aionz allocate error in catcom_alc, istat = ', istat
            call wexit( 'catcom', trim( cmsg ) )
          endif

! initial set
          aionz(1:nzsmx) = 0.0_8
          catmz(1:nzsmx) = ' '
          dratm(1:nzsmx) = ' '
          ip_cxr(1:nzsmx) = 0
          ip_ion(1:nzsmx) = 0
          ip_plt(1:nzsmx) = 0
          ip_prb(1:nzsmx) = 0
          ip_prc(1:nzsmx) = 0
          ip_rad(1:nzsmx) = 0
          ip_rec(1:nzsmx) = 0
          natmz(1:nzsmx) = 0
          nxs(1:nzsmx) = 0
          nyp(1:nzsmx) = 0
          pin(1:ndxrw,1:ndxkn,1:nzsmx) = 0
          plg(1:ndxrw,1:ndxkn,1:nzsmx) = 0
          prk(1:ndxkn,1:nzsmx) = -1
          pxa(1:2,1:ndxkn,1:nzsmx) = 0.0_8
          pxd(1:2,1:ndxkn,1:nzsmx) = 0.0_8
          pxi(1:2,1:ndxkn,1:nzsmx) = 0.0_8
          pxs(1:2,1:ndxkn,1:nzsmx) = 0.0_8
          xdaty(1:ndxdt,1:nzsmx) = 0.0_8
          xdatz(1:ndmis,1:ndxkn,1:nzsmx) = 0.0_8
          xdevl(1:ndxkn,1:nzsmx) = ' '
          xdnam(1:ndxkn,1:nzsmx) = ' '
          xdrnk(1:ndxkn,1:nzsmx) = 0
          xdsn(1:ndxkn,1:nzsmx) = ' '
          xdtx1(1:ndxx1,1:ndxkn,1:nzsmx) = 0.0_8
          xdtx2(1:ndxx2,1:ndxkn,1:nzsmx) = 0.0_8
          xjend(1:ndxkn,1:nzsmx) = 0
          xjnum(1:ndxkn,1:nzsmx) = 0
          xjsta(1:ndxkn,1:nzsmx) = 0
          xwmax(1:ndxrw,1:ndxkn,1:nzsmx) = 0.0_8
          xwmin(1:ndxrw,1:ndxkn,1:nzsmx) = 0.0_8
          xwmlt(1:ndxrw,1:ndxkn,1:nzsmx) = 0.0_8
          xwnum(1:ndxrw,1:ndxkn,1:nzsmx) = 0
          xwspc(1:ndxrw,1:ndxkn,1:nzsmx) = ' '
          xwunt(1:ndxrw,1:ndxkn,1:nzsmx) = ' '
          xwvar(1:ndxrw,1:ndxkn,1:nzsmx) = ' '
        endif
      case( 2 )
! deallocate
        if( allocated( aionz ) ) then
          deallocate( aionz, catmz, dratm, ip_cxr, ip_ion, ip_plt
     >        , ip_prb, ip_prc, ip_rad, ip_rec, natmz, nxs, nyp, pin
     >        , plg, prk, pxa, pxd, pxi, pxs, xdaty, xdatz, xdevl, xdnam
     >        , xdrnk, xdsn, xdtx1, xdtx2, xjend, xjnum, xjsta, xwmax
     >        , xwmin, xwmlt, xwnum, xwspc, xwunt, xwvar
     >        , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'aionz deallocate error in catcom_alc, istat = ', istat
            call wexit( 'catcom', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
