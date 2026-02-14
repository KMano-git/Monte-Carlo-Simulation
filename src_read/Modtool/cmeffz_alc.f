! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cmeffz
! nzsmx for SOLDOR
      subroutine cmeffz_alc( kk )
      use cmeffz, only : dzan, vdnz, vdnz0, vsez, vsez0, vsezg, vsezz
     >    , wai, wfab, wfeb, wfrab, wkab, wma, wza, xzflz, xzwtm, xzwtp
      use csize,  only : mdsp, ndmis, ndx, ndy, nzsmx
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
        if( .not. allocated( dzan ) ) then
          allocate( dzan(ndx,ndy), vdnz(ndx,ndy,0:ndmis,nzsmx)
     >      , vdnz0(ndx,ndy,0:ndmis,nzsmx), vsez(ndx,ndy)
     >      , vsez0(ndx,ndy), vsezg(ndx,ndy), vsezz(ndx,ndy), wai(mdsp)
     >      , wfab(mdsp,mdsp), wfeb(mdsp), wfrab(mdsp,mdsp)
     >      , wkab(mdsp,mdsp), wma(mdsp), wza(mdsp), xzflz(ndx,ndy)
     >      , xzwtm(ndx,ndy), xzwtp(ndx,ndy)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'dzan allocate error in cmeffz_alc, istat = ', istat
            call wexit( 'cmeffz', trim( cmsg ) )
          endif

! initial set
          dzan(1:ndx,1:ndy) = 0.0_8
          vdnz(1:ndx,1:ndy,0:ndmis,1:nzsmx) = 0.0_8
          vdnz0(1:ndx,1:ndy,0:ndmis,1:nzsmx) = 0.0_8
          vsez(1:ndx,1:ndy) = 0.0_8
          vsez0(1:ndx,1:ndy) = 0.0_8
          vsezg(1:ndx,1:ndy) = 0.0_8
          vsezz(1:ndx,1:ndy) = 0.0_8
          wai(1:mdsp) = 0.0_8
          wfab(1:mdsp,1:mdsp) = 0.0_8
          wfeb(1:mdsp) = 0.0_8
          wfrab(1:mdsp,1:mdsp) = 0.0_8
          wkab(1:mdsp,1:mdsp) = 0.0_8
          wma(1:mdsp) = 0.0_8
          wza(1:mdsp) = 0.0_8
          xzflz(1:ndx,1:ndy) = 0.0_8
          xzwtm(1:ndx,1:ndy) = 0.0_8
          xzwtp(1:ndx,1:ndy) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( dzan ) ) then
          deallocate( dzan, vdnz, vdnz0, vsez, vsez0, vsezg, vsezz, wai
     >      , wfab, wfeb, wfrab, wkab, wma, wza, xzflz, xzwtm, xzwtp
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'dzan deallocate error in cmeffz_alc, istat = ', istat
            call wexit( 'cmeffz', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
