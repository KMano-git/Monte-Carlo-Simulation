! added dynamic allocation of arrays by kamata 2022/05/29
! allocate variables in cplmet
      subroutine cplmet_alc( kk )
      use cplmet, only : gare, gdsv, gwtm, gwtp, hare, hdsp, hdsv, hdxm
     >    , hdxp, hgdx, hgdy, hpit, hvol, hvsb, hwtm, hwtp, icel, icmax
     >    , icmin, jcel, jnxm, jnxp, jtmax, jtmin, kgdx, kgdy, kreg
     >    , romn, vlmn
      use csize,  only : ndx, ndy
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
        if( .not. allocated( gare) ) then
          allocate( gare(ndx,ndy), gdsv(ndx,ndy,4), gwtm(ndx,ndy)
     >      , gwtp(ndx,ndy), hare(ndx,ndy), hdsp(ndx,ndy,4)
     >      , hdsv(ndx,ndy,4), hdxm(ndx,ndy), hdxp(ndx,ndy)
     >      , hgdx(ndx,ndy), hgdy(ndx,ndy), hpit(ndx,ndy), hvol(ndx,ndy)
     >      , hvsb(ndx,ndy), hwtm(ndx,ndy), hwtp(ndx,ndy), icel(ndx,ndy)
     >      , icmax(ndx), icmin(ndx), jcel(ndx,ndy), jnxm(ndx,ndy)
     >      , jnxp(ndx,ndy), jtmax(ndy), jtmin(ndy), kgdx(ndx,ndy,4)
     >      , kgdy(ndx,ndy,4), kreg(ndx,ndy), romn(ndy), vlmn(ndx)
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'gare allocate error in cplmet_alc, istat = ', istat
            call wexit( 'cplmet', trim( cmsg ) )
          endif

! initial set
          gare(1:ndx,1:ndy) = 0.0_8
          gdsv(1:ndx,1:ndy,1:4) = 0.0_8
          gwtm(1:ndx,1:ndy) = 0.0_8
          gwtp(1:ndx,1:ndy) = 0.0_8
          hare(1:ndx,1:ndy) = 0.0_8
          hdsp(1:ndx,1:ndy,1:4) = 0.0_8
          hdsv(1:ndx,1:ndy,1:4) = 0.0_8
          hdxm(1:ndx,1:ndy) = 0.0_8
          hdxp(1:ndx,1:ndy) = 0.0_8
          hgdx(1:ndx,1:ndy) = 0.0_8
          hgdy(1:ndx,1:ndy) = 0.0_8
          hpit(1:ndx,1:ndy) = 0.0_8
          hvol(1:ndx,1:ndy) = 0.0_8
          hvsb(1:ndx,1:ndy) = 0.0_8
          hwtm(1:ndx,1:ndy) = 0.0_8
          hwtp(1:ndx,1:ndy) = 0.0_8
          icel(1:ndx,1:ndy) = 0
          icmax(1:ndx) = 0
          icmin(1:ndx) = 0
          jcel(1:ndx,1:ndy) = 0
          jnxm(1:ndx,1:ndy) = 0
          jnxp(1:ndx,1:ndy) = 0
          jtmax(1:ndy) = 0
          jtmin(1:ndy) = 0
          kgdx(1:ndx,1:ndy,1:4) = 0
          kgdy(1:ndx,1:ndy,1:4) = 0
          kreg(1:ndx,1:ndy) = 0
          romn(1:ndy) = 0.0_8
          vlmn(1:ndx) = 0.0_8
        endif
      case( 2 )
! deallocate
        if( allocated( gare ) ) then
          deallocate( gare, gdsv, gwtm, gwtp, hare, hdsp, hdsv, hdxm
     >      , hdxp, hgdx, hgdy, hpit, hvol, hvsb, hwtm, hwtp, icel
     >      , icmax, icmin, jcel, jnxm, jnxp, jtmax, jtmin, kgdx, kgdy
     >      , kreg, romn, vlmn
     >      , stat = istat )
          if( istat /= 0 ) then
            write(cmsg,'(a,i6)')
     >        'gare deallocate error in cplmet_alc, istat = ', istat
            call wexit( 'cplmet', trim( cmsg ) )
          endif
        endif
      end select

      return
      end
