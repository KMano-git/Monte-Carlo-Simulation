!**********************************************************************
      subroutine set_edrfl
!**********************************************************************
!
!     set reflection point with aro in the main plasma
!
!----------------------------------------------------------------------
      use cimcom, only : frf, ipcm, lpcm, ndrf, nrf, rf1_c, rf1_r
     >    , rf1_ro, rf1_z, rf2_c, rf2_r, rf2_ro, rf2_z
      use cplmet, only : icaxs, icspx
      use cpmpls, only : romp
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer i, i1
      real*8  x1a, y1a, x1b, y1b, x2a, y2a, x2b, y2b
      real*8  wkx(ndrf), wky(ndrf); integer wkc(ndrf)
      integer  ker, ic1, ic2, ko1, ko2, md
! function
      integer    imox, imoy
!
      write(n6,'(/2x,"*** set_edrfl ***")')
      write(n6,'(2x,"icspx, icaxs =",2i5)') icspx, icaxs
      write(n6,'(2x,"romp =",1p10e12.3)') (romp(i),i=1,icaxs)
!
!::avoid warning message for undefined variables lpcm, ipcm
      lpcm = -1
      ipcm = -1
!
      write(n6,'(2x,"nrf =",i4,"  ro =",2f10.5)') nrf, rf1_ro, rf2_ro
      if( nrf.gt.ndrf ) then
      write(n6,'(2x,"dimension error  nrf+1.gt.ndrf ",2i6)')
     >   nrf+1, ndrf
      call wexit( 'set_edrfl', 'array size over' )
      stop
      endif
!
      call posman(rf1_ro, nrf, wkx, wky, wkc)
      do i = 1, nrf
      rf1_r(i) = wkx(i)
      rf1_z(i) = wky(i)
      rf1_c(i) = wkc(i)
      enddo
      rf1_r(nrf+1) = rf1_r(1)
      rf1_z(nrf+1) = rf1_z(1)
      rf1_c(nrf+1) = rf1_c(1)
!
      call posman(rf2_ro, nrf, wkx, wky, wkc)
      do i = 1, nrf
      rf2_r(i) = wkx(i)
      rf2_z(i) = wky(i)
      rf2_c(i) = wkc(i)
      enddo
      rf2_r(nrf+1) = rf2_r(1)
      rf2_z(nrf+1) = rf2_z(1)
      rf2_c(nrf+1) = rf2_c(1)
!
      frf = nrf
!
!::debug write
      md = nrf/5
      if( md.eq.0 ) md = 1
      write(n6,'(/2x,"*** set_edrfl ***")')
      ker = 0
      write(n6,'(2x,3x,"i",3x,"xp1",7x,"yp1",7x,"ln1",6x,"ix1",2x,"iy1",
     > 2x,"k1",3x,"xp2",7x,"yp2",7x,"ln2",6x,"ix2",2x,"iy2",2x,"k2")')
      do i = 1, nrf+1
      i1 = i + 1
      if( i.eq.nrf+1 ) i1 = 2
      x1a = rf1_r(i); y1a = rf1_z(i); x1b = rf1_r(i1); y1b = rf1_z(i1)
      x2a = rf2_r(i); y2a = rf2_z(i); x2b = rf2_r(i1); y2b = rf2_z(i1)
      ic1 = rf1_c(i)
      ic2 = rf2_c(i)
      call mchkin(x1a,y1a,ic1,ko1)
      call mchkin(x2a,y2a,ic2,ko2)
      ker = ker + ko1 + ko2
!
      if( mod(i,md).eq.0 .or. ko1.ne.0 .or. ko2.ne.0 ) then
      write(n6,'(2x,i4,3f10.5,2i5,i4,3f10.5,2i5,i4)') i
     > ,x1a,y1a,sqrt((x1b-x1a)**2+(y1b-y1a)**2),imox(ic1),imoy(ic1),ko1
     > ,x2a,y2a,sqrt((x2b-x2a)**2+(y2b-y2a)**2),imox(ic2),imoy(ic2),ko2
      endif
      enddo
      write(n6,'(2x,"ker =",i3)') ker
!
      if( ker.ne.0 ) then
        write(n6,'(2x,"Stop at sub. set_edrfl")')
        call wexit( 'set_edrfl', 'ker > 0' )
      endif
!
      return
      end
