!**********************************************************************
      subroutine ntwlst
!**********************************************************************
      use cntcom, only : chwl, cwal, e0wl, icwl, igtw, ihwl, ikwl, inwl
     >    , iplx, iply, ipmp, ipwl, isxw, isyw, ixyw, mknd, mseg, ncmax2
     >    , next, nhwl, nopm, nowl2, npew, npmp, npsw, npwl2, nwpm, rcwl
     >    , tywl, xpnt, ypnt
      use cplmet, only : itpve, itpvs, itsle, itsls, jtmax
      use csize,  only : ndpm
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  nw, iws, iwe, iw, ic, i, k, ii
      integer  ih
      real*8   zrcymx, zrcymn
! function
      integer  imox, imoy
!
      write(n6,'(/2x,"*** ntwlst ***")')
!
!::wall type
      do nw = 1, nowl2
      iws = npsw(nw)
      iwe = npew(nw)
      zrcymx = -1.0d20
      zrcymn = +1.0d20
      do iw = iws, iwe
      if( (iw.eq.iwe) .and. (iws.ne.iwe) ) cycle
      zrcymx = dmax1( zrcymx, rcwl(iw) )
      zrcymn = dmin1( zrcymn, rcwl(iw) )
      enddo
      write(n6,'(2x,"nw =",i2,2x,a,2x,2i5,"  rcy =",2f9.4)')
     >   nw,cwal(nw),iws,iwe,zrcymn,zrcymx
      enddo
!
!::wall type
      write(n6,'(2x,"jmax(sol) =",i5,"  (prv) =",i5,"  itsls =",i5,
     >   "  itpve =",i5)') jtmax(itsle),jtmax(itpvs),itsls,itpve
      write(n6,'(2x)')
      do nw = 1, nowl2
      iws = npsw(nw)
      iwe = npew(nw)
      write(n6,'(2x,"nw =",i2,2x,a,2x,2i5,"  inwl =",2i4)')
     >  nw,cwal(nw),iws,iwe,inwl(iws),inwl(iwe)
      write(n6,'(6x,"iw",4x,"knd",3x,"ip",4x,"ic",4x,"igx",3x,
     >   "igy",3x,"ipx",3x,"ipy",3x,"xp",6x,"yp",6x,"jdp",3x,"idp",
     >   3x,"rcy",3x,"e0w",5x,"igtw",2x,"inwl",2x,"ihwl")')
      do iw = iws, iwe
      ic = icwl(iw)
      write(n6,'(2x,8i6,2f8.3,2i6,2f8.3,2i6,2x,a,2x,a)')
     > iw,ikwl(iw),ipwl(iw),icwl(iw),imox(ic),imoy(ic),
     > iplx(icwl(iw)),iply(icwl(iw)),xpnt(ipwl(iw)),ypnt(ipwl(iw)),
     > isxw(iw),isyw(iw),rcwl(iw),e0wl(iw),igtw(iw)
     >,inwl(iw),          tywl(iw),chwl(ihwl(iw))
      enddo
      enddo
!
!::albedo wall
      ii = 0
      do ih = 1, nhwl
      if( chwl(ih)(1:1).eq."a" ) then
      ii = ii + 1
      if( ii.gt.ndpm ) goto 920
      nwpm(ii) = ih
      endif
      enddo
      nopm = ii
!
      write(n6,'(/2x,"albedo wall  nopm =",i3)') nopm
      do i = 1, nopm
      ih = nwpm(i)
      zrcymx = 0.0d0
      zrcymn = 1.0d0
      iws = +999
      iwe = -999
      do iw = 1, npwl2
      if( ihwl(iw).ne.ih ) cycle
      iws = min0( iws, iw )
      iwe = max0( iws, iw )
      if( rcwl(iw).ge.1.0d0 ) cycle
      if( rcwl(iw).le.0.0d0 ) cycle
      zrcymx = dmax1( zrcymx, rcwl(iw) )
      zrcymn = dmin1( zrcymn, rcwl(iw) )
      enddo
      write(n6,'(2i5,2x,a,2x,2i5,"  rcy =",2f9.4)')
     > i, ih, chwl(ih), iws, iwe, zrcymn, zrcymx
      enddo
!
!::albedo tile
      write(n6,'(/2x,"albedo wall segment   npmp =",i3)') npmp
      do i = 1, npmp
      iw = ipmp(i)
      ic = icwl(iw)
      write(n6,'(2x,8i6,2f8.3,2i6,2f8.3,2i6,2x,a,2x,a)')
     > iw,ikwl(iw),ipwl(iw),icwl(iw),imox(ic),imoy(ic),
     > iplx(icwl(iw)),iply(icwl(iw)),xpnt(ipwl(iw)),ypnt(ipwl(iw)),
     > isxw(iw),isyw(iw),rcwl(iw),e0wl(iw),ixyw(iw),igtw(iw)
     >,tywl(iw),chwl(ihwl(iw))
      enddo
!
!::cell conection
      write(n6,'(/2x,"cell conection for gate or albedo or shebron")')
      write(n6,'(2x,5x,"next(ic,ik)  when k=ikwl")')
      write(n6,'(2x,1x,"nw",2x,"tywl",5x,"iw",3x,"xpnt",5x,"ypnt",2x,
     >  "igtw",5x,"ic",3x,"mx",3x,"my",2x,"msg",3x,"ik",4x,"mx1",2x,
     >  "my1",4x,"mx2",2x,"my2",4x,"mx3",2x,"my3",4x,"mx4",2x,"my4")')
      do nw = 1, nowl2
      iws = npsw(nw)
      iwe = npew(nw)
      do iw = iws, iwe
      ic = icwl(iw)
      if( igtw(iw).eq.0 ) cycle
      if( ic.gt.0 ) then
      write(n6,'(2x,i3,2x,a,2x,i5,2f9.4,i4,i7,2i5,2i5,4(2x,"/",3i5))')
     > nw, tywl(iw), iw, xpnt(ipwl(iw)), ypnt(ipwl(iw)), igtw(iw),
     > ic, imox(ic), imoy(ic), mseg(ic), ikwl(iw),
     > (mknd(ic,k),imox(next(ic,k)),imoy(next(ic,k)),k=1,4)
      else
      write(n6,'(2x,i3,2x,a,2x,i5,2f9.4,i4,i7,2i5,2i5,4(2x,"/",3i5))')
     > nw, tywl(iw), iw, xpnt(ipwl(iw)), ypnt(ipwl(iw)), igtw(iw),
     > ic, imox(ic), imoy(ic), 0, ikwl(iw)
      endif
      enddo
      enddo
!
!::cell connection
      write(n6,'(/2x,4x,"ic",4x,"ix",3x,"iy",3x,"k",4x,"iw",4x,
     >  "tywl",3x,"icwl",3x,"igtw",3x,"rcwl",5x,"e0wl",5x,"ip",4x
     >  "xp",7x,"yp")')
      do ic = 1, ncmax2
      do k = 1, mseg(ic)
      if( mknd(ic,k).ge.0 ) cycle
      iw = -mknd(ic,k)
      if( index("pgas",tywl(iw)(1:1)).gt.0 .or.
     >    index("pgas",tywl(iw)(3:3)).gt.0 ) then
      write(n6,'(2x,i7,2i5,i5,i7,2x,a,2x,i7,i4,2f9.4,i7,2f9.4)')
     >  ic, imox(ic),imoy(ic),
     >    k, iw, tywl(iw), icwl(iw), igtw(iw), rcwl(iw), e0wl(iw),
     >    ipwl(iw), xpnt(ipwl(iw)), ypnt(ipwl(iw))
      endif
      enddo
      enddo
!
      return
!
 920  continue
      call wexit("ntwlst","too many albedo   ii.gt.ndpm")
      end
!
!**********************************************************************
      integer function imox(ic)
!**********************************************************************
      use cntcom, only : migx
      implicit none

!::local variables
      integer, intent(in) :: ic
!
      if( ic.le.0 ) then
        imox = 0
      else
        imox = migx(ic)
      endif
!
      return
      end
!
!**********************************************************************
      integer function imoy(ic)
!**********************************************************************
      use cntcom, only : migy
      implicit none

!::local variables
      integer, intent(in) :: ic
!
      if( ic.le.0 ) then
        imoy = 0
      else
        imoy = migy(ic)
      endif
!
      return
      end
