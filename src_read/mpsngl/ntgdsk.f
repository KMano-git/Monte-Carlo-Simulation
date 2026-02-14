!**********************************************************************
      subroutine ntgdsk(nft,cact)
!**********************************************************************
      use cntcom,      only : bvx, bvy, bvz, chwl, cswl, cwal, e0wl
     >    , grax, grxp, gzax, gzxp, iaxs, icgt, icwl, igtw, ihwl, ikwl
     >    , impl, inwl, ipgt, iplx, iply, ipmp, ipwl, ispx, isxw, isyw
     >    , ivb1, ivb2, ivs1, ivs2, iwgt, iwl1, iwl2, ixyw, jdp1, jdp2
     >    , jsl1, jsl2, jvb1, jvb2, jxp1, jxp2, mcel, mclx, mcly, mcly2
     >    , mgrd, migx, migy, mknd, mrgn, mseg, mxjw, mxjw2, myit, ncmax
     >    , ncmax2, negt, next, ngmax, ngmax2, nhwl, nocl, nogd
     >    , nogdv, nogt, nowl, nowl2, npew, npmp, npsw, npwl, npwl2
     >    , nsgt, rcwl, rcywl, snwl, temwl, tywl, volm, xegt, xpnt, xsgt
     >    , yegt, ypnt, ysgt
      use csize,       only : ndad, ndgt, ndgtp, ndmc, ndmg, ndms, ndpm
     >    , ndwh, ndwl, ndwp, ndx, ndy, nvxd, nvyd
      use cunit,       only : lmspe, lmype, n6
      use mpi!,         only : mpi_bcast
      implicit none
!
!::argument
      integer,   intent(in) :: nft
      character, intent(in) :: cact*(*)
!
!::local variables
      integer   nndmg, nndmc, nndms, nndx, nndy, nnvxd, nnvyd
      integer   nndwl, nndwh, nndwp
      integer   nndgt, nndgtp, nndad, nndpm
      integer   iw, ic, mmigx, mmigy, ks, kg, ngs, nge, ng
      integer   ip, icx, icy, i
      real*8    bdir
      character cver*80

! functions
      integer   lenx
!
!::dummy
      integer noclv(0:nvxd,0:nvyd)
!
      noclv = 0
!
!::version
      cver = "/home/shimizu/Gmesh/soc/ntgdsk.f   04/04/30"
!
      write(n6,'(/2x,"*** ntgdsk ***  ",a)') cact(1:lenx(cact))
      write(n6,'(2x,"version : ",a)') cver(1:lenx(cver))
!
      if( lmype.ne.lmspe ) goto 1000
!
!::master
      if( nft.le.0 ) then
      write(n6,'(2x,"nft =",i3)') nft
      call wexit("ntgdsk","nft.le.0")
      endif
!
!-----------------------------------------------------------------------
!::write
!-----------------------------------------------------------------------
      if( cact(1:1).eq."w" ) then
      write(nft) cver
!-----
      write(nft) ndmg, ndmc, ndms, ndx, ndy, nvxd, nvyd
      write(nft)
     >    xpnt,ypnt,volm,bvx,bvy,bvz,grax,gzax,grxp,gzxp
     >   ,mclx,mcly,mcly2,mxjw,mxjw2,myit
     >   ,nogd,nocl,mcel,nogdv,noclv,mseg,mgrd,mknd
     >   ,mrgn,migx,migy,next,iplx,iply
     >   ,ngmax,ngmax2,ncmax,ncmax2
     >   ,jdp1,jxp1,jsl1,jsl2,jxp2,jdp2
     >   ,iwl1,ispx,iwl2,impl,iaxs,ivs1,ivs2
     >   ,jvb1, jvb2, ivb1, ivb2
!----
      write(nft) ndwl, ndwh, ndwp
      write(nft)
     >   rcywl, temwl, rcwl, e0wl, cswl, snwl
     >  ,cwal, chwl, tywl
     >  ,npsw, npew, inwl, ikwl, ipwl
     >  ,icwl, isxw, isyw, ixyw, igtw
     >  ,ihwl, ipmp
     >  ,nowl, npwl, npmp, nowl2, npwl2, nhwl
!-----
      write(nft) ndgt, ndgtp, ndad, ndpm
      write(nft)
     >   xsgt,xegt,ysgt,yegt
     >  ,nsgt,negt,iwgt,ipgt,icgt,nogt
      endif
!
!-----------------------------------------------------------------------
!::read
!-----------------------------------------------------------------------
      if( cact(1:1).eq."r" ) then
      read(nft) cver
!
      write(n6,'(2x,"version : ",a)') cver(1:lenx(cver))
!-----
      read(nft) nndmg, nndmc, nndms, nndx, nndy, nnvxd, nnvyd
!
      write(n6,'(2x,"dimension size (grid)")')
      write(n6,'(2x,"ndmg/dsk =",2i5,"  ndmc/dsk =",2i5,
     > "  ndms/dsk =",2i5)') ndmg, nndmg, ndmc, nndmc, ndms, nndms
      write(n6,'(2x,"ndx/dsk  =",2i5,"  ndy/dsk  =",2i5,
     > "  nvxd/dsk =",2i5,"  nvyd/dsk =",2i5)')
     > ndx,nndx, ndy,nndy, nvxd,nnvxd,nvyd,nnvyd
      if( ndmg.ne.nndmg .or. ndmc.ne.nndmc .or. ndms.ne.nndms .or.
     >    ndx.ne.nndx   .or. ndy.ne.nndy   .or.
     >    nvxd.ne.nnvxd .or. nvyd.ne.nnvyd ) goto 910
!
      read(nft)
     >    xpnt,ypnt,volm,bvx,bvy,bvz,grax,gzax,grxp,gzxp
     >   ,mclx,mcly,mcly2,mxjw,mxjw2,myit
     >   ,nogd,nocl,mcel,nogdv,noclv,mseg,mgrd,mknd
     >   ,mrgn,migx,migy,next,iplx,iply
     >   ,ngmax,ngmax2,ncmax,ncmax2
     >   ,jdp1,jxp1,jsl1,jsl2,jxp2,jdp2
     >   ,iwl1,ispx,iwl2,impl,iaxs,ivs1,ivs2
     >   ,jvb1, jvb2, ivb1, ivb2
!----
      read(nft) nndwl, nndwh, nndwp
!
      write(n6,'(2x,"dimension size (wall)")')
      write(n6,'(2x,"ndwl =",2i5,"  ndwh =",2i5,"  ndwp =",2i5)')
     >   ndwl,nndwl,ndwh,nndwh,ndwp,nndwp
      if( ndwl.ne.nndwl .or. ndwh.ne.nndwh .or. ndwp.ne.nndwp ) goto 910
!
      read(nft)
     >   rcywl, temwl, rcwl, e0wl, cswl, snwl
     >  ,cwal, chwl, tywl
     >  ,npsw, npew, inwl, ikwl, ipwl
     >  ,icwl, isxw, isyw, ixyw, igtw
     >  ,ihwl, ipmp
     >  ,nowl, npwl, npmp, nowl2, npwl2, nhwl
!-----
      read(nft) nndgt, nndgtp, nndad, nndpm
!
      write(n6,'(2x,"dimension size (gate)")')
      write(n6,'(2x,"ndgt =",2i5,"  ndgtp =",2i5,"  ndad =",2i5,
     > "  ndpm =",2i5)') ndgt,nndgt,ndgtp,nndgtp,ndad,nndad,ndpm,nndpm
      if( ndgt.ne.nndgt .or. ndgtp.ne.nndgtp .or. ndad.ne.nndad .or.
     >    ndpm.ne.nndpm ) goto 910

      read(nft)
     >   xsgt,xegt,ysgt,yegt
     >  ,nsgt,negt,iwgt,ipgt,icgt,nogt
!
!::KSFUJI  confirm    mcel(j,i) = 0 for dummy cell

      endif
!
!-----------------------------------------------------------------------
!::correction of b-direction
!-----------------------------------------------------------------------
      write(n6,'(2x,"[direction of b-field]")')
!-----
      call ntbdir(bdir)
!-----
      write(n6,'(2x)')
      if( bdir.lt.0.0 ) then
      write(n6,'(2x,"[change b-direction so as to b*e_sita > 0")')
      do i = 0, ndmc
      bvx(i) = -bvx(i)
      bvy(i) = -bvy(i)
      enddo
      endif
!
!::KSFUJI correction of gate
      call edtgate
      call lstgate
!
!----------------------------------------------------------------------
!::send ntgdsk
!----------------------------------------------------------------------
 1000 continue

!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      if( cact(1:1).eq."w" ) return
!
!--grid
      write(n6,'(5x,"jdp1 =",i3,"  jxp1 =",i3,"  jsl1 =",i3,
     >  "  jsl2 =",i3,"  jxp2 =",i3,"  jdp2 =",i3)')
     >   jdp1,jxp1,jsl1,jsl2,jxp2,jdp2
      write(n6,'(5x,"ivs1 =",i3,"  iwl1 =",i3,"  ispx =",i3,
     >  "  iwl2 =",i3,"  ivs2 =",i3,"  impl =",i3,"  iaxs =",i3)')
     >   ivs1,iwl1,ispx,iwl2,ivs2,impl,iaxs
      write(n6,'(5x,"jvb1 =",i3,"  jvb2 =",i3,"  ivb1 =",i3,
     >  "  ivb2 =",i3)') jvb1,jvb2,ivb1,ivb2
      write(n6,'(5x,"axis =",2f8.4,"  xpnt =",2f8.4)')
     >   grax,gzax,grxp,gzxp
      write(n6,'(5x,"nowl =",3i6,"  npwl =",3i6,"  npmp =",i5,
     > "  nhwl =",i3)') nowl,nowl2,ndwl,npwl,npwl2,ndwp,npmp,nhwl
!--wall
      write(n6,'(6x,"iw",4x,"knd",3x,"ip",4x,"ic",4x,"igx",3x,
     >   "igy",3x,"ipx",3x,"ipy",3x,"xp",6x,"yp",6x,"jdp",3x,"idp",
     >   3x,"rcy",3x,"ivid",2x,"igtw")')
      do iw = 1, npwl2, 50
      ic = icwl(iw)
      mmigx = 0
      mmigy = 0
      if( ic.ne.0 ) then
        mmigx = migx(ic)
        mmigy = migy(ic)
      endif
      write(n6,'(2x,8i6,2f8.3,2i6,f8.3,i6,2x,a,2x,a)')
     > iw,ikwl(iw),ipwl(iw),icwl(iw),mmigx,mmigy,
     > iplx(icwl(iw)),iply(icwl(iw)),xpnt(ipwl(iw)),ypnt(ipwl(iw)),
     > isxw(iw),isyw(iw),rcwl(iw),igtw(iw),tywl(iw)
     >,chwl(ihwl(iw))
      enddo
!--gate
      write(n6,'(2x,"nogt =",i3)') nogt
      do ks = 1, 2
      if( ks.eq.1 ) write(n6,'(2x,a)') "gate in plasma wall"
      if( ks.eq.2 ) write(n6,'(2x,a)') "gate in vacume wall"
      do kg = 1, nogt
      ngs = nsgt(kg,ks)
      nge = negt(kg,ks)
      write(n6,'(2x,i2,2x,2i4,2x,2f9.4,2x,2f9.4)') kg
     >  ,ngs,nge,xsgt(kg,ks),ysgt(kg,ks),xegt(kg,ks),yegt(kg,ks)
      if( ngs.eq.0 ) cycle
      do ng = ngs, nge
      iw = iwgt(ng)
      ip = ipgt(ng)
      ic = icgt(ng)
      icx = 0; icy = 0
      if( ic.gt.0 ) then
      icx = migx(ic); icy = migy(ic)
      endif
      write(n6,'(2x,3i6,2f10.4,3i6)')
     >   ikwl(iw),iw,ip,xpnt(ip),ypnt(ip),ic,icx,icy
      enddo
      enddo
      enddo
!
      return
!
!::error
 910  continue
      call wexit("ntgdsk","dimension size error")
      end
