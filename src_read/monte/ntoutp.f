!**********************************************************************
      subroutine ntoutp(nft,ctyp)
!**********************************************************************
      use cntcom, only : bvx, bvy, bvz, cstyp, cwal, den0, deng, eng0
     >    , engg, flbr, flxin, icbr, icwl, ikbr, iplx, iply, mcel, migx
     >    , migy, nbr, npew, npsw, prbr, tfbr, vlp0
      use cntmnt, only : dotn, dotn2, sdotn, sdotn2, sn0, ssn, ssp
     >    , sumsn, sumsp, sumwe, sumwi, swe, swi, totsn, totsp, totwe
     >    , totwi
      use cntpls, only : dene, teme, temi, vflw
      use cntwcn, only : wabs, wend, werr, wion, wnrm, wpmp, wreg, wsbr
     >    , wsum, wssn, wssp, wswe, wswi, wtot
      use cplcom, only : nion, trad
      use cplmet, only : hvol, icel, itmax, itpve, itsle, jcel, jtmax 
      use csize,  only : ndsp
      use csonic, only : itim
      use cunit,  only : n6
      implicit none

!::arguments
      integer,   intent(in) :: nft
      character, intent(in) :: ctyp*(*)
!
!::local varaibles
      integer  i, ig, j, jout, ic, icx, icy, ia, it, jt, jmax
      integer  nw, iws, iwe, iw, ll, jw
      real*8   zsn(ndsp),zsp(ndsp),zwi,zwe
      real*8   zsum, zflx, zprb
! function
      integer   lenx
!
      write(nft,'(/2x,"*** ntoutp ***  (",a,")","  itim =",i8,
     > "   cstyp =",a,"  dotn =",1pe10.3)')
     >  ctyp,itim,cstyp(1:lenx(cstyp)),dotn
!
!----------------------------------------------------------------------
!::conservation of weit
!----------------------------------------------------------------------
      if( index(ctyp,'conv').gt.0 ) then
      write(nft,'(2x,"wtot =",1p2e14.6,"  werr =",1pe14.6,"  wion =",
     > 1p2e14.6,"  wabs =",1pe14.6,"  wpmp =",1pe14.6)')
     >  wtot,wsum, werr, wion,wion, wabs, wpmp
      write(nft,'(2x,"wreg =",1pe14.6,"  od/sl/id =",1p3e12.4,
     > "  op/ip =",1p2e12.4,"  ed/hc =",1p2e12.4,"  vc/?? =",1p2e12.4)')
     >  wion, wreg(1), wreg(2), wreg(3), wreg(4), wreg(5), wreg(6),
     >  wreg(7), wreg(8)+wreg(9), wreg(0)
      write(nft,'(2x,"wend =",1p5e14.6)') (wend(i),i=1,5)
!
      endif
!
!----------------------------------------------------------------------
!::den0,eng0
!----------------------------------------------------------------------
      if( index(ctyp,'den0').gt.0 ) then
      ig = 1
      do it = 1, itmax
      if( it.ne.25 ) cycle
      jmax = jtmax(it)
      write(nft,'(3x,2x,"i",3x,"j",4x,"ic",2x,"icx",1x,"icy",4x,
     >   "dene",8x,"teme",8x,"temi",8x,"den0",8x,"eng0",8x,"vlp0"
     >   ,8x,"totn0",7x,"bvx",9x,"bvy",9x,"bvz")')
      do jw = 1, jmax
      j = jcel(jw,it)
      i = icel(jw,it)
      jout = 0
      ic = mcel(j,i)
      if( ic.le.0 ) cycle
      icx = iplx(ic)
      icy = iply(ic)
      write(nft,'(2x,2i4,i6,2i4,1p7e12.3,1p3e11.2)')
     >  i,j,ic,icx,icy,dene(ic),teme(ic),temi(ic)
     > ,den0(ic,ig),eng0(ic,ig),vlp0(ic,ig),sn0(icx,icy,ig)
     > ,bvx(ic),bvy(ic),bvz(ic)
      if( dabs(vlp0(ic,ig)).gt.1.0d6 ) then
        write(nft,'(2x,"Error vlp0 > 1.0d6 ")')
      endif
      enddo
      enddo
      endif
!
!----------------------------------------------------------------------
!::Sp momentum source
!----------------------------------------------------------------------
      if( index(ctyp,'smon').gt.0 ) then
      ig = 1
      ia = 1
      do ll = 1, 2
      if( ll.eq.1 ) it = itsle
      if( ll.eq.2 ) it = 17
      jmax = jtmax(it)
!
      write(nft,'(3x," jt",3x,"ic",2x,"icx",1x,"icy",4x,
     >   "dene",8x,"teme",8x,"temi",8x,"vflw",8x,"den0",8x,"deng",8x
     >   ,"tN0",9x,"tSn",9x,"tSp")')

      do jt = 2, jmax-1
      icx = jcel(jt,it)
      icy = icel(jt,it)
      jout = 0
      if( jt.le.30 .or. jt.gt.jmax-30 ) jout = 1
      if( jout.eq.0 ) cycle
      ic = mcel(icx,icy)
      if( ic.le.0 ) cycle
      write(nft,'(2x,i4,i6,2i4,1p9e12.3)')
     >  jt,ic,icx,icy,dene(ic),teme(ic),temi(ic),vflw(ic,ig)
     > ,den0(ic,ig),deng(ic,1)
     > ,sn0(icx,icy,ia),ssn(icx,icy,ia),ssp(icx,icy,ia)
      enddo
!
      enddo
      endif
!
!----------------------------------------------------------------------
!::total source terms
!----------------------------------------------------------------------
      if( index(ctyp,'tots').gt.0 ) then
      ia = 1
      write(n6,'(2x,"NSourc dotn =",1pe15.7,"  sn,sp,wi,we,wr =",
     >  1p5e15.7)') dotn2,sumsn(10,ia),sumsp(10,ia),sumwi(10),sumwe(10)
      write(n6,'(2x,"NTotal dotn =",1pe15.7,"  sn,sp,wi,we,wr =",
     >  1p5e15.7)') sdotn2,totsn(10,ia),totsp(10,ia),totwi(10),
     >  totwe(10),trad
      return
!
      do 310 ia = 1, nion
      zsn(ia) = 0.0d0
      zsp(ia) = 0.0d0
 310  continue
      zwi = 0.0d0
      zwe = 0.0d0
!
      do 320 it = 2, itmax-1
      if( it.eq.itpve ) goto 320
      jmax = jtmax(it)
      do 330 jt = 2, jmax-1
      icx = jcel(jt,it)
      icy = icel(jt,it)
      zwi = zwi + swi(icx,icy)*hvol(icx,icy)
      zwe = zwe + swe(icx,icy)*hvol(icx,icy)
      do 340 ia = 1, nion
      zsn(ia) = zsn(ia) + ssn(icx,icy,ia)*hvol(icx,icy)
      zsp(ia) = zsp(ia) + ssp(icx,icy,ia)*hvol(icx,icy)
 340  continue
 330  continue
 320  continue
!
      write(n6,'(2x,"Integ dotn =",1pe15.7,"  sn,sp,wi,we,wr =",
     >  1p5e15.7)') sdotn,zsn(ia),zsp(ia),zwi,zwe,trad
      endif
!
!----------------------------------------------------------------------
!::denisity at wall
!----------------------------------------------------------------------
      if( index(ctyp,'wln0').gt.0 ) then
      ig = 1
      do nw = 1, 4
      iws = npsw(nw)
      iwe = npew(nw)
      write(nft,'(/2x,"*** den0 at wall ***  ",a,2x,a)') cwal(nw),cstyp
      write(nft,'(3x,2x,"nw",2x,"iw",4x,"ic",2x,"icx",1x,"icy",4x,
     >   "dene",8x,"teme",8x,"temi",8x,"den0",8x,"eng0",8x,"vlp0"
     >   ,8x,"totn0",7x,"deng",7x,"engg")')
      do iw = iws, iwe
      ic  = icwl(iw)
      if( ic.le.0 ) cycle
      icx = iplx(ic)
      icy = iply(ic)
      write(nft,'(2x,2i4,i6,2i4,1p7e12.3,1p3e11.2)')
     >  nw,iw,ic,icx,icy,dene(ic),teme(ic),temi(ic)
     > ,den0(ic,ig),eng0(ic,ig),vlp0(ic,ig),sn0(icx,icy,ig)
     > ,deng(ic,1), engg(ic,1)
      enddo
      enddo
      endif
!
!----------------------------------------------------------------------
!::integral Sp in inner divertor
!----------------------------------------------------------------------
      if( index(ctyp,'dbsp').gt.0 ) then
      call ntdbsp(nft)
      endif
!
!----------------------------------------------------------------------
!::source terms at a certain cell
!----------------------------------------------------------------------
      if( index(ctyp,'sorc').gt.0 ) then
      ig = 1
      write(nft,'(3x,2x,"ic",2x,"icx",1x,"icy",4x,
     >   "dene",8x,"teme",8x,"temi",8x,"den0",8x,"eng0",8x,"vlp0"
     >   ,8x,"totn0",7x,"ssn",9x,"ssp",9x,"swe",9x,"swi")')
      ic = 5249
      icx = iplx(ic)
      icy = iply(ic)
      write(nft,'(2x,i6,2i4,1p7e12.3,1p4e11.2)')
     >  ic,icx,icy,dene(ic),teme(ic),temi(ic)
     > ,den0(ic,ig),eng0(ic,ig),vlp0(ic,ig),sn0(icx,icy,ig)
     > ,wssn(ic,ig),wssp(ic,ig),wswe(ic),wswi(ic)
      endif
!
!----------------------------------------------------------------------
!::emitted neutral flux
!----------------------------------------------------------------------
      if( index(ctyp,'emfl').gt.0 ) then
      write(nft,'(2x,"cstyp =",a,"  tflx =",1p2e12.3,"  wtot =",
     >  1pe12.3)') cstyp(1:lenx(cstyp)),tfbr,flxin,wtot
      ia = 1
      zsum = 0.0d0
      do i = 1, nbr
      ic = icbr(i)
      zflx = wsbr(ic,ia)/wnrm*flxin
      zsum = zsum + wsbr(ic,ia)
      zprb = zsum/wnrm
      write(nft,'(2x,i4,i3,i6,2i4,1p4e12.3)') i,ikbr(i),icbr(i)
     > ,migx(icbr(i)),migy(icbr(i)),flbr(i),zflx,prbr(i),zprb
      enddo
      endif
!
      return
      end
