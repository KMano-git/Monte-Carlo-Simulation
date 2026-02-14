!***********************************************************************
      subroutine pst_wrad(kk)
!***********************************************************************
      use cntcom, only : cfbr, chwl, flbr, flxin, icbr, icwl, ihwl
     >    , iplx, iply, iptl, ipwl, ltrc, migx, migy, mrgn, nbr, ncmax
     >    , nhwl, npew, npsw, prbr, tfbr, volm, xpnt, ypnt
      use cntwcn, only : whta, wtot, wwal
      use cntwfl, only : xare, xdeg
      use cplcom, only : nsmp_wrad, wrad, vne, vte
      use cplpst, only : dwrad, nphtn
      use cplwrd, only : wime
      use csize,  only : ndbr, ndmc, ndwp, ndxy
      use csonic, only : limp
      use cunit,  only : n6
      use mpi!,    only : mpi_wtime
      use mod_externalgrid
      implicit none
!
!::argument
      integer, intent(in) :: kk
!
!::local variables
      real*8  hwrd(ndmc), drad(ndmc), wlrad(ndwp)
      real*8  wbrw(ndmc), whtw(ndmc,8)
      real*8  fclim, zmax, zsum, zrad
      integer ii, j, i, nsiz, ic, imd, n6w
      integer icx, icy, nsamp, icon
      integer  nw, iw, iws, iwe, ip, ip1
      real*8  dradmx, zvl, ate, ane, awr
      real*8  sum, dw
      real*8  xst, yst, xen, yen, dlen
      real*8  x1, y1, x2, y2, x0, y0, th
      real*8  smrad, sum1, sum2, sum3
      real*8  cptm0, cptm1
      integer icst, icen, iwen, ih, ir
      real*8  tot1, tot2
!
      write(n6,'(/2x,"*** pst_wrad ***  kk =",i2)') kk
!
!::KSFUJI
      if( kk.eq.0 ) then
      dwrad(1:ndxy,1:4) = 0.0d0
      return
      endif
!
!::output file
      n6w = 21
      open(unit=n6w,file="zdpwrad")
!
      write(n6w,'(/2x,"*** pst_wrad ***  radiation to the plate")')
!
!::input data
      fclim = 1.0d-3
      nsamp = nsmp_wrad                      !KH 160420 see plinpt
      write(n6,'(2x,"nsamp = ",i9)') nsamp
!
      ltrc  = 0
      nphtn = nsamp
!
      cptm0 = MPI_WTIME()
!
!::wrad
      hwrd(1:ndmc) = 0.0d0
      if( limp.eq.0 ) then
        ii = 0
        do ic = 1, ncmax
          j = iplx(ic)
          i = iply(ic)
          if( j.le.0 .or. i.le.0 ) cycle
          ii = ii + 1
          hwrd(ic) = wrad(j,i)
        enddo
        write(n6w,'(2x,"set hwrd from wrad(j,i) in soldor  ",2i7)')
     >   ncmax, ii
!
      else
        ii = 0
        do ic = 1, ncmax
          hwrd(ic) = wime(ic)
          ii = ii + 1
        enddo
        write(n6w,'(2x,"set hwrd from wime(ic) in impmc  ",2i7)')
     >   ncmax, ii
      endif
!
!----------------------------------------------------------------------
!::max value
!----------------------------------------------------------------------
      ii = 0
      zmax = -1.0d20
      zsum = 0.0d0
      do ic = 1, ncmax
        zrad = -hwrd(ic)*volm(ic)
        zsum = zsum + zrad
        zmax = dmax1( zmax, zrad )
        drad(ic) = zrad
        ii = ii + 1
      enddo
      dradmx = zmax
!
      write(n6w,'(2x,"No of cell =",2i6,"  tot-rad =",1pe12.3,
     >    "  max-rad =",1pe12.3)') ncmax, ii, zsum, zmax
!
      if( zsum.le.0.0d0 ) then
        write(n6w,*)"pst_wrad is skipped: sum_wrad =0.0"
        return
      endif
!
!----------------------------------------------------------------------
!::limit
!----------------------------------------------------------------------
      zsum = 0.0d0
      zmax = -1.0d20
      ii = 0
      do ic = 1, ncmax
        if( drad(ic)/dradmx .le. fclim ) then
          drad(ic) = 0.0d0
          cycle
        endif
        ii = ii + 1
        zsum = zsum + drad(ic)
        zmax = dmax1(zmax,drad(ic))
      enddo
      smrad = zsum
!
      write(n6w,'(2x,"fclim =",1pe12.3)') fclim
      write(n6w,'(2x,"No of cell =",2i6,"  tot-rad =",1pe12.3,
     >    "  max-rad =",1pe12.3)') ncmax, ii, zsum, zmax
!
!----------------------------------------------------------------------
!::clear /cntbir/
!----------------------------------------------------------------------
      cfbr = "phton"
      nbr  = 0
      tfbr = 0.0d0
      nsiz = ndbr
      call setd( prbr, nsiz, 0.0d0 )
      call setd( flbr, nsiz, 0.0d0 )
!
!----------------------------------------------------------------------
!::birth profile
!----------------------------------------------------------------------
      zsum = 0.0d0
      ii   = 0
      do ic = 1, ncmax
        zvl = drad(ic)
        if( zvl.le.0.0d0 ) cycle
        zsum = zsum + zvl
        ii = ii + 1
        if( ii.gt.ndbr ) goto 910
        flbr(ii) = zvl
        prbr(ii) = zsum
        icbr(ii) = ic
      enddo
!
!----------------------------------------------------------------------
!::normalization & max value
!----------------------------------------------------------------------
      nbr = ii
      tfbr = zsum
      zmax = -1.0d20
      do i = 1, nbr
        zmax = dmax1( zmax, flbr(i) )
        prbr(i) = prbr(i)/zsum
      enddo
      prbr(nbr) = 1.00001d0
!
      write(n6w,'(2x,"nbr =",i5,"  tfbr =",1pe12.3)') nbr, tfbr
!
!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
      write(n6w,'(2x)')
      write(n6w,'(5x,"i",3x,"ic",3x,"ix",2x,"iy",3x,"pr",10x,"sn",
     >  10x,"Ne",10x,"Te",10x,"Wr")')
      imd = max0( nbr/20, 1 )
      do i = 1, nbr
        if( i.eq.1 .or. i.eq.nbr .or. mod(i-1,imd).eq.0 ) then
          ic = icbr(i)
          icx = iplx(ic)
          icy = iply(ic)
          ate = vte(icx,icy)
          ane = vne(icx,icy)
          awr = drad(ic)/volm(ic)
          write(n6w,'(2x,i4,i6,2i4,1p6e12.3)')
     >      i, ic, migx(ic), migy(ic), prbr(i), flbr(i),
     >      ane, ate, awr
        endif
      enddo
!
!----------------------------------------------------------------------
!::clear scoreing variables
!----------------------------------------------------------------------
      flxin = 0.0d0
      call ntcler   ! set dotn
!
      wbrw(1:ndmc) = 0.0d0
      whtw(1:ndmc,1:8) = 0.0d0
!
!----------------------------------------------------------------------
!::monte calculation
!----------------------------------------------------------------------
      write(n6w,'(2x)')
!
      write(n6w,'(2x,"   including toloidal effect")')
      do i = 1, nsamp
        iptl = i
        call ntsta_teff
        call ntfolw_teff(xst,yst,xen,yen,icst,icen,iwen,dlen)
        wbrw(icst) = wbrw(icst) + 1.0d0
        ih = ihwl(iwen)
        if( ih.le.0 .or. ih.gt.8 ) cycle
        whtw(icst,ih) = whtw(icst,ih) + 1.0d0
      enddo
!
!----------------------------------------------------------------------
!::heat flux onto the wall
!----------------------------------------------------------------------
      write(n6w,'(2x)')
      write(n6w,'(2x,"nw",5x,"iw",5x,"ic",3x,"ix",3x,"iy",3x,
     >  "ipx",2x,"ipy",2x,"x0",8x,"y0",9x,"th",7x,"dare",8x,
     >  "dsmp",8x,"dwrd",8x,"are",8x,"smp",9x,"wrd")')
!
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
!
      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
        tot1 = 0.0d0
        tot2 = 0.0d0
!
        do iw = iws, iwe
          ip  = ipwl(iw)
          ic  = icwl(iw)
          if(use_exdata) then
            ip1 = ipwl2(iw)
            call get_position_byic(ic,ip,.true.,x1)
            call get_position_byic(ic,ip,.false.,y1)
            call get_position_byic(ic,ip1,.true.,x2)
            call get_position_byic(ic,ip1,.false.,y2)
          else
            if(iw==iwe) cycle
            ip1 = ipwl(iw+1)
            x1  = xpnt(ip)
            y1  = ypnt(ip)
            x2  = xpnt(ip1)
            y2  = ypnt(ip1)
          endif
          x0  = 0.5d0*(x1+x2)
          y0  = 0.5d0*(y1+y2)
          th  = xdeg(iw)
!
          wlrad(iw) = smrad/dfloat(nsamp)*wwal(iw)/xare(iw)
          sum1 = sum1 + xare(iw)
          sum2 = sum2 + wwal(iw)
          sum3 = sum3 + wlrad(iw)*xare(iw)
!
          tot1 = tot1 + wwal(iw)
          tot2 = tot2 + wlrad(iw)*xare(iw)
!
          write(n6w,'(2x,i2,2x,i5,i7,2i5,2i5,2f10.5,f10.3,1p6e12.3)')
     >      nw, iw, ic, migx(ic), migy(ic), iplx(ic), iply(ic),
     >      x0, y0, th,xare(iw), wwal(iw), wlrad(iw), sum1, sum2, sum3
        enddo
!
        write(n6w,'(2x,"sum-wwal =",1pe12.4,"  power =",1pe12.4,
     >    "  integ(wlrad) =",1pe12.4)') tot1, 
     >    smrad/dfloat(nsamp)*tot1, tot2
!
      enddo
      write(n6w,'(/2x,"are =",1pe12.3,"  wtot =",1p2e12.3,
     >  "  wrad =",1p2e12.3)') sum1, wtot, sum2, smrad, sum3
!
      write(n6w,'(/4x,"ic",3x,"ipx",2x,"ipy",2x,"irg",3x,"wbrw",
     >  2x,8(3x,a4,5x))') (trim(chwl(ih)),ih=1,8)
      do ic = 1, ncmax
      ir = mrgn(ic)
      if( ir.ne.1 .and. ir.ne.3 ) cycle
      if( wbrw(ic).le.200.0d0 ) cycle
      write(n6w,'(2x,i7,2i5,i4,2x,1pe12.3,2x,1p8e12.3)')
     >  ic,iplx(ic),iply(ic),mrgn(ic),wbrw(ic),(whtw(ic,ih),ih=1,8)
      enddo
!
!::whta
      write(n6w,'(/2x,"whta(in), whta(rfl), whta(pas/abs)")')
      do ih = 1, nhwl
      write(n6w,'(2x,i3,2x,a,2x,1p3e12.3)')
     >  ih, chwl(ih), whta(ih,1), whta(ih,2), whta(ih,3)
      enddo
!
      cptm1 = MPI_WTIME()
      write(n6w,'(/2x," total cpu = ", f10.3, " sec.","  nsamp =",i10)')
     >  cptm1-cptm0, nsamp
!
!----------------------------------------------------------------------
!::wlrad for monte ==> wl
!----------------------------------------------------------------------
      nphtn = nsamp
      dwrad(1:ndxy,1:4) = 0.0d0
!
      write(n6w,'(2x,"dwrad = wlrad*xare  Monte ==> Soldor")')
      write(n6w,'(2x,"nphtn =",i10)') nphtn
      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
        sum = 0.0d0
        do iw = iws, iwe
          if(.not.use_exdata .and. iw==iwe) cycle
          dw = wlrad(iw)*xare(iw)
          sum = sum + dw
          ic = icwl(iw)
          j = iplx(ic)
          i = iply(ic)
          if( j.eq.0 .or. i.eq.0 ) cycle
          if( nw.eq.2 .or. nw.eq.4 ) then
            dwrad(i,nw) = dw
            write(n6w,'(2x,"check  ",i3,i6,i7,2i5,1p2e12.3)')
     >        nw, iw, ic, j, i, wlrad(iw), dwrad(i,nw)
          endif
        enddo
      enddo
!
      close(n6w)
!
      return
!
!----------------------------------------------------------------------
!::error
!----------------------------------------------------------------------
 910  continue
      write(n6,'(2x,"ii.gt.ndbr ",2i7)') ii,ndbr
      call wexit("pst_wrad","dimension error in recomb.")
      end