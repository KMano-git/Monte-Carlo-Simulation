!***********************************************************************
      subroutine imgphys(cspt)
!***********************************************************************
      use cimcom, only : aimas, azmas, arbz, csbz, flbz, gfbz, icbz
     >    , iebz, ikbz, isbz, iwbz, nbz, ndbz, prbz, sflux, snbz, tfbz
     >    , yfomt, yfphy
      use cimctl, only : cdirz
      use cntcom, only : cswl, icwl, ikwl, ipwl, isxw, isyw, npew, npsw
     >    , snwl, xpnt, ypnt, ncmax
      use cntwfl, only : xare, xdeg, xfhti
      use cphcns, only : cmp
      use csize,  only : ndwp
      use csonic, only : itim, time
      use cunit,  only : lmspe, lmype, n6
      use mod_externalgrid
      implicit none 
!
!::argument
      character, intent(in) :: cspt*(*)
!
!::local variables
      integer  kspphy
      real*8   flxion(ndwp), flxphy(ndwp), flxspt(ndwp)
      real*8   totion, totphy, totspt
      integer  nw, iw, iws, iwe, i6, iwmx, j, i
      integer  ipmg, ipmg1
      integer  ii, ic
      real*8   x1, y1, x2, y2, x0, y0, th, zfi, zne, zni, zte, zti
      real*8   zcs, zed, yd, zflx(ndwp), zmax, fmin, zsum
      real*8   zfn
! function
      real(8)    sputphy
      integer    imox, imoy, nsptin
!
      write(n6,'(2x,"*** imgphys ***")')
!
!::flag
      kspphy = 0
      if( dabs(aimas-1.0d0).le.0.01d0  ) kspphy = nsptin("H")
      if( dabs(aimas-2.0d0).le.0.01d0  ) kspphy = nsptin("D")
      if( dabs(aimas-3.0d0).le.0.01d0  ) kspphy = nsptin("T")
      write(n6,'(2x,"aimas =",f7.3,"  kspphy =",i2,"  yfphy =",f7.3)')
     >   aimas, kspphy, yfphy
      if( kspphy.eq.0 ) then
        call wexit("imgphys","kspphy = 0")
      endif
!
!::clear
      call setd( flxion, ndwp, 0.0d0 )  ! D+ ion
      call setd( flxphy, ndwp, 0.0d0 )  ! physical sputtering
      call setd( flxspt, ndwp, 0.0d0 )  ! sputtering used in impmc
!
!::debug write
      i6 = 21
      if( lmype.eq.lmspe ) then  ! KSOPEN null
        open(unit=i6, file=trim(cdirZ)//"/"//cspt//".txt")
      else
        open(unit=i6, file="/dev/null")
      endif
      write(i6,'(2x,"physical sputtering  itim =",i7,"  time =",
     >   1pe14.6)') itim, time
      write(i6,'(2x,"nw",3x,"iw",3x,"x0",9x,"y0",8x,"th",8x,"Are",
     >  8x,"j",4x,"i",4x,"Ne",9x,"Te",9x,"Ti",9x,"Ed",9x,
     >  "Yd",9x,"Fion",7x,"Fphy",7x,"Fspt")')
!
!::flux
      totion = 0.0d0
      totphy = 0.0d0
      totspt = 0.0d0
!
      do nw = 2, 4, 2 ! only for inner diverter(nw=2) or outer diverter(nw=4)
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe
          ipmg  = ipwl(iw)
          if(use_exdata) then
            ipmg1 = ipwl2(iw)
          else
            if(iw == iwe) cycle
            ipmg1 = ipwl(iw+1)
          endif
          ic = icwl(iw)
          if(ic .le. ncmax .or. .not.use_exdata) then
            x1 = xpnt(ipmg)
            y1 = ypnt(ipmg)
            x2 = xpnt(ipmg1)
            y2 = ypnt(ipmg1)
          elseif(ic .le. ncmax+vac_ele_size)then
            x1 = vac_grid_x(ipmg)
            y1 = vac_grid_y(ipmg)
            x2 = vac_grid_x(ipmg1)
            y2 = vac_grid_y(ipmg1)
          elseif(ic .le. ncmax+vac_ele_size+pri_ele_size)then
            x1 = pri_grid_x(ipmg)
            y1 = pri_grid_y(ipmg)
            x2 = pri_grid_x(ipmg1)
            y2 = pri_grid_y(ipmg1)
          else
            x1 = vgx_EX(ipmg)
            y1 = vgy_EX(ipmg)
            x2 = vgx_EX(ipmg1)
            y2 = vgy_EX(ipmg1)
          endif
          x0  = 0.5d0*(x1+x2)
          y0  = 0.5d0*(y1+y2)
          th  = xdeg(iw)
          iwmx = iw
!
!::Note   xfhti = Ni*Csd*h*cosw   xare = 2*pi*R*dl
          zfi = xfhti(iw,nw)
!
!::impact energy (Eirene code)
!       Ed = 5/2*Ti + 1/2*Md*Cs**2 + fai          fai = 3*Te
!       Ec = 5/2*Ti + 1/2*Mc*Cs**2 + <z>*fai      <z> = 2
!
          j = isxw(iw)
          i = isyw(iw)
          if( j.le.0 .or. i.le.0 ) cycle
          call gdpls(j,i,zne,zni,zte,zti,zcs)
!
          zed = 2.5d0*zti + 0.5d0*aimas*cmp*zcs**2 + 3.0d0*zte
!
!::sputtering of normal incidence
          yd = sputphy( zed, 0.0d0, kspphy )
!
!::angular effect of incidence
          yd = yd*yfphy
!
!::sputterd flux
          flxion(iw) = zfi
          flxphy(iw) = zfi*yd
!
!::emitted neutral impurity
          flxspt(iw) = flxphy(iw)
!
!::total flux
          totion  = totion + flxion(iw)*xare(iw)
          totphy  = totphy + flxphy(iw)*xare(iw)
          totspt  = totspt + flxspt(iw)*xare(iw)
!
!::debug write
          write(i6,'(2x,i2,i5,2f11.6,f9.3,1pe11.3,
     >       2i5,1p5e11.3,1p3e11.3)')
     >      nw, iw, x0, y0, th, xare(iw),
     >      j, i, zne, zte, zti, zed,  yd,
     >      flxion(iw), flxphy(iw), flxspt(iw)
        enddo !iw
      enddo !nw
!
      write(i6,'(2x)')
      write(i6,'(2x,"Total flux  Fion =",1pe11.3,"  Fphy =",1pe11.3,
     >  "  Fspt =",1pe11.3)') totion, totphy,
     >   totspt
!
!::no sputtering
      if( totspt.le.0.0d0 ) then
        write(i6,'(2x,"No physical sputtering   Fspt =",1pe11.3)')
     >    totspt
        write(n6,'(2x,"No physical sputtering   Fspt =",1pe11.3)')
     >    totspt
        nbz = 0
        tfbz = 0.0d0
        gfbz = 0.0d0
        return
      endif
!
      close(i6)
!
!::physical sputtering
      call setd( zflx, ndwp, 0.0d0 )
      do iw = 1, iwmx
        zflx(iw) = flxspt(iw)
      enddo
!
!::maximum
      zmax = -1.0d20
      do iw = 1, iwmx
        zmax = dmax1( zmax, zflx(iw)*xare(iw) )
      enddo
      write(n6,'(2x,"max-flux =",1pe12.3)') zmax
      if( zmax.le.0.0d0 ) then
        call wexit("imgphys","flx-max < 0.0")
      endif
!
!::birth profile
      fmin = zmax*yfomt
      zsum = 0.0d0
      ii = 0
      do iw = 1, iwmx
        zfn = zflx(iw)*xare(iw)
        if( zfn.le.fmin ) cycle
        zsum = zsum + zfn
        ii = ii + 1
        if( ii.gt.ndbz ) call wexit("imgphys","ii.gt.ndbz")
        prbz(ii) = zsum
        flbz(ii) = zfn
        arbz(ii) = xare(iw)
        icbz(ii) = icwl(iw)
        ikbz(ii) = ikwl(iw)
        isbz(ii) = ipwl(iw)
        if(use_exdata) then
          iebz(ii) = ipwl2(iw)
        else
          iebz(ii) = ipwl(iw+1)
        endif
        iwbz(ii) = iw
        csbz(ii) = cswl(iw)
        snbz(ii) = snwl(iw)
      enddo
!
!::normalization
      nbz = ii
      tfbz = zsum
      gfbz = zmax
      do i = 1, nbz
        prbz(i) = prbz(i)/zsum
      enddo
      prbz(nbz) = 1.0001d0
!
!::output
      write(n6,'(4x,"i",4x,"iw",4x,"x0",7x,"y0",7x,"th",7x,"prb",9x,
     >  "spt",9x,"are",9x,"ix",3x,"iy")')
      zsum = 0.0d0
      do i = 1, nbz
        iw = iwbz(i)
        ic = icbz(i)
        if(use_exdata) then
          call get_position_byic(ic,isbz(i),.true.,x1)
          call get_position_byic(ic,iebz(i),.true.,x2)
          call get_position_byic(ic,isbz(i),.false.,y1)
          call get_position_byic(ic,iebz(i),.false.,y2)
        else
          x1 = xpnt(isbz(i))
          x2 = xpnt(iebz(i))
          y1 = ypnt(isbz(i))
          y2 = ypnt(iebz(i))
        endif
        x0 = 0.5d0*(x1+x2)
        y0 = 0.5d0*(y1+y2)
        zsum = zsum + flbz(i)
        if( y0.gt.0.0d0 ) cycle
        if( i.le.5 .or. i.ge.nbz-2 ) then
          write(n6,'(2x,2i5,2f9.4,f9.3,1p3e12.3,2i5)') i,iw,x0,y0,
     >      xdeg(iw),prbz(i), flbz(i)/arbz(i), arbz(i), imox(ic),
     >      imoy(ic)
        endif
      enddo
!
      sflux = zsum
      write(n6,'(2x,10x,18x,5x,"total",2x,1p2e12.3)') tfbz, sflux
!
      write(n6,'(2x,"azmas =",f10.3)') azmas
      return
      end
