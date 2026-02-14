!***********************************************************************
      subroutine imgzspt( cspt, kis )
!***********************************************************************
      use cimcom, only : arbz, azmas, csbz, flbz, gfbz, icbz, iebz, ikbz
     >    , isbz, iwbz, nbz, ndbz, prbz, sflux, snbz, tfbz, yfomt
     >    , yfphy
      use cimctl, only : cdirz
      use cntcom, only : cswl, icwl, ikwl, ipwl, isxw, isyw, npew, npsw
     >    , snwl, xpnt, ypnt, ncmax
      use cntwfl, only : xare, xdeg
      use cphcns, only : cev
      use csize,  only : ndwp
      use csonic, only : itim, time
      use cunit,  only : lmspe, lmype, n6
      use czwflx, only : zw_catmz, zw_eflx, zw_ismax, zw_pflx
      use mod_externalgrid
      implicit none

! arguments
      character, intent(in) :: cspt*(*)
      integer,   intent(in) :: kis      ! impurity number. This intends a source impurity(Ar), NOT a targert impurity(Carbon).

! local variables
      integer  kspphy
      real*8   flxz(ndwp), flxspt(ndwp)
      real*8   totz, totspt, totez
      integer  nw, iw, iws, iwe, i6, iwmx, j, i
      integer  ipmg, ipmg1
      integer  ii, ic, iz
      real*8   x1, y1, x2, y2, x0, y0, th, zfz, ztz,zne, zni, zte, zti
      real*8   zcs, zez, yd, zflx(ndwp), zmax, fmin, zsum
      real*8   zfn
! function
      real(8)    sputphy
      integer    imox, imoy, nsptin

      write(n6,'(2x,"*** imgzspt kis = ", i3, " ***")') kis

      if(zw_ismax(kis) .le. 0) then
        write(n6,*) zw_catmz(kis), "skip imgzspt due to zw_ismax = 0"
        nbz = 0
        tfbz = 0.0d0
        gfbz = 0.0d0
        return ! Yamoto
      endif

!::flag
      kspphy = nsptin(zw_catmz(kis))
      write(n6,'(2x,"catmz =",a3,"  kspphy =",i2,"  yfphy =",f7.3)')
     >   zw_catmz(kis), kspphy, yfphy
      if( kspphy.eq.0 ) then
      call wexit("imgphys","kspphy = 0")
      endif

!::clear
      call setd( flxz, ndwp, 0.0d0 )  ! ion
      call setd( flxspt, ndwp, 0.0d0 )  ! sputtering used in impmc

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
     >  8x,"j",4x,"i",4x,"Ne",9x,"Te",9x,"Ti",9x,"aveEz",6x,
     >  "Yz_eff",5x,"Flxz",7x,"Fspt")')

!::flux
      totz = 0.0d0
      totspt = 0.0d0

      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
        zne = 0.0d0
        zni = 0.0d0
        zte = 0.0d0
        zti = 0.0d0
        zcs = 0.0d0
        i = 0
        j = 0
        do iw = iws, iwe
          ipmg  = ipwl(iw)
          if(use_exdata) then
            ipmg1 = ipwl2(iw)
          else
            if(iw == iwe) cycle
            ipmg1 = ipwl(iw+1)
          endif
          ic = icwl(iw)
          if(ic.le.ncmax .or. .not.use_exdata) then
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

          ii = 0
          totez = 0.0d0

          do iz = 0, zw_ismax(kis)
            zfz = zw_pflx(iz,iw,kis)
            if(zfz.le.0.0d0) cycle
            ztz = zw_eflx(iz,iw,kis) / zfz / cev

!::impact energy
            if(nw==2 .or. nw==4) then
              j = isxw(iw)
              i = isyw(iw)
              if( j.le.0 .or. i.le.0 ) cycle
              call gdpls(j,i,zne,zni,zte,zti,zcs)
              zez = ztz + 3.0d0*zte*iz
            else
              zez = ztz
            endif
            ii = ii + 1
            totez = totez + ztz

!::sputtering of normal incidence
            yd = sputphy( zez, 0.0d0, kspphy )

!::angular effect of incidence
            yd = yd*yfphy

!::sputterd flux
            flxz(iw) = flxz(iw) + zfz

!::emitted neutral impurity
            flxspt(iw) = flxspt(iw) + zfz*yd

!::total flux
            totz  = totz + flxz(iw)
            totspt  = totspt + flxspt(iw)
          enddo ! iz

          if(ii.gt.0) totez = totez / ii
          yd = 0.0d0
          if(flxz(iw).gt.0.0d0) yd = flxspt(iw)/flxz(iw)

!::debug write
          write(i6,'(2x,i2,i5,2f11.6,f9.3,1pe11.3,
     >     2i5,1p6e11.3,1p4e11.3)')
     >    nw, iw, x0, y0, th, xare(iw),
     >    j, i, zne, zte, zti, totez, yd,
     >    flxz(iw), flxspt(iw)
        enddo ! iw
      enddo ! nw

      write(i6,'(2x)')
      write(i6,'(2x,"Total flux  Fion =",1pe11.3,"  Fspt =",1pe11.3)')
     >      totz, totspt

!::no sputtering
      if( totspt.le.0.0d0 ) then
      write(i6,'(2x,"No physical sputtering   Fspt =",1pe11.3)') totspt
      write(n6,'(2x,"No physical sputtering   Fspt =",1pe11.3)') totspt
      nbz = 0
      tfbz = 0.0d0
      gfbz = 0.0d0
      return
      endif

      close(i6)

!::physical sputtering
      call setd( zflx, ndwp, 0.0d0 )
      do iw = 1, iwmx
      zflx(iw) = flxspt(iw)
      enddo

!::maximum
      zmax = -1.0d20
      do iw = 1, iwmx
      zmax = dmax1( zmax, zflx(iw)*xare(iw) )
      enddo
      write(n6,'(2x,"max-flux =",1pe12.3)') zmax
      if( zmax.le.0.0d0 ) then
        call wexit("imgzspt","flx-max < 0.0")
      endif

!::birth profile
      fmin = zmax*yfomt
      zsum = 0.0d0
      ii = 0
      do iw = 1, iwmx
      zfn = zflx(iw)*xare(iw)
      if( zfn.le.fmin ) cycle
      zsum = zsum + zfn
      ii = ii + 1
      if( ii.gt.ndbz ) goto 910
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

!::normalization
      nbz = ii
      tfbz = zsum
      gfbz = zmax
      do i = 1, nbz
      prbz(i) = prbz(i)/zsum
      enddo
      prbz(nbz) = 1.0001d0

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
      write(n6,'(2x,2i5,2f9.4,f9.3,1p3e12.3,2i5)') i,iw,x0,y0,xdeg(iw),
     >  prbz(i), flbz(i)/arbz(i), arbz(i), imox(ic), imoy(ic)
      endif
      enddo

      sflux = zsum
      write(n6,'(2x,10x,18x,5x,"total",2x,1p2e12.3)') tfbz, sflux

!xxx  azmas = svmas
      write(n6,'(2x,"azmas =",f10.3)') azmas
      return

 910  continue
      call wexit("imgphys","ii.gt.ndbz")
      end
