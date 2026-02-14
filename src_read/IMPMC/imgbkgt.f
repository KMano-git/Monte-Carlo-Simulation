!***********************************************************************
      subroutine imgbkgt
!***********************************************************************
      use cimcom, only : arbz, csbz, flbz, gfbz, icbz, iebz, ikbz, isbz
     >    , iwbz, nbz, ndbz, prbz, snbz, prbz, sflux, tfbz
      use cimpuf, only : bk_are, bk_deg, bk_flx, bk_iiw, bk_imx, bk_ity
     >    , bk_mag, bk_nty
      use cntcom, only : cswl, icwl, ikwl, ipwl, snwl, xpnt, ypnt
      use cunit,  only : n6
      use mod_externalgrid, only : use_exdata, ipwl2, get_position_byic
      implicit none
!
!::local variables
      integer :: ii, i, m, iw, ic
      real(8) :: zmax, zsmn, zfn, zsum, x0, y0, x1, x2, y1, y2
      real(8) :: dgbz(ndbz)   ! degree
! function
      integer :: imox, imoy
!
      write(n6,'(/2x,"*** imgbkgt ***  bk_nty =",i2,"  bk_mag =",
     >  1p2e12.4)') bk_nty, (bk_mag(i),i=1,bk_nty)
!
      if( bk_nty.le.0 ) then
        nbz = 0
        tfbz = 0.0d0
        gfbz = 0.0d0
        return
      endif
!
      ii = 0
      zmax = -1.0d20
      zsmn = 0.0d0
      do i = 1, bk_imx
        m    = bk_ity(i)
        iw   = bk_iiw(i)
        zfn  = bk_flx(i)*bk_mag(m)
        zsmn = zsmn + zfn
        if( zfn.le.0.0d0 ) cycle
!
        zmax = dmax1( zmax, zfn )
        ii = ii + 1
        if( ii.gt.ndbz ) call wexit("imgbkgt","ii.gt.ndbz")
        prbz(ii) = zsmn
        flbz(ii) = zfn
        arbz(ii) = bk_are(i)
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
        dgbz(ii) = bk_deg(i)
      enddo
!
!::return
      nbz = ii
      tfbz = zsmn
      gfbz = zmax
      if(nbz .eq. 0) then
        sflux = 0.0d0
        prbz  = 0.0d0
        return
      else
!
!::normalization
        do i = 1, nbz
          prbz(i) = prbz(i)/zsmn
        enddo
        prbz(nbz) = 1.0001d0
!
!::output
        write(n6,'(4x,"i",4x,"iw",4x,"x0",7x,"y0",7x,"th",7x,"prb",9x,
     >    "spt",9x,"are",9x,"ix",3x,"iy")')
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
          write(n6,'(2x,2i5,2f9.4,f9.3,1p3e12.3,2i5)') i,iw,x0,y0
     >     , dgbz(i),prbz(i), flbz(i)/arbz(i)
     >     , arbz(i), imox(ic), imoy(ic)
        enddo
!
        sflux = zsum
        write(n6,'(2x,10x,18x,2x,"total",2x,1p2e12.3)') tfbz, sflux
        return
      endif
!
      return
      end
