!***********************************************************************
      subroutine imgchem(cspt)
!***********************************************************************
      use cimcom, only : arbz, csbz, flbz, gfbz, icbz, iebz, ikbz, iwbz
     >    , isbz, nbz, ndbz, prbz, sflux, snbz, tfbz, ychem, yfomt
     >    , mdl_roth, rothfit_flx, rothfit_eps
      use cimctl, only : cdirz
      use cntcom, only : cswl, icwl, ikwl, ipwl, npew, npsw, snwl, xpnt
     >    , ypnt, ncmax
      use cntwfl, only : xare, xcfl, xdeg, xfhta, xfhti, xfhtm, xopa
     >    , xopi, xopm, xsno
      use csize,  only : ndwp
      use csonic, only : itim, time
      use cunit,  only : lmspe, lmype, n6, wh_dir
      use mod_externalgrid
      implicit none
!
!::argument
      character, intent(in) :: cspt*(*)
!
!::local variables
      integer  iw, iwmx, nw, iws, iwe
      integer  ipmg, ipmg1
      integer  ii, i, ic, kop
      real*8  zfn, totfi, totfa, totfm
      real*8  zmax, fmin, zsum, xc, yc, x1, y1, x2, y2
      real*8  x0, y0, th, zfi, zfa, zfm
      real*8  zflxt(ndwp), zflxi(ndwp), zflxa(ndwp), zflxm(ndwp)
      real*8  zfspt(ndwp)
      real*8  yspt, yspt_iw
      integer  i6, i_roth
! function
      integer    imox, imoy
!
      write(n6,'(2x,"*** imgchem ***  ",a)') cspt
!
      yspt = ychem
!
      write(n6,'(2x,"ychem  =",f8.4,"  yfomt =",1pe12.3)') ychem, yfomt
      write(n6,'(2x,"  cspt =","[",a,"]","  yspt =",f8.4)')
     >   trim(cspt), yspt
!
      do i = 1, xsno
        kop = 0
        if(  index(cspt,xcfl(i)).gt.0 ) kop = 1
        xopi(i) = kop
        xopa(i) = kop
        xopm(i) = kop
        if( xcfl(i).ne."idp" .and. xcfl(i).ne."odp" ) then
          xopi(i) = 0
        endif
      enddo
!
!::flag (input data)
      write(n6,'(2x,"xsno = " i2)') xsno
      write(n6,'(2x,"xcfl (loc)= ",10(a3,2x))')    (xcfl(i),i=1,xsno)
      write(n6,'(2x,"xopi (ion)= ",10(2x,i1,2x))') (xopi(i),i=1,xsno)
      write(n6,'(2x,"xopa (atm)= ",10(2x,i1,2x))') (xopa(i),i=1,xsno)
      write(n6,'(2x,"xopm (mol)= ",10(2x,i1,2x))') (xopm(i),i=1,xsno)
!
!::clear
      call setd( zflxt, ndwp, 0.0d0 )
      call setd( zflxa, ndwp, 0.0d0 )
      call setd( zflxm, ndwp, 0.0d0 )
      call setd( zfspt, ndwp, 0.0d0 )
!
!::debug write
      i6 = 21
      if( lmype.eq.lmspe ) then  ! KSOPEN null
        open(unit=i6, file=trim(cdirZ)//"/"//cspt//".txt")
      else
        open(unit=i6, file="/dev/null")
      endif
      write(i6,'(2x,"chemical sputtering_",a,"  itim =",i7,"  time =",
     >  1pe14.6)') cspt, itim, time
      write(i6,'(2x,"nw",3x,"iw",3x,"Xc",9x,"Yc",8x,"th",8x,"are",
     >  9x,"tot",9x,"itot",9x,"atot",9x,"mtot",9x,"Fchm")')
!
!:: sputtering output
      i_roth = 121
      if(mdl_roth .and. lmype.eq.lmspe) then
        open(unit=i_roth, file="../"//trim(wh_dir)//"/zyspt.txt"
     >    ,status='replace')
        write(i_roth,'("#########################################")')
        write(i_roth,'("# chemical sputtering rate for each wall")')
        write(i_roth,'("# (J. Roth et al NF (2004) 44 L21-25)")')
        write(i_roth,'("# Y = ychem/(1+(f/f_0)^e)")')
        write(i_roth,'("# f_0:",1pe12.4," e:",1pe12.4)')
     >   rothfit_flx,rothfit_eps
        write(i_roth,'("#########################################")')
        write(i_roth,'(3x,"iw",4x,"yspt",10x,"flux")')
      else
        open(unit=i_roth, file="/dev/null")
      endif
!
!::flux
      totfi = 0.0d0
      totfa = 0.0d0
      totfm = 0.0d0
!
      do nw = 1, 4
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
          if(ic.le.ncmax .or. .not.use_exdata) then
            x1  = xpnt(ipmg)
            y1  = ypnt(ipmg)
            x2  = xpnt(ipmg1)
            y2  = ypnt(ipmg1)
          elseif(ic .le. ncmax+vac_ele_size) then
            x1  = vac_grid_x(ipmg)
            y1  = vac_grid_y(ipmg)
            x2  = vac_grid_x(ipmg1)
            y2  = vac_grid_y(ipmg1)
          elseif(ic .le. ncmax+vac_ele_size+pri_grid_size) then
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
          zfi = 0.0d0
          zfa = 0.0d0
          zfm = 0.0d0
          do i = 1, xsno
            if( xopi(i).eq.1 ) zfi = zfi + xfhti(iw,i)
            if( xopa(i).eq.1 ) zfa = zfa + xfhta(iw,i)
            if( xopm(i).eq.1 ) zfm = zfm + xfhtm(iw,i)
          enddo
!
          zflxt(iw) = zfi + zfa + zfm
          zflxi(iw) = zfi
          zflxa(iw) = zfa
          zflxm(iw) = zfm
!
          if(mdl_roth) then
            yspt_iw = yspt/(1.0d0+(zflxt(iw)/rothfit_flx)**rothfit_eps) ! use J.Roth model.
            zfspt(iw) = zflxt(iw)*yspt_iw
            ! file output
            write(i_roth,'(i5,2x,1pe12.4,2x,1pe12.4)') 
     >        iw, yspt_iw, zfspt(iw)
          else
            zfspt(iw) = zflxt(iw)*yspt ! default
          endif
!
          totfi     = totfi + zfi*xare(iw)
          totfa     = totfa + zfa*xare(iw)
          totfm     = totfm + zfm*xare(iw)
!
          write(i6,'(2x,i2,i5,2f11.6,f9.3,1p6e12.4)')
     >      nw,iw, x0, y0, th, xare(iw), zflxt(iw),
     >      zflxi(iw), zflxa(iw), zflxm(iw), zfspt(iw)
        enddo
      enddo
!
      close(i_roth)
!
      write(i6,'(2x)')
      write(i6,'(2x,"  d-flux     tot/ion/atm/mol = ",1p4e12.4)')
     >  totfi+totfa+totfm, totfi, totfa, totfm
      write(i6,'(2x,"  cd4-flux   tot/ion/atm/mol = ",1p4e12.4)')
     > yspt*(totfi+totfa+totfm), yspt*totfi, yspt*totfa, yspt*totfm
      close(i6)
!
!::maximum
      zmax = -1.0d20
      do iw = 1, iwmx
        zmax = dmax1( zmax, zfspt(iw)*xare(iw) )
      enddo
      write(n6,'(2x,"max-flux =",1pe12.3)') zmax
      if( zmax.lt.0.0d0 ) then
        call wexit("imgchem","flx-max < 0.0")
      endif
!
!::birth profile  yomt : flux omit (minium flux)
      fmin = zmax*yfomt
      zsum = 0.0d0
      ii = 0
      do iw = 1, iwmx
        zfn = zfspt(iw)*xare(iw)
        if( zfn.le.fmin ) cycle
        zsum = zsum + zfn
        ii = ii + 1
        if( ii.gt.ndbz ) call wexit("imgchem","ii.gt.ndbz")
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
      if(ii .eq. 0) then
        call wexit("imgchem","no birth particle")
      endif
      nbz = ii
      tfbz = zsum
      gfbz = zmax
      do i = 1, nbz
        prbz(i) = prbz(i)/zsum
      enddo
      prbz(nbz) = 1.0001d0
!
!::output
      write(n6,'(4x,"i",4x,"iw",4x,"xc",7x,"yc",7x,"th",7x,"prb",9x,
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
        xc = 0.5d0*(x1+x2)
        yc = 0.5d0*(y1+y2)
        zsum = zsum + flbz(i)
        if( yc.gt.0.0d0 ) cycle
        if( i.le.5 .or. i.ge.nbz-2 ) then
          write(n6,'(2x,2i5,2f9.4,f9.3,1p3e12.3,2i5)') i,iw,xc,yc
     >     , xdeg(iw), prbz(i), flbz(i)/arbz(i), arbz(i)
     >     , imox(ic), imoy(ic)
        endif
      enddo
!
      sflux = zsum
      write(n6,'(2x,10x,18x,5x,"total",2x,1p2e12.3)') tfbz, sflux
!
      return
      end
