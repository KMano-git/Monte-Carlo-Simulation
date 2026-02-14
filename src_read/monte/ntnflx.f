!**********************************************************************
      subroutine ntnflx(kflx,tfion,tfntl)
!**********************************************************************
!
!      neutral birth profile
!
!  argument
! ------------------------------------------------------------
!   kflx       : 1:(sol)/2:(idp)/3:(prv)/4:(odp)/5:(man)
!              : 6:(puf)/7:(vol)
!   tfion      : total ion flux
!   tfntl      : total neutral flux
!
! common variables
! ------------------------------------------------------------
!   icbr(ndbr) : cell number in monte carlo
!
!   ikbr(ndbr) : type of cell boundary  1/2/3/4 when s-source
!   isbr(ndbr) : start point number in monte carlo
!   iebr(ndbr) : end point number in monte carlo
!   iwbr(ndbr) : flag of wall(1) / plasma surface(0)
!   iwpbr(ndbr): index of plasma surface isxp(iwp), isyp(iwp)
!   csbr(ndbr) : cosine of cell boundary
!   snbr(ndbr) : sine of cell boundary
!   eibr(ndbr) : injection energy of H+ ions
!   e0br(ndbr) : injection energy of H2 neutrals (wall temperature)
!
!   nbr        : number of cell
!   prbr(ndbr) : cumulative probabirity
!   flbr(ndbr) : local neutral flux
!   fmbr(ndbr) : local neutral flux in monte cal.
!   tfbr       : total neutral flux
!   cfbr       : name of flux  sol/idp/prv/odp/man/puf/vol  see ntmset
!
!----------------------------------------------------------------------
      use cntcom, only : cfbr, csbr, csps, e0br, e0pfa, e0pfm, e0ps
     >    , eibr, eips, flbr, flxin, icbr, icps, iebr, ikbr, ikps
     >    , ipfiw, ipfmp, ipfmx, iplx, iply, ipps, isbr, iwbr, iwpbr
     >    , lemd2, lempf, lemwl, nbr, ncmax, npep, npf, npsp, pfflx
     >    , pfmag, prbr, snbr, snps, tfbr
     >    , use_gaspufftime
     >    , iwpbr_wall, imps
      use cntmnt, only : mcflx
      use cplcom, only : flps, flvl
      use csize,  only : ndbr
      use csonic, only : itim
      use cunit,  only : cdgr, mygrp, n6
      use cntpfctl, only : fpfctl
      implicit none
!
!::argument
      integer, intent(in)  :: kflx
      real(8), intent(out) :: tfion, tfntl
!
!::local variables
      integer ia, ii, nw, ifgwl, iws, iwe, iw, i, m, ic, j
     > , ii_loop 
      integer ldbg/1/; save ldbg
      real*8  zsmi, zsmn, zfi, zfn, asn, zmax
!
!::clear /cntbir/
      cfbr = mcflx(kflx)
      nbr  = 0
      tfbr = 0.0d0
      call setd( prbr, ndbr, 0.0d0 )
      call setd( flbr, ndbr, 0.0d0 )
!
!::ion and neutral flux
      tfion = 0.0d0
      tfntl = 0.0d0
!
!::plasma ion species
      ia = 1
!
!::model
      ifgwl = -1
!
!----------------------------------------------------------------------
!::surface source (sol/idp/prv/odp)
!----------------------------------------------------------------------
!::plasma surface
      if( kflx.ge.1 .and. kflx.le.4 ) then
        nw  = kflx
        if( nw.eq.1 .or. nw.eq.3 ) then
          if( lemwl.eq.0 ) ifgwl = 1  !  D0, 3eV, inward
          if( lemwl.eq.1 ) ifgwl = 2  !  D0, Ti,  ouward
        elseif( nw.eq.2 .or. nw.eq.4 ) then
          ifgwl = 3
        endif
!
        ii = 0
        zsmi = 0.0d0
        zsmn = 0.0d0
        iws = npsp(nw)
        iwe = npep(nw)
        do iw = iws, iwe-1
          zfi = flps(iw,ia)
          zfn = dmax1(zfi,0.0d0)
!
          zsmi = zsmi + zfi
          zsmn = zsmn + zfn
          if( zfn.le.0.0d0 ) cycle
!
          if( mygrp.ne.cdgr(2) ) cycle  ! KSFUJI
          ii = ii + 1
          flbr(ii) = zfn
          prbr(ii) = zsmn
          icbr(ii) = icps(iw)
          ikbr(ii) = ikps(iw)
          isbr(ii) = ipps(iw)
          iebr(ii) = ipps(iw+1)
          iwbr(ii) = ifgwl
          iwpbr(ii)= iw
          csbr(ii) = csps(iw)
          snbr(ii) = snps(iw)
          eibr(ii) = eips(iw)
          if(lemd2.eq.1)eibr(ii)=0.0d0
          e0br(ii) = e0ps(iw)
        enddo
        ! use imps to get wall boundary index on diverter
        do ii_loop = 1, ii
          iwpbr_wall(ii_loop) = imps(iwpbr(ii_loop))
        enddo
        !
!
!----------------------------------------------------------------------
!::surface (main)   Note  clock wise
!----------------------------------------------------------------------
      elseif( kflx.eq.5 ) then
        ifgwl = 0
        zsmi = 0.0d0
        zsmn = 0.0d0
        nw = 5
        do iw = npsp(nw), npep(nw)-1
          zfi  = flps(iw,ia)
          zsmi = zsmi + zfi
        enddo
!
!----------------------------------------------------------------------
!::puff source
!----------------------------------------------------------------------
      elseif( kflx.eq.6 ) then
        !::time dependance gaspuff
        if(use_gaspufftime)then
          call gaspuff_update
        endif
        if( lempf.eq.0 ) ifgwl = 4
        if( lempf.eq.1 ) ifgwl = 5
        ii = 0
        zsmi = 0.0d0
        zsmn = 0.0d0
        if( npf.gt.0 ) then
          do i = 1, ipfmx
            m  = ipfmp(i)
            iw = ipfiw(i)
            zfn  = pfflx(i)*pfmag(m)*fpfctl(m)
            zsmn = zsmn + zfn   ! tfntl
            if( zfn.le.0.0d0 ) cycle
!
            if( mygrp.ne.cdgr(2) ) cycle  ! KSFUJI
            ii = ii + 1
            flbr(ii) = zfn
            prbr(ii) = zsmn
            icbr(ii) = icps(iw)
            ikbr(ii) = ikps(iw)
            isbr(ii) = ipps(iw)
            iebr(ii) = ipps(iw+1)
            iwbr(ii) = ifgwl
            iwpbr(ii)= iw
            csbr(ii) = csps(iw)
            snbr(ii) = snps(iw)
            eibr(ii) = eips(iw)
            if( ifgwl.eq.4 ) e0br(ii) = e0pfa
            if( ifgwl.eq.5 ) e0br(ii) = e0pfm
          enddo
        endif
!
!----------------------------------------------------------------------
!::volume source (recombination H0 ==> H+ )
!----------------------------------------------------------------------
      elseif( kflx.eq.7 ) then
        ifgwl = 0
        zsmi = 0.0d0
        zsmn = 0.0d0
        ii   = 0
        do ic = 1, ncmax
          j = iplx(ic)
          i = iply(ic)
          if( j.le.0 .or. i.le.0 ) cycle
          asn = flvl(j,i,ia)
          if( asn.le.0.0d0 ) cycle
          zsmn = zsmn + asn
          ii = ii + 1
          if( ii.gt.ndbr ) goto 910
          flbr(ii) = asn
          prbr(ii) = zsmn
          icbr(ii) = ic
        enddo
!
!----------------------------------------------------------------------
!::volume source (recombination H+ ==> H0 )
!----------------------------------------------------------------------
      elseif( kflx.eq.8 ) then
        ifgwl = 0
        zsmi = 0.0d0
        zsmn = 0.0d0
        ii   = 0
        do ic = 1, ncmax
          j = iplx(ic)
          i = iply(ic)
          if( j.le.0 .or. i.le.0 ) cycle
          asn = flvl(j,i,ia)
          if( asn.le.0.0d0 ) cycle
          zsmi = zsmi + asn
        enddo
!
!----------------------------------------------------------------------
!::error
!----------------------------------------------------------------------
      else
        write(n6,'(/2x,"*** ntnflx ***  kflx =",i3)') kflx
        call wexit("ntnflx","1<= kflx <= 7")
      endif
!
!----------------------------------------------------------------------
!::normalization & max value
!----------------------------------------------------------------------
      if( zsmn.gt.0.0d0 ) then
        nbr = ii
        tfbr = zsmn
        zmax = -1.0d20
        do i = 1, nbr
          zmax = dmax1( zmax, flbr(i) )
          prbr(i) = prbr(i)/zsmn
        enddo
      endif
!
!----------------------------------------------------------------------
!::total flux
!----------------------------------------------------------------------
      tfion = zsmi
      tfntl = zsmn
!
!::common variables
      flxin = tfbr
!
!----------------------------------------------------------------------
!::debug write
!----------------------------------------------------------------------
      if( mygrp.eq.cdgr(1) ) return
      if( itim.gt.0 ) return
      if( ldbg.eq.0 ) return
      call debg_ntnflx(kflx,ifgwl)
      return
!
!----------------------------------------------------------------------
!::error
!----------------------------------------------------------------------
 910  continue
      write(n6,'(2x,"ii.gt.ndbr ",2i7)') ii,ndbr
      call wexit("ntnflx","dimension error in recomb.")
!
!**********************************************************************
      contains
      subroutine gaspuff_update
!**********************************************************************
      use cntcom, only : pfmag, ndpf
      implicit none
!
!::local variables
      integer m
      real*8 interpolate_value

      do m = 1, ndpf
        call timevalue_inter(interpolate_value,m)
        pfmag(m) = interpolate_value
      enddo
      end subroutine gaspuff_update

!**********************************************************************
      subroutine timevalue_inter(output_val,target_m)
!**********************************************************************
      use csonic, only:time
      use cplcom, only:qtim_ini
      use cntcom, only:pufftime_size, pufftime, pfmag_time
      implicit none
!
!::arguments
      integer, intent(in) :: target_m
      real*8, intent(out) :: output_val
!
!::local variables
      integer i_time
      real*8 t_next,t_preb,val_next,val_preb,time_accu

      if(time.ne.0.0d0) then
        time_accu = time
      else
        time_accu = qtim_ini ! first step for continue calculation
      endif

      ! time less than time-profile
      if(time_accu .le. pufftime(1)) then
        output_val = pfmag_time(target_m,1)
        return
      endif
      ! time over than time-profile
      if(time_accu .ge. pufftime(pufftime_size)) then
        output_val = pfmag_time(target_m,pufftime_size)
        return
      endif

      do i_time = 2, pufftime_size
        t_next = pufftime(i_time)
        if(time_accu .le. t_next) then
          t_preb = pufftime(i_time-1)
          val_next = pfmag_time(target_m,i_time)
          val_preb = pfmag_time(target_m,i_time-1)
          ! liner interpolate of value
          output_val = val_preb +
     >      (val_next-val_preb)/(t_next-t_preb) * (time_accu-t_preb)
          return
        endif
      enddo
      end subroutine timevalue_inter
!::
      end subroutine ntnflx
!
!**********************************************************************
      subroutine debg_ntnflx(kflx,ifgwl)
!**********************************************************************
      use cntcom, only : cfbr, csbr, e0br, eibr, flbr, icbr, iebr, ikbr
     >    , iplx, iply, isbr, iwbr, lemdp, lempf, lemwl, migx, migy, nbr
     >    , ncmax, prbr, snbr, temin_rec, tfbr, volm
      use cntpls, only : dene, deni, teme
      use cplcom, only : flvl, gfvl
      use cunit,  only : n6
      implicit none
!
!::argument
      integer, intent(in) :: kflx, ifgwl
!
!::local variables
      real*8  snup(10)
      integer ncup(10)
      integer ia, ic, j, i,imd
      real*8  asn, asv, zrt
      real*8  ate, ane
      character  cmsg(-1:5)*21
!                   12345678901234567890
      data cmsg(-1)/"-1: error invalid"/
      data cmsg(0)/" 0: no effective"/
      data cmsg(1)/" 1: wl D0 3eV inward"/
      data cmsg(2)/" 2: wl D0 Ti  outward"/
      data cmsg(3)/" 3: dp"/
      data cmsg(4)/" 4: pf D0 3eV  inward"/
      data cmsg(5)/" 5: pf D2 e0pf inward"/
! function
      integer  lenx
      real(8)  svhrcm
!
!----------------------------------------------------------------------
!::plasma surface or wall
!----------------------------------------------------------------------
      if( kflx.ne.7 ) then
      write(n6,'(/2x,"*** ntnflx(",a,") ***  kflx =",i2,"  nbr =",i5,
     >  "  tfbr =",1pe14.6)') cfbr(1:lenx(cfbr)), kflx, nbr, tfbr
!
      if( nbr.gt.0 ) then
      write(n6,'(2x,"ifgwl =",2i3,2x,a,2x,"lemwl,lemdp,lempf =",3i3)')
     >   ifgwl, iwbr(1), cmsg(ifgwl), lemwl, lemdp, lempf
      write(n6,'(3x,"i",2x,"pr",11x,"ic",5x,"ix",3x,"iy",1x,"ik",2x,
     >  "is",4x,"ie",4x,"wbr",2x,"cs",6x,"sn",6x,"e0",6x,"ei")')
      do i = 1, nbr
      if( i.ge.10 .and. i.le.nbr-10 ) cycle
      write(n6,'(i4,f13.9,i6,2i5,i3,2i6,i5,3f8.4,f8.2,i3)')
     > i, prbr(i), icbr(i), migx(icbr(i)), migy(icbr(i)),
     > ikbr(i),isbr(i),iebr(i),iwbr(i),csbr(i),snbr(i),e0br(i),eibr(i)
      enddo
      endif
      return
      endif
!
!----------------------------------------------------------------------
!::volume source
!----------------------------------------------------------------------
      call setd( snup, 10, 0.0d0 )
      call seti( ncup, 10, 0 )
!
      ia = 1
!
      do ic = 1, ncmax
        j = iplx(ic)
        i = iply(ic)
        if( j.le.0 .or. i.le.0 ) cycle
        asn = flvl(j,i,ia)
        zrt = asn/gfvl(ia)
        snup(5) = snup(5) + asn
        ncup(5) = ncup(5) + 1
        if( zrt.gt.1.0d-4 ) then
          snup(4) = snup(4) + asn
          ncup(4) = ncup(4) + 1
        endif
        if( zrt.gt.1.0d-3 ) then
          snup(3) = snup(3) + asn
          ncup(3) = ncup(3) + 1
        endif
        if( zrt.gt.1.0d-2 ) then
          snup(2) = snup(2) + asn
          ncup(2) = ncup(2) + 1
        endif
        if( zrt.gt.1.0d-1 ) then
          snup(1) = snup(1) + asn
          ncup(1) = ncup(1) + 1
        endif
      enddo
!
      write(n6,'(/2x,"*** ntnflx(",a,") ***  kflx =",i2,"  nbr =",i5,
     >  "  tfbr =",1pe14.6)') cfbr(1:lenx(cfbr)), kflx, nbr, tfbr
!
      do i = 1, 5
      write(n6,'(2x,"ratio >",1pe9.2,"  num. cell =",i5,"  tot-sn =",
     >  1pe12.3)') 1.0d0/10.0d0**i, ncup(i), snup(i)
      enddo
!
      write(n6,'(2x)')
      write(n6,'(3x,"i",3x,"ic",3x,"ix",2x,"iy",3x,"pr",10x,"sn",
     >  10x,"Ne",10x,"Te",10x,"sigv",8x,"Sn")')
      imd = max0( nbr/20, 1 )
      do i = 1, nbr, imd
        ic = icbr(i)
        if( teme(ic).le.0.0 .or. dene(ic).le.0.0 ) cycle
        ate = teme(ic)
        ane = dene(ic)
        ate = dmax1(ate,temin_rec)  ! svhrcm
        asv = svhrcm(ate,ane)
        asn = dene(ic)*deni(ic,ia)*asv*volm(ic)
        write(n6,'(i4,i6,2i4,1p6e12.3)')
     >    i, ic, migx(ic), migy(ic), prbr(i), flbr(i),
     >    ane, ate, asv, asn
      enddo
!
      return
      end