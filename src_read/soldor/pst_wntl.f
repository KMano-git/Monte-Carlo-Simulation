!***********************************************************************
      subroutine pst_wntl(kk)
!***********************************************************************
      use cntcom, only : npew, npsw
      use cntmnt, only : mfmax, pfl_ntl, vcsrc, vflux, visrc
      use cntwcn, only : swabs, swion, swnrm, swsum, wabs, wcsrc, wehta
     >    , wehtm, wfhta, wfhtm, wion, wisrc, wnrm, wsum
      use cntwfl, only : dwntl, iywl, jxwl, xare, xcfl
      use cphcns, only : cev
      use csize,  only : ndwp, ndxy
      use csonic, only : itim, time
      use cunit,  only : n6
      use mod_externalgrid, only : use_exdata
      implicit none
!
!::argument
      integer, intent(in) :: kk
!
!::local variables
      integer  nw, iw, iws, iwe, ixy
      integer  i, iflx, isrc, i6
      real*8   zfa, zfm, zea, zem, tfntl, znrm, zmx
      real*8   sfhta(ndwp), sfhtm(ndwp), sehta(ndwp), sehtm(ndwp)
!
      write(n6,'(/2x,"*** pst_wntl ***  kk =",i2)') kk
!
!::clear
      sfhta(1:ndwp) = 0.0d0
      sfhtm(1:ndwp) = 0.0d0
      sehta(1:ndwp) = 0.0d0
      sehtm(1:ndwp) = 0.0d0
      dwntl(1:ndxy,1:4) = 0.0d0
      if( kk.eq.0 ) return
!
!::sumup for source
      do iflx = 1, mfmax
        isrc  = visrc(iflx)
        if( isrc.le.0 ) cycle
        tfntl = pfl_ntl(iflx)
!
        xcfl(isrc) = vcsrc(isrc)
        write(n6,'(2x,"iflx/isrc =",i2,"/",i2,2x,a,2x,1pe14.6," =>",
     >     1pe14.6)')
     >    iflx, isrc, vcsrc(isrc), vflux(isrc), tfntl
!
        wisrc = isrc
        wcsrc = vcsrc(isrc)
        call setvmwork( 1, isrc )
        call ntwcon
!
        write(n6,'(2x,"wsum =",1p2e12.4,"  wnrm =",1p2e12.4,
     >    "  wion =",1p2e12.4,"  wabs =",1p2e12.4)')
     >    wsum, swsum, wnrm, swnrm, wion, swion, wabs, swabs
!
        znrm = tfntl/swnrm
        zfa = 0.0d0
        zfm = 0.0d0
        zea = 0.0d0
        zem = 0.0d0
!
        do nw = 1, 4
          iws = npsw(nw)
          iwe = npew(nw)
          do iw = iws, iwe
            zfa = zfa + wfhta(iw)
            zfm = zfm + wfhtm(iw)
            zea = zea + wehta(iw)
            zem = zem + wehtm(iw)
            sfhta(iw) = sfhta(iw) + wfhta(iw)*znrm
            sfhtm(iw) = sfhtm(iw) + wfhtm(iw)*znrm
            sehta(iw) = sehta(iw) + wehta(iw)*znrm*cev
            sehtm(iw) = sehtm(iw) + wehtm(iw)*znrm*cev
          enddo  ! loop(iw)
        enddo  ! loop(nw)
!
        zfa = zfa*znrm
        zfm = zfm*znrm
        zea = zea*znrm
        zem = zem*znrm
!
        if(zfa.ne.0.0d0 .and. zfm.ne.0.0d0)then
          write(n6,'(2x,"flux =",1p2e12.4,"  avE0 =",1p2e12.4)')
     >     zfa, zfm, zea/zfa, zem/zfm
        endif
      enddo
!
!::power density
      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe
          if(.not.use_exdata .and. iw==iwe) cycle
          sfhta(iw) = sfhta(iw)/xare(iw)
          sfhtm(iw) = sfhtm(iw)/xare(iw)
          sehta(iw) = sehta(iw)/xare(iw)
          sehtm(iw) = sehtm(iw)/xare(iw)
        enddo
      enddo
!
!::debug write
      i6 = 21
      open(unit=i6,file="Wntl.txt")
      write(i6,'(2x,"neutral flux to wall  itim =",i8,"  time =",
     >   1pe14.6)') itim, time
      do nw = 1, 4
        write(i6,'(2x,3x,"nw",3x,"iw",2x,"ixy",2x,"xare",8x,
     >    "Fa",10x,"Fm",10x,"Qa",10x,"Qm",10x,"Qt")')
        iws = npsw(nw)
        iwe = npew(nw)
        zfa = 0.0d0
        zfm = 0.0d0
        zea = 0.0d0
        zem = 0.0d0
        zmx = 0.0d0
        do iw = iws, iwe
          if(.not.use_exdata .and. iw==iwe) cycle
          zfa = zfa + sfhta(iw)*xare(iw)
          zfm = zfm + sfhtm(iw)*xare(iw)
          zea = zea + sehta(iw)*xare(iw)
          zem = zem + sehtm(iw)*xare(iw)
          zmx = dmax1( zmx, sehta(iw)+sehtm(iw) )
!
          if( nw.eq.1 .or. nw.eq.3 ) ixy = jxwl(iw)
          if( nw.eq.2 .or. nw.eq.4 ) ixy = iywl(iw)
          if( ixy.gt.0 ) then
            dwntl(ixy,nw) = (sehta(iw)+sehtm(iw))*xare(iw)
          endif
!
          write(i6,'(2x,3i5,1p6e12.4)')
     >      nw, iw, ixy, xare(iw), sfhta(iw), sfhtm(iw),
     >       sehta(iw), sehtm(iw), sehta(iw)+sehtm(iw)
        enddo
        write(i6,'(/2x,10x,"total",1p5e12.4)')
     >     zfa, zfm, zea, zem, zea+zem
        write(i6,'(2x,10x, "max  ",1pe12.4)') zmx
        write(i6,'(2x)')
      enddo
      close(i6)
!
      call ntwflx
!
      return
      end
