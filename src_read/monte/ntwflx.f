!**********************************************************************
      subroutine ntwflx
!**********************************************************************
!
!    neutral flux onto the wall
!       normalization factor : swnrm    see sub. ntwcon
!
!----------------------------------------------------------------------
      use cntcom,  only : npew, npsw, npsp, npep
      use cntmnt,  only : mfmax, pfl_ntl, vcsrc, visrc
      use cntwcn,  only : swnrm, wcsrc, weema, weemm, wees, wehta, wehtm
     >    , wfhta, wfhtm, wisrc
      use cntwedf, only : wedf, xwedf, wedf_line, xwedf_line
      use cntwfl,  only : xare, xcfl, xeema, xeemm, xees
     >    , xehta, xehtm, xfhta, xfhtm, xsno, xtea, xtees, xtem, xtfa
     >    , xtfm
      use cphcns,  only : cev
      use csize,   only : ndmsr, ndwp
      use csonic,  only : itim, time
      use cunit,   only : n6
      use mod_externalgrid, only : use_exdata
      implicit none
!
!::local variables
      integer  nsiz, nw, iw, iws, iwe, iflx, isrc, i, is, iw_p
      real*8   area
      real*8 zfa, zfm, zea, zem, zeea, zeem, znrm, tfntl
!
      write(n6,'(/2x,"*** ntwflx ***  itim =",i7,"  time =",1pe14.6)')
     >    itim, time
!
!::clear
      nsiz = ndwp*(ndmsr+1)
      call setd( xfhta, nsiz, 0.0d0 )
      call setd( xfhtm, nsiz, 0.0d0 )
      call setd( xehta, nsiz, 0.0d0 )
      call setd( xehtm, nsiz, 0.0d0 )
      call setd( xeema, nsiz, 0.0d0 )
      call setd( xeemm, nsiz, 0.0d0 )
      call setd( xees,  nsiz, 0.0d0 )
      xwedf = 0.0d0
      xwedf_line = 0.0d0
!
!::sumup for source
      do iflx = 1, mfmax
        isrc  = visrc(iflx)
        if( isrc.le.0 ) cycle
        tfntl = pfl_ntl(iflx)
        xcfl(isrc) = vcsrc(isrc)
        wisrc = isrc
        wcsrc = vcsrc(isrc)
        call setvmwork( 1, isrc )
        call ntwcon
!
        zfa = 0.0d0
        zfm = 0.0d0
        zea = 0.0d0
        zem = 0.0d0
        znrm = 1.0d0/swnrm*tfntl

        do nw = 1, 4
          iws = npsw(nw)
          iwe = npew(nw)
          !:: wall boundary
          do iw = iws, iwe
            xfhta(iw,0) = xfhta(iw,0) + wfhta(iw)*znrm
            xfhtm(iw,0) = xfhtm(iw,0) + wfhtm(iw)*znrm
            xehta(iw,0) = xehta(iw,0) + wehta(iw)*znrm*cev
            xehtm(iw,0) = xehtm(iw,0) + wehtm(iw)*znrm*cev
            xeema(iw,0) = xeema(iw,0) + weema(iw)*znrm*cev
            xeemm(iw,0) = xeemm(iw,0) + weemm(iw)*znrm*cev
!
            xfhta(iw,isrc) = xfhta(iw,isrc) + wfhta(iw)*znrm
            xfhtm(iw,isrc) = xfhtm(iw,isrc) + wfhtm(iw)*znrm
            xehta(iw,isrc) = xehta(iw,isrc) + wehta(iw)*znrm*cev
            xehtm(iw,isrc) = xehtm(iw,isrc) + wehtm(iw)*znrm*cev
            xeema(iw,isrc) = xeema(iw,isrc) + weema(iw)*znrm*cev
            xeemm(iw,isrc) = xeemm(iw,isrc) + weemm(iw)*znrm*cev
          enddo  ! loop(iw)
!
          !:: SOL boundary
          iws = npsp(nw)
          iwe = npep(nw)
          do iw_p = iws, iwe
            xees(iw_p,0)    = xees(iw_p,0)    + wees(iw_p)*znrm*cev
            xees(iw_p,isrc) = xees(iw_p,isrc) + wees(iw_p)*znrm*cev
          enddo ! loop(iw_p)
!
        enddo  ! loop(nw)
        xwedf(:,:,:) = xwedf(:,:,:) + wedf(:,:,:)*znrm
        xwedf_line(:,:,:) = xwedf_line(:,:,:) + wedf_line(:,:,:)*znrm
        xsno = isrc
      enddo  ! loop(iflx)
      xcfl(0) = "tot"
!
!::flux density
      do isrc = 0, xsno
        do nw = 1, 4
          iws = npsw(nw)
          iwe = npew(nw)
          do iw = iws, iwe
            zfa = 0.0d0
            zfm = 0.0d0
            zea = 0.0d0
            zem = 0.0d0
            zeea = 0.0d0
            zeem = 0.0d0
            if( xare(iw) > 0.0d0 ) then
              area = 1.d0/xare(iw)
              zfa = xfhta(iw,isrc)*area
              zfm = xfhtm(iw,isrc)*area
              zea = xehta(iw,isrc)*area
              zem = xehtm(iw,isrc)*area
              zeea = xeema(iw,isrc)*area
              zeem = xeemm(iw,isrc)*area
            endif
            xfhta(iw,isrc) = zfa
            xfhtm(iw,isrc) = zfm
            xehta(iw,isrc) = zea
            xehtm(iw,isrc) = zem
            xeema(iw,isrc) = zeea
            xeemm(iw,isrc) = zeem
          enddo
        enddo
      enddo
!
!::total flux
      nsiz = ndmsr+1
      call setd( xtfa, nsiz, 0.0d0 )
      call setd( xtfm, nsiz, 0.0d0 )
      xtea = 0.0d0
      xtem = 0.0d0
      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
        do iw = iws, iwe
          if(.not.use_exdata .and. iw==iwe) cycle
          do is = 0, xsno
            area = xare(iw)
            xtfa(is) = xtfa(is) + xfhta(iw,is)*area
            xtfm(is) = xtfm(is) + xfhtm(iw,is)*area
            xtea(is,1) = xtea(is,1) + xehta(iw,is)*area
            xtem(is,1) = xtem(is,1) + xehtm(iw,is)*area
            xtea(is,2) = xtea(is,2) + xeema(iw,is)*area
            xtem(is,2) = xtem(is,2) + xeemm(iw,is)*area
          enddo
        enddo
      enddo
!
      xtees = 0.0d0
      do nw = 1,3,2 ! only for sol/prv region
        iws = npsp(nw)
        iwe = npep(nw)
        do iw_p = iws, iwe
          xtees(nw) = xtees(nw) + xees(iw_p,0)
        enddo
      enddo

      return
      end