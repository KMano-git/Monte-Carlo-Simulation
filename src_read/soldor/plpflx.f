!***********************************************************************
      subroutine plpflx
!***********************************************************************
      use cntcom, only : isxp, isyp, lrcmd, npep, npsp, temin_rec
      use cntmnt, only : mcflx, mfmax, n6_src, pfl_ion, pfl_ntl
      use cplcom, only : ama, flps, flvl, gfvl, nion, nlp, tfps, tfvl
     >    , vna, vne, vte
      use cplmet, only : hvol, icel, itmax, itpve, jcel, jtmax
      use cplqcn, only : qfx_cv, qfx_df, qfy_df
      use csize,  only : ndwl, ndwp, ndx, ndy
      use csonic, only : itim, lstop, time
      implicit none
!
!::local variables
      integer  k1, k2, k3
      integer  it, nt, ia, m1a, nw, iw, j, i, iflx
      integer  jmax, jt
      real*8   zsum, zfi, tfion, tfntl
      real*8   ane, ani, ate, asv, asn, svhrcm, zmax, rclim
!
      rclim = 1.0d-3
!
!-----------------------------------------------------------------------
!::clear
!-----------------------------------------------------------------------
      call plqcon_cl
!
      do k2=1,nion
        do k1=1,ndwl
          tfps(k1,k2) = 0.0d0
        enddo
      enddo
!
      do k2=1,nion
        do k1=1,ndwp
          flps(k1,k2) = 0.0d0
        enddo
      enddo
!
      do k1=1,nion
        tfvl(k1) = 0.0d0
        gfvl(k1) = 0.0d0
      enddo
!
      do k3=1,nion
        do k2=1,ndy
          do k1=1,ndx
            flvl(k1,k2,k3) = 0.0d0
          enddo
        enddo
      enddo
!
!-----------------------------------------------------------------------
!::ion flux
!-----------------------------------------------------------------------
      do nt = 2, itmax-1
        if( nt.eq.itpve ) cycle
!
!::clear  ! KSFUJI
        call pxcler(nt)
!
!::convection term
        call pxconv(nt)
        if( lstop.eq.1 ) return
!
!::viscosity term
        call pxdiff(nt)
        call pxvisc(nt)
!
!::residual of F^psi -- gg
        call pxrdyv(nt)
        call pxrdvp(nt)  ! vpinch
!
      enddo  ! loop(nt)
!
      do ia = 1, nion
        m1a = 2*ia - 1
!
!::flux (sol)
        nw = 1
        zsum = 0.0d0
        do iw = npsp(nw), npep(nw)-1
          j = isxp(iw)
          i = isyp(iw)
          zfi = -qfy_df(j,i,m1a)/ama(ia)
          flps(iw,ia) = zfi
          zsum = zsum + zfi
        enddo
        tfps(nw,ia) = zsum
!
!::flux (idp)
        nw = 2
        zsum = 0.0d0
        do iw = npsp(nw), npep(nw)-1
          j = isxp(iw)
          i = isyp(iw)
          zfi = (qfx_cv(j,i,m1a)+qfx_df(j,i,m1a))/ama(ia)
          flps(iw,ia) = zfi
          zsum = zsum + zfi
        enddo
        tfps(nw,ia) = zsum
!
!::flux (prv)
        nw = 3
        zsum = 0.0d0
        do iw = npsp(nw), npep(nw)-1
          j = isxp(iw)
          i = isyp(iw)
          zfi = qfy_df(j,i,m1a)/ama(ia)
          flps(iw,ia) = zfi
          zsum = zsum + zfi
        enddo
        tfps(nw,ia) = zsum
!
!::flux (odp)
        nw = 4
        zsum = 0.0d0
        do iw = npsp(nw), npep(nw)-1
          j = isxp(iw)
          i = isyp(iw)
          zfi = -(qfx_cv(j,i,m1a)+qfx_df(j,i,m1a))/ama(ia)
          flps(iw,ia) = zfi
          zsum = zsum + zfi
        enddo
        tfps(nw,ia) = zsum
!
!::flux (man)
        nw = 5
        zsum = 0.0d0
        do iw = npsp(nw), npep(nw)-1
          j = isxp(iw)
          i = isyp(iw)
          zfi = qfy_df(j,i,m1a)/ama(ia)
          flps(iw,ia) = zfi
          zsum = zsum + zfi
        enddo
        tfps(nw,ia) = zsum
!
      enddo  !  loop(ia)
!
!-----------------------------------------------------------------------
!::recombination ( H+ + e ==> H0 )
!-----------------------------------------------------------------------
      if( lrcmd.eq.1 ) then
        do ia = 1, nion
          zmax = -1.0d20
          do it = 2, itmax-1
            if( it.eq.itpve ) cycle
            jmax = jtmax(it)
            do jt = 2, jmax-1
              j = jcel(jt,it)
              i = icel(jt,it)
              ane = vne(j,i)
              ani = vna(j,i,ia)
              ate = vte(j,i)
              if( ane.le.0.0d0 .or. ate.le.0.0d0 ) cycle
              if( ate.ge.10.0d0 ) cycle  !  2.0 ==> 10.0   2008/07/25
              ate = dmax1(ate,temin_rec)  ! svhrcm
              asv = svhrcm(ate,ane)
              asn = ane*ani*asv*hvol(j,i)
              zmax = dmax1( zmax, asn )
              flvl(j,i,ia) = asn
            enddo  ! loop(jt)
          enddo  ! loop(it)
          gfvl(ia) = zmax
!
          zsum = 0.0d0
          do it = 2, itmax-1
            if( it.eq.itpve ) cycle
            jmax = jtmax(it)
            do jt = 2, jmax-1
              j = jcel(jt,it)
              i = icel(jt,it)
              if( flvl(j,i,ia)/gfvl(ia) .le. rclim ) then
                flvl(j,i,ia) = 0.0d0
              endif
              zsum = zsum + flvl(j,i,ia)
            enddo
          enddo
          tfvl(ia) = zsum
!
        enddo  ! loop(ia)
      endif
!
!-----------------------------------------------------------------------
!::birth of neutral particle
!-----------------------------------------------------------------------
      do iflx = 1, mfmax
        call ntnflx(iflx,tfion,tfntl)
        pfl_ion(iflx) = tfion
        pfl_ntl(iflx) = tfntl
      enddo
!
      return
      end