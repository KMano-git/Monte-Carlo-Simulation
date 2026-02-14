!*********************************************************************
      subroutine imp_pack
!*********************************************************************
      use cimcom, only : fi, ien, igi, igs, igtem, igtno, il, iml
     >    , ir, is, itg, ivd, ivs, ndmp, npmax, ntags, pemt, rr, tt, v
     >    , vr, vv, vvr, vvz, vz, wght, wstp, zz
      use cimntl, only : sint, srn, stb, stm, svx, svy, svz, sx0, sy0
     >    , sz0
      use cunit,  only : n6
      implicit none

!::local variables
      integer :: ii, jj, i, im, imax0, ip, npstp

      write(n6,'(/2x,"*** imp_pack ***")')
!::number of particle with wght = 0.0 (stop particle)
      npstp = 0
      do ip = 1, npmax
        if( wght(ip) == 0.0d0 ) npstp = npstp + 1
!        if( ien(ip).ge.6 ) npstp = npstp + 1 !SY this could be better
        if( ien(ip).ge.6 ) wght(ip)=0.0d0
      enddo

      write(n6,'(4x,"imp_pack   npmax =",i6,"  npstp =",i6)')
     >    npmax, npstp
      if( npstp == 0 ) return

      imax0 = npmax

      ii = 0
      jj = npmax + 1

 110  continue
      ii = ii + 1
      if( ii > npmax ) goto 200
      if( wght(ii) .ne. 0.0d0 ) goto 110

 120  continue
      jj = jj - 1
      if( jj <= ii ) goto 200
      if( wght(jj) .eq. 0.0d0 ) goto 120

!-----
!::rr(ii) (w=0)  <==  rr(jj) (w>0)
      ir    (ii) = ir    (jj)
      is    (ii) = is    (jj)
      il    (ii) = il    (jj)
      ivs   (ii) = ivs   (jj)
      ivd   (ii) = ivd   (jj)
      itg   (ii) = itg   (jj)
      igs   (ii) = igs   (jj)
      igi   (ii) = igi   (jj)
      iml   (ii) = iml   (jj)
      ien   (ii) = ien   (jj)
      igtno (ii) = igtno (jj)
      igtem (ii) = igtem (jj)
      tt    (ii) = tt    (jj)
      rr    (ii) = rr    (jj)
      zz    (ii) = zz    (jj)
      fi    (ii) = fi    (jj)
      vr    (ii) = vr    (jj)
      vz    (ii) = vz    (jj)
      vv    (ii) = vv    (jj)
      ntags(ii) = ntags(jj) ! SY Source tag
      pemt  (ii) = pemt  (jj)
      wght  (ii) = wght  (jj)
      wstp  (ii) = wstp  (jj)
      vvr   (ii) = vvr   (jj)
      vvz   (ii) = vvz   (jj)
      v     (ii) = v     (jj)

      sx0   (ii) = sx0   (jj)
      sy0   (ii) = sy0   (jj)
      sz0   (ii) = sz0   (jj)
      svx   (ii) = svx   (jj)
      svy   (ii) = svy   (jj)
      svz   (ii) = svz   (jj)
      srn   (ii) = srn   (jj)
      sint  (ii) = sint  (jj)
      stm   (ii) = stm   (jj)
      stb   (ii) = stb   (jj)
!zeroset jj
!-----
      ir    (jj) = 0
      is    (jj) = 0
      il    (jj) = 0
      ivs   (jj) = 0
      ivd   (jj) = 0
      itg   (jj) = 0
      igs   (jj) = 0
      igi   (jj) = 0
      iml   (jj) = 0
      ien   (jj) = 0
      igtno (jj) = 0
      igtem (jj) = -9
      ntags (jj) = 0
      tt    (jj) = 0.0d0
      rr    (jj) = 0.0d0
      zz    (jj) = 0.0d0
      fi    (jj) = 0.0d0
      vr    (jj) = 0.0d0
      vz    (jj) = 0.0d0
      vv    (jj) = 0.0d0
      pemt  (jj) = 0.0d0
      wght  (jj) = 0.0d0
      wstp  (jj) = 0.0d0
      vvr   (jj) = 0.0d0
      vvz   (jj) = 0.0d0
      v     (jj) = 0.0d0

      sx0   (jj) = 0.0d0
      sy0   (jj) = 0.0d0
      sz0   (jj) = 0.0d0
      svx   (jj) = 0.0d0
      svy   (jj) = 0.0d0
      svz   (jj) = 0.0d0
      srn   (jj) = 0.0d0
      sint  (jj) = 0.0d0
      stm   (jj) = 0.0d0
      stb   (jj) = 0.0d0

!      wght(jj) = wghtsv
!      ien(jj) = wghtsv ! SY just zeroset jj is much better
      ! other variables(jj) will be overwritten
!-----
      goto 110

 200  continue
      do i = npmax, 1, -1
        im = i
        if( wght(i) .ne. 0.0d0 ) goto 210
      enddo
      im = 0
 210  continue
      npmax = im
      write(n6,'(4x,"imp_pack   npmax =",i6," =>",i6)')
     >   imax0, npmax
      write(n6,'(2x,"wght =",20f8.5)') (wght(ip),ip=1,20)

      return
      end
