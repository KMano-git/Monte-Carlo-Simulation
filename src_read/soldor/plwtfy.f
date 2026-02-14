!***********************************************************************
      subroutine plwtfy
!***********************************************************************
!
!   weit for y-flux (upwind scheme)
!
!      q(j,i+1/2) = xwtm(j,i)*q(j,i) + xwtp(j,i)*q(j,i+1)
!
!   index check  O.K.
!
!      do it = 2, itmax-1     it = itpvs (22)
!      jmax = jtmax(it)          j/jp  = 1-31, 114-144
!      do jw = 1, jmax-1      it = itmps (32)
!        j  = jcel(jw,it)        j/jp  = 32-113
!        jp = jcel(jw+1,it)
!          ------
!      enddo               No dubly defined index (j,i)
!
!   metric
!
!     zf11 = 0.0 because hdsv(j,i,kce/w/n/s) = 0.0  ! mistake
!           No found hdsv in this routine.
!
!-----------------------------------------------------------------------
      use cplcom, only : ama, iupwd, nion, q1a, vdda, vnag, xfdfa, xwtam
     >    , xwtap
      use cplmet, only : gdsv, gwtm, gwtp, icel, itmax, itpve, jcel
     >    , jtmax, kce, kcn, kcs, kcw
      implicit none
!
!::local variables
      integer  it, i, ip, jmax, jw, j, jm, ia
      real*8   zdaan, Ade, Adw, Adn, Ads, zrae, zraw
      real*8   zf11, zwtm, zwtp
!
      do it = 1, itmax-1
        if( it.eq.itpve ) cycle
!
        i  = icel(1,it)
        ip = i + 1
        jmax = jtmax(it)
!
        do jw = 2, jmax-1
          j  = jcel(jw,it)
          jm = jcel(jw-1,it)
!
          do ia = 1, nion
!
!::d-coef (Da) at cell boundary (j,i+1/2) & (j,i-1/2)
            zdaan =
     >       1.0d0/(gwtm(j,i)/vdda(j,i,ia)+gwtp(j,i)/vdda(j,ip,ia))
            Ade = gdsv(j,i,kce)*zdaan
            Adw = gdsv(j,i,kcw)*zdaan
            Adn = gdsv(j,i,kcn)*zdaan
            Ads = gdsv(j,i,kcs)*zdaan
!
!::parameter (roa,va) at E(j+1/2,i+-1/2), W(j-1/2,i+-1/2)
            zrae = ama(ia)*vnag(j, i,ia)
            zraw = ama(ia)*vnag(jm,i,ia)
!
!::velocity across cell boundary
            zf11 = -Ade*zrae+Adw*zraw-Adn*q1a(j,ip,ia)+Ads*q1a(j,i,ia)
!
!::weit
            zwtm = gwtm(j,i)
            zwtp = gwtp(j,i)
!
!::up-wind
            if( iupwd(j,i).eq.1 ) then
              zwtm = 1.0d0
              zwtp = 0.0d0
              if( zf11.lt.0.0d0 ) then
                zwtm = 0.0d0
                zwtp = 1.0d0
              endif
            endif
!
            xwtam(j,i,ia) = zwtm
            xwtap(j,i,ia) = zwtp
            xfdfa(j,i,ia) = zf11
!
          enddo  ! loop(ia)
        enddo  ! loop(jw)
      enddo  ! loop(it)
!
      end
