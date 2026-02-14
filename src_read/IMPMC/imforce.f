!***********************************************************************
      subroutine imforce(ip,tauz)
!***********************************************************************
!
!   kinetic thermal force
!          D. Reiser at al., Nucl. Fusion 38 (1998) 165
!
!   comt   ip  ic   ix iy is   tt          Ez        Ti
!   gdti      glti      grdti     gdte      glte      grdte     ftot
!   ftot1     x         bet0      bet       vz0      dvz
!
!-----------------------------------------------------------------------
      use cimcom, only : aimas, ami, amz, azmas, cfgte, cfgti, fforc
     >    , gdte, gdti, glte, glti, ir, is, lfgti, lfgtzf, tforc, vvr
     >    , vz
      use cntpls, only : temi, vflw, zefm
      use cphcns, only : cev
      implicit none
!
!::argument
      integer, intent(in) :: ip
      real(8), intent(in) :: tauz
!
!::local variables
      integer  ic, iz
      real*8   grdti, grdte
!
!::local variables (ftot)
      real(8) :: za, zfca, zcti, zcte
      real(8) :: zti, zvf, zvi2, xx2, xxp2, zbet, ckp0
      real(8) :: zcti1
      integer ia
!
      ic  = ir(ip)
      iz  = is(ip)
      if( iz.le.0 ) return
      if( temi(ic).le.0.0d0 ) return

!
!---------------------------------------------------------------------
!::Fz = electric field + e-thermal + i-thermal force
      ckp0 = 1.0d0/0.5657d0
      ia  = 1
      zvf = vflw(ic,ia)
      za  = iz
      zti = temi(ic)
      zvi2 = 2.0d0*cev*zti/ami
      xx2  = (vvr(ip)+(vz(ip)-zvf)**2)/zvi2
      xxp2 = (vz(ip)-zvf)**2/zvi2
      zbet = 1.5d0*ckp0*(1.0d0-2.0d0*xxp2)*exp(-xx2)
      zcti1 = (1.0d0+aimas/azmas)*za**2*zbet
      zfca  = 1.0d0
      if( lfgtzf.eq.1 ) zfca  = 1.0d0/za**2*(za/zefm(ic)-1.0d0)*za
      zcti  = cfgti(iz)*zfca
      zcte  = cfgte(iz)*zfca
      if( lfgti.eq.1 ) zcti  = zcti1*zfca
!
!::new limit of temperature gradient
      grdti = dmin1( glti(ic), dabs(gdti(ic)) )
      grdti = dsign( grdti, gdti(ic) )
      grdte = dmin1( glte(ic), dabs(gdte(ic)) )
      grdte = dsign( grdte, gdte(ic) )
      tforc(ip) = (zcti*grdti+zcte*grdte)*cev ! Yamoto: Thermal force plot
      fforc(ip) = amz*(zvf-vz(ip))/tauz

      return
      end
