!**********************************************************************
      subroutine set_gcoef
!**********************************************************************
      use cimcom, only : aimas, azmas, cfez, cfgte, cfgti, ishez, ishte
     >    , ishti, ismax, lfgtzf, ndis
      use cunit,  only : n6
      implicit none
!
!::local variables
      real*8   zmu, sq2, sqmu
      real(8) :: za, zef, zfca
      real(8), dimension(0:ndis) :: hfgte, hfgti, hfca
      integer  i
!
!----------------------------------------------------------------------
!::alfz betz
!----------------------------------------------------------------------
      zmu  = azmas/(azmas+aimas)
      sq2  = sqrt(2.0)
      sqmu = sqrt(zmu)
      do i = 1, ismax
      za   = i
      cfez(i)  = za
      cfgte(i) = 0.71*za*za
!<==  cfgti(i) ~ 2.20*za*za
      cfgti(i) =-3.0*( 1.0-zmu
     >      -5.0*sq2*za*za*(1.1*zmu*zmu*sqmu-0.35*zmu*sqmu) )
     >      / ( 2.6-2.0*zmu+5.4*zmu*zmu )
!
      if( ishez.eq.0 ) cfez(i)  = 0.0d0
      if( ishte.eq.0 ) cfgte(i) = 0.0d0
      if( ishti.eq.0 ) cfgti(i) = 0.0d0
!
      zef  = 1.0d0
      zfca = 1.0d0
      if( lfgtzf.eq.1 ) zef = 2.5d0
      if( lfgtzf.eq.1 )   zfca = 1.0d0/za**2*(za/zef-1.0d0)*za
      hfgte(i) = cfgte(i)*zfca
      hfgti(i) = cfgti(i)*zfca
      hfca(i)  = zfca
      enddo
!
      write(n6,'(/2x,"*** set_gcoef ***")')
!
      write(n6,602) aimas,azmas,zmu
  602 format(3x,'aimas =',f9.3,'  azmas =',f9.3,'  zmu =',f9.3)
      write(n6,604) 'cfez  =',(cfez(i),i=1,ismax)
      write(n6,604) 'cfgte =',(cfgte(i),i=1,ismax)
      write(n6,604) 'cfgti =',(cfgti(i),i=1,ismax)
!
      write(n6,'(/2x,"lfgtzf =",i2,"  zef =",f8.4)') lfgtzf, zef
      if( lfgtzf.eq.1 ) then
      write(n6,'(2x,"cfgte = 0.71*za**2 ==> 0.71*(za/zef-1)*za")')
      write(n6,'(2x,"cfgti = 2.20*za**2 ==> 2.20*(za/zef-1)*za")')
      endif
      write(n6,604) 'hfca  =',(hfca(i),i=1,ismax)
      write(n6,604) 'hfgte =',(hfgte(i),i=1,ismax)
      write(n6,604) 'hfgti =',(hfgti(i),i=1,ismax)
      write(n6,604) 'rtzfI =',(hfgti(i)/cfgti(i),i=1,ismax)
      write(n6,604) 'rtzfE =',(hfgte(i)/cfgte(i),i=1,ismax)
  604 format(3x,a7,1p10e12.3/(3x,7x,1p10e12.3))
!
!----------------------------------------------------------------------
!::total force
!----------------------------------------------------------------------
!  Fz(ic,iz) = cfez(iz)*gdez(ic)+cfgte(iz)*gdte(ic)+cfgti(iz)*gdti(ic)
!
      return
      end
