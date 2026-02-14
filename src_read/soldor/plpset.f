!**********************************************************************
      subroutine plpset
!**********************************************************************
      use cplcom, only : aion, ama, aza, cfab, cfeb, cfrab, ckab, clpa
     >    , crod, nion
      use cphcns, only : cme, cmp
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  ia, ib
      real*8   zaab
! 
      write(n6,'(/2x,"*** plpset ***")')
!
!-----------------------------------------------------------------------
!::physical constant
!-----------------------------------------------------------------------
      do ia = 1, nion
      ama(ia)  = aion(ia)*cmp
      clpa(ia) = 0.0d0
      enddo
      clpa(nion) = 1.0d0
!
      write(n6,'(2x," cme  = ",1pe11.3)') cme
      write(n6,'(2x," aion = ",1p6e11.3)') (aion(ia),ia=1,nion)     
      write(n6,'(2x," ama  = ",1p6e11.3)') (ama (ia),ia=1,nion)
      write(n6,'(2x," clpa = ",1p6e11.3)') (clpa(ia),ia=1,nion)
!
!-----------------------------------------------------------------------
!::coulomb logatithm      see  Braginskii or BALDUR code
!-----------------------------------------------------------------------
!   clge = ln( 10.9e10*Te/sqrt(ne)/zeff )            Te < 50 eV
!        = ln( 1.54e10*Te*sqrt(Te)/sqrt(ne)/zeff )   Te > 50 eV
!   clgi = ln( 1.54e10*Ti*sqrt(Te)/sqrt(ne)/zeff ) 
!                Note  Te [eV], ne [1/cm**3]
!                          
!-----------------------------------------------------------------------
!::collision frequency    see  Bramms code
!-----------------------------------------------------------------------
!  1/tau_eb = cfeb(ib)   *clge*nb/Te**1.5
!  1/tau_ab = cfab(ia,ib)*clgi*nb/Ti**1.5
!                Note.  Te, Ti [eV],  nb [1/m**3]
!
!    eta_a  = 0.96*na*Ti / sum_b{1/tau_ab}
!    ke     = 3.16*ne*Te/me / sum_b{1/tau_eb}
!    ki     = sum_a{ 3.90*na*Ti/mp / sum_b{sqrt(Aa*Ab)*(1/tau_ab)} }
!    feq    = sum_b{ 3*me/mp*(1/Ab*1/tau_eb) }
!
      do ia = 1, nion
      cfeb(ia) = aza(ia)**2/3.4411d11
      do ib = 1, nion     
      cfab(ia,ib) = aza(ia)**2*aza(ib)**2
     >  /(sqrt((aion(ia)+aion(ib))/2.0d0)*2.0853d13)
      ckab(ia,ib) = cfab(ia,ib)*sqrt(aion(ia)*aion(ib))
      enddo
      enddo
!
!-----------------------------------------------------------------------
!::friction force
!-----------------------------------------------------------------------
!   fr-force = cfrab(ia,ib)*clg*na*nb/Ti**1.5
      do ia = 1, nion
      do ib = 1, nion
      zaab = aion(ia)*aion(ib)/(aion(ia)+aion(ib))
      cfrab(ia,ib) = 6.7819d-14*cmp*sqrt(zaab)*aza(ia)**2*aza(ib)**2
      if(ia.eq.ib) cfrab(ia,ib) = 0.0d0
      enddo
      enddo
!
!-----------------------------------------------------------------------
!::cronecker delta
!-----------------------------------------------------------------------
      do ia = 1, nion
      do ib = 1, nion
      crod(ia,ib) = 0.0d0
      enddo
      crod(ia,ia) = 1.0d0
      enddo
!
      return
      end
