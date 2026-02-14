!**********************************************************************
      subroutine ntscor_cl(krct)
!**********************************************************************
!
!      collision estimator
!
!    source term   N0 * Ni*<sig*v> * dlt_(moment,energy)
!
!       N0 = Norm * weit
!       Ni*<sig*v> =  frequency of calling this routine
!
!    Note
!   --irct = 9      "DS5"   [ e  + H2+ =>  H+  + H+ + 2e ]
!    ntscor_fl : zwi = zdn*(1.5d0*ze0*cev+derc(2)*2.0d0*cev)
!    ntscor_tr : zwi = zdn*(0.5d0*tmas*vel2+derc(2)*2.0d0*cev)
!    ntscor_cl : zwi = zdn*(0.5d0*tmas*vel2*2.0d0+derc(2)*2.0d0*cev)
!          because tmas = mass(D2) in fl & tr tmas = mass(D) in cl
!
!----------------------------------------------------------------------
      use celcom, only : dcl_ipl, pcl_eng, pcl_mvp
      use cntcom, only : derc, elrc, evelb, icpo, igas, igasb, mrgn
     >    , noia, rmas, scfc, tmas, tmasb, vel2, velpb, weit, weitb
      use cntpls, only : temi, vflw
      use cntscg, only : sn_mc, sp_mc, tsn_mc, tsp_mc, twe_mc, twi_mc
     >    , we_mc, wi_mc
      use cphcns, only : cev
      implicit none
!
!--argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer  krct
      integer, intent(in) :: krct
!----------------------------------------------------------------------
!xc--common variables
!x      real*8   vflw(1,1),temi(1)                            ! cntpls
!x      integer  icpo, igas, igasb, noia(1);  real*8 rmas(1)  ! cntcom
!x      real*8   scfc(2),elrc(2),derc(1)                      ! cntcom
!x      real*8   weit,weitb,tmas,tmasb,velpb,evelb,vel2       ! cntcom
!x      real*8   sn_mc(1,1),sp_mc(1,1),we_mc(1),wi_mc(1)      ! cntscg
!x      real*8   pcl_mvp,pcl_eng,pcl_ipl                      ! celcom
!x      integer   dcl_ipl
!x      real*8   cev                                          ! cphcns
!----------------------------------------------------------------------
!--local variables
      integer  ic, ig, ia, ig2, ia2
      integer  irc, irg
      real*8   zdn,zsn,zsp,zwe,zwi,zsp2
!----------------------------------------------------------------------
!
      ic   = icpo
      irg  = mrgn(ic)
!
!::electron ionziation  H0(ig) + e ==> H+(ig) + 2e
!
      if( krct.eq.1 ) then
      irc = 1
      ig  = igasb
      ia  = noia(ig)
      zdn = (weitb-weit) * scfc(1)
      zsn =  zdn
      zsp =  zdn*tmasb*velpb
      zwe = -zdn*elrc(1)*cev
      zwi =  zdn*evelb*cev
      sn_mc(ic,ia) = sn_mc(ic,ia) + zsn
      sp_mc(ic,ia) = sp_mc(ic,ia) + zsp
      we_mc(ic)    = we_mc(ic)    + zwe
      wi_mc(ic)    = wi_mc(ic)    + zwi
!-----
      tsn_mc(irc,irg) = tsn_mc(irc,irg) + zsn
      tsp_mc(irc,irg) = tsp_mc(irc,irg) + zsp
      twe_mc(irc,irg) = twe_mc(irc,irg) + zwe
      twi_mc(irc,irg) = twi_mc(irc,irg) + zwi
!
!::charge exchange  H0(ig) + H+(ig2) ==> H+(ig) + H0(ig2)
!                   after collision      igasb     igas
!
      elseif( krct.eq.2 ) then
      irc = 2
      ig  = igasb
      ig2 = igas
      ia  = noia(ig)
      ia2 = noia(ig2)
!
      zdn = weit * scfc(2)
      zsn = zdn
      zsp = zdn*tmasb*velpb
      zsp2= zdn*rmas(ig2)*vflw(ic,ig2)
      zwi = zdn*evelb*cev
     >     -zdn*(1.5d0*temi(ic)*cev
     >                      +0.5d0*rmas(ig2)*vflw(ic,ig2)**2)
      sn_mc(ic,ia ) = sn_mc(ic,ia ) + zsn
      sn_mc(ic,ia2) = sn_mc(ic,ia2) - zsn   !%% if check CX source
      sp_mc(ic,ia ) = sp_mc(ic,ia ) + zsp
      sp_mc(ic,ia2) = sp_mc(ic,ia2) - zsp2
      wi_mc(ic)     = wi_mc(ic)     + zwi
!-----
      tsn_mc(irc,irg) = tsn_mc(irc,irg) + zsn
      tsn_mc(irc,irg) = tsn_mc(irc,irg) - zsn   !%% if check CX source
      tsp_mc(irc,irg) = tsp_mc(irc,irg) + zsp
      tsp_mc(irc,irg) = tsp_mc(irc,irg) - zsp2
      twi_mc(irc,irg) = twi_mc(irc,irg) + zwi
!
!--irct = 3  "EL1"    [ H0(ig) + H+(ig2) => H0(ig) + H+(ig2) ]
      elseif( krct.eq.3 ) then
      irc = 3
      ia2 = dcl_ipl
      zdn = weit * scfc(3)
      zsp = -zdn*pcl_mvp
      zwi = -zdn*pcl_eng
      sp_mc(ic,ia2) = sp_mc(ic,ia2) + zsp
      wi_mc(ic)     = wi_mc(ic)     + zwi
!-----
      tsp_mc(irc,irg) = tsp_mc(irc,irg) + zsp
      twi_mc(irc,irg) = twi_mc(irc,irg) + zwi
!
!--irct = 4      "DS1"   [ e  + H2  =>  H0  + H0 + e ]
      elseif( krct.eq.4 ) then
      irc = 4
      zdn = weit * scfc(4)
      zwe = -zdn*elrc(1)*cev
      we_mc(ic) = we_mc(ic) + zwe
!-----
      twe_mc(irc,irg) = twe_mc(irc,irg) + zwe
!
!---irct = 5      "DS2"   [ e  + H2  =>  H2+ + 2e ]
      elseif( krct.eq.5 ) then
      irc = 5
      zdn = weit * scfc(5)
      zwe = -zdn*elrc(2)*cev
      we_mc(ic) = we_mc(ic) + zwe
!-----
      twe_mc(irc,irg) = twe_mc(irc,irg) + zwe
!
!---irct = 6      "DS3"   [ e  + H2  =>  H+  + H0 + e ]
      elseif( krct.eq.6 ) then
      irc = 6
      ia  = noia(igas)
      zdn = weit * scfc(6)
      zsn = zdn
      zwi = zdn*(0.5d0*tmas*vel2+derc(3)*cev)   ! tmas H+
      zwe = -zdn*elrc(3)*cev
      sn_mc(ic,ia) = sn_mc(ic,ia) + zsn
      wi_mc(ic)    = wi_mc(ic)    + zwi
      we_mc(ic)    = we_mc(ic)    + zwe
!-----
      tsn_mc(irc,irg) = tsn_mc(irc,irg) + zsn
      twi_mc(irc,irg) = twi_mc(irc,irg) + zwi
      twe_mc(irc,irg) = twe_mc(irc,irg) + zwe
!
!--irct = 7  "EL2"    [ H2(ig) + H+(ig2) => H2(ig) + H+(ig2) ]
      elseif( krct.eq.7 ) then
      irc = 7
      ia2 = dcl_ipl
      zdn = weit * scfc(7)
      zsp = -zdn*pcl_mvp
      zwi = -zdn*pcl_eng
      sp_mc(ic,ia2) = sp_mc(ic,ia2) + zsp
      wi_mc(ic)     = wi_mc(ic)     + zwi
!-----
      tsp_mc(irc,irg) = tsp_mc(irc,irg) + zsp
      twi_mc(irc,irg) = twi_mc(irc,irg) + zwi
!
!--irct = 8      "DS4"   [ e  + H2+ =>  H+  + H0 + e ]
      elseif( krct.eq.8 ) then
      irc = 8
      ia  = noia(igas)
      zdn = weit * scfc(8)
      zsn = zdn
      zwi = zdn*(0.5d0*tmas*vel2+derc(1)*cev)
      zwe = -zdn*elrc(1)*cev
      sn_mc(ic,ia) = sn_mc(ic,ia) + zsn
      wi_mc(ic)    = wi_mc(ic)    + zwi
      we_mc(ic)    = we_mc(ic)    + zwe
!-----
      tsn_mc(irc,irg) = tsn_mc(irc,irg) + zsn
      twi_mc(irc,irg) = twi_mc(irc,irg) + zwi
      twe_mc(irc,irg) = twe_mc(irc,irg) + zwe
!
!--irct = 9      "DS5"   [ e  + H2+ =>  H+  + H+ + 2e ]
      elseif( krct.eq.9 ) then
      irc = 9
      ia  = noia(igas)
      zdn = weit * scfc(9)
      zsn = zdn*2.0d0
      zwi = zdn*(0.5d0*tmas*vel2*2.0d0+derc(2)*2.0d0*cev)
      zwe = -zdn*elrc(2)*cev
      sn_mc(ic,ia) = sn_mc(ic,ia) + zsn
      wi_mc(ic)    = wi_mc(ic)    + zwi
      we_mc(ic)    = we_mc(ic)    + zwe
!-----
      tsn_mc(irc,irg) = tsn_mc(irc,irg) + zsn
      twi_mc(irc,irg) = twi_mc(irc,irg) + zwi
      twe_mc(irc,irg) = twe_mc(irc,irg) + zwe
!
!--irct = 10      "DS6"   [ e  + H2+ =>  H0  + H0 ]
      elseif( krct.eq.10 ) then
      irc = 10
      zdn = weit * scfc(10)
      zwe = -zdn*elrc(3)*cev
      we_mc(ic) = we_mc(ic) + zwe
!-----
      twe_mc(irc,irg) = twe_mc(irc,irg) + zwe
!
!--irct = 11      "RC"    [ e  + H+  =>  H0 ]
      elseif( krct.eq.11 ) then
      irc = 11
!          calculate source terms due to recombination at ntrecm
      call wexit("ntscor_cl","wrong krct = 11")

! Collision estimator for n-n collision has not been implemented yet---- toku
      elseif( krct.eq.12 ) then
      irc = 12
      elseif( krct.eq.13 ) then
      irc = 13
!
!--irct = ?
      else
      call wexit("ntscor_cl","wrong krct ?")
      endif
!
      return
      end
