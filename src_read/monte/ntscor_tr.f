!**********************************************************************
      subroutine ntscor_tr(dltm)
!**********************************************************************
!
!    trak length estimator
!
!     mrct = 1 (D0),  = 2 (D2),  = 3 (D2+)
!     dltm  :  test particle stayed in a cell (icpo) during dltm
!
!      Note.  when check CX-event, comment this statement
!          ssnt(ic,ia2) = ssnt(ic,ia2) - zsn  !%% if check CX source
!
!----------------------------------------------------------------------
      use BGProf, only : den0BG, dengBG
      use celcom, only : dtr_eng, dtr_mvp
      use cntcom, only : derc, elrc, evel, icpo, igas, mrct, mrgn, ngas
     >    , noia, nrct, rmas, scfc, srct, tmas, vel2, velp, weit
      use cntpls, only : deni, vflw, temi
      use cntscg, only : sn_mt, sp_mt, tsn_mt, tsp_mt, twe_mt, twi_mt
     >    , we_mt, wi_mt
      use cphcns, only : cev
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   real*8 dltm
      real(8), intent(in) :: dltm
!
!::local varaiables
! modified 2/2 lines organize local variables and include files by kamata 2021/06/16
!ik   integer llp, ic, ig, ia, ig2, ia2, irc, irg
!ik   real*8  zfrq, zels, zdn, zsn, zsp, zsp2, zwe, zwi
      integer ic, ig, ia, ig2, ia2, irc, irg
      real*8  zdn, zsn, zsp, zsp2, zwe, zwi
      real*8  zn0, zni
!
!----------------------------------------------------------------------
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   data llp/0/
!
      ic   = icpo
      ig   = igas
      irg  = mrgn(ic)
!
!::vacume
      if( nrct.eq.0 ) return
!
!----------------------------------------------------------------------
!::D0
!----------------------------------------------------------------------
      if( mrct.eq.1 ) then
      zn0 = dltm*weit
!
!::electron ionziation [1] H0(ig) + e ==> H+(ig) + 2e
!
      irc = 1
      ia  = noia(ig)
      zdn = srct(1)*zn0 * scfc(1)
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   zsn = zdn
      zsp = zdn*tmas*velp
      zwe =-zdn*elrc(1)*cev
      zwi = zdn*evel*cev
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   sn_mt(ic,ia) = sn_mt(ic,ia) + zsn
      sn_mt(ic,ia) = sn_mt(ic,ia) + zdn
      sp_mt(ic,ia) = sp_mt(ic,ia) + zsp
      we_mt(ic)    = we_mt(ic)    + zwe
      wi_mt(ic)    = wi_mt(ic)    + zwi
!-----
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   tsn_mt(irc,irg) = tsn_mt(irc,irg) + zsn
      tsn_mt(irc,irg) = tsn_mt(irc,irg) + zdn
      tsp_mt(irc,irg) = tsp_mt(irc,irg) + zsp
      twe_mt(irc,irg) = twe_mt(irc,irg) + zwe
      twi_mt(irc,irg) = twi_mt(irc,irg) + zwi
!
!::charge exchange  H0(ig) + H+(ig2) ==> H+(ig) + H0(ig2)
!
      irc = 2
      ia  = noia(ig)
      do ig2 = 1, ngas
      ia2 = noia(ig2)
      zdn = srct(1+ig2)*zn0 * scfc(2)
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   zsn = zdn
      zsp = zdn*tmas*velp
      zsp2= zdn*rmas(ig2)*vflw(ic,ig2)
      zwi = zdn*evel*cev
     >     -zdn*(1.5d0*temi(ic)*cev
     >                      +0.5d0*rmas(ig2)*vflw(ic,ig2)**2)
! modified 2/2 lines organize local variables and include files by kamata 2021/06/16
!ik   sn_mt(ic,ia ) = sn_mt(ic,ia ) + zsn
!ik   sn_mt(ic,ia2) = sn_mt(ic,ia2) - zsn  !%% if check CX source
      sn_mt(ic,ia ) = sn_mt(ic,ia ) + zdn
      sn_mt(ic,ia2) = sn_mt(ic,ia2) - zdn  !%% if check CX source
      sp_mt(ic,ia ) = sp_mt(ic,ia ) + zsp
      sp_mt(ic,ia2) = sp_mt(ic,ia2) - zsp2
      wi_mt(ic)     = wi_mt(ic)     + zwi
!-----
! modified 2/2 lines organize local variables and include files by kamata 2021/06/16
!ik   tsn_mt(irc,irg) = tsn_mt(irc,irg) + zsn
!ik   tsn_mt(irc,irg) = tsn_mt(irc,irg) - zsn  !%% if check CX source
      tsn_mt(irc,irg) = tsn_mt(irc,irg) + zdn
      tsn_mt(irc,irg) = tsn_mt(irc,irg) - zdn  !%% if check CX source
      tsp_mt(irc,irg) = tsp_mt(irc,irg) + zsp
      tsp_mt(irc,irg) = tsp_mt(irc,irg) - zsp2
      twi_mt(irc,irg) = twi_mt(irc,irg) + zwi
      enddo   !  loop (ig2)
!
!::elactic collision  [3]  "EL1"  irct = 3
!                         [ H0(ig) + H+(ig2) => H0(ig) + H+(ig2) ]
      irc = 3
      ia  = noia(ig)
      do ig2 = 1, ngas
      ia2 = noia(ig2)
      call ntsctr_el(irc,ia2)
      zni = deni(ic,ig2)
      zdn = zni*zn0 * scfc(3)
      zsp = -zdn*dtr_mvp
      zwi = -zdn*dtr_eng
      sp_mt(ic,ia2) = sp_mt(ic,ia2) + zsp
      wi_mt(ic)     = wi_mt(ic)     + zwi
!-----
      tsp_mt(irc,irg) = tsp_mt(irc,irg) + zsp
      twi_mt(irc,irg) = twi_mt(irc,irg) + zwi
      enddo   !  loop (ig2)

!::elastic collision  [12]  "EL3"  irct = 12
!                         [ H0(ig) + H0(ig) => H0(ig) + H0(ig) ]
      irc = 12
      ia  = noia(ig)
      do ig2 = 1, ngas
      ia2 = noia(ig2)
      call ntsctr_el(irc,1)
      zni = den0BG(ic,1)
      zdn = zni*zn0 * scfc(irc)
      zsp = -zdn*dtr_mvp
      zwi = -zdn*dtr_eng
!      sp_mt(ic,1) = sp_mt(ic,1) + zsp
!      wi_mt(ic)     = wi_mt(ic)     + zwi
!-----
!      tsp_mt(irc,irg) = tsp_mt(irc,irg) + zsp
!      twi_mt(irc,irg) = twi_mt(irc,irg) + zwi
      enddo   !  loop (ig2)
!      write(n6,*) 'EL3 has been called. irc/scfc(irc)=',
!     &     irc,scfc(irc); call flush(n6)
!
!::elastic collision  [13]  "EL4"  irct = 13
!                         [ H0(ig) + H2(ig) => H0(ig) + H2(ig) ]
      irc = 13
      ia  = noia(ig)
      do ig2 = 1, ngas
      ia2 = noia(ig2)
      call ntsctr_el(irc,ia2)
      zni = dengBG(ic,ig2)
      zdn = zni*zn0 * scfc(irc)
      zsp = -zdn*dtr_mvp
      zwi = -zdn*dtr_eng
!      sp_mt(ic,ia2) = sp_mt(ic,ia2) + zsp
!      wi_mt(ic)     = wi_mt(ic)     + zwi
!-----
!      tsp_mt(irc,irg) = tsp_mt(irc,irg) + zsp
!      twi_mt(irc,irg) = twi_mt(irc,irg) + zwi
      enddo   !  loop (ig2)
!      write(n6,*) 'EL4 has been called. irc/scfc(irc)=',
!     &     irc,scfc(irc); call flush(n6)
!
!----------------------------------------------------------------------
!::D2
!----------------------------------------------------------------------
      elseif( mrct.eq.2 ) then
      ia  = noia(ig)
      zn0 = dltm*weit
!
!--irct = 4      "DS1"   [ e  + H2  =>  H0  + H0 + e ]
      irc = 4
      zdn = srct(1)*zn0 * scfc(4)
      zwe = -zdn*elrc(1)*cev
      we_mt(ic) = we_mt(ic) + zwe
!-----
      twe_mt(irc,irg) = twe_mt(irc,irg) + zwe
!
!---irct = 5      "DS2"   [ e  + H2  =>  H2+ + 2e ]
      irc = 5
      zdn = srct(2)*zn0 * scfc(5)
      zwe = -zdn*elrc(2)*cev
      we_mt(ic) = we_mt(ic) + zwe
!-----
      twe_mt(irc,irg) = twe_mt(irc,irg) + zwe
!
!---irct = 6      "DS3"   [ e  + H2  =>  H+  + H0 + e ]
      irc = 6
      zdn = srct(3)*zn0 * scfc(6)
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   zsn = zdn
      zwi = zdn*(0.5d0*tmas*vel2*0.5d0+derc(3)*cev)  ! tmas,vel2 H2
      zwe = -zdn*elrc(3)*cev
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   sn_mt(ic,ia) = sn_mt(ic,ia) + zsn
      sn_mt(ic,ia) = sn_mt(ic,ia) + zdn
      wi_mt(ic)    = wi_mt(ic)    + zwi
      we_mt(ic)    = we_mt(ic)    + zwe
!-----
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   tsn_mt(irc,irg) = tsn_mt(irc,irg) + zsn
      tsn_mt(irc,irg) = tsn_mt(irc,irg) + zdn
      twi_mt(irc,irg) = twi_mt(irc,irg) + zwi
      twe_mt(irc,irg) = twe_mt(irc,irg) + zwe
!
!::elactic collision  [7]  "EL2"  irct = 7
!                         [ H2(ig) + H+(ig2) => H2(ig) + H+(ig2) ]
      irc = 7
      ia  = noia(ig)
      do ig2 = 1, ngas
      ia2 = noia(ig2)
      call ntsctr_el(irc,ia2)
      zni = deni(ic,ig2)
      zdn = zni*zn0 * scfc(7)
      zsp = -zdn*dtr_mvp
      zwi = -zdn*dtr_eng
      sp_mt(ic,ia2) = sp_mt(ic,ia2) + zsp
      wi_mt(ic)     = wi_mt(ic)     + zwi
!-----
      tsp_mt(irc,irg) = tsp_mt(irc,irg) + zsp
      twi_mt(irc,irg) = twi_mt(irc,irg) + zwi
      enddo   !  loop (ig2)

!::elastic collision  [14]  "EL5"  irct = 14
!                         [ D2(ig) + D0(ig) => D2(ig) + D0(ig) ]
      irc = 14
      ia  = noia(ig)
      do ig2 = 1, ngas
      ia2 = noia(ig2)
      call ntsctr_el(irc,ia2)
      zni = den0BG(ic,ig2)
      zdn = zni*zn0 * scfc(irc)
      zsp = -zdn*dtr_mvp
      zwi = -zdn*dtr_eng
!      sp_mt(ic,ia2) = sp_mt(ic,ia2) + zsp
!      wi_mt(ic)     = wi_mt(ic)     + zwi
!-----
!      tsp_mt(irc,irg) = tsp_mt(irc,irg) + zsp
!      twi_mt(irc,irg) = twi_mt(irc,irg) + zwi
      enddo   !  loop (ig2)
!      write(n6,*) 'EL5 has been called. irc/scfc(irc)=',
!     &     irc,scfc(irc) !; call flush(n6)
!
!::elastic collision  [15]  "EL6"  irct = 15
!                         [ H2(ig) + H2(ig) => H2(ig) + H2(ig) ]
      irc = 15
      ia  = noia(ig)
      do ig2 = 1, ngas
      ia2 = noia(ig2)
      call ntsctr_el(irc,ia2)
      zni = dengBG(ic,ig2)
!      if(zni>1.d0) write(6,'(A6,4E12.4)') 'H2H2',xpos,ypos,scfc(irc),
!     &     dengBG(ic,ia2)
      zdn = zni*zn0 * scfc(irc)
      zsp = -zdn*dtr_mvp
      zwi = -zdn*dtr_eng
!      sp_mt(ic,ia2) = sp_mt(ic,ia2) + zsp
!      wi_mt(ic)     = wi_mt(ic)     + zwi
!-----
!      tsp_mt(irc,irg) = tsp_mt(irc,irg) + zsp
!      twi_mt(irc,irg) = twi_mt(irc,irg) + zwi
      enddo   !  loop (ig2)
!      write(n6,*) 'EL6 has been called. irc/scfc(irc)=',
!     &     irc,scfc(irc) !; call flush(n6)
!
!----------------------------------------------------------------------
!::D2+
!----------------------------------------------------------------------
      elseif( mrct.eq.3 ) then
      ia  = noia(ig)
      zn0 = dltm*weit
!
!--irct = 8      "DS4"   [ e  + H2+ =>  H+  + H0 + e ]
      irc = 8
      zdn = srct(1)*zn0 * scfc(8)
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   zsn = zdn
      zwi = zdn*(0.5d0*tmas*vel2*0.5d0+derc(1)*cev)
      zwe = -zdn*elrc(1)*cev
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   sn_mt(ic,ia) = sn_mt(ic,ia) + zsn
      sn_mt(ic,ia) = sn_mt(ic,ia) + zdn
      wi_mt(ic)    = wi_mt(ic)    + zwi
      we_mt(ic)    = we_mt(ic)    + zwe
!-----
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   tsn_mt(irc,irg) = tsn_mt(irc,irg) + zsn
      tsn_mt(irc,irg) = tsn_mt(irc,irg) + zdn
      twi_mt(irc,irg) = twi_mt(irc,irg) + zwi
      twe_mt(irc,irg) = twe_mt(irc,irg) + zwe
!
!--irct = 9      "DS5"   [ e  + H2+ =>  H+  + H+ + 2e ]
      irc = 9
      zdn = srct(2)*zn0 * scfc(9)
      zsn = zdn*2.0d0
      zwi = zdn*(0.5d0*tmas*vel2+derc(2)*2.0d0*cev)
      zwe = -zdn*elrc(2)*cev
      sn_mt(ic,ia) = sn_mt(ic,ia) + zsn
      wi_mt(ic)    = wi_mt(ic)    + zwi
      we_mt(ic)    = we_mt(ic)    + zwe
!-----
      tsn_mt(irc,irg) = tsn_mt(irc,irg) + zsn
      twi_mt(irc,irg) = twi_mt(irc,irg) + zwi
      twe_mt(irc,irg) = twe_mt(irc,irg) + zwe
!
!--irct = 10      "DS6"   [ e  + H2+ =>  H0  + H0 ]
      irc = 10
      zdn = srct(3)*zn0 * scfc(10)
      zwe = -zdn*elrc(3)*cev
      we_mt(ic) = we_mt(ic) + zwe
!-----
      twe_mt(irc,irg) = twe_mt(irc,irg) + zwe
!
!----------------------------------------------------------------------
!::D+  ???   mrct = 4   irct = 11
!----------------------------------------------------------------------
      elseif( mrct.eq.4 ) then
!
!--irct = 11      "RC"    [ e  + H+  =>  H0 ]
!          calculate source terms due to recombination at ntrecm
      irc = 11
      call wexit("ntscor_tr","wrong ktyp = 4")
!
      else
      call wexit("ntscor_tr","wrong ktyp ?")
      endif
!
      return
      end
