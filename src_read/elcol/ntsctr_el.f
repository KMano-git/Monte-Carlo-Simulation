!***********************************************************************
      subroutine ntsctr_el(irc,ipl)
!***********************************************************************
!
!     calculate scoring data
!           for elastic collision by track-length estimator
!
!    Note   momentum & energy change of test partilce
!
!           source terms in soldor
!               Sp// = - N0*Ni*dtr_mvp     = - N0*Ni*dcl_mvp
!               Swi  = - N0*Ni*dtr_eng     = - N0*Ni*dcl_eng
!
!               dtr_mvp = <dlt[m0*vp]*sig*v>   track estimator
!               dcl_mvp =  dlt[m0*vp]*<sig*v>  collision-estimator
!
!       irc :  index of reaction    see sub. ntcros
!                 = 3  (D0 + D+)  = 7  (D2 + D+)
!       ipl :  plasma ion
!
!     argument of I(Etm,Tim)
!           Etm = 1/2*mu*vt**2 = mu/mt*Et  [eV/amu]
!           Tim = 1/2*mu*vi**2 = mu/mi*Ti  [eV/amu]
!
!           Etm*[eV]/mu = 1/2*v0**2
!           Tim*[eV]/mu = Ti*[eV]/mi  in degas2-code
!
!        t : test particle (neutral)
!        i : ion  particle
!
!    Note.   Etm = 1/2*mu*{(velx-vfx)^2+(vely-vfy)^2+(velz-vfz)^2}/cev
!                                                      Mar/23/02
!-----------------------------------------------------------------------
      use celcom,  only : dtr_eng, dtr_mvp, dtr_mvx, dtr_mvy, dtr_mvz
     >    , el_i_10, el_i_11, el_i_12, el_rate, ptr_eng, ptr_mvp
     >    , ptr_mvx, ptr_mvy, ptr_mvz, ptr_sgv
      use cntcom,  only : bvx, bvy, bvz, icpo, rmas, tlmt_el, tmas, velx
     >    , vely, velz
      use cntpls,  only : vflw, temi
      use cphcns,  only : cev, cmp
      use cxdcom,  only : xwmlt
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer  irc,ipl
      integer, intent(in) :: irc, ipl
!
!::common varaibles
!x      real*8   tmas,rmas(2),evel,temi(2)
!x      real*8   velx,vely,velz,vflw(1,2),bvx(2),bvy(2),bvz(2)
!x      integer  icpo
!x      data  tmas/1.0/,evel/1.0/,icpo/1/,temi/2*1.0/
!x      data  velx,vely,velz/3*1.0/,vflw,bvx,bvy,bvz/8*1.0/
!
!::local variables
      integer  kel,i_10,i_11,i_12,i_rate
! modified 2/2 lines organize local variables and include files by kamata 2021/06/16
!ik   real*8   amt,ami,amu,amr,ame,etm,tim,tim2,tim3
!ik   real*8   vtx,vty,vtz,vfx,vfy,vfz,va2,va,utx,uty,utz,ut2,ut
      real*8   amt,ami,amu,amr,ame,etm,tim,tim2
      real*8   vfx,vfy,vfz,va2,utx,uty,utz,ut2,ut
      real*8   zuv,sp0
! modified 1/3 lines organize local variables and include files by kamata 2021/06/16
!ik   real*8   xdt_eval2, intg_10, intg_11, intg_12
      real*8   intg_10, intg_11, intg_12
! function
      real(8)    xdt_eval2
!
!-----------------------------------------------------------------------
!
!::set index
      if( irc.eq.3 ) then
         kel = 1
      elseif( irc.eq.7 ) then
         kel = 2
      elseif(( irc.eq.12) .OR. ( irc.eq.14) .OR.
     >       ( irc.eq.13) .OR. ( irc.eq.15) ) then
         dtr_mvx = 0.d0 ! -utx*amr*sp0
         dtr_mvy = 0.d0 ! -uty*amr*sp0
         dtr_mvz = 0.d0 ! -utz*amr*sp0
         dtr_mvp = 0.d0 ! dtr_mvx*bvx(icpo)+dtr_mvy*bvy(icpo)+dtr_mvz*bvz(icpo)
         dtr_eng = 0.d0 ! -ame*( 0.5d0*(amt+ami)*sp0*zuv-0.5d0*ami*intg_12 )
         return
      else
         call wexit("ntsctr_el","incorrect irc")
      endif
      i_10 = el_I_10(kel)
      i_11 = el_I_11(kel)
      i_12 = el_I_12(kel)
      i_rate = el_rate(kel)
!
!::mass
      amt = tmas
      ami = rmas(ipl)
      amu = cmp
      amr = amt*ami/(amt+ami)
      ame = 2.0d0*amt*ami/(amt+ami)**2
!
!::velocity
! deleted 3 lines organize local variables and include files by kamata 2021/06/16
!ik   vtx = velx
!ik   vty = vely
!ik   vtz = velz
      vfx = vflw(icpo,ipl)*bvx(icpo)
      vfy = vflw(icpo,ipl)*bvy(icpo)
      vfz = vflw(icpo,ipl)*bvz(icpo)
      va2 = 2.0d0*temi(icpo)*cev/ami
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   va  = dsqrt(va2)
!
!::velocity in the frame moving with flow velocity
! modified 3/3 lines organize local variables and include files by kamata 2021/06/16
!ik   utx = vtx - vfx
!ik   uty = vty - vfy
!ik   utz = vtz - vfz
      utx = velx - vfx
      uty = vely - vfy
      utz = velz - vfz
      ut2 = utx**2+uty**2+utz**2
      ut  = dsqrt(ut2)
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   zuv = utx*vtx+uty*vty+utz*vtz
      zuv = utx*velx+uty*vely+utz*velz
!
!::ena,tmb
!x    etm = amu/amt*evel        !  error
      etm = 0.5d0*amu*ut2/cev
      tim = amu/ami*temi(icpo)
      tim2 = dmax1( tim, tlmt_el )    !  Sp
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   tim3 = tim                      !  Wi heating at low Ti
!
!::integration  (cgs => MKS unit)
      intg_10 = xdt_eval2(i_10,etm,tim2)*xwmlt(1,i_10)
      intg_11 = xdt_eval2(i_11,etm,tim2)*xwmlt(1,i_11)
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   intg_12 = xdt_eval2(i_12,etm,tim3)*xwmlt(1,i_12)
      intg_12 = xdt_eval2(i_12,etm,tim )*xwmlt(1,i_12)
!
!::momentum & energy change of test particle
      sp0 = intg_11/ut - 0.5d0*va2/ut2*intg_10
!
      dtr_mvx = -utx*amr*sp0
      dtr_mvy = -uty*amr*sp0
      dtr_mvz = -utz*amr*sp0
      dtr_mvp = dtr_mvx*bvx(icpo)+dtr_mvy*bvy(icpo)+dtr_mvz*bvz(icpo)
      dtr_eng = -ame*( 0.5d0*(amt+ami)*sp0*zuv-0.5d0*ami*intg_12 )
!
      return
!
!-----------------------------------------------------------------------
!::debug check
      ptr_sgv = xdt_eval2(i_rate,etm,tim)*xwmlt(1,i_rate)
      ptr_mvx = dtr_mvx/ptr_sgv
      ptr_mvy = dtr_mvy/ptr_sgv
      ptr_mvz = dtr_mvz/ptr_sgv
      ptr_mvp = dtr_mvp/ptr_sgv
      ptr_eng = dtr_eng/ptr_sgv
!
!-----< test_elsrc >
!x      write(n6,'(a4,i8,3i5,1p13e11.2)')
!x     >  "trac",iptl,istp,iplx(icpo),iply(icpo)
!x     > ,amt,ptr_sgv,vtx,ptr_mvx/amt,vty,ptr_mvy/amt,vtz,ptr_mvz/amt
!x     > ,vtx*bvx(icpo)+vty*bvy(icpo)+vtz*bvz(icpo),ptr_mvp/amt
!x     > ,0.5d0*(vtx**2+vty**2+vtz**2),ptr_eng/amt
!
      return
      end
